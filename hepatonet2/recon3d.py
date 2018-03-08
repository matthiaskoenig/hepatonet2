"""
Importing data from RECON3D in the database.

Loads reactions, species, genes and gene-reactions association in the
internal data base.

"""
import os
import json
import scipy.io
from pprint import pprint
from collections import defaultdict
import numpy as np


def loadmat(filename):
    """
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    """
    data = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def _check_keys(dict):
    """
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    """
    for key in dict:
        if isinstance(dict[key], scipy.io.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict


def _todict(matobj):
    """
    A recursive function which constructs from matobjects nested dictionaries
    """
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, scipy.io.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

# Read the Matlab files and store in the repository

import re

def _value_by_idx(array, idx):
    """ Returns the value or None. """
    value = array[idx]
    if isinstance(value, np.ndarray) and len(value) == 0:
        return None
    if isinstance(value, (np.int, np.uint, np.int64, np.int16, np.uint8)):
        return int(value)
    return value


def parse_species(model, repo_dir):
    """ Creates the species files from the given mat model.

    :param model:
    :return:
    """
    substances = defaultdict(dict)
    species = defaultdict(dict)
    compartments = set()
    pattern = r'(.*)\[(.*)\]'

    metCharges = model['metCharges']
    metCHEBIID = model['metCHEBIID']
    metFormulas = model['metFormulas']
    metHMDBID = model['metHMDBID']
    metInChIString = model['metInChIString']
    metKEGGID = model['metKEGGID']
    metNames = model['metNames']
    metPdMap = model['metPdMap']
    metPubChemID = model['metPubChemID']
    metReconMap = model['metReconMap']
    metSmiles = model['metReconMap']

    for idx, met in enumerate(model["mets"]):
        # parse substance and compartment
        match = re.match(pattern, met)
        substance_id = match.group(1)
        compartment_id = match.group(2)
        compartments.add(compartment_id)

        # create substance
        if not substance_id in substances:
            s = dict()
            s['reconid'] = substance_id
            s['charge'] = _value_by_idx(metCharges, idx)
            s['formula'] = _value_by_idx(metFormulas, idx)
            s['chebi'] = _value_by_idx(metCHEBIID, idx)
            s['hmdb'] = _value_by_idx(metHMDBID, idx)
            s['inchi'] = _value_by_idx(metInChIString, idx)
            s['kegg'] = _value_by_idx(metKEGGID, idx)
            s['reconname'] = _value_by_idx(metNames, idx)
            s['pubchem'] = _value_by_idx(metPubChemID, idx)
            s['smilex'] = _value_by_idx(metSmiles, idx)
            s = {k: v for k, v in s.items() if v is not None}
            substances[substance_id] = s

        # create species
        m = dict()
        m['reconid'] = met
        m['substance'] = substance_id
        m['compartment'] = compartment_id
        m['pdmap'] = _value_by_idx(metPdMap, idx)
        m['reconmap'] = _value_by_idx(metReconMap, idx)
        m = {k: v for k, v in m.items() if v is not None}
        species[met] = m

    return substances, species, compartments


def parse_reactions(model):
    """ Creates the species files from the given mat model.

    :param model:
    :return:
    """
    reactions = defaultdict(dict)
    bounds = defaultdict(dict)
    exchanges = None
    sinks = None
    compartments = set()
    # pattern = r'([A-Z|0-9|_]*)(.*)'

    ub = model['ub']
    lb = model['lb']
    rxnCOG = model['rxnCOG']
    rxnConfidenceScores = model['rxnConfidenceScores']
    rxnECNumbers = model['rxnECNumbers']
    rxnKEGGID = model['rxnKEGGID']
    rxnKeggOrthology = model['rxnKeggOrthology']
    rxnNames = model['rxnNames']
    rxnNotes = model['rxnNotes']
    rxnReconMap = model['rxnReconMap']
    rxnReferences = model['rxnReferences']
    subSystems = model['subSystems']

    mets = model['mets']
    S = model['S']

    for idx, rxn in enumerate(model["rxns"]):
        # parse substance and compartment
        # match = re.match(pattern, rxn)
        # reaction_id = match.group(1)
        # compartment_id = match.group(2)
        # print(rxn, reaction_id, compartment_id)

        # store bounds
        bounds[rxn] = {
            'ub': _value_by_idx(ub, idx),
            'lb': _value_by_idx(lb, idx),
        }

        # create species
        r = dict()
        r['recon_id'] = rxn
        r['cog'] = _value_by_idx(rxnCOG, idx)
        r['recon_confidence'] = _value_by_idx(rxnConfidenceScores, idx)
        r['ec'] = _value_by_idx(rxnECNumbers, idx)
        r['kegg'] = _value_by_idx(rxnKEGGID, idx)
        r['kegg_orthology'] = _value_by_idx(rxnKeggOrthology, idx)
        r['recon_name'] = _value_by_idx(rxnNames, idx)
        r['recon_notes'] = _value_by_idx(rxnNotes, idx)
        r['recon_map'] = _value_by_idx(rxnReconMap, idx)
        r['recon_references'] = _value_by_idx(rxnReferences, idx)
        r['subsystems'] = _value_by_idx(subSystems, idx)

        # create reaction equation (irreversibility via bounds)
        col = S[:, idx]
        s_indices = np.nonzero(col)[0]
        left, right = [], []

        for s_idx in s_indices:
            stoichiometry = float(col[s_idx].data[0])
            met = mets[s_idx]
            if stoichiometry < 0:
                left.append((stoichiometry, met))
            else:
                right.append((stoichiometry, met))

        reversibility = "<=>"
        if bounds[rxn]['ub'] == 0:
            reversibility = "<="
        if bounds[rxn]['lb'] == 0:
            reversibility = "=>"
        if bounds[rxn]['lb'] == 0 and bounds[rxn]['ub'] == 0:
            reversibility = "|"

        # print(left, reversibility, right)
        r['rxn_left'] = left
        r['rxn_right'] = right
        r['rxn_reversibility'] = reversibility
        r = {k: v for k, v in r.items() if v is not None}
        reactions[rxn] = r

    return reactions, bounds


def _value_by_key(row, key):
    """ Returns the value or None from given DataFrame row """
    value = row[key].iloc[0]
    if pd.isna(value):
        return None
    return value


def parse_genes(model, info_df):
    # rules
    recon_genes = model['genes']
    genes = {}
    for gid in recon_genes:
        if gid == "0":
            print("invalid gene:", gid)
            continue
        entrez_id, entrez_version = [int(t) for t in gid.split('.')]

        row = info_df[info_df["EntrezGene ID"] == entrez_id]

        g = dict()
        g["entrez"] = entrez_id
        g["entrez_version"] = entrez_version

        if len(row) > 0:
            g["ensemble"] = _value_by_key(row, "Ensembl Gene ID")
            g["hgnc"] = _value_by_key(row, "HGNC symbol")
            g["description"] = _value_by_key(row, "Description")
            g["mim"] = _value_by_key(row, "MIM Gene Accession")
            g["uniprot"] = _value_by_key(row, "UniProt/SwissProt ID")
            g["hpa"] = _value_by_key(row, "Human Protein Atlas Antibody ID")
            g["go"] = _value_by_key(row, "GO Term Accession")
            g["go_name"] = _value_by_key(row, "GO Term Name")
            g["go_definition"] = _value_by_key(row, "GO Term Definition")

        g = {k: v for k, v in g.items() if v is not None}
        genes[gid] = g

    return genes


def parse_gene_associations(model, outdir):
    """ Parses the gene association rules.

    :param model:
    :param outdir:
    :return:
    """
    rules = model['rules']
    genes = model['genes']
    rxns = model['rxns']

    associations = dict()
    pattern = re.compile(r'x\([0-9]+\)')

    for k in range(len(rules)):
        rule = _value_by_idx(rules, k)
        rxn_id = rxns[k]
        if rule:
            # replace gene ids in rule
            matches = re.findall(pattern, rule)
            for match in matches:
                idx = int(match[2:(len(match)-1)]) - 1  # zero indexed python, one indexed matlab
                rule = rule.replace(match, genes[idx])

            associations[rxn_id] = rule

    # rules
    return associations


if __name__ == "__main__":
    dir_cur = os.path.dirname(os.path.realpath(__file__))
    dir_recon = os.path.join(dir_cur, '..', 'input', 'models', 'recon3d')
    RECON3D_MODEL_MAT = os.path.join(dir_recon, "Recon3DModel_301.mat")
    RECON3D_MAT = os.path.join(dir_recon, "Recon3D_301.mat")
    RECON3D_S4 = os.path.join(dir_recon, "nbt.4072-S4.xlsx")
    print(RECON3D_S4)

    # Create a parts repository from given model
    repo_dir = os.path.join(dir_cur, "..", "repository")

    mat = loadmat(RECON3D_MODEL_MAT)
    model = mat["Recon3DModel"]

    if False:
        print("*** SPECIES ***")
        substances, species, compartments = parse_species(model, repo_dir)

        print("substances:", len(substances))
        with open(os.path.join(repo_dir, "substances.json"), "w") as f:
            json.dump(substances, f, sort_keys=True, indent=2)
        print("species:", len(species))
        with open(os.path.join(repo_dir, "species.json"), "w") as f:
            json.dump(species, f, sort_keys=True, indent=2)

        print("*** REACTIONS ***")
        reactions, bounds = parse_reactions(model, repo_dir)
        print("reactions:", len(reactions))
        with open(os.path.join(repo_dir, "reactions.json"), "w") as f:
            json.dump(reactions, f, sort_keys=True, indent=2)

        print("*** GENES ***")
        # additional gene information
        # "nbt.4072-S4.xlsx" "Supplement Data File 8"
        # biomart services
        import pandas as pd
        xls = pd.ExcelFile(RECON3D_S4)

        info_df = xls.parse("Supplement Data File 8", skiprows=1)
        # print(info_df.head())
        genes = parse_genes(model, info_df)
        print("genes:", len(genes))
        with open(os.path.join(repo_dir, "genes", "human_genes.json"), "w") as f:
            json.dump(genes, f, sort_keys=True, indent=2)

    print("*** GENE ASSOCIATIONS ***")
    associations = parse_gene_associations(model, repo_dir)
    print("associations:", len(associations))
    with open(os.path.join(repo_dir, "genes", "human_gene_associations.json"), "w") as f:
        json.dump(associations, f, sort_keys=True, indent=2)







