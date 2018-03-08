"""
Importing data from RECON3D in the database.

Loads reactions, species, genes and gene-reactions association in the
internal data base.

"""
import os
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
    if isinstance(value, np.int64):
        return int(value)
    return value


def create_species(model, repo_dir):
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


def create_reactions(model, repo_dir):
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
    print(type(S))
    print(S.shape)

    for idx, rxn in enumerate(model["rxns"]):
        # parse substance and compartment
        # match = re.match(pattern, rxn)
        # reaction_id = match.group(1)
        # compartment_id = match.group(2)
        # print(rxn, reaction_id, compartment_id)

        # store bounds
        bounds[rxn] = {
            'ub': ub[idx],
            'lb': lb[idx]
        }

        # create species
        reactions[rxn] = {
            'recon_id': rxn,
            'cog': rxnCOG[idx],
            'recon_confidence': rxnConfidenceScores[idx],
            'ec': rxnECNumbers[idx],
            'kegg': rxnKEGGID[idx],
            'kegg_orthology': rxnKeggOrthology[idx],
            'recon_name': rxnNames[idx],
            'recon_notes': rxnNotes[idx],
            'recon_map': rxnReconMap[idx],
            'recon_references': rxnReferences[idx],
            'subsystems': subSystems[idx],
        }

        # create reaction equation (irreversibility via bounds)
        col = S[:, idx]
        s_indices = np.nonzero(col)[0]
        left, right = [], []

        for s_idx in s_indices:
            stoichiometry = col[s_idx].data[0]
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
        reactions['rxn_left'] = left
        reactions['rxn_right'] = right
        reactions['rxn_reversibility'] = reversibility

    return reactions, bounds


def create_genes(model, out_dir):
    # rules
    pass


def create_gene_associations(model, outdir):
    """ Parses the gene association rules.

    :param model:
    :param outdir:
    :return:
    """
    # rules
    pass


if __name__ == "__main__":
    dir_cur = os.path.dirname(os.path.realpath(__file__))
    dir_recon = os.path.join(dir_cur, '..', 'input', 'models', 'recon3d')
    RECON3D_MODEL_MAT = os.path.join(dir_recon, "Recon3DModel_301.mat")
    RECON3D_MAT = os.path.join(dir_recon, "Recon3D_301.mat")

    # Create a parts repository from given model
    repo_dir = os.path.join(dir_cur, "..", "repository")

    mat = loadmat(RECON3D_MODEL_MAT)
    model = mat["Recon3DModel"]

    print("*** SPECIES ***")
    substances, species, compartments = create_species(model, repo_dir)
    # pprint(substances)
    # pprint(species)
    # pprint(compartments)

    # TODO: store information as JSON
    import json
    print("substances:", len(substances))
    with open(os.path.join(repo_dir, "substances.json"), "w") as f:
        json.dump(substances, f, sort_keys=True, indent=2)
    print("species:", len(species))
    with open(os.path.join(repo_dir, "species.json"), "w") as f:
        json.dump(species, f, sort_keys=True, indent=2)


    exit()

    print("*** REACTIONS ***")
    reactions, bounds = create_reactions(model, repo_dir)
    # pprint(bounds)
    # pprint(reactions)

    print("*** GENES ***")
    reactions, bounds = create_reactions(model, repo_dir)

    print("*** GENE ASSOCIATIONS ***")
    reactions, bounds = create_reactions(model, repo_dir)

    # additional gene information
    # "nbt.4072 - S4.xlsx" "Supplement Data File 8"
    # biomart services




