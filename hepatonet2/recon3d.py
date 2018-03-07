"""
Importing data from RECON3D in the database.

Loads reactions, species, genes and gene-reactions association in the
internal data base.

"""
import os
import scipy.io
from pprint import pprint


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


def create_species(model, out_dir):
    """ Creates the species files from the given mat model.

    :param model:
    :return:
    """

    # metCharges
    # metCHEBIID
    # metFormulas
    # metHMDBID
    # metInCHIString
    # metKEGGId
    # metNames
    # metPdMap
    # metPubChemID
    # metReconMap
    # metSmiles
    # mets

    # compartmentsrule


def create_reactions(model, out_dir):
    """ Creates the species files from the given mat model.

    :param model:
    :return:
    """
    # rxns

    # parse stoichiometry (S)
    # ub, lb

    # rxnCOG
    # rxnConfidenceScores
    # rxnECNumbers
    # rxnKEGGID
    # rxnKeggOrthology
    # rxnNames
    # rxnNotes
    # rxnReconMap
    # rxnReferences
    # subsystems

    # compartments


def create_genes(model, out_dir):
    # rules

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

    mat = loadmat(RECON3D_MODEL_MAT)
    model = mat["Recon3DModel"]
    pprint(list(model.keys()))

    # additional gene information
    "nbt.4072 - S4.xlsx" "Supplement Data File 8"
    # biomart services




