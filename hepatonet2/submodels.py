"""
Create SBML models for the available submodels
"""
import os
import json
from collections import defaultdict

import libsbml

dir_cur = os.path.dirname(os.path.realpath(__file__))
repo_dir = os.path.join(dir_cur, "..", "repository")
models_dir = os.path.join(dir_cur, "..", "models")
input_dir = os.path.join(dir_cur, "..", "input")


def parse_submodels():
    """ Parses the subsystem information.

    :return:
    """
    subsystems = defaultdict(list)
    with open(os.path.join(repo_dir, "reactions.json"), "r") as f_reactions:
        reactions = json.load(f_reactions)
        for rid, r in reactions.items():
            subsystem = r.get("subsystems")
            subsystems[subsystem].append(rid)

        with open(os.path.join(repo_dir, "subsystems.json"), "w") as f_subsystems:
            json.dump(subsystems, f_subsystems, sort_keys=True, indent=2)
    return subsystems

def _normalize_to_sid(text):
    """ Process text to be SID conform. """
    # FIXME: use regular expressions for multiple replacements
    text = text.replace(" ", '_')
    text = text.replace("/", '_')
    text = text.replace("-", "_")
    text = text.lower()
    return text


def create_sbml_for_subsystem(subsystem, organism):
    """ Creates the SBML model for the given subsystem.

    :param subsystem:
    :param organism:
    :return:
    """
    # load RECON3D reference
    path_recon = os.path.join(input_dir, "models", "recon3d", "Recon3D.xml")
    print(path_recon)
    doc_recon = libsbml.readSBMLFromFile(path_recon)  # type: libsbml.SBMLDocument
    model_recon = doc_recon.getModel()  # type: libsbml.Model


    doc = libsbml.SBMLDocument(3, 2)  # type: libsbml.SBMLDocument
    # TODO: add fbc package
    model = doc.createModel()  # type: libsbml.Model
    model.setId(_normalize_to_sid(subsystem))

    # reaction ids for subsystem
    reaction_ids = set()
    compartment_ids = set()
    species_ids = set()
    geneproduct_ids = set()

    with open(os.path.join(repo_dir, "reactions.json"), "r") as f_reactions:
        reactions = json.load(f_reactions)

        for rid, rinfo in reactions.items():
            if rinfo.get("subsystems") == subsystem:
                reaction_ids.add(rid)
                r_recon = model_recon.getReaction(rid)  # type: libsbml.Reaction
                for s in r_recon.getListOfReactants():  # type: libsbml.SpeciesReference
                    species_ids.add(s.getSpecies())
                for s in r_recon.getListOfProducts():
                    species_ids.add(s.getSpecies())
                for s in r_recon.getListOfModifiers():
                    species_ids.add(s.getSpecies())

                # parse the geneProduct information

    # compartments
    for sid in species_ids:
        c_recon = model_recon.getSpecies(sid)
        compartment_ids.add(c_recon.getId())

    # add compartments in model
    for cid in compartment_ids:
        c_recon = model_recon.getCompartment(cid)  # type: libsbml.Compartment
        c = model.createCompartment()  # type: libsbml.Compartment
        c.setId(c_recon.getId())
        c.setConstant(c_recon.getConstant())
        c.setMetaId(c_recon.getMetaId())
        c.setName(c_recon.getName())
        c.setSBOTerm(c_recon.getSBOTerm())
        # c.setUnits(c_recon.getUnits())  # FIXME (no units so far)
        c.setAnnotation(c_recon.getAnnotation())

    # add species in model
    for sid in species_ids:
        s_recon = model_recon.getSpecies(sid)  # type: libsbml.Species
        s = model.createSpecies()  # type: libsbml.Species
        s.setId(s_recon.getId())
        s.setBoundaryCondition(s_recon.getBoundaryCondition())
        s.setCompartment(s_recon.getCompartment())
        s.setConstant(s_recon.getConstant())
        s.setHasOnlySubstanceUnits(s_recon.getHasOnlySubstanceUnits())
        s.setName(s_recon.getName())

        # fbc
        s_recon_fbc = s_recon.getPlugin("fbc")  # type: libsbml.FbcSpeciesPlugin
        s_fbc = s.getPlugin("fbc")  # type: libsbml.FbcSpeciesPlugin
        s_fbc.setCharge(s_recon_fbc.getCharge())
        s_fbc.setChemicalFormula(s_recon_fbc.getChemicalFormula())






    with open(os.path.join(repo_dir, "reactions.json"), "r") as f_reactions, \
         open(os.path.join(repo_dir, "species.json"), "r") as f_species, \
         open(os.path.join(repo_dir, "substances.json"), "r") as f_substances:
        reactions = json.load(f_reactions)
        species = json.load(f_species)
        substances = json.load(f_substances)

        for rid, rinfo in reactions.items():
            if rinfo.get("subsystems") == subsystem:
                # create reaction
                print(rid)
                r = model.createReaction()  # type: libsbml.Reaction
                recon_id = rinfo.get('recon_id')
                r.setId("R_{}".format(recon_id))
                recon_name = rinfo.get('recon_name')
                if recon_name:
                    r.setName(recon_name)
                reversible = True
                reversibility = rinfo.get('reversibility')
                if reversibility and reversibility == "=>":
                    reversible = False
                r.setReversible(reversible)

                # add species for the reactions
                for stoichiometry, species_id in rinfo['rxn_left']:
                    print(stoichiometry, species_id)

                    specie = species[species_id]
                    compartment_id = specie['compartment']

                for stoichiometry, species_id in rinfo['rxn_right']:
                    print(stoichiometry, species_id)

                # add compartments for the species

                # add reactions which belong to the subsystem

                # add gene associations for the subsystem

                # add annotations


    return doc


if __name__ == "__main__":
    subsystems_dict = parse_submodels()

    # Create submodels
    subsystems = ['Glycolysis/gluconeogenesis']
    # subsystems = subsystems_dict.keys()
    organism = "human"
    for subsystem in subsystems:
        doc = create_sbml_for_subsystem(subsystem=subsystem, organism=organism)
        sbml_path = os.path.join(models_dir, "subsystems", "{}_{}.xml".format(organism, _normalize_to_sid(subsystem)))
        libsbml.writeSBMLToFile(doc, sbml_path)
        print(sbml_path)


