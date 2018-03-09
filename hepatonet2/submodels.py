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

def parse_submodels():

    subsystems = defaultdict(list)
    with open(os.path.join(repo_dir, "reactions.json"), "r") as f_reactions:
        reactions = json.load(f_reactions)
        for rid, r in reactions.items():
            subsystem = r.get("subsystems")
            subsystems[subsystem].append(rid)

        with open(os.path.join(repo_dir, "subsystems.json"), "w") as f_subsystems:
            json.dump(subsystems, f_subsystems, sort_keys=True, indent=2)


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
    doc = libsbml.SBMLDocument(3, 2)  # type: libsbml.SBMLDocument
    # TODO: add fbc package
    model = doc.createModel()
    model.setId(_normalize_to_sid(subsystem))

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


                for stoichiometry, species_id in rinfo['rxn_left']:
                    print(stoichiometry, species_id)

                # add compartments for the species

                # add reactions which belong to the subsystem

                # add gene associations for the subsystem

    return doc


if __name__ == "__main__":
    parse_submodels()

    # Create submodels
    subsystems = ['Glycolysis/gluconeogenesis']
    organism = "human"
    for subsystem in subsystems:
        doc = create_sbml_for_subsystem(subsystem=subsystem, organism=organism)
        sbml_path = os.path.join(models_dir, "subsystems", "{}_{}.xml".format(organism, _normalize_to_sid(subsystem)))
        libsbml.writeSBMLToFile(doc, sbml_path)
        print(sbml_path)


