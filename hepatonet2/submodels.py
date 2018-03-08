"""
Create SBML models for the available submodels
"""
import os
import json
from collections import defaultdict

import libsbml

dir_cur = os.path.dirname(os.path.realpath(__file__))
repo_dir = os.path.join(dir_cur, "..", "repository")

def parse_submodels():

    subsystems = defaultdict(list)
    with open(os.path.join(repo_dir, "reactions.json"), "r") as f_reactions:
        reactions = json.load(f_reactions)
        for rid, r in reactions.items():
            subsystem = r.get("subsystems")
            subsystems[subsystem].append(rid)

        with open(os.path.join(repo_dir, "subsystems.json"), "w") as f_subsystems:
            json.dump(subsystems, f_subsystems, sort_keys=True, indent=2)


def create_sbml_submodels(subsystem, organism):
    doc = libsbml.SBMLDocument(3, 1)  # type: libsbml.SBMLDocument
    model = doc.createModel()
    mid = subsystem.replace(" ", '_')
    mid = mid.replace("-", "_")

    model.setId(mid)

    with open(os.path.join(repo_dir, "subsystems.json"), "r") as f_reactions:
        reactions = json.load(f_reactions)

        for rid, r in reactions.items():
            if r.get("subsystems") == subsystem:

            # add reactions which belong to the subsystem
            


if __name__ == "__main__":
    parse_submodels()