"""
Create SBML models for the available submodels
"""
import os
import json
from collections import defaultdict

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

if __name__ == "__main__":
    parse_submodels()