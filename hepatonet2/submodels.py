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


def parse_submodels_recon3d():
    """ Parses the subsystem information from SBML file

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

    # parse subsystems
    subsystems = defaultdict(list)
    model_recon_groups = model_recon.getPlugin("groups")  # type: libsbml.GroupsModelPlugin
    for group in model_recon_groups.getListOfGroups():  # type: libsbml.Group
        key = group.getName()
        members = []
        for member in group.getListOfMembers():  # type: libsbml.Member
            members.append(member.getIdRef())
        subsystems[key] = members

    sbmlns = libsbml.SBMLNamespaces(3, 2, 'fbc', 2)
    doc = libsbml.SBMLDocument(sbmlns)  # type: libsbml.SBMLDocument
    model = doc.createModel()  # type: libsbml.Model
    model.setId(_normalize_to_sid(subsystem))

    # reaction ids for subsystem
    reaction_ids = set(subsystems[subsystem])
    compartment_ids = set()
    species_ids = set()
    geneproducts = dict()
    print(reaction_ids)


    # GPR rules
    def process_association(association, gene_products=None):
        """ Recursively get gene products from GPR. """
        if gene_products is None:
            gene_products = []

        if association.isFbcOr() or association.isFbcAnd():
            for c in association.getListOfAssociations():
                process_association(c, gene_products)
        elif association.isGeneProductRef():
            gid = association.getGeneProduct()
            gene_products.append(gid)

        return gene_products

    # species & gene products from reaction
    for rid in reaction_ids:
        r_recon = model_recon.getReaction(rid)  # type: libsbml.Reaction
        for s in r_recon.getListOfReactants():  # type: libsbml.SpeciesReference
            species_ids.add(s.getSpecies())
        for s in r_recon.getListOfProducts():
            species_ids.add(s.getSpecies())
        for s in r_recon.getListOfModifiers():
            species_ids.add(s.getSpecies())

        # parse the geneProduct information
        r_recon_fbc = r_recon.getPlugin("fbc")  # type: libsbml.FbcReactionPlugin
        if r_recon_fbc:
            gpa = r_recon_fbc.getGeneProductAssociation()  # type: libsbml.GeneProductAssociation
            if gpa is not None:
                association = gpa.getAssociation()  # type: libsbml.FbcAssociation
                # list of gene products
                geneproducts[rid] = process_association(association)

    print(species_ids)

    # compartments
    for sid in species_ids:
        s_recon = model_recon.getSpecies(sid)  # type: libsbml.Species
        c_recon = model_recon.getCompartment(s_recon.getCompartment())
        compartment_ids.add(c_recon.getId())
    print(compartment_ids)

    # add compartments
    for cid in compartment_ids:
        c_recon = model_recon.getCompartment(cid)  # type: libsbml.Compartment
        c = model.createCompartment()  # type: libsbml.Compartment
        c.setId(c_recon.getId())
        c.setConstant(c_recon.getConstant())
        c.setMetaId(c_recon.getMetaId())
        c.setName(c_recon.getName())
        c.setSBOTerm(c_recon.getSBOTerm())
        # c.setUnits(c_recon.getUnits())  # FIXME (no units so far), register in listOfUnitDefinitions
        c.setAnnotation(c_recon.getAnnotation())

    # add species
    for sid in species_ids:
        s_recon = model_recon.getSpecies(sid)  # type: libsbml.Species
        s = model.createSpecies()  # type: libsbml.Species
        s.setId(s_recon.getId())
        s.setMetaId(s_recon.getMetaId())
        s.setBoundaryCondition(s_recon.getBoundaryCondition())
        s.setCompartment(s_recon.getCompartment())
        s.setConstant(s_recon.getConstant())
        s.setHasOnlySubstanceUnits(s_recon.getHasOnlySubstanceUnits())
        s.setName(s_recon.getName())
        s.setSBOTerm(s.getSBOTerm())
        s.setAnnotation(s_recon.getAnnotation())

        # fbc
        s_recon_fbc = s_recon.getPlugin("fbc")  # type: libsbml.FbcSpeciesPlugin
        s_fbc = s.getPlugin("fbc")  # type: libsbml.FbcSpeciesPlugin
        s_fbc.setCharge(s_recon_fbc.getCharge())
        s_fbc.setChemicalFormula(s_recon_fbc.getChemicalFormula())

    # gene products
    model_fbc = model.getPlugin("fbc")  # type: libsbml.FbcModelPlugin
    model_recon_fbc = model_recon.getPlugin("fbc")  # type: libsbml.FbcModelPlugin
    for rid, gpids in geneproducts.items():
        for gid in gpids:
            gp_recon = model_recon_fbc.getGeneProduct(gid)  # type: libsbml.GeneProduct
            gp = model_fbc.createGeneProduct()  # type: libsbml.GeneProduct
            gp.setId(gp_recon.getId())
            gp.setLabel(gp_recon.getLabel())
            gp.setMetaId(gp_recon.getMetaId())
            gp.setAnnotation(gp_recon.getAnnotation())

            # TODO: add the ensg information (from external file)

    # add reactions
    for rid in reaction_ids:
        r_recon = model_recon.getReaction(rid)  # type: libsbml.Reaction
        r = model.createReaction()  # type: libsbml.Reaction
        r.setId(r_recon.getId())
        r.setMetaId(r_recon.getMetaId())
        r.setAnnotation(r_recon.getAnnotation())
        r.setName(r_recon.getName())
        r.setSBOTerm(r_recon.getSBOTerm())
        r.setFast(r_recon.getFast())
        r.setReversible(r_recon.getReversible())

        # fbc
        r_recon_fbc = r_recon.getPlugin("fbc")  # type: libsbml.FbcReactionPlugin
        r_fbc = r.getPlugin("fbc")  # type: libsbml.FbcReactionPlugin
        r_fbc.setLowerFluxBound(r_recon_fbc.getLowerFluxBound())
        r_fbc.setUpperFluxBound(r_recon_fbc.getUpperFluxBound())

        for sref_recon in r_recon.getListOfReactants():  # type: libsbml.SpeciesReference
            sid = sref_recon.getSpecies()
            s = model.getSpecies(sid)
            r.addReactant(s)  # type: libsbml.SpeciesReference
            sref = r.getReactant(sid)
            sref.setConstant(sref_recon.getConstant())
            sref.setStoichiometry(sref_recon.getStoichiometry())
            sref.setSBOTerm(sref_recon.getSBOTerm())

        for sref_recon in r_recon.getListOfProducts():  # type: libsbml.SpeciesReference
            sid = sref_recon.getSpecies()
            s = model.getSpecies(sid)
            r.addProduct(s)  # type: libsbml.SpeciesReference
            sref = r.getProduct(sid)
            sref.setConstant(sref_recon.getConstant())
            sref.setStoichiometry(sref_recon.getStoichiometry())
            sref.setSBOTerm(sref_recon.getSBOTerm())

        # FIXME: implement modifiers on reactions
        '''          
        for sref_recon in r_recon.getListOfModifiers():  # type: libsbml.SpeciesReference
            s = model.getSpecies(sref_recon.getSpecies())
            r.addModifier(s)  # type: libsbml.SpeciesReference
            sref = r.getModifier(s.getId())
            sref.setSBOTerm(sref_recon.getSBOTerm())
        '''

        # add GPA information
        # TODO: implement


    # mapping to mouse
    # TODO: implement

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


