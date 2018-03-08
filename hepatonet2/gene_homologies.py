"""
Query homologs from ensmble for given genes.

http://ensemblgenomes.org/info/access/rest

searches the homolog genes
https://rest.ensembl.org/homology/id/ENSG00000106633?content-type=application/json&target_species=mouse
"""
import os
import json
import requests
import json
import time
from pprint import pprint


class EnsemblRestClient(object):
    """ Client to query ensemble. """

    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint):

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        url = self.server + endpoint
        print(url)
        response = requests.get(url)
        self.req_count += 1

        return response

    def get_homologs(self, symbol, species="human"):

        endpoint = '/homology/id/{}/?content-type=application/json&target_species={}'.format(symbol, species)

        response = self.perform_rest_action(endpoint)
        if response.status_code == 200:
            json = response.json()
            data = json['data']
            homologies = data[0]['homologies']
            return homologies
        else:
            print(endpoint)
            print(response)
            return None

    def get_xrefs(self, symbol):
        """ Cross references. """
        # https://rest.ensembl.org/xrefs/id/ENSG00000106633?content-type=application/json
        # https://rest.ensembl.org/xrefs/id/ENSG00000106633?content-type=application/json;external_db=Uniprot_gn
        endpoint = '/xrefs/id/{}?content-type=application/json'.format(symbol)

        response = self.perform_rest_action(endpoint)
        if response.status_code == 200:
            json = response.json()
            return json
        else:
            print(endpoint)
            print(response)
            return None


def get_homologies(symbols, out_dir, species=['mouse', 'rat']):
    """ Get all homologies for given gene symbols in the species

    :param symbols: List of symbols
    :param species:
    :return:
    """
    client = EnsemblRestClient()

    Nsymbols = len(symbols)
    for s in species:
        homologies = {}

        for k, symbol in enumerate(symbols):
            print("[{}|{}] {}".format(k, Nsymbols, symbol))
            homologies[symbol] = client.get_homologs(symbol=symbol, species=s)

        with open(os.path.join(out_dir, "{}_homologies.json".format(s)), "w") as f:
            json.dump(homologies, f, sort_keys=True, indent=2)



if __name__ == '__main__':

    import warnings

    # Get all the homolog information from ensemble for genes
    dir_cur = os.path.dirname(os.path.realpath(__file__))
    repo_dir = os.path.join(dir_cur, "..", "repository")
    print(repo_dir)

    retrieve_homologies = False

    species = ["mouse", "rat"]

    # -------------------------------------------
    # Calculate homologies, this takes very long
    # -------------------------------------------
    if retrieve_homologies:

        symbols = []
        with open(os.path.join(repo_dir, "genes", "human_genes.json"), "r") as f_human_genes:
            human_genes = json.load(f_human_genes)
            symbols = [g['ensemble'] for g in human_genes.values() if 'ensemble' in g]

        get_homologies(symbols=symbols, out_dir=os.path.join(repo_dir, "genes"), species=species)

    # -------------------------------------------
    # Perform gene mapping
    # -------------------------------------------
    client = EnsemblRestClient()
    s = "mouse"
    with open(os.path.join(repo_dir, "genes", "human_genes.json"), "r") as f_human_genes:
        human_genes = json.load(f_human_genes)

        entrez_map = {}

        with open(os.path.join(repo_dir, "genes", "{}_homologies.json".format(s)), "r") as f_homo:
            homologies = json.load(f_homo)

            for k, gid in enumerate(human_genes):
                gene = human_genes[gid]
                entrez = gid.split('.')[0]

                ensemble = gene.get('ensemble', None)
                if ensemble:
                    # look for the homology information
                    homologs = homologies.get(ensemble)
                    print(homologs)

                    if len(homologs) == 0:
                        warnings.warn("No homolog for: {}".format(ensemble))
                        continue
                    if len(homologs) > 1:
                        warnings.warn("More than one homolog for: {}".format(ensemble))
                    homolog = homologs[0]

                    target = homolog.get('target')
                    if target:
                        symbol = target['id']
                        g = {
                            'ensemble': symbol,
                            'protein_id': target.get('protein_id'),
                            'homology_type': homolog.get('type'),
                        }
                        xrefs = client.get_xrefs(symbol)
                        #pprint(xrefs)
                        for xref in xrefs:
                            dbname = xref['dbname']
                            dbid = xref['primary_id']
                            if dbname == "EntrezGene":
                                g['entrez'] = dbid
                            elif dbname == "HGNC":
                                g['hgnc'] = dbid
                            elif dbname == "MIM_GENE":
                                g['mim'] == dbid
                            elif dbname == "Uniprot_gn":
                                if "uniprot" in g:
                                    g['uniprot'].append(dbid)
                                else:
                                    g['uniprot'] = [dbid]
                        print('*' * 40)
                        pprint(g)
                        print('*' * 40)

                if k == 4:
                    break



    # some problems with gene transcripts which have to be solved (alternative versions from gene to protein)



    # process the information
    '''
    for s in species:

        # load the homologies
        with open(os.path.join(repo_dir, 'genes', "mouse_homologies.json"), "r") as f:
            mouse_homologies = json.load(f)
        pprint(mouse_homologies)

        # generate genes and gene associations
    '''

