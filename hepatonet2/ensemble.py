"""
Query homologs from ensmble for given genes.
"""
import os
import json
import requests
import json
import time
from pprint import pprint

# Mapping based on ENSG
# https://rest.ensembl.org/homology/id/ENSG00000106633?content-type=application/json&target_species=mouse
# => searches the homolog protein_ids

# TODO: for all ENSG ids with mouse and rat !

# Than we have ENSG ids
# https://rest.ensembl.org/xrefs/id/ENSG00000106633?content-type=application/json
# https://rest.ensembl.org/xrefs/id/ENSG00000106633?content-type=application/json;external_db=Uniprot_gn

# TODO: perform mapping to uniprot identifiers & than double check the UniProtIDs based on uniprot evidence


class EnsemblRestClient(object):
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

    # Get all the homolog information from ensemble for genes
    dir_cur = os.path.dirname(os.path.realpath(__file__))
    repo_dir = os.path.join(dir_cur, "..", "repository")

    if False:
        # Calculate homologies, this takes very long
        symbols = []
        with open(os.path.join(repo_dir, "genes", "human_genes.json"), "r") as f_human_genes:
            human_genes = json.load(f_human_genes)
            symbols = [g['ensemble'] for g in human_genes.values() if 'ensemble' in g]

        get_homologies(symbols=symbols, out_dir=os.path.join(repo_dir, "genes"), species=["mouse", "rat"])

