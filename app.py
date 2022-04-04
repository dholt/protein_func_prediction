#!/usr/bin/env python3

# Install required Python libraries:
# pip3 install -U requests progress

import sys
import requests
import json
import pprint
import itertools
from progress.bar import IncrementalBar

uniprot_base_url = "https://www.uniprot.org/uniprot"   # /uniprot_id.fasta
deepgo_api_url = "https://deepgo.cbrc.kaust.edu.sa/deepgo/api/create"
deepgo_group_name = "Molecular Function"
deepgo_batch_size = 10
uniprot_ident = []
fasta_data = {}
results = {}

# Data structure sent to DeepGo API
deepgo_json = {
    "version": "1.0.3",
    "data_format": "fasta",
    "data": "",
    "threshold": "0.3"
}

# Function for iterating over dictionary in batches
# https://stackoverflow.com/questions/28022223/how-to-iterate-over-a-dictionary-n-key-value-pairs-at-a-time
def chunked(it, size):
    it = iter(it)
    while True:
        p = tuple(itertools.islice(it, size))
        if not p:
            break
        yield p

# Read input file and parse Uniprot identifiers
with open(sys.argv[1], 'r') as data_file:
    # skip the first line in the file containing headers
    lines = data_file.readlines()[1:]
    # read each line in input file and parse out identifier
    for line in lines:
        uniprot_ident.append((line.split('|')[0].split(' ')[-1]))

# For each identifier, get Fasta file
print("Fetching Fasta data for", len(uniprot_ident), "identifiers")
bar = IncrementalBar('Progress', max = len(uniprot_ident))
for ident in uniprot_ident:
    f = requests.get(uniprot_base_url + '/' + ident + '.fasta')
    fasta_data[ident] = f.text
    bar.next()
bar.finish()

# For each identifier/fasta combo, submit batches to the DeepGo API for prediction results
#for ident, fasta in fasta_data.items():
for chunk in chunked(fasta_data.items(), deepgo_batch_size):
    # gather batches of fasta data to submit to the API
    for ident, fasta in chunk:
        deepgo_json['data'] += fasta
    #print(deepgo_json)
    # submit data to API
    r = requests.post(deepgo_api_url, data=deepgo_json)
    #print(ident)
    # iterate through predictions for each protein
    for protein in r.json()['predictions']:
        #pprint.pprint(protein)
        # parse out group of prediction results we're looking for
        for group in protein['functions']:
            if group['name'] == deepgo_group_name:
                #pprint.pprint(group['functions'])
                # add predictions to final data structure
                results[protein['protein_info'].split('|')[1]] = group['functions']

pprint.pprint(results)
