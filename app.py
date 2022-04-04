#!/usr/bin/env python3
#
# Takes as input a file of UniProt identifiers and returns protein function predictions from the DeepGo web API
# Example usage: ./app.py <filename>
# Example file data:
#   SeqID Pdomain_length Pdomain_start Pdomain_end UniProtID
#   19   309   326   634  Q9ULL5|PRR12_HUMAN Proline-rich protein 12 OS=Homo sapiens OX=9606
#   73   115   437   551  Q8NDF8|PAPD5_HUMAN Terminal nucleotidyltransferase 4B OS=Homo sapi
#
# Results are output to a CSV file named: <input_filename>_results.csv
#
# Dependencies:
#   Install required Python libraries:
#       pip3 install -U requests progress

import sys
import requests
import json
import pprint
import itertools
from progress.bar import IncrementalBar
import csv

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

#pprint.pprint(fasta_data)

# For each identifier/fasta combo, submit batches to the DeepGo API for prediction results
#for ident, fasta in fasta_data.items():
print("Fetching DeepGO predictions in batches of", deepgo_batch_size)
bar = IncrementalBar('Progress', max = len(fasta_data)/deepgo_batch_size)
for chunk in chunked(fasta_data.items(), deepgo_batch_size):
    # gather batches of fasta data to submit to the API
    for ident, fasta in chunk:
        deepgo_json['data'] += fasta
    #print(deepgo_json)
    # submit data to API
    r = requests.post(deepgo_api_url, data=deepgo_json)
    #print(ident)
    # iterate through predictions for each protein
    try:
        for protein in r.json()['predictions']:
            #pprint.pprint(protein)
            # parse out group of prediction results we're looking for
            for group in protein['functions']:
                if group['name'] == deepgo_group_name:
                    #pprint.pprint(group['functions'])
                    # add predictions to final data structure
                    results[protein['protein_info'].split('|')[1]] = group['functions']
    # Some of the fasta data seems to have invalid characters? Example API response:
    #   { "detail": "JSON parse error - Invalid control character at: line 4 column 106 (char 159)" }
    except KeyError as e:
        #print("Failure gathering predictions:", e)
        pass

    bar.next()

bar.finish()
#pprint.pprint(results)

# Write data to CSV file
# create header row
header = ['ident']
for prediction_group in results.values():
    for prediction in prediction_group:
        if prediction[1] not in header:
            header.append(prediction[1])

# modify data to fit header
data = []
for ident, prediction_group in results.items():
    row = []
    row.append(ident)
    for column in header[1:]:
        found = False
        for prediction in prediction_group:
            if prediction[1] == column:
                row.append(prediction[2])
                found = True
                continue
        if not found:
            row.append(0)
    data.append(row)

with open(sys.argv[1]+'_results.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    writer.writerow(header)

    # write multiple rows
    writer.writerows(data)

print("Finished, results file:", sys.argv[1]+'_results.csv')
