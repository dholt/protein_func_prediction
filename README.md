# Protein function predictor using DeepGOWeb

This is a helper application for converting a list of UniProt identifiers to a CSV file
of protein function predictions using the (DeepGOWeb)[https://deepgo.cbrc.kaust.edu.sa/deepgo/] API

Takes as input a file of UniProt identifiers and returns protein function predictions from the DeepGo web API

Example usage: `./app.py <filename>`

```
Example file data:
  SeqID Pdomain_length Pdomain_start Pdomain_end UniProtID
  19   309   326   634  Q9ULL5|PRR12_HUMAN Proline-rich protein 12 OS=Homo sapiens OX=9606
  73   115   437   551  Q8NDF8|PAPD5_HUMAN Terminal nucleotidyltransferase 4B OS=Homo sapi
```

Results are output to a CSV file named: `<input_filename>_results.csv`

Dependencies:
  Install required Python libraries:
      `pip3 install -U requests progress`
