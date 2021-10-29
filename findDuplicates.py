#!/usr/bin/python

import json
import pandas as pd
import argparse
import numpy as np, scipy.stats as st
import os
import re
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("blobdir", help="The blob directory containing read coverage in .json format")
parser.add_argument("buscoTable", help="The full table output of BUSCO (usually called: full_table.tsv for BUSCO5 or full_table_{outprefix}.tsv for BUSCO3)")
parser.add_argument("fasta", help=".fasta BUSCO was run on in -geno mode")
args = parser.parse_args()

path = args.blobdir

for filename in os.listdir(path):
    if re.match(r".*reads_cov\.json", filename):
         with open(os.path.join(path, filename), 'r') as f:
                readCov = json.loads(f.read())
                
coverage = pd.json_normalize(readCov, record_path=['values'])

df = pd.read_csv(args.buscoTable, sep='\t', comment='#', usecols=[0, 1, 2], header=None)
df.columns = ["busco id", "Status", "record"]
df["record"] = df["record"].apply(str)

df_contig_lengths = pd.DataFrame(columns=("record", "length"))
for record in SeqIO.parse(args.fasta, "fasta"):
    df_contig_lengths = df_contig_lengths.append({"record" : record.name, "length" : len(record.seq)}, ignore_index = True)

df_contig_lengths["record"] = df_contig_lengths["record"].apply(str)
    
def findDuplicates(df1, df2, cov):
	df1 = pd.merge(left = df1, right = df2)
	df1 = pd.merge(left = df1, right = cov, left_index=True, right_index=True)
	df1 = df1.rename(columns={0: "coverage"})    
	df1["length"] = df1["length"].apply(pd.to_numeric)
	complete = df1[df1["Status"] == "Complete"]["record"].unique()
	df2 = df1[~df1["record"].isin(complete)].dropna().sort_values("busco id")
	numberdupes = df2["record"].value_counts().to_frame("Number BUSCOs")
	numberdupes["record"] = numberdupes.index
	df3 = df2[df2["Status"] == "Duplicated"].drop_duplicates(subset=["record"]).sort_values(["busco id", "length"])
	buscos = df3["busco id"]
	chooseLongest = df3[buscos.isin(buscos[buscos.duplicated()])].sort_values(["busco id"])
	idx = chooseLongest.groupby(['busco id'])['length'].transform(max) == chooseLongest['length']
	keepthese = chooseLongest[idx]
	df4 = df3[~df3["record"].isin(keepthese["record"])].dropna().drop_duplicates("record")
	df5 = pd.merge(left = df4[["record", "length", "coverage"]], right = numberdupes, how = "left", left_on = "record", right_on = "record")
	low,high = st.t.interval(0.95, len(df5["coverage"])-1, loc=np.mean(df5["coverage"]), scale=st.sem(df5["coverage"]))
	lessthanhalflow = low/2
	morethantwicehigh = high*2
	conditions = [
        (df5['coverage'] <= lessthanhalflow),
        (df5['coverage'] > lessthanhalflow) & (df5['coverage'] <= morethantwicehigh),
        (df5['coverage'] > morethantwicehigh)
        ]
	values = ['low', 'medium', 'high']
	df5['coverage class'] = np.select(conditions, values)
	return df5

dupdf = findDuplicates(df, df_contig_lengths, coverage)

prefix = args.fasta.split('/')[-1].split('.fasta')[0]

outfilename = prefix + "_duplicates.txt"

dupdf.to_csv(outfilename, index=False, sep='\t', header=True)

print("Find Duplicates From BUSCO is finished.")
