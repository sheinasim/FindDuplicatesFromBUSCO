#!/usr/bin/python

import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("buscoTable", help="The full table output of BUSCO (usually called: full_table.tsv for BUSCO5 or full_table_{outprefix}.tsv for BUSCO3)")
parser.add_argument("fasta", help=".fasta BUSCO was run on in -geno mode")
args = parser.parse_args()

df = pd.read_csv(args.buscoTable, sep='\t', comment='#')
df.columns = ["Busco id", "Status", "Sequence", "Gene Start", "Gene End", "Strand", "Score", "Length", "OrthoDB url", "Description"]
df = df[["Busco id", "Status", "Sequence"]]

df_contig_lengths = pd.DataFrame(columns=("Sequence", "Length"))
for record in SeqIO.parse(args.fasta, "fasta"):
    df_contig_lengths = df_contig_lengths.append({"Sequence" : record.name, "Length" : len(record.seq)}, ignore_index = True)

def findDuplicates(df1, df2):
	df1 = pd.merge(left = df1, right = df2)
	df1["Length"] = df1["Length"].apply(pd.to_numeric)
	complete = df1[df1["Status"] == "Complete"]["Sequence"].unique()
	df2 = df1[~df1["Sequence"].isin(complete)].dropna().sort_values("Busco id")
	numberdupes = df2["Sequence"].value_counts().to_frame("Number BUSCOs")
	numberdupes["Sequence"] = numberdupes.index
	df3 = df2[df2["Status"] == "Duplicated"].drop_duplicates(subset=["Sequence"]).sort_values(["Busco id", "Length"])
	buscos = df3["Busco id"]
	chooseLongest = df3[buscos.isin(buscos[buscos.duplicated()])].sort_values(["Busco id"])
	idx = chooseLongest.groupby(['Busco id'])['Length'].transform(max) == chooseLongest['Length']
	keepthese = chooseLongest[idx]
	df4 = df3[~df3["Sequence"].isin(keepthese["Sequence"])].dropna().drop_duplicates("Sequence")
	df5 = pd.merge(left = df4[["Sequence", "Length"]], right = numberdupes, how = "left", left_on = "Sequence", right_on = "Sequence")
	return df5

dupdf = findDuplicates(df, df_contig_lengths)

outfilename = args.fasta.split(".fasta")[0] + "_duplicates.txt"

dupdf.to_csv(outfilename, index=False, sep='\t')

print("Find Duplicates From BUSCO is finished.")
