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
df_all = df[["Busco id", "Status", "Sequence"]]
duplicated = df[df["Status"] == "Duplicated"]['Sequence'].unique()

df_contig_lengths = pd.DataFrame(columns=('Sequence', 'Length'))
for record in SeqIO.parse(args.fasta, "fasta"):
    df_contig_lengths = df_contig_lengths.append({'Sequence' : record.name, 'Length' : len(record.seq)}, ignore_index = True)

def findDuplicates(df1, df2, duptigs):
	df1 = pd.merge(left = df1, right = df2)
	df1['Length'] = df1['Length'].apply(pd.to_numeric)
	df3 = df1.sort_values('Length', ascending=False).drop_duplicates('Busco id').sort_index()
	keepers = df3['Sequence'].unique()
	duplicate_contigs = np.setdiff1d(duptigs, keepers)
	return duplicate_contigs

dupArray = findDuplicates(df_all, df_contig_lengths, duplicated)

np.savetxt("duplicates.txt", dupArray, delimiter='\n', fmt='%s')
