import pandas as pd
import json


# List all chromosomes
chrom_lst = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M"]

# Import transcription factors dataframe
tf = pd.read_csv("epic_TF.csv.gz", index_col=0).drop('nb_TF', axis=1).rename(columns={'TF_name': 'Transcription_Factors'})
tf.to_csv("tf.csv")

# Import manifest dataframe
fields = ["IlmnID", "CHR", "MAPINFO"]
dtype_chr = {'CHR': 'str'}
manifest = pd.read_csv("hg38_manifest.csv", usecols=fields, index_col=0, dtype=dtype_chr)

# Merge the two dataframes
test = pd.merge(manifest, tf, left_index=True, right_index=True)
test.to_csv("test.csv")

# For each chromosome, save a gzipped csv file and the indices of the positions in the chromosome
for i in chrom_lst:
	test_chrom = test.loc[test["CHR"] == i]
	test_chrom.to_csv(f"Transcription_Factors/tf_chr_{i}.csv.gz", compression="gzip")
	index_chrom = test_chrom.MAPINFO.to_list()
	with open(f"Transcription_Factors/tf_chr_{i}.json", "w") as f:
		json.dump(index_chrom, f)
