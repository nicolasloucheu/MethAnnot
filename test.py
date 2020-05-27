import os
import json
import numpy as np
import pandas as pd
from functools import reduce


chrom = 2
start = 102200815
end = 102300815

directory = "static/data/Annotations"
annots_names = next(os.walk(directory))[1]
indices_annots = []
annotations = []
total_options = {}


for i in range(len(annots_names)):
	for file in os.listdir(f"static/data/Annotations/{annots_names[i]}"):
		if file.endswith(".json") and f"chr_{chrom}." in file:
			print(file)
			with open(f"static/data/Annotations/{annots_names[i]}/{file}", 'rb') as json_file:
				index_annot = json.load(json_file)
			current_ind = [j for j, x in enumerate(index_annot) if (x >= start and x <= end)]
			indices_annots.append(current_ind)
	for file in os.listdir(f"static/data/Annotations/{annots_names[i]}"):
		if file.endswith(".csv.gz") and f"chr_{chrom}." in file:
			try:
				annot_file = pd.read_csv(f"static/data/Annotations/{annots_names[i]}/{file}", compression='gzip', index_col=0, skiprows = range(1, indices_annots[i][0]+1), nrows = (indices_annots[i][-1]-indices_annots[i][0]+1))
			except:
				annot_file = pd.DataFrame(columns=["CHR", "MAPINFO"])
				annot_file[f"{annots_names[i]}"] = np.nan
			annotations.append(annot_file)


annots = reduce(lambda x, y: pd.merge(x, y, left_index=True, right_index=True, on = ['CHR', 'MAPINFO'], how='outer'), annotations)

for x in annots_names:
	annots_options = []
	list_annots = annots.loc[:,x].dropna().to_list()
	for i in list_annots:
		each_cpg_list = i.split(',')
		for j in each_cpg_list:
			annots_options.append(j)
	total_options[x] = pd.Series(annots_options).unique().tolist()