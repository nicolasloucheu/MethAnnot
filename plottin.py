import numpy as np
import pandas as pd
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly
import json
import pickle
import plotly.offline as pyo
import os
from functools import reduce

pd.set_option('display.max_rows', None)

# Create the dataframes from region information and sample paths
def create_dfs(chrom, start, end, sample_name):
	start = int(start)
	end = int(end)

	# Getting indices to import only the data in the wanted region
	with open(f"static/data/chrom_lists/chrom_{chrom}_lst.txt", 'rb') as fp:
		chrom_list_controls = pickle.load(fp)
	ind_controls = [i for i, x in enumerate(chrom_list_controls) if (x >= start and x <= end)]
	with open(f"static/data/genes_index/genes_{chrom}_pos.txt", 'rb') as fp:
		genes_list = pickle.load(fp)
	ind_genes = [i for i, x in enumerate(genes_list) if (x >= start and x <= end)]

	# Prop_len will be a dynamic variable defining the length of the current x-axis so that object sizes will vary depending on that length. 
	# They will not be too big nor too small (example: triangles showing the direction of genes have to vary depending on the length of the x-axis)
	prop_len = (end-start)/60

	# If there are beta-values of controls, import those values, else import an empty dataframe
	if len(ind_controls) > 0:
		bv_means_controls = pd.read_csv(f"static/data/chrom_means/chrom_{chrom}_means.csv.gz", 
										compression="gzip", index_col=0, skiprows = range(1, ind_controls[0]+1), nrows = (ind_controls[-1]-ind_controls[0]+1))
	else:
		bv_means_controls = pd.DataFrame()


	# For each sample, if there are beta-values, import them, else import an empty dataframe
	# Same for z-scores of samples
	bv_sample = []
	z_scores = []
	for sample_value in sample_name:
		with open(f"instance/uploads/{sample_value}/{sample_value}_chrom_{chrom}_lst.txt", 'rb') as fp:
			chrom_list_sample = pickle.load(fp)
		ind_sample = [i for i, x in enumerate(chrom_list_sample) if (x >= start and x <= end)]

		with open(f"instance/uploads/{sample_value}/{sample_value}_index_{chrom}.txt", 'rb') as fp:
			z_list_sample = pickle.load(fp)
		ind_z = [i for i, x in enumerate(z_list_sample) if (x >= start and x <= end)]

		if len(ind_sample) > 0:
			bv_sample.append(pd.read_csv(f"instance/uploads/{sample_value}/{sample_value}_chrom_{chrom}.csv.gz", compression='gzip', index_col=0, skiprows = range(1, ind_sample[0]+1), nrows = (ind_sample[-1]-ind_sample[0]+1)))
		else:
			bv_sample.append(pd.DataFrame())
		if len(ind_z) > 0:
			z_scores.append(pd.read_csv(f"instance/uploads/{sample_value}/{sample_value}_z_{chrom}.csv.gz", compression='gzip', index_col=0, skiprows = range(1, ind_z[0]+1), nrows = (ind_z[-1]-ind_z[0]+1)))
		else:
			z_scores.append(pd.DataFrame())

	# Genes import and parsing so the genes don't overlap (depending on prop_len)
	if len(ind_genes) > 0:
		df_genes = pd.read_csv(f"static/data/genes/genes_{chrom}.csv.gz", compression='gzip', skiprows = range(1, ind_genes[0]+1), nrows = (ind_genes[-1]-ind_genes[0]+1), usecols=['external_gene_name', 'START', 'END', 'STRAND'])
	
		genes_df_lst = []
		gene_start = []
		gene_end = []
		gene_strand = []
		gene_name = []
		for i in df_genes.external_gene_name.unique():
			sub_df = df_genes.loc[df_genes.external_gene_name == i].reset_index()
			gene_start.append(min(sub_df.START))
			gene_end.append(max(sub_df.END))
			gene_strand.append(sub_df.STRAND.loc[0])
			gene_name.append(i)
		sub_genes = pd.DataFrame({'Gene_name': gene_name, 'Gene_start':gene_start, 'Gene_end': gene_end, 'Gene_strand': gene_strand})
		sub_genes.loc[0, 'track'] = 0
		for i in range(1, len(sub_genes)):
			b1 = sub_genes.loc[i, 'Gene_start']
			track = 0
			e0_index = max(sub_genes.loc[sub_genes['track'] == track].index.to_list())
			e0 = sub_genes.loc[e0_index, 'Gene_end']
			while b1 <= e0+(prop_len*5):
				track += 1
				try:
					e0_index = max(sub_genes.loc[sub_genes['track'] == track].index.to_list())
					e0 = sub_genes.loc[e0_index, 'Gene_end']
				except:
					break
			sub_genes.loc[i, 'track'] = track

	else:
		df_genes = pd.DataFrame()
		sub_genes = pd.DataFrame()

	directory = "static/data/Annotations"
	annots_names = next(os.walk(directory))[1]
	indices_annots = []
	annotations = []
	total_options = {}


	for i in range(len(annots_names)):
		for file in os.listdir(f"static/data/Annotations/{annots_names[i]}"):
			if file.endswith(".json") and f"chr_{chrom}." in file:
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

	print(annots)

	for x in annots_names:
		annots_options = []
		list_annots = annots.loc[:,x].dropna().to_list()
		for i in list_annots:
			each_cpg_list = i.split(',')
			for j in each_cpg_list:
				annots_options.append(j)
		total_options[x] = pd.Series(annots_options).unique().tolist()

	return bv_means_controls, bv_sample, z_scores, sub_genes, annots, annots_names, total_options



def create_plot(bv_means_controls, bv_sample, z_scores, sub_genes, start, end, chrom, x_range, y_range, col_sample, annots, annots_values):

	start = int(start)
	end = int(end)
	shapes = []

	# If there are no controls, no genes or no samples in the region, the figure will be empty and range from start to end
	if bv_means_controls.empty and sub_genes.empty and len(bv_sample) == 0:
		fig = go.Figure(
			layout = dict(
				xaxis = dict(
					range = [start, end]
				)
			)
		)
		len_drops = 0

	# If at least one of them is not empty, it will generate a figure
	else:
		# Defining the range of the plot if the user zoomd/moved or not
		range_plot = [0, 0]
		if x_range != None and y_range != None:
			range_plot = x_range
		else:
			range_plot[0] = start
			range_plot[1] = end

		# Redefining prop_len depending on the zoom
		prop_len = (range_plot[1] - range_plot[0])/60

		# Create the figure
		fig = make_subplots(
			rows=2,
			cols=1,
			shared_xaxes=True,
			vertical_spacing=0.02,
			row_heights=[0.3, 0.7]
		)

		# Adding trace for percentile 1% of controls
		fig.add_trace(
			go.Scatter(
				x = bv_means_controls['MAPINFO'],
				y = bv_means_controls['PER01'],
				mode='none',
				line = dict(width = 0),
				legendgroup="controls",
				showlegend=False,
				hoverinfo='none'
			),
			row=2,
			col=1
		)

		#Adding trace for percentile 5% of controls
		fig.add_trace(
			go.Scatter(
				x = bv_means_controls['MAPINFO'],
				y = bv_means_controls['PER05'],
				fill = 'tonexty',
				mode = 'none',
				fillcolor = 'rgba(255, 0, 0, .05)',
				legendgroup="controls",
				showlegend=False,
				hoverinfo='none'
			),
			row=2,
			col=1
		)

		#Adding trace for percentile 10% of controls
		fig.add_trace(
			go.Scatter(
				x = bv_means_controls['MAPINFO'],
				y = bv_means_controls['PER10'],
				fill = 'tonexty',
				mode = 'none',
				fillcolor = 'rgba(255, 0, 0, .1)',
				legendgroup="controls",
				showlegend=False,
				hoverinfo='none'
			),
			row=2,
			col=1
		)

		#Adding trace for percentile 25% of controls
		fig.add_trace(
			go.Scatter(
				x = bv_means_controls['MAPINFO'],
				y = bv_means_controls['PER25'],
				fill = 'tonexty',
				mode = 'none',
				fillcolor = 'rgba(255, 0, 0, .2)',
				legendgroup="controls",
				showlegend=False,
				hoverinfo='none'
			),
			row=2,
			col=1
		)

		#Adding trace for mean of controls
		fig.add_trace(
			go.Scatter(
				x = bv_means_controls['MAPINFO'],
				y = bv_means_controls['MEAN'],
				fill = 'tonexty',
				mode = 'markers+lines',
				fillcolor = 'rgba(255, 0, 0, .4)',
				line = dict(color = 'rgba(255, 0, 0, 1)'),
				legendgroup="controls",
				name = "Controls",
				customdata = bv_means_controls.index,
				hovertemplate = 'Beta-value: %{y:.2f}<br>Position: %{x}<br>CpG name: %{customdata}'
			),
			row=2,
			col=1
		)

		#Adding trace for percentile 75% of controls
		fig.add_trace(
			go.Scatter(
				x = bv_means_controls['MAPINFO'],
				y = bv_means_controls['PER75'],
				fill = 'tonexty',
				mode = 'none',
				fillcolor = 'rgba(255, 0, 0, .4)',
				legendgroup="controls",
				showlegend=False,
				hoverinfo='none'
			),
			row=2,
			col=1
		)

		#Adding trace for percentile 90% of controls
		fig.add_trace(
			go.Scatter(
				x = bv_means_controls['MAPINFO'],
				y = bv_means_controls['PER90'],
				fill = 'tonexty',
				mode = 'none',
				fillcolor = 'rgba(255, 0, 0, .2)',
				legendgroup="controls",
				showlegend=False,
				hoverinfo='none'
			),
			row=2,
			col=1
		)

		#Adding trace for percentile 95% of controls
		fig.add_trace(
			go.Scatter(
				x = bv_means_controls['MAPINFO'],
				y = bv_means_controls['PER95'],
				fill = 'tonexty',
				mode = 'none',
				fillcolor = 'rgba(255, 0, 0, .1)',
				legendgroup="controls",
				showlegend=False,
				hoverinfo='none'
			),
			row=2,
			col=1
		)

		#Adding trace for percentile 99% of controls
		fig.add_trace(
			go.Scatter(
				x = bv_means_controls['MAPINFO'],
				y = bv_means_controls['PER99'],
				fill = 'tonexty',
				mode = 'none',
				fillcolor = 'rgba(255, 0, 0, .05)',
				legendgroup="controls",
				showlegend=False,
				hoverinfo='none'
			),
			row=2,
			col=1
		)

		# For each sample plot the beta-values
		for i in range(len(bv_sample)):
			try:
				fig.add_trace(
					go.Scatter(
						x = bv_sample[i]['MAPINFO'], 
						y = bv_sample[i].iloc[:,0],
						line = dict(
							color = col_sample[bv_sample[i].columns[0]],
							width = 2
						),
						name = f"{bv_sample[i].columns[0]}",
						legendgroup = f"group_{bv_sample[i].columns[0]}",
						mode = 'markers+lines',
						customdata = bv_sample[i].index,
						hovertemplate = 'Beta-value: %{y:.2f}<br>Position: %{x}<br>CpG name: %{customdata}'
					),
				row=2,
				col=1
				)
			except:
				pass

		shapes.append(
			dict(
			# Horizontal Line
				type="line",
				x0=start,
				y0=0,
				x1=end,
				y1=0,
				line=dict(
					color="black",
					dash="dot",
				),
				xref="x2",
				yref="y2"
			)
		)

		# Show genes as rectangle
		for i in range(len(sub_genes)):
			shapes.append(
				dict(
					type='rect',
					x0 = sub_genes.loc[i]['Gene_start'],
					y0 = -0.15-(0.13*sub_genes.loc[i]['track']),
					x1 = sub_genes.loc[i]['Gene_end'],
					y1 = -0.1-(0.13*sub_genes.loc[i]['track']),
					line = dict(
						width = 1
					),
					fillcolor = 'rgba(230, 154, 89, 0.7)',
					xref="x2",
					yref="y2"
				)
			)

			# Compute where to plot the Gene name
			a = 0
			if (sub_genes.loc[i]['Gene_start'] < range_plot[0] and sub_genes.loc[i]['Gene_end'] < range_plot[1] and \
				sub_genes.loc[i]['Gene_end'] > range_plot[0]):
				a = (range_plot[0] + sub_genes.loc[i]['Gene_end'])/2
			elif (sub_genes.loc[i]['Gene_end'] > range_plot[1] and sub_genes.loc[i]['Gene_start'] > range_plot[0] and \
				  sub_genes.loc[i]['Gene_start'] < range_plot[1]):
				a = (range_plot[1] + sub_genes.loc[i]['Gene_start'])/2
			elif (sub_genes.loc[i]['Gene_start'] > range_plot[0] and sub_genes.loc[i]['Gene_start'] < range_plot[1] and \
				  sub_genes.loc[i]['Gene_end'] < range_plot[1] and\
				  sub_genes.loc[i]['Gene_end'] > range_plot[0]):
				a = (sub_genes.loc[i]['Gene_start'] + sub_genes.loc[i]['Gene_end'])/2
			elif (sub_genes.loc[i]['Gene_start'] < range_plot[0] and sub_genes.loc[i]['Gene_end'] > range_plot[1]):
				a = (range_plot[1] + range_plot[0])/2
			fig.add_trace(
				go.Scatter(
					x = [a],
					y = [-0.09-(0.13*sub_genes.loc[i]['track'])],
					text = [sub_genes.loc[i]['Gene_name']],
					mode = 'text',
					showlegend=False,
					hoverinfo='none',
					textfont = dict(
						size = 13
					),
					textposition="top center"
				),
				row=2,
				col=1
			)

			# Compute where to plot the triangle showing the gene strand
			xplace = sub_genes.loc[i]['Gene_end'] if sub_genes.loc[i]['Gene_strand'] == 1 else sub_genes.loc[i]['Gene_start']
			genes_dir = xplace+(0.5*prop_len) if sub_genes.loc[i]['Gene_strand'] == 1 else xplace-(0.5*prop_len)
			yhaut = -0.095-(0.13*sub_genes.loc[i]['track'])
			ybas = -0.155-(0.13*sub_genes.loc[i]['track'])
			ymid = -0.125-(0.13*sub_genes.loc[i]['track'])
			fig.add_trace(
				go.Scatter(
					x = [xplace, xplace, genes_dir, xplace],
					y = [yhaut, ybas, ymid, yhaut],
					fill = 'toself',
					mode = 'none',
					fillcolor = 'black',
					showlegend=False,
					hoverinfo='none'
				),
				row=2,
				col=1
			)

		###Plotting z-scores
		max_tmp_y2 = []
		for i in range(len(z_scores)):
			try:
				fig.add_trace(
					go.Scatter(
						x = z_scores[i]['MAPINFO'], 
						y = z_scores[i].iloc[:,0],
						line = dict(
							color = col_sample[z_scores[i].columns[0]]
						),
						name = f"{z_scores[i].columns[0]}",
						legendgroup = f"group_{z_scores[i].columns[0]}",
						showlegend = False,
						mode = 'markers',
						customdata = z_scores[i].index,
						hovertemplate = 'Z-score: %{y:.2f}<br>Position: %{x}<br>CpG name: %{customdata}'
					),
				row=1,
				col=1
				)
				max_tmp_y2.append(max(z_scores[i].iloc[:,0]))
			except:
				pass
		try:
			max_y2 = max(6, max(max_tmp_y2))
		except:
			max_y2 = 6

		#Threshold
		shapes.append(
			dict(
			# Horizontal Line
				type="line",
				x0=start,
				y0=5,
				x1=end,
				y1=5,
				line=dict(
					color="red"
				),
				xref="x1",
				yref="y1"
			)
		)

		#Plot the annotations the user selected

		color_palette = ['rgba(231,148,26,1)', 'rgba(194,16,12,1)', 'rgba(55,186,77)', 'rgba(34,86,177,1)', 'rgba(239,138,76,1)', 'rgba(37,124,112,1)', 'rgba(34,3,44,1)']

		to_del = []
		annots_to_plot = {}
		annots_merge = pd.DataFrame(columns=['CHR', 'MAPINFO'])

		for annot_name in annots_values:
			if annots_values[annot_name] == None:
				to_del.append(annot_name)

		for i in to_del:
			del annots_values[i]

		for i in annots_values:
			annots_to_plot.setdefault(i, {})
			for j in annots_values[i]:
				tmp = annots.loc[annots[i].str.contains(j, na=False), ["CHR", "MAPINFO", i]].rename(columns={i: j})
				tmp[j] = j
				annots_merge = pd.merge(tmp, annots_merge, left_index=True, right_index=True, on=['CHR', 'MAPINFO'], how='outer')
				annots_to_plot[i][j] = tmp.MAPINFO.to_list()

		list_len = [0]
		for i in annots_to_plot:
			list_len.append(len(annots_to_plot[i]))
		list_len = np.cumsum(list_len).tolist()
		len_drops = list_len[-1]

		for i, j in enumerate(annots_to_plot):
			shapes.append(
				dict(
					type="line",
					x0=start,
					y0=-0.2-(max(sub_genes.track)*0.13)-0.1*i-list_len[i]*0.1,
					x1=end,
					y1=-0.2-(max(sub_genes.track)*0.13)-0.1*i-list_len[i]*0.1,
					line=dict(
						color="black",
						dash="dot",
					),
					xref="x2",
					yref="y2"
				)
			)
			for k, l in enumerate(annots_to_plot[j]):
				annots_x = annots_to_plot[j][l]
				annots_y = [-0.3-(max(sub_genes.track)*0.13)-0.1*i-list_len[i]*0.1-0.1*k] * len(annots_x)
				sub_df = annots_merge.loc[annots_merge['MAPINFO'].isin(annots_to_plot[j][l]), ['MAPINFO', l]]
				group_annots = [j for x in range(len(annots_x))]
				pos_map = sub_df.MAPINFO.to_list()
				index_annots = sub_df.index.to_list()
				name_annots = sub_df[l].to_list()
				my_text = [f'Position: {i}<br>CpG site: {j}<br>Group of annotation: {k}<br>Annotation: {l}' for i, j, k, l in zip(pos_map, index_annots, group_annots, name_annots)] 
				fig.add_trace(
					go.Scatter(
						x = annots_x,
						y = annots_y,
						mode = 'markers',
						line = dict(
							color = color_palette[i],
							width = 2
						),
						showlegend=False,
						name = '',
						text = my_text,
						hoverinfo = "text"
					),
					row=2,
					col=1
				)
				fig.add_trace(
					go.Scatter(
						x = [range_plot[0]],
						y = [annots_y[0]],
						text = l,
						mode = 'text',
						showlegend=False,
						hoverinfo='none',
						textfont = dict(
							size = 13
						),
						textposition = 'middle right'
					),
					row=2,
					col=1
				)


		#--------------------------------------------

		# Define the y ticks to show
		# fig.update_yaxes(
		# 	tickvals=[0, 0.2, 0.4, 0.6, 0.8, 1]
		# )

		# Defines the layout of the plot
		fig.update_layout(
			shapes = shapes,
			yaxis1 =  dict(
				title_text = "Z-score",
				title_font = dict(
					color = "white"
				),
				tickfont = dict(
					color = "white"
				),
				range = [0, (max_y2 + (max_y2/10))]
			),
			xaxis2 = dict(
				range = range_plot,
				title_text = 'Position',
				title_font = dict(
					color = 'white'
				),
				tickfont = dict(
					color = 'white'
				)
			),
			yaxis2 = dict(
				title_text = "Beta-value",
				title_font = dict(
					color = 'white'
				),
				tickfont = dict(
					color = 'white'
				),
				tickvals = [0, 0.2, 0.4, 0.6, 0.8, 1]
			),
			height=750+(len_drops*40),
			paper_bgcolor = '#393833',
			legend = dict(
				font = dict(
					color = 'white'
				)
			)
		)

	# Output some important configuration variables
	config = {
		'config': {
			"start": start,
			"end": end
		}
	}

	# Create the JSON object with all the information
	graph = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
	graph_comp = json.loads(graph)
	graph_comp.update(config)
	graphJSON = json.dumps(graph_comp)

	return graphJSON