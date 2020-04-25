import numpy as np
import pandas as pd
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly
import json
import pickle
import plotly.offline as pyo


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
	with open(f"static/data/TF_index/TF_{chrom}_lst.txt", 'rb') as fp:
		TF_index = pickle.load(fp)
	ind_TF = [i for i, x in enumerate(TF_index) if (x >= start and x <= end)]
	with open(f"static/data/annots_index/annots_chr_{chrom}_lst.txt", 'rb') as fp:
		annots_index = pickle.load(fp)
	ind_annots = [i for i, x in enumerate(annots_index) if (x >= start and x <= end)]

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
	

	# TFBS import
	if len(ind_TF) > 0:
		df_TF = pd.read_csv(f"static/data/TF_pos/TF_{chrom}.csv.gz", compression='gzip', skiprows = range(1, ind_TF[0]+1), nrows = (ind_TF[-1]-ind_TF[0]+1), index_col=0, keep_default_na=False)
	else:
		df_TF = pd.DataFrame()
	
	# Other annotations import
	if len(ind_annots) > 0:
		df_annots = pd.read_csv(f"static/data/annots_pos/annots_chr_{chrom}.csv.gz", compression='gzip', skiprows = range(1, ind_annots[0]+1), nrows = (ind_annots[-1]-ind_annots[0]+1), index_col=0, keep_default_na=False)
	else:
		df_annots = pd.DataFrame()
	
	if df_TF.empty:
		TF_options = []
	else:
		TF_options = sorted(list(df_TF['TF_name'].unique()))
	if df_annots.empty:
		cpg_options = []
		hmm_options = []
		enh_dis = True
	else:
		cpg_options = sorted(list(df_annots['CpG_Annotations'].unique()))
		hmm_options = sorted(list(df_annots['ChromHMM'].unique()))
		if df_annots.Enhancers_Annotations.unique().all() != 'NA':
			enh_dis = False
		else:
			enh_dis = True


	# z-scores import


	return TF_options, cpg_options, hmm_options, enh_dis, bv_means_controls, bv_sample, z_scores, sub_genes, df_TF, df_annots



def create_plot(bv_means_controls, bv_sample, z_scores, sub_genes, df_TF, df_annots, start, end, chrom, TF_drop, cpg_annots, chromhmm, enhancers, x_range, y_range, col_sample):

	start = int(start)
	end = int(end)
	shapes = []

	# If there are no controls, no genes or no samples in the region, the figure will be empty and range from start to end
	if bv_means_controls.empty or sub_genes.empty or len(bv_sample) == 0:
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

		# Plot the TFBS that the user selected
		if TF_drop != None:
			if len(TF_drop) > 0:
				if TF_drop[0] != None:
					#hline TF
					shapes.append(
						dict(
							type="line",
							x0=start,
							y0=-0.2-(max(sub_genes.track)*0.13),
							x1=end,
							y1=-0.2-(max(sub_genes.track)*0.13),
							line=dict(
								color="black",
								dash="dot",
							),
							xref="x2",
							yref="y2"
						)
					)
					for i in range(len(TF_drop)):
						TF_fil = df_TF.loc[df_TF['TF_name'] == TF_drop[i]]
						TF_x = df_TF.MAPINFO.loc[df_TF['TF_name'] == TF_drop[i]].reset_index(drop=True).to_list()
						TF_y = [-0.28-(max(sub_genes.track)*0.13)-0.15*i] * len(TF_x)
						pos_map_TF = TF_fil.MAPINFO.to_list()
						index_TF = TF_fil.index.to_list()
						name_TF = TF_fil.TF_name.to_list()
						my_text = [f'Position: {i}<br>CpG site: {j}<br>Transcription Factor: {k}' for i, j, k in zip(pos_map_TF, index_TF, name_TF)] 
						fig.add_trace(
							go.Scatter(
								x = TF_x,
								y = TF_y,
								mode = 'markers',
								line = dict(
									color = 'rgba(220, 69, 31, 1)',
									width = 2
								),
								legendgroup="TF",
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
								y = [TF_y[0]],
								text = TF_drop[i],
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
				if TF_drop[0] == None:
					TF_lim = 0
				else:
					TF_lim = len(TF_drop)*0.15
			else:
				TF_lim = 0
		else:
			TF_lim = 0


		# Plot the CpG locations (islands, shelves, shores, inter) that the user selected
		if cpg_annots != None and cpg_annots != []:
			if len(cpg_annots) != 0:
				#hline cpg
				shapes.append(
					dict(
						type="line",
						x0=start,
						y0=-0.2-(max(sub_genes.track)*0.13)-TF_lim,
						x1=end,
						y1=-0.2-(max(sub_genes.track)*0.13)-TF_lim,
						line=dict(
							color="black",
							dash="dot",
						),
						xref="x2",
						yref="y2"
					)
				)
				for i in range(len(cpg_annots)):
					cpg_fil = df_annots.loc[df_annots['CpG_Annotations'] == cpg_annots[i]]
					cpg_x = df_annots.MAPINFO.loc[df_annots['CpG_Annotations'] == cpg_annots[i]].reset_index(drop=True).to_list()
					cpg_y = [-0.28-(max(sub_genes.track)*0.13)-TF_lim-(0.15*i)] * len(cpg_x)
					pos_map_cpg = cpg_fil.MAPINFO.to_list()
					index_cpg = cpg_fil.index.to_list()
					name_cpg = cpg_fil.CpG_Annotations.to_list()
					my_text = [f'Position: {i}<br>CpG site: {j}<br>CpG annotation: {k}' for i, j, k in zip(pos_map_cpg, index_cpg, name_cpg)] 
					fig.add_trace(
						go.Scatter(
							x = cpg_x,
							y = cpg_y,
							mode = 'markers',
							line = dict(
								color = 'rgba(108, 60, 51, 1)',
								width = 2
							),
							legendgroup="cpg",
							showlegend=False,
							text = my_text,
							hoverinfo = 'text'
						),
						row=2,
						col=1
					)  
					fig.add_trace(
						go.Scatter(
							x = [range_plot[0]],
							y = [cpg_y[0]],
							text = cpg_annots[i],
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
			cpg_lim = len(cpg_annots)*0.15
		else:
			cpg_lim = 0


		# Plot the chromatin conformation that the user selected
		if chromhmm != None and chromhmm != []:
			if len(chromhmm) != 0:
				#hline hmm
				shapes.append(
					dict(
						type="line",
						x0=start,
						y0=-0.2-(max(sub_genes.track)*0.13)-TF_lim-cpg_lim,
						x1=end,
						y1=-0.2-(max(sub_genes.track)*0.13)-TF_lim-cpg_lim,
						line=dict(
							color="black",
							dash="dot",
						),
						xref="x2",
						yref="y2"
					)
				)
				for i in range(len(chromhmm)):
					hmm_fil = df_annots.loc[df_annots['ChromHMM'] == chromhmm[i]]
					hmm_x = df_annots.MAPINFO.loc[df_annots['ChromHMM'] == chromhmm[i]].reset_index(drop=True).to_list()
					hmm_y = [-0.28-(max(sub_genes.track)*0.13)-TF_lim-cpg_lim-(0.15*i)] * len(hmm_x)
					pos_map_hmm = hmm_fil.MAPINFO.to_list()
					index_hmm = hmm_fil.index.to_list()
					name_hmm = hmm_fil.ChromHMM.to_list()
					my_text = [f'Position: {i}<br>CpG site: {j}<br>Chromatin State: {k}' for i, j, k in zip(pos_map_hmm, index_hmm, name_hmm)] 
					fig.add_trace(
						go.Scatter(
							x = hmm_x,
							y = hmm_y,
							mode = 'markers',
							line = dict(
								color = 'rgba(198, 39, 21, 1)',
								width = 2
							),
							legendgroup="hmm",
							showlegend=False,
							text = my_text,
							hoverinfo = 'text'
						),
						row=2,
						col=1
					)  
					fig.add_trace(
						go.Scatter(
							x = [range_plot[0]],
							y = [hmm_y[0]],
							text = chromhmm[i],
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
			hmm_lim = len(chromhmm)*0.15
		else:
			hmm_lim = 0

		# Plot the enhancers if the user selected the checkbox
		if enhancers != None:
			if enhancers == 'true' and len(df_annots.loc[df_annots['Enhancers_Annotations'] == 'enhancers_fantom']) > 0:
				enh_lim = 0.15
				#hline enhancers
				shapes.append(
					dict(
						type="line",
						x0=start,
						y0=-0.2-(max(sub_genes.track)*0.13)-TF_lim-cpg_lim-hmm_lim,
						x1=end,
						y1=-0.2-(max(sub_genes.track)*0.13)-TF_lim-cpg_lim-hmm_lim,
						line=dict(
							color="black",
							dash="dot",
						),
						xref="x2",
						yref="y2"
					)
				)
				enh_fil = df_annots.loc[df_annots['Enhancers_Annotations'] == 'enhancers_fantom']
				enh_x = df_annots.MAPINFO.loc[df_annots['Enhancers_Annotations'] == 'enhancers_fantom'].reset_index(drop=True).to_list()
				enh_y = [-0.28-(max(sub_genes.track)*0.13)-TF_lim-cpg_lim-hmm_lim] * len(enh_x)
				pos_map_enh = enh_fil.MAPINFO.to_list()
				index_enh = enh_fil.index.to_list()
				my_text = [f'Position: {i}<br>CpG site: {j}<br>Fantom Enhancer' for i, j in zip(pos_map_enh, index_enh)] 
				fig.add_trace(
					go.Scatter(
						x = enh_x,
						y = enh_y,
						mode = 'markers',
						line = dict(
							color = 'rgba(111, 82, 151, 1)',
							width = 2
						),
						name = 'Enhancer',
						legendgroup="enh",
						showlegend=False,
						text = my_text,
						hoverinfo = 'text'
					),
					row=2,
					col=1
				)
				fig.add_trace(
					go.Scatter(
						x = [range_plot[0]],
						y = [enh_y[0]],
						text = 'Enhancers',
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
			else:
				enh_lim = 0
		else:
			enh_lim = 0

		# len_drops will be the same as prop_len but for the y-axis. It will allows to plot a bigger graph if there are a lot of annotations selected.
		# Thanks to that, the beta values graph is not flattened
		len_drops = TF_lim + cpg_lim + hmm_lim + enh_lim




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
		max_y2 = max(6, max(max_tmp_y2))

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
			height=750+(len_drops*500),
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
			"end": end,
			"len_drops": len_drops
		}
	}

	# Create the JSON object with all the information
	graph = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
	graph_comp = json.loads(graph)
	graph_comp.update(config)
	graphJSON = json.dumps(graph_comp)

	return graphJSON