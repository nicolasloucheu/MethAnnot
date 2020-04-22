import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly
import json
import pickle
import plotly.offline as pyo



def create_dfs(chrom, start, end, sample_name):
	start = int(start)
	end = int(end)
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

	prop_len = (end-start)/60
	if len(ind_controls) > 0:
		bv_means_controls = pd.read_csv(f"static/data/chrom_means/chrom_{chrom}_means.csv.gz", 
										compression="gzip", index_col=0, skiprows = range(1, ind_controls[0]+1), nrows = (ind_controls[-1]-ind_controls[0]+1))
	else:
		bv_means_controls = pd.DataFrame()


	bv_sample = []
	for sample_value in sample_name:
		with open(f"instance/uploads/{sample_value}/{sample_value}_chrom_{chrom}_lst.txt", 'rb') as fp:
			chrom_list_sample = pickle.load(fp)
		ind_sample = [i for i, x in enumerate(chrom_list_sample) if (x >= start and x <= end)]
		if len(ind_sample) > 0:
			bv_sample.append(pd.read_csv(f"instance/uploads/{sample_value}/{sample_value}_chrom_{chrom}.csv.gz", compression='gzip', index_col=0, skiprows = range(1, ind_sample[0]+1), nrows = (ind_sample[-1]-ind_sample[0]+1)))

	#genes import
	
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
	
	#TF import
	if len(ind_TF) > 0:
		df_TF = pd.read_csv(f"static/data/TF_pos/TF_{chrom}.csv.gz", compression='gzip', skiprows = range(1, ind_TF[0]+1), nrows = (ind_TF[-1]-ind_TF[0]+1), index_col=0, keep_default_na=False)
	else:
		df_TF = pd.DataFrame()
	
	#Annots import
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

	return TF_options, cpg_options, hmm_options, enh_dis, bv_means_controls, bv_sample, sub_genes, df_TF, df_annots



def create_plot(bv_means_controls, bv_sample, sub_genes, df_TF, df_annots, start, end, chrom, TF_drop, cpg_annots, chromhmm, enhancers, x_range, y_range, col_sample):

	start = int(start)
	end = int(end)
	if bv_means_controls.empty or sub_genes.empty or len(bv_sample) == 0:
		fig = go.Figure(
			layout = dict(
				xaxis = dict(
					range = [start, end]
				)
			)
		)
	else:
		range_plot = [0, 0]
		if x_range != None and y_range != None:
			range_plot = x_range
		else:
			range_plot[0] = start
			range_plot[1] = end
		prop_len = (range_plot[1] - range_plot[0])/60
		fig = go.Figure()
		fig.add_trace(
			go.Scatter(
				x = bv_means_controls['MAPINFO'],
				y = bv_means_controls['PER01'],
				mode='none',
				line = dict(width = 0),
				legendgroup="controls",
				showlegend=False,
				hoverinfo='none'
			)
		)

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
			)
		)

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
			)
		)

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
			)
		)

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
			)
		)
		print(bv_means_controls)

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
			)
		)

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
			)
		)

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
			)
		)

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
			)
		)

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
						mode = 'markers+lines',
						customdata = bv_sample[i].index,
						hovertemplate = 'Beta-value: %{y:.2f}<br>Position: %{x}<br>CpG name: %{customdata}'
					)
				)
			except:
				pass

		fig.add_shape(
				# Line Horizontal
					type="line",
					x0=start,
					y0=0,
					x1=end,
					y1=0,
					line=dict(
						color="black",
						dash="dot",
					)
		)

		for i in range(len(sub_genes)):
			fig.add_shape(
				type='rect',
				x0 = sub_genes.loc[i]['Gene_start'],
				y0 = -0.15-(0.13*sub_genes.loc[i]['track']),
				x1 = sub_genes.loc[i]['Gene_end'],
				y1 = -0.1-(0.13*sub_genes.loc[i]['track']),
				line = dict(
					width = 1
				),
				fillcolor = 'rgba(230, 154, 89, 0.7)'
			)

			sub_genes.loc[(sub_genes['Gene_end'] > 60500000) & (sub_genes['Gene_start'] < 60700000)]
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
				)
			)
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
				)
			)

		if TF_drop != None:
			if len(TF_drop) > 0:
				if TF_drop[0] != None:
					#hline TF
					fig.add_shape(
						type="line",
						x0=start,
						y0=-0.2-(max(sub_genes.track)*0.13),
						x1=end,
						y1=-0.2-(max(sub_genes.track)*0.13),
						line=dict(
							color="black",
							dash="dot",
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
							)
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
							)
						)
				if TF_drop[0] == None:
					TF_lim = 0
				else:
					TF_lim = len(TF_drop)*0.15
			else:
				TF_lim = 0
		else:
			TF_lim = 0

		if cpg_annots != None and cpg_annots != []:
			if len(cpg_annots) != 0:
				#hline cpg
				fig.add_shape(
					type="line",
					x0=start,
					y0=-0.2-(max(sub_genes.track)*0.13)-TF_lim,
					x1=end,
					y1=-0.2-(max(sub_genes.track)*0.13)-TF_lim,
					line=dict(
						color="black",
						dash="dot",
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
						)
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
						)
					)
			cpg_lim = len(cpg_annots)*0.15
		else:
			cpg_lim = 0

		if chromhmm != None and chromhmm != []:
			if len(chromhmm) != 0:
				#hline hmm
				fig.add_shape(
					type="line",
					x0=start,
					y0=-0.2-(max(sub_genes.track)*0.13)-TF_lim-cpg_lim,
					x1=end,
					y1=-0.2-(max(sub_genes.track)*0.13)-TF_lim-cpg_lim,
					line=dict(
						color="black",
						dash="dot",
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
						)
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
						)
					)
			hmm_lim = len(chromhmm)*0.15
		else:
			hmm_lim = 0

		if enhancers != None:
			if enhancers == 'true' and len(df_annots.loc[df_annots['Enhancers_Annotations'] == 'enhancers_fantom']) > 0:
				enh_lim = 0.15
				#hline enhancers
				fig.add_shape(
					type="line",
					x0=start,
					y0=-0.2-(max(sub_genes.track)*0.13)-TF_lim-cpg_lim-hmm_lim,
					x1=end,
					y1=-0.2-(max(sub_genes.track)*0.13)-TF_lim-cpg_lim-hmm_lim,
					line=dict(
						color="black",
						dash="dot",
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
					)
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
					)
				)
			else:
				enh_lim = 0
		else:
			enh_lim = 0

		len_drops = TF_lim + cpg_lim + hmm_lim + enh_lim
		#--------------------------------------------


		fig.update_yaxes(
			tickvals=[0, 0.2, 0.4, 0.6, 0.8, 1]
		)

		fig.update_layout(
			xaxis = dict(
				range = range_plot,
				title_text = 'Position',
				title_font = dict(
					color = 'white'
				),
				tickfont = dict(
					color = 'white'
				)
			),
			yaxis = dict(
				title_text = "Beta-value",
				title_font = dict(
					color = 'white'
				),
				tickfont = dict(
					color = 'white'
				)
			),
			height=750+(len_drops*500),
			paper_bgcolor = '#393833',
			legend = dict(
				font = dict(
					color = 'white'
				)
			)
	#         hovermode='x unified'
		)

	config = {
		'config': {
			"start": start,
			"end": end
		}
	}

	graph = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
	graph_comp = json.loads(graph)
	graph_comp.update(config)
	graphJSON = json.dumps(graph_comp)

	return graphJSON