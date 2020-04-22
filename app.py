from flask import Flask, render_template, request, session, jsonify
import json
import plotly
import plotly.offline as py
import plotly.graph_objs as go
import os
import re
from plottin import create_dfs, create_plot
import pandas as pd
import pickle
 
app = Flask(__name__)



chrom_lst = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M"]
samples = []
uploads_dir = os.path.join(app.instance_path, 'uploads')
if not os.path.exists(uploads_dir):
	os.makedirs(uploads_dir)
 

# Main page
@app.route('/', methods = ['POST', 'GET'])
def home():
	# Defining return variables
	TF_options = []
	cpg_options = []
	hmm_options = []
	enh_dis = True
	region = ''
	error_region = ''
	region_mes = ''
	ready_to_plot = False
	start = 0
	end = 0
	col_sample = {}
	col_check_val = 'off'

	# If a file is given by the user
	if request.method == "POST" and 'file' in request.files:
		file = request.files.getlist('file')
		if len(file) > 1:
			init_dir = file[0].filename.split('/')[0]
			dir_name = os.path.join(uploads_dir, init_dir)
			samples.append(init_dir)
			if not os.path.exists(dir_name):
				os.makedirs(dir_name)
			for i in file:
				i.save(os.path.join(uploads_dir, i.filename))

	# For each sample already added, check if it is deleted or the color is changed
	for sample in samples:
		if request.method == "POST" and f'del_{sample}' in request.form:
			samples.remove(sample)

		if request.method == "POST" and f'onoffswitch_{sample}' in request.form:
			col_sample[sample] = "#00f"
		else:
			col_sample[sample] = '#297822'

	# If a region has been given by the user
	if request.method == "POST" and 'region-input' in request.form:
		region = request.form['region-input']
		if region != '':
			# Try to parse the region and return appropriate message
			try:
				chrom = re.search('chr(.*):', region).group(1)
				start = re.search(':(.*)-', region).group(1).replace(',', '')
				end = region.split('-')[-1].replace(',', '')
				if chrom not in chrom_lst:
					error_region = 'error'
					region_mes = f'{chrom} is not a chromosome. Chromosome must be between 1 and 22, or X or Y.'
				elif start > end:
					error_region = 'error'
					region_mes = 'Start must be lower than end'
				else:
					if len(samples) > 0:
						ready_to_plot = True
					else:
						region_mes = 'You must provide at least one sample to plot'
			except:
				error_region = 'error'
				region_mes = 'You have to provide the region as chr{chr}:{start}-{end}. Example: chr8:60,670,000-60,700,000'

	# If at least one sample and a correct region have been given by the user, do all the computation and return the graph
	if ready_to_plot:
		# Create dataframes from initial information (chromosome, start, end and samples paths)
		TF_options, cpg_options, hmm_options, enh_dis, bv_means_controls, bv_sample, sub_genes, df_TF, df_annots = create_dfs(chrom, start, end, samples)
		#Store all these informations in session for the variables
		session['chrom'] = chrom
		session['start'] = start
		session['end'] = end
		session['samples'] = samples
		session['col_sample'] = col_sample
		# Store in a csv file for dataframes
		bv_means_controls.to_csv('instance/tmp/bv_means_controls.csv.gz', compression='gzip')
		bv_json = [i.to_json(orient='split') for i in bv_sample]
		pickle.dump(bv_json, open("instance/tmp/bv_json.p", "wb"))
		sub_genes.to_csv('instance/tmp/sub_genes.csv.gz', compression='gzip')
		df_TF.to_csv('instance/tmp/df_TF.csv.gz', compression='gzip')
		df_annots.to_csv('instance/tmp/df_annots.csv.gz', compression='gzip')

		# Create the plot from that information
		plotly_plot = create_plot(bv_means_controls, bv_sample, sub_genes, df_TF, df_annots, start, end, chrom, None, None, None, None, None, None, col_sample)
	else:
		plotly_plot = ""

	#Rener everything in the html file
	return render_template("index.html", samples=samples, region=region, error_region=error_region, region_mes=region_mes, plotly_plot=plotly_plot, ready_to_plot=ready_to_plot, TF_options=TF_options, cpg_options=cpg_options, hmm_options=hmm_options, enh_dis=enh_dis, start=start, end=end, col_sample=col_sample)


# Hidden route to compute the resulting graph when some annotations are added. It will return only the graph, not the html path.
@app.route('/update_graph', methods=['GET', 'POST'])
def change_features():
	# Get all the information from session or csv files
	chrom = session.get('chrom', None)
	start = session.get('start', None)
	end = session.get('end', None)
	samples = session.get('samples', None)
	col_sample = session.get('col_sample', None)
	bv_means_controls = pd.read_csv('instance/tmp/bv_means_controls.csv.gz', index_col=0, compression='gzip')
	bv_sample = [pd.read_json(i, orient='split') for i in pickle.load(open("instance/tmp/bv_json.p", "rb"))]
	sub_genes = pd.read_csv('instance/tmp/sub_genes.csv.gz', index_col=0, compression='gzip')
	df_TF = pd.read_csv('instance/tmp/df_TF.csv.gz', index_col=0, compression='gzip')
	df_annots = pd.read_csv('instance/tmp/df_annots.csv.gz', index_col=0, compression='gzip')

	# Get zoom/move information from the plot
	x_range = [float(i) for i in request.args.getlist('xrange[]')]
	y_range = [float(i) for i in request.args.getlist('yrange[]')]

	# Get dropdowns and checkbox annotations information
	TF_value = request.args.getlist('TF_value[]')
	cpg_value = request.args.getlist('cpg_value[]')
	hmm_value = request.args.getlist('hmm_value[]')
	enh_val = request.args['enh_val']
	# And store that information
	session['TF_value'] = TF_value
	session['cpg_value'] = cpg_value
	session['hmm_value'] = hmm_value
	session['enh_val'] = enh_val
	graphJSON= create_plot(bv_means_controls, bv_sample, sub_genes, df_TF, df_annots, start, end, chrom, TF_value, cpg_value, hmm_value, enh_val, x_range, y_range, col_sample)
	return graphJSON


# Hidden route to compute the resulting graph when the user modifies the graph (zoom, move, ...). It will return only the graph, not the html path.
@app.route('/update_zoom', methods=['GET', 'POST'])
def change_zoom():

	# Get all the information from session or csv files
	chrom = session.get('chrom', None)
	start = session.get('start', None)
	end = session.get('end', None)
	samples = session.get('samples', None)
	col_sample = session.get('col_sample', None)
	bv_means_controls = pd.read_csv('instance/tmp/bv_means_controls.csv.gz', index_col=0, compression='gzip')
	bv_sample = [pd.read_json(i, orient='split') for i in pickle.load(open("instance/tmp/bv_json.p", "rb"))]
	sub_genes = pd.read_csv('instance/tmp/sub_genes.csv.gz', index_col=0, compression='gzip')
	df_TF = pd.read_csv('instance/tmp/df_TF.csv.gz', index_col=0, compression='gzip')
	df_annots = pd.read_csv('instance/tmp/df_annots.csv.gz', index_col=0, compression='gzip')
	TF_value = session.get('TF_value', None)
	cpg_value = session.get('cpg_value', None)
	hmm_value = session.get('hmm_value', None)
	enh_val = session.get('enh_val', None)

	# Get zoom/move information from the plot
	x_range = [float(i) for i in request.args.getlist('xrange[]')]
	y_range = [float(i) for i in request.args.getlist('yrange[]')]

	# Get dropdowns and checkbox annotations information and store it again
	TF_value = request.args.getlist('TF_value[]')
	cpg_value = request.args.getlist('cpg_value[]')
	hmm_value = request.args.getlist('hmm_value[]')
	enh_val = request.args['enh_val']
	session['TF_value'] = TF_value
	session['cpg_value'] = cpg_value
	session['hmm_value'] = hmm_value
	session['enh_val'] = enh_val
	graphJSON= create_plot(bv_means_controls, bv_sample, sub_genes, df_TF, df_annots, start, end, chrom, TF_value, cpg_value, hmm_value, enh_val, x_range, y_range, col_sample)
	return graphJSON


# Hidden route to compute the resulting graph when a sample color is changed. It will return only the graph, not the html path.
@app.route('/update_color', methods=['GET', 'POST'])
def change_color():

	# Get all the information from session or csv files
	chrom = session.get('chrom', None)
	start = session.get('start', None)
	end = session.get('end', None)
	samples = session.get('samples', None)
	col_sample = session.get('col_sample', None)
	bv_means_controls = pd.read_csv('instance/tmp/bv_means_controls.csv.gz', index_col=0, compression='gzip')
	bv_sample = [pd.read_json(i, orient='split') for i in pickle.load(open("instance/tmp/bv_json.p", "rb"))]
	sub_genes = pd.read_csv('instance/tmp/sub_genes.csv.gz', index_col=0, compression='gzip')
	df_TF = pd.read_csv('instance/tmp/df_TF.csv.gz', index_col=0, compression='gzip')
	df_annots = pd.read_csv('instance/tmp/df_annots.csv.gz', index_col=0, compression='gzip')
	TF_value = session.get('TF_value', None)
	cpg_value = session.get('cpg_value', None)
	hmm_value = session.get('hmm_value', None)
	enh_val = session.get('enh_val', None)

	# Get zoom/move information from the plot
	x_range = [float(i) for i in request.args.getlist('xrange[]')]
	y_range = [float(i) for i in request.args.getlist('yrange[]')]

	# Get color information of samples
	color_sample = request.args['color']
	id_checkbox = request.args['id_checkbox'].split('switch_')[-1]
	col_sample[id_checkbox] = color_sample
	session['col_sample'] = col_sample

	graphJSON = create_plot(bv_means_controls, bv_sample, sub_genes, df_TF, df_annots, start, end, chrom, TF_value, cpg_value, hmm_value, enh_val, x_range, y_range, col_sample)
	return graphJSON


# Hidden route to compute the resulting graph when the region is shifted (left or right). It will return only the graph, not the html path.
@app.route('/update_region', methods = ['GET', 'POST'])
def change_region():
	# Get all the information from session
	new_start = request.args['start']
	new_end = request.args['end']
	chrom = session.get('chrom', None)
	samples = session.get('samples', None)
	col_sample = session.get('col_sample', None)
	session['TF_value'] = None
	session['cpg_value'] = None
	session['hmm_value'] = None
	session['enh_val'] = None

	# Modifying the start and end session information
	session['start'] = new_start
	session['end'] = new_end

	# Creating the dataframes for new region
	TF_options, cpg_options, hmm_options, enh_dis, bv_means_controls, bv_sample, sub_genes, df_TF, df_annots = create_dfs(chrom, new_start, new_end, samples)

	# Put possible annotation for this new region in a dict to append to JSON object
	options_dict = {"TF_options": TF_options, "cpg_options": cpg_options, "hmm_options": hmm_options, "enh_dis": enh_dis}

	# Saving dataframes to csv files
	bv_means_controls.to_csv('instance/tmp/bv_means_controls.csv.gz', index_col=0, compression='gzip')
	bv_json = [i.to_json(orient='split') for i in bv_sample]
	pickle.dump(bv_json, open("instance/tmp/bv_json.p", "wb"))
	sub_genes.to_csv('instance/tmp/sub_genes.csv.gz', index_col=0, compression='gzip')
	df_TF.to_csv('instance/tmp/df_TF.csv.gz', index_col=0, compression='gzip')
	df_annots.to_csv('instance/tmp/df_annots.csv.gz', index_col=0, compression='gzip')

	# Create the plot
	plotly_plot = create_plot(bv_means_controls, bv_sample, sub_genes, df_TF, df_annots, new_start, new_end, chrom, None, None, None, None, None, None, col_sample)

	# Adding the new possible annotations information
	graphs = json.loads(plotly_plot)
	graphs.update(options_dict)
	graphJSON = json.dumps(graphs)

	return graphJSON



if __name__ == "__main__":
	app.secret_key='secret123'
	app.run(debug = True)