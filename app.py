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
tmp_dir = os.path.join("instance", 'tmp')
if not os.path.exists(tmp_dir):
	os.makedirs(tmp_dir)

# Main page
@app.route('/', methods = ['POST', 'GET'])
def home():
	# Defining return variables
	region = ''
	error_region = ''
	region_mes = ''
	ready_to_plot = False
	start = 0
	end = 0
	col_sample = session.get('col_sample', {})
	col_check_val = 'off'
	annots_values = {}



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
			col_sample[init_dir] = '#297822'

	# For each sample already added, check if it is deleted or the color is changed
	for sample in samples:
		if request.method == "POST" and f'del_{sample}' in request.form:
			samples.remove(sample)


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
		bv_means_controls, bv_sample, z_scores, sub_genes, annots, annots_names, total_options = create_dfs(chrom, start, end, samples)
		#Store all these informations in session for the variables
		session['chrom'] = chrom
		session['start'] = start
		session['end'] = end
		session['samples'] = samples
		session['col_sample'] = col_sample
		# Store in a csv file for dataframes
		bv_means_controls.to_csv('instance/tmp/bv_means_controls.csv.gz', compression='gzip')
		annots.to_csv('instance/tmp/annots.csv.gz', compression='gzip')
		bv_json = [i.to_json(orient='split') for i in bv_sample]
		pickle.dump(bv_json, open("instance/tmp/bv_json.p", "wb"))
		z_json = [i.to_json(orient='split') for i in z_scores]
		pickle.dump(z_json, open("instance/tmp/z_json.p", "wb"))
		sub_genes.to_csv('instance/tmp/sub_genes.csv.gz', compression='gzip')

		# Create the plot from that information
		plotly_plot = create_plot(bv_means_controls, bv_sample, z_scores, sub_genes, start, end, chrom, None, None, col_sample, annots, annots_values)
	else:
		plotly_plot = ""
		annots_names = []
		total_options = []
		annots_json = ""
		names_corrected = []

	top_z = []
	mean_lst = {}
	for sample_value in samples:
		top_z.append(pd.read_csv(f"instance/uploads/{sample_value}/top_z_scores.csv.gz", compression="gzip", index_col=0, nrows=100))
		with open(f"instance/uploads/{sample_value}/z_score_mean.pickle", 'rb') as fp:
			z_mean = pickle.load(fp)
		mean_lst[sample_value] = z_mean

	#Rener everything in the html file
	return render_template("index.html", samples=samples, region=region, error_region=error_region, region_mes=region_mes, plotly_plot=plotly_plot, ready_to_plot=ready_to_plot, start=start, end=end, col_sample=col_sample, top_z=top_z, mean_lst=mean_lst, annots_names=annots_names, total_options=total_options)


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
	z_scores = [pd.read_json(i, orient='split') for i in pickle.load(open("instance/tmp/z_json.p", "rb"))]
	sub_genes = pd.read_csv('instance/tmp/sub_genes.csv.gz', index_col=0, compression='gzip')
	annots = pd.read_csv('instance/tmp/annots.csv.gz', compression='gzip', index_col=0)

	# Get zoom/move information from the plot
	x_range = [float(i) for i in request.args.getlist('xrange[]')]
	y_range = [float(i) for i in request.args.getlist('yrange[]')]

	# Get dropdowns and checkbox annotations information
	annots_values = json.loads(request.args.getlist('annots')[0])
	# And store that information
	session['annots_values'] = annots_values
	graphJSON= create_plot(bv_means_controls, bv_sample, z_scores, sub_genes, start, end, chrom, x_range, y_range, col_sample, annots, annots_values)
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
	z_scores = [pd.read_json(i, orient='split') for i in pickle.load(open("instance/tmp/z_json.p", "rb"))]
	sub_genes = pd.read_csv('instance/tmp/sub_genes.csv.gz', index_col=0, compression='gzip')
	annots = pd.read_csv('instance/tmp/annots.csv.gz', compression='gzip', index_col=0)


	# Get zoom/move information from the plot
	x_range = [float(i) for i in request.args.getlist('xrange[]')]
	y_range = [float(i) for i in request.args.getlist('yrange[]')]

	# Get dropdowns and checkbox annotations information and store it again
	annots_values = json.loads(request.args.getlist('annots')[0])
	session['annots_values'] = annots_values
	graphJSON= create_plot(bv_means_controls, bv_sample, z_scores, sub_genes, start, end, chrom, x_range, y_range, col_sample, annots, annots_values)
	
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
	z_scores = [pd.read_json(i, orient='split') for i in pickle.load(open("instance/tmp/z_json.p", "rb"))]
	sub_genes = pd.read_csv('instance/tmp/sub_genes.csv.gz', index_col=0, compression='gzip')
	annots = pd.read_csv('instance/tmp/annots.csv.gz', compression='gzip', index_col=0)

	# Get zoom/move information from the plot
	x_range = [float(i) for i in request.args.getlist('xrange[]')]
	y_range = [float(i) for i in request.args.getlist('yrange[]')]

	# Get color information of samples
	color_sample = request.args['color']
	id_checkbox = request.args['id_checkbox'].split('switch_')[-1]
	col_sample[id_checkbox] = color_sample
	session['col_sample'] = col_sample

	annots_values = json.loads(request.args.getlist('annots')[0])
	session['annots_values'] = annots_values

	graphJSON= create_plot(bv_means_controls, bv_sample, z_scores, sub_genes, start, end, chrom, x_range, y_range, col_sample, annots, annots_values)
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
	session['annots_values'] = None

	# Modifying the start and end session information
	session['start'] = new_start
	session['end'] = new_end

	# Creating the dataframes for new region
	bv_means_controls, bv_sample, z_scores, sub_genes, annots, annots_names, total_options = create_dfs(chrom, new_start, new_end, samples)

	# Saving dataframes to csv files
	bv_means_controls.to_csv('instance/tmp/bv_means_controls.csv.gz', compression='gzip')
	bv_json = [i.to_json(orient='split') for i in bv_sample]
	pickle.dump(bv_json, open("instance/tmp/bv_json.p", "wb"))
	z_json = [i.to_json(orient='split') for i in z_scores]
	pickle.dump(z_json, open("instance/tmp/z_json.p", "wb"))
	sub_genes.to_csv('instance/tmp/sub_genes.csv.gz', compression='gzip')
	annots.to_csv('instance/tmp/annots.csv.gz', compression='gzip')

	# Create the plot
	plotly_plot = create_plot(bv_means_controls, bv_sample, z_scores, sub_genes, new_start, new_end, chrom, None, None, col_sample, annots, {})

	# Adding the new possible annotations information
	graphs = json.loads(plotly_plot)
	graphs.update(total_options)
	graphJSON = json.dumps(graphs)

	return graphJSON



# Hidden route to compute the resulting graph when the region is selected from the top z-scores
@app.route('/show_z_region', methods = ['GET', 'POST'])
def z_region():
	# Get all the information from session
	new_chrom = request.args['chrom']
	new_start = request.args['start']
	new_end = request.args['end']
	samples = session.get('samples', None)
	col_sample = session.get('col_sample', None)
	session['annots_values'] = None

	# Modifying the start and end session information
	session['start'] = new_start
	session['end'] = new_end
	session['chrom'] = new_chrom

	# Creating the dataframes for new region
	bv_means_controls, bv_sample, z_scores, sub_genes, annots, annots_names, total_options = create_dfs(new_chrom, new_start, new_end, samples)

	# Saving dataframes to csv files
	bv_means_controls.to_csv('instance/tmp/bv_means_controls.csv.gz', compression='gzip')
	bv_json = [i.to_json(orient='split') for i in bv_sample]
	pickle.dump(bv_json, open("instance/tmp/bv_json.p", "wb"))
	z_json = [i.to_json(orient='split') for i in z_scores]
	pickle.dump(z_json, open("instance/tmp/z_json.p", "wb"))
	sub_genes.to_csv('instance/tmp/sub_genes.csv.gz', compression='gzip')
	annots.to_csv('instance/tmp/annots.csv.gz', compression='gzip')

	# Create the plot
	plotly_plot = create_plot(bv_means_controls, bv_sample, z_scores, sub_genes, new_start, new_end, new_chrom, None, None, col_sample, annots, {})

	# Adding the new possible annotations information
	graphs = json.loads(plotly_plot)
	graphs.update(total_options)
	graphJSON = json.dumps(graphs)

	return graphJSON




if __name__ == "__main__":
	app.secret_key='secret123'
	app.run(debug = True)