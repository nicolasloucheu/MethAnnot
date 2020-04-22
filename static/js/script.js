$(document).ready(function(){
    $('[data-toggle="tooltip"]').tooltip();   
});

$(function() {
	// Sidebar toggle behavior
	$('#toggle-btn').on('click', function() {
		$('#sidebar, #content').toggleClass('active');
	});
});

$('#update-btn').on('click',function(){
	var gd = document.getElementById('chart')
	var xRange = gd.layout.xaxis.range;
	var yRange = gd.layout.yaxis.range;
	var TF_multiselect = $('#TF_drop').val();
	var cpg_multiselect = $('#cpg_annots').val();
	var hmm_multiselect = $('#chromhmm').val();

	if ($('.enhancers').prop('checked') == null) {
		var enh_val = false
	}
	else {
		var enh_val = $('.enhancers').prop('checked')
	}

	$.ajax({
		url: "/update_graph",
		type: "GET",
		contentType: 'application/json;charset=UTF-8',
		data: {
			'TF_value': TF_multiselect,
			'cpg_value': cpg_multiselect,
			'hmm_value': hmm_multiselect,
			'enh_val': enh_val,
			'xrange': xRange,
			'yrange': yRange

		},
		dataType:"json",
		success: function (data) {
			var len_drops = 750 + data.config.len_drops*500
			var config = {
				responsive: true,
				showTips: false,
				modeBarButtons: [
					[
						'toImage',
					],
					[
						'zoomIn2d',
						'zoomOut2d',
					],
					[
						{
							name: 'myResetScale2d',
							title: 'Reset axes',
							icon: Plotly.Icons.home,
							click: function(gd) {
								Plotly.relayout(gd, 'xaxis.range', [data.config.start, data.config.end])
							}
						}	
					]
				],
				toImageButtonOptions: {width: 1920, height: len_drops}
			}
			Plotly.newPlot('chart', data.data, data.layout, config);
		}
	});
})


$('#chart').on('plotly_relayout',function(){
	var gd = document.getElementById('chart')
	var xRange = gd.layout.xaxis.range;
	var yRange = gd.layout.yaxis.range;
	var TF_multiselect = $('#TF_drop').val();
	var cpg_multiselect = $('#cpg_annots').val();
	var hmm_multiselect = $('#chromhmm').val();
	
	if ($('#enhancers').prop('checked') == null) {
		var enh_val = false
	}
	else {
		var enh_val = $('#enhancers').prop('checked')
	}

	$.ajax({
		url: "/update_zoom",
		type: "GET",
		contentType: 'application/json;charset=UTF-8',
		data: {
			'xrange': xRange,
			'yrange': yRange,
			'TF_value': TF_multiselect,
			'cpg_value': cpg_multiselect,
			'hmm_value': hmm_multiselect,
			'enh_val': enh_val
		},
		dataType:"json",
		success: function (data) {
			var len_drops = 750 + data.config.len_drops*500
			var config = {
				responsive: true,
				showTips: false,
				modeBarButtons: [
					[
						'toImage',
					],
					[
						'zoomIn2d',
						'zoomOut2d',
					],
					[
						{
							name: 'myResetScale2d',
							title: 'Reset axes',
							icon: Plotly.Icons.home,
							click: function(gd) {
								Plotly.relayout(gd, 'xaxis.range', [data.config.start, data.config.end])
							}
						}	
					]
				],
				toImageButtonOptions: {width: 1920, height: len_drops}
			}
			Plotly.newPlot('chart', data.data, data.layout, config);
		}
	});
})


$('.onoffswitch-checkbox:checkbox').on('click', function(){
	var gd = document.getElementById('chart')
	var xRange = gd.layout.xaxis.range;
	var yRange = gd.layout.yaxis.range;
	var $this = $(this)

	if ($this.is(":checked")) {
		var id_checkbox = $this.attr("id");
		var color = "#00f"
	} 
	else {
		var id_checkbox = $this.attr("id");
		var color = "#297822"
	}

	$.ajax({
		url: "update_color",
		type:'GET',
		contentType: 'application/json;charset=UTF-8',
		data: {
			"id_checkbox": id_checkbox,
			"color": color,
			'xrange': xRange,
			'yrange': yRange
		},
		dataType:"json",
		success: function (data) {
			var len_drops = 750 + data.config.len_drops*500
			var config = {
				responsive: true,
				showTips: false,
				modeBarButtons: [
					[
						'toImage',
					],
					[
						'zoomIn2d',
						'zoomOut2d',
					],
					[
						{
							name: 'myResetScale2d',
							title: 'Reset axes',
							icon: Plotly.Icons.home,
							click: function(gd) {
								Plotly.relayout(gd, 'xaxis.range', [data.config.start, data.config.end])
							}
						}	
					]
				],
				toImageButtonOptions: {width: 1920, height: len_drops}
			}
			Plotly.newPlot('chart', data.data, data.layout, config);
		}
	});

});



$('#shift-left').on('click', function() {
	var region = document.getElementById('region-input').value
	var shift = document.getElementById('region-shift').value


	if ( region && shift ) {
		var chrom = region.split(":")[0]
		var start = region.split(":")[1].split("-")[0]
		var end = region.split(":")[1].split("-")[1]

		var new_start = parseInt(start) - parseInt(shift)
		var new_end = parseInt(end) - parseInt(shift)

		$("#region-input").val(chrom + ":" + new_start + "-" + new_end)

		$(".selectpicker").selectpicker();
		$('.selectpicker').val(null)
		$('.selectpicker').selectpicker('refresh');

		$.ajax({
			url: "update_region",
			type:'GET',
			contentType: 'application/json;charset=UTF-8',
			data: {
				"chrom": chrom,
				"start": new_start,
				'end': new_end
			},
			dataType:"json",
			success: function (data) {
				var len_drops = 750 + data.config.len_drops*500
				var config = {
					responsive: true,
					showTips: false,
					modeBarButtons: [
						[
							'toImage',
						],
						[
							'zoomIn2d',
							'zoomOut2d',
						],
						[
							{
								name: 'myResetScale2d',
								title: 'Reset axes',
								icon: Plotly.Icons.home,
								click: function(gd) {
									Plotly.relayout(gd, 'xaxis.range', [data.config.start, data.config.end])
								}
							}	
						]
					],
					toImageButtonOptions: {width: 1920, height: len_drops}
				}
				Plotly.newPlot('chart', data.data, data.layout, config);
				var TF_options = data.TF_options
				var cpg_options = data.cpg_options
				var hmm_options = data.hmm_options
				var enh_dis = data.enh_dis
				var $TF_drop = $('#TF_drop');
				var $cpg_annots = $('#cpg_annots');
				var $chromhmm = $('#chromhmm');
				$(".selectpicker").selectpicker();
				$TF_drop.html('');
				$cpg_annots.html('');
				$chromhmm.html('');
				$.each(TF_options, function(value, item) {
					$TF_drop.append('<option>' + item + '</option>');
				});
				$.each(cpg_options, function(value, item) {
					$cpg_annots.append('<option>' + item + '</option>');
				});
				$.each(hmm_options, function(value, item) {
					$chromhmm.append('<option>' + item + '</option>');
				});
				$('.selectpicker').selectpicker('refresh');

				if (enh_dis) {
					$(".enhancers").prop("disabled", true);
				} else {
					$(".enhancers").prop("disabled", false);
				}
			}
		});
	}
});




$('#shift-right').on('click', function() {
	var region = document.getElementById('region-input').value
	var shift = document.getElementById('region-shift').value


	if ( region && shift ) {
		var chrom = region.split(":")[0]
		var start = region.split(":")[1].split("-")[0]
		var end = region.split(":")[1].split("-")[1]

		var new_start = parseInt(start) + parseInt(shift)
		var new_end = parseInt(end) + parseInt(shift)

		$("#region-input").val(chrom + ":" + new_start + "-" + new_end)

		$.ajax({
			url: "update_region",
			type:'GET',
			contentType: 'application/json;charset=UTF-8',
			data: {
				"chrom": chrom,
				"start": new_start,
				'end': new_end
			},
			dataType:"json",
			success: function (data) {
				var len_drops = 750 + data.config.len_drops*500
				var config = {
					responsive: true,
					showTips: false,
					modeBarButtons: [
						[
							'toImage',
						],
						[
							'zoomIn2d',
							'zoomOut2d',
						],
						[
							{
								name: 'myResetScale2d',
								title: 'Reset axes',
								icon: Plotly.Icons.home,
								click: function(gd) {
									Plotly.relayout(gd, 'xaxis.range', [data.config.start, data.config.end])
								}
							}	
						]
					],
					toImageButtonOptions: {width: 1920, height: len_drops}
				}
				Plotly.newPlot('chart', data.data, data.layout, config);
				var TF_options = data.TF_options
				var cpg_options = data.cpg_options
				var hmm_options = data.hmm_options
				var enh_dis = data.enh_dis
				var $TF_drop = $('#TF_drop');
				var $cpg_annots = $('#cpg_annots');
				var $chromhmm = $('#chromhmm');
				$(".selectpicker").selectpicker();
				$TF_drop.html('');
				$cpg_annots.html('');
				$chromhmm.html('');
				$.each(TF_options, function(value, item) {
					$TF_drop.append('<option>' + item + '</option>');
				});
				$.each(cpg_options, function(value, item) {
					$cpg_annots.append('<option>' + item + '</option>');
				});
				$.each(hmm_options, function(value, item) {
					$chromhmm.append('<option>' + item + '</option>');
				});
				$('.selectpicker').selectpicker('refresh');

				if (enh_dis) {
					$(".enhancers").prop("disabled", true);
				} else {
					$(".enhancers").prop("disabled", false);
				}
			}
		});
	}
});