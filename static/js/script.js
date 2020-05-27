// TOOLTIPS
$(document).ready(function(){
    $('[data-toggle="tooltip"]').tooltip();   
});

// Sidebar toggle behavior
$(function() {
	$('#toggle-btn').on('click', function() {
		$('#sidebar, #content').toggleClass('active');
	});
});


//Graph update when changing annotations
$('#update-btn').on('click',function(){
	var gd = document.getElementById('chart')
	var xRange = gd.layout.xaxis2.range;
	var yRange = gd.layout.yaxis2.range;
	var annots_names = document.getElementsByClassName('selectpicker');
	var annots = {};

	for (i = 0; i < annots_names.length; i++) {
		annots[annots_names[i].id] = ($("#" + annots_names[i].id).val());
	}

	$.ajax({
		url: "/update_graph",
		type: "GET",
		contentType: 'application/json;charset=UTF-8',
		data: {
			'annots': JSON.stringify(annots),
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
								Plotly.relayout(gd, 'xaxis2.range', [data.config.start, data.config.end])
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


//Graph update when touching the graph (changing zoom, moving)
$('#chart').on('plotly_relayout',function(){
	var gd = document.getElementById('chart')
	var xRange = gd.layout.xaxis2.range;
	var yRange = gd.layout.yaxis2.range;
	var annots_names = document.getElementsByClassName('selectpicker');
	var annots = {};

	for (i = 0; i < annots_names.length; i++) {
		annots[annots_names[i].id] = ($("#" + annots_names[i].id).val());
	}

	$.ajax({
		url: "/update_zoom",
		type: "GET",
		contentType: 'application/json;charset=UTF-8',
		data: {
			'xrange': xRange,
			'yrange': yRange,
			'annots': JSON.stringify(annots)
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
								Plotly.relayout(gd, 'xaxis2.range', [data.config.start, data.config.end])
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


//Graph update when changing color of a sample
$('.onoffswitch-checkbox:checkbox').on('click', function(){
	var gd = document.getElementById('chart')
	var xRange = gd.layout.xaxis2.range;
	var yRange = gd.layout.yaxis2.range;
	var annots_names = document.getElementsByClassName('selectpicker');
	var annots = {};

	for (i = 0; i < annots_names.length; i++) {
		annots[annots_names[i].id] = ($("#" + annots_names[i].id).val());
	}


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
			'yrange': yRange,
			'annots': JSON.stringify(annots)
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



//Graph update when shifting to the left
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
									Plotly.relayout(gd, 'xaxis2.range', [data.config.start, data.config.end])
								}
							}	
						]
					],
					toImageButtonOptions: {width: 1920, height: len_drops}
				}
				Plotly.newPlot('chart', data.data, data.layout, config);
				var annots_names = document.getElementsByClassName('selectpicker');
				var annots_values = {}

				for (i = 0; i < annots_names.length; i++) {
					annots_values[annots_names[i].id] = data[annots_names[i].id];
				};
				$(".selectpicker").selectpicker();
				for (i = 0; i < annots_names.length; i++) {
					$("#" + annots_names[i].id).html('');
					$.each(annots_values[annots_names[i].id], function(value, item){
						$("#" + annots_names[i].id).append('<option>' + item + '</option>');
					});
				};
				$('.selectpicker').selectpicker('refresh');

				document.getElementById("region-show").innerHTML = "REGION: " + chrom + ":" + new_start + "-" + new_end
			}
		});
	}
});



//Graph update when shifting to the right
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
									Plotly.relayout(gd, 'xaxis2.range', [data.config.start, data.config.end])
								}
							}	
						]
					],
					toImageButtonOptions: {width: 1920, height: len_drops}
				}
				Plotly.newPlot('chart', data.data, data.layout, config);
				var annots_names = document.getElementsByClassName('selectpicker');
				var annots_values = {}

				for (i = 0; i < annots_names.length; i++) {
					annots_values[annots_names[i].id] = data[annots_names[i].id];
				};
				$(".selectpicker").selectpicker();
				for (i = 0; i < annots_names.length; i++) {
					$("#" + annots_names[i].id).html('');
					$.each(annots_values[annots_names[i].id], function(value, item){
						$("#" + annots_names[i].id).append('<option>' + item + '</option>');
					});
				};
				$('.selectpicker').selectpicker('refresh');
				document.getElementById("region-show").innerHTML = "REGION: " + chrom + ":" + new_start + "-" + new_end
			}
		});
	}
});

$('.top_z').change(function(event){
	var z_region = event.target.value;
	var new_chrom_tmp = z_region.split("CHR: ")[1].split(" -")[0]
	if (isNaN(new_chrom_tmp)) {
		var new_chrom = new_chrom_tmp
	}else{
		var new_chrom = Math.round(new_chrom_tmp)
	}
	var new_start = Number(z_region.split("Position: ")[1].split(".")[0]) - 50000
	var new_end = Number(z_region.split("Position: ")[1].split(".")[0]) + 50000

	$("#region-input").val('chr' + new_chrom + ":" + new_start + "-" + new_end)

	$.ajax({
		url: "show_z_region",
		type:'GET',
		contentType: 'application/json;charset=UTF-8',
		data: {
			"chrom": new_chrom,
			"start": new_start,
			'end': new_end
		},
		dataType:"json",
		success: function (data) {
			var len_drops = 750 + 0*500
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
								Plotly.relayout(gd, 'xaxis2.range', [data.config.start, data.config.end])
							}
						}	
					]
				],
				toImageButtonOptions: {width: 1920, height: len_drops}
			}
			Plotly.newPlot('chart', data.data, data.layout, config);
			var annots_names = document.getElementsByClassName('selectpicker');
			var annots_values = {}

			for (i = 0; i < annots_names.length; i++) {
				annots_values[annots_names[i].id] = data[annots_names[i].id];
			};
			$(".selectpicker").selectpicker();
			for (i = 0; i < annots_names.length; i++) {
				$("#" + annots_names[i].id).html('');
				$.each(annots_values[annots_names[i].id], function(value, item){
					$("#" + annots_names[i].id).append('<option>' + item + '</option>');
				});
			};
			$('.selectpicker').selectpicker('refresh');
			document.getElementById("region-show").innerHTML = "REGION: chr" + new_chrom + ":" + new_start + "-" + new_end
		}
	});
});