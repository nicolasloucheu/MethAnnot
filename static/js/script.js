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
	var enh_val = $('#enhancers').prop('checked'); /* it will return true or false */

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
            Plotly.newPlot('chart', data);
        }
    });
    $('#TF_drop').selectpicker('refresh');
})


$('#chart').on('plotly_relayout',function(){
	var gd = document.getElementById('chart')
	var xRange = gd.layout.xaxis.range;
	var yRange = gd.layout.yaxis.range;
    var TF_multiselect = $('#TF_drop').val();
    var cpg_multiselect = $('#cpg_annots').val();
    var hmm_multiselect = $('#chromhmm').val();
    var enh_val = $('#enhancers').prop('checked'); /* it will return true or false */

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
            Plotly.newPlot('chart', data);
        }
    });
})


$('input:checkbox').on('click', function(){
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
            Plotly.newPlot('chart', data);
        }
    });

});