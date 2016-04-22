var params =  {
  // 'x': 30,
  // 'y': 50,
  'line-width': 3,
  'line-length': 40,
  'text-margin': 10,
  'font-size': 14,
  'font': 'normal',
  'font-family': 'Helvetica',
  'font-weight': 'normal',
  'font-color': 'black',
  'line-color': 'grey',
  'element-color': 'black',
  'fill': 'white',
  'yes-text': 'yes',
  'no-text': 'no',
  'arrow-end': 'block',
  'scale': 1,
  'symbols': {
    'inputoutput': {
      'font-color': 'red',
      'element-color': 'red',
      'fill': 'yellow'
    },
    'operation':{
      'fill': 'steelblue'
    }
  },
  'flowstate' : {
    'analysis' : { 'fill' : '#C45879', 'font-size' : 12, 'yes-text' : 'visualize', 'no-text' : 'differential expression' },
    'degs' : { 'fill' : '#C45879', 'font-size' : 12, 'yes-text' : 'enrich against\ngene sets', 'no-text' : 'search for\nsmall molecules' },
    'visualization' : { 'fill' : '#C45879', 'font-size' : 12, 'yes-text' : 'static', 'no-text' : 'interactive' }
  }
};

$(document).ready(function() {
    $.ajax({
        url : "flowchart.txt",
        dataType : "text",
        success : function(data) {
            var chart = flowchart.parse(data);
            chart.drawSVG('canvas', params);
        }
    });
});

