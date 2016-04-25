var params =  {
  // 'x': 30,
  // 'y': 50,
  'line-width': 2,
  'line-length': 0.01,
  'text-margin': 3,
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
  'scale': 0.9,
  'symbols': {
    'inputoutput': {
      'font-color': 'red',
      'element-color': 'red',
      'fill': 'yellow'
    },
    'operation':{
      'fill': 'lightblue'
    },
    'condition':{
      'line-length': 10,
    }
  },
  'flowstate' : {
    'analysis' : { 'fill' : '#C45879', 'font-size' : 12, 'yes-text' : '\nvisualize', 'no-text' : '\ndifferential expression' },
    'degs' : {'line-length':60,  'fill' : '#C45879', 'font-size' : 12, 'yes-text' : '\nenrich against\ngene sets', 'no-text' : 'search for\nsmall molecules' },
    'visualization' : {'line-length':30, 'fill' : '#C45879', 'font-size' : 12, 'yes-text' : '\nstatic', 'no-text' : '\ninteractive' }
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

