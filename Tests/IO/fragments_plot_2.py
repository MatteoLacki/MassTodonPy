from highcharts import Highchart
H = Highchart(width=750, height=1000)


categories = range(1,len(fasta)+1)
L = len(fasta)
labels = [ '<tr>c'+str(c)+'</tr><tr>z'+str(L-c)+'</tr>' for c in categories]

options = {
   'chart':     { 'type': 'bar' },
   'title':     { 'text': 'Intensities of Fragments' },
   'subtitle':  { 'text': 'as estimated by the MassTodon' },
   'xAxis': {
       'categories': labels,
       'title': {'text': 'null'},
       'labels':{
           'useHTML':True
       }
   },
   'yAxis': {
       'min': 0,
       'title': {
           'text': 'Population (millions)',
           'align': 'high'
       },
       'labels': {
           'overflow': 'justify'
       }
   },
   'tooltip': {
       'valueSuffix': ' millions'
   },
   'plotOptions': {
       'bar': {
           'dataLabels': {
               'enabled': True
           }
       }
   },
   'legend': {
       'layout': 'vertical',
       'align': 'right',
       'verticalAlign': 'top',
       'x': -40,
       'y': 80,
       'floating': True,
       'borderWidth': 1
   },
   'credits': {
       'enabled': False
   }
 }

H.set_dict_options(options)
H.add_data_set( data = [ fragments['c'][i-1] for i in categories ],
                type = 'bar',
                name = 'c'     )


categories.reverse()
H.add_data_set( data = [ -fragments['z'][i-1] for i in categories ],
                type = 'bar',
                name = 'z'     )

#    series: [{
#        name: 'Year 1800',
#        data: [107, 31, 635, 203, 2]
#    }, {
#        name: 'Year 1900',
#        data: [133, 156, 947, 408, 6]
#    }, {
#        name: 'Year 2012',
#        data: [1052, 954, 4250, 740, 38]
#    }]
# }

H.save_file('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/fragments')
