from highcharts import Highchart

def make_fragment_pyramid_plot(fasta, fragments):
    H = Highchart(width=750, height=1000)
    fasta_list = [f for f in fasta]
    categories = range(1,len(fasta)+1)
    options = {
        'chart': {'type': 'bar'},
        'title': {'text': 'Intensities of Fragments'},
        'subtitle': {'text': 'as estimated by the MassTodon'},
        'xAxis': [{
          'categories': categories,
          'reversed': False,
          'labels': {
            'step': 1,
            'format': 'z{value}'
          }
        }, {
          'opposite': True,
          'reversed': True,
          'linkedTo': 0,
          'labels': {
              'step': 1,
              'format': 'c{value}'
          }
        }],
        'yAxis': {
          'title': {'text': None},
          'labels': {
                'formatter': 'function () {\
                    return Math.abs(this.value);\
                }'
            }
        },
        'plotOptions': {
          'series': { 'stacking': 'normal' }
        },
        'tooltip': {
            'valueDecimals': 0,
            'headerFormat': '<table>',
            'pointFormat': "<tr><td style='color:{series.color};padding:0'>{series.name}{point.x}: <b>{point.y:.0f}</b> </td></tr>",
            'footerFormat': '</table>',
            'shared':     True,
            'useHTML':  True
        }
        # 'tooltip': {
        #     'formatter': "function () {\
        #         var s = '';\
        #         var d;\
        #         $.each(this.points, function () {\
        #             if (this.x < 0){\
        #                 d = this.x;\
        #                 s += '<br/>' + this.series.name + x +  ': ' +\
        #                     this.y;\
        #             } else {\
        #                 d = -int(this.x);\
        #                 s += '<br/>' + this.series.name + x +  ': ' +\
        #                     this.y;\
        #             }\
        #         });\
        #         return s;\
        #     }",
        #     'useHTML':  True,
        #     'shared':   True
        # }
     }
    H.set_dict_options(options)
    H.add_data_set( data = [ fragments['c'][i-1] for i in categories ],
                    type = 'bar',
                    name = 'c'     )


    categories.reverse()
    H.add_data_set( data = [ -fragments['z'][i-1] for i in categories ],
                    type = 'bar',
                    name = 'z'     )
    return H
H.save_file('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/fragment_pyramid')
