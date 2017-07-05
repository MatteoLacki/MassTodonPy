from highcharts import Highchart
from collections import Counter, defaultdict

def make_precursor_intensity_plot(precursors_intensitites):
    H = Highchart(width=750, height=600)
    charges = precursors_intensitites.keys()
    categories = range(1,max(charges)+1)
    values = [ int(precursors_intensitites[q]) for q in categories ]
    options = {
    	'title': {'text': 'Total Intensity of Precursors'},
        'subtitle':{'text': 'for different charge states'},
        'xAxis': {
            'categories': categories,
            'title': { 'text': 'Charge State'},
            'labels':{ 'format': 'q = {value}'}
        },
        'yAxis': {
            'min': 0,
            'title': {
                'text': 'Intensity',
                'align': 'high'
            },
            'labels': {'overflow': 'justify'}
        },
        'plotOptions': {
            'bar': {'dataLabels': { 'enabled': True }}
        }
    }
    H.set_dict_options(options)
    H.add_data_set(data=values,type='bar', name='Intensity')
    return H

# H = make_precursor_intensity_plot(precursors_intensitites)
# H.save_file('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/precursorCharge')
# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/options.json','w') as h:
#     json.dump(H.option, h)
#
# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/data.json','w') as h:
#     json.dump(H.data, h)


def make_etnod_ptr_probability_plot(algos):
    categories = ('PTR', 'ETnoD')
    H = Highchart(width=750, height=600)
    for cat in categories:
        H.add_data_set(
            data = [ algos[a][0][cat]*100 for a in algos],
            name = cat,
            type = 'bar'
        )
    options = {
        'chart': { 'type': 'bar' },
        'title': { 'text': 'Probabilities of PTR and ETnoD'},
        'subtitle':{'text': 'as estimated using three pairing algorithms: Basic, Intermediate, and Advanced'},
        'xAxis': { 'categories': [ a.split('_')[0].title() for a in algos]},
        'yAxis':[ {
            'min': 0,
            'max': 100,
            'title': {'text': 'Probability of ETnoD'},
            'labels':{ 'format': '{value:.0f} %'},
            'opposite': True,
            },
            {
            'min': 0,
            'max': 100,
            'title': {'text': 'Probability of PTR'},
            'labels':{ 'format': '{value:.0f} %'},
            'reversed': True
        }],
        'legend': {
            'reversed': True
        },
        'plotOptions': {
            'series': {
                'stacking': 'normal'
            }
        }}
    H.set_dict_options(options)
    return H

# import json
# H = make_etnod_ptr_probability_plot(algos)
# H.save_file('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/etnod_ptr_plots')
# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/etnodPtr_options.json','w') as h:
#     json.dump(H.option, h)
#
# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/etnodPtr_data.json','w') as h:
#     json.dump(H.data, h)

def make_branching_ratio_plot(branching_ratios, branching_ratio, pk):
    H = Highchart(width=750, height=600)
    options = {
        'chart': {
            'type': 'line'
        },
        'title': {
            'text': 'Branching Ratios'
        },
        'subtitle': {
            'text': 'as calculated for different number of reactions'
        },
        'xAxis': {
            'categories': pk,
            'labels':{
                'format': '{value} reaction(s)'
            },
        },
        'yAxis': {
            'title': {
                'text': 'Branching Ratio = PTR/ETnoD'
            }
        },
        'plotOptions': {
            'line': {
                'dataLabels': {
                    'enabled': True
                },
                'enableMouseTracking': True
            }
        },
        'tooltip': {
            'valueDecimals': 4
        }
        }
    H.set_dict_options(options)
    H.add_data_set( data        = branching_ratios,
                    dataLabels  = {'format': '{y:.2f}'},
                    name        = 'branching ratios for different number of reactions'     )
    H.add_data_set( data        = [(0,branching_ratio), (max(pk)-1,branching_ratio) ],
                    type        = 'line',
                    dataLabels  = {'format': '{y:.2f}'},
                    name        = 'branching ratio aggregated over all reactions'     )
    return H

# import json
# H = make_branching_ratio_plot(branching_ratios, branching_ratio, pk)
# H.save_file('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/branching_ratio')
# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/lines_option.json','w') as h:
#     json.dump(H.option, h)
# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/lines_data.json','w') as h:
#     json.dump(H.data, h)


def make_fragmentation_prob_plot(fasta, algos):
    H = Highchart(width=2000, height=500)
    fasta_list = [f for f in fasta]
    options = {
        'chart': {'type': 'column'},
        'title': {'text': 'Probabilities of Fragmentation'},
        'subtitle': {'text': 'as obtained by pairing ions by three different algorithms: basic, intermediate, and advanced.'},
        'xAxis': {
            'categories': fasta_list,
            'crosshair': True,
            'title':{'text':'Amino Acidic Sequence: N terminus to the left, C terminus to the right'}
        },
        'yAxis': {
            'min': 0,
            'title': {
                'text': 'Probability of Fragmentation'
            },
            'labels':{
                'format': '{value:.1f} %'
            },
        },
        'tooltip': {
            'valueDecimals': 1,
            'headerFormat': '<table>',
            'pointFormat': "<tr><td style='color:{series.color};padding:0'>{series.name}:</td>\
                <td style='padding:0'><b>{point.y:.3f} %</b></td></tr>",
            'footerFormat': '</table>',
            'shared':     True,
            'useHTML':  True
        },
        'plotOptions': {
            'column': {
                'pointPadding': 0.2,
                'borderWidth': 0
            }
        }
     }
    H.set_dict_options(options)
    for algo in algos:
        probs, counts = algos[algo]
        H.add_data_set(
            data = [probs[i]*100 for i in xrange(len(fasta))],
            name = algo.split('_')[0].title() ,
            type = 'column'
        )
    return H

# H = make_fragmentation_prob_plot(fasta, algos)
# H.save_file('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/fragmentation_probs')
# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/fragmentation_probs_options.json','w') as h:
#     json.dump(H.option, h)
#
# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/fragmentation_probs_data.json','w') as h:
#     json.dump(H.data, h)

#TODO: fix this mother-fucking function. I hate these Highcharts....
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
# H.save_file('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/fragment_pyramid')



def make_highcharts(fasta, Q, raw_estimates, algos):
    '''Prepare the outputs of MassTodon for highcharts.'''
    precursors  = defaultdict(Counter)
    precursors_intensitites = Counter()
    fragments   = defaultdict(Counter)
    for r in raw_estimates:
        for M in r['alphas']:
            if M['molType'] == 'precursor':
                ETnoDs  = M['g']
                PTRs    = Q-M['q']-M['g']
                ReactionsNo = ETnoDs + PTRs
                precursors_intensitites[M['q']] += M['estimate']
                if ReactionsNo>0:
                    precursors[ReactionsNo]['PTR'] += M['estimate']*PTRs
                    precursors[ReactionsNo]['ETnoD'] += M['estimate']*ETnoDs
            else:
                frag_type   = M['molType'][0]
                frag_number = int(M['molType'][1:])
                fragments[frag_type][frag_number] += M['estimate']
    pk = precursors.keys()
    pk.sort()
    branching_ratios = []
    for i in pk:
        PTRs = precursors[i]['PTR']
        ETnoDs = precursors[i]['ETnoD']
        try:
            branching_ratio = PTRs/float(ETnoDs)
        except ZeroDivisionError:
            branching_ratio = None
        branching_ratios.append( branching_ratio )
    probs, counts = algos['basic_analysis']
    try:
        branching_ratio = probs['PTR']/probs['ETnoD']
    except ZeroDivisionError:
        branching_ratio = None
    precursor_intensity_plot    = make_precursor_intensity_plot(precursors_intensitites)
    etnod_ptr_probability_plot  = make_etnod_ptr_probability_plot(algos)
    branching_ratio_plot        = make_branching_ratio_plot(branching_ratios, branching_ratio, pk)
    fragmentation_prob_plot     = make_fragmentation_prob_plot(fasta, algos)
    fragment_pyramid_plot       = make_fragment_pyramid_plot(fasta, fragments)

    return precursor_intensity_plot, etnod_ptr_probability_plot, branching_ratio_plot, fragmentation_prob_plot, fragment_pyramid_plot
