from highcharts import Highchart
from collections import Counter, defaultdict

def make_highchart(width, height, options, series):
    H = Highchart(width=width, height=height)
    H.set_dict_options(options)
    for s in series:
        H.add_data_set(**s)
    return H


def prepare_precursor_intensities(precursors_intensitites):
    charges    = precursors_intensitites.keys()
    categories = range(1,max(charges)+1)
    series = [
        {
            'name': 'Intensity',
            'data': [ int(precursors_intensitites[q]) for q in categories ],
            'type': 'bar'
        }
    ]
    options= {
        'title': {
            'text': 'Total Intensity of Precursors'
        },
        'subtitle':{
            'text': 'for different charge states'
        },
        'xAxis': {
            'categories':   categories,
            'title': {
                'text': 'Charge State'
            },
            'labels':{
                'format': 'q = {value}'
            }
        },
        'yAxis': {
            'min':  0,
            'title': {
                'text': 'Intensity',
                'align': 'high'
            },
            'labels': {
                'overflow': 'justify'
            }
        },
        'plotOptions': {
            'bar': {
                'dataLabels':{
                    'enabled': True
                }
            }
        }
    }
    return {'options': options,
            'series' : series   }


def prepare_etnod_ptr_probabilities(algos):
    categories = ('PTR', 'ETnoD')
    series = [  {   'name': cat,
                    'data': [ algos[a][0][cat]*100 for a in algos ],
                    'type': 'bar'
                } for cat in categories ]

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
        }
    }
    return {'options': options,
            'series' : series   }


def prepare_branching_ratios(branching_ratios, branching_ratio, pk):
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
    series = [
        {   'name': 'branching ratios for different number of reactions',
            'data':        branching_ratios,
            'dataLabels':  {'format': '{y:.2f}'}     },

        {   'name': 'branching ratio aggregated over all reactions',
            'data': [(0,branching_ratio), (max(pk)-1,branching_ratio) ],
            'type': 'line',
            'dataLabels':  {'format': '{y:.2f}'}     }
    ]
    return {'options': options,
            'series':  series}


def prepare_fragmentation_probs(fasta, algos):
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
    series = []
    for algo in algos:
        probs, counts = algos[algo]
        series.append(
            {
                'data': [probs[i]*100 for i in xrange(len(fasta))],
                'name': algo.split('_')[0].title(),
                'type': 'column'
            }
        )
    return {'options': options,
            'series' : series  }

#TODO: fix this function. check if fragments are correctly matched
def prepare_fragment_pyramid(fasta, fragments):
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
    series = [
        {
            'name': 'c',
            'data': [ fragments['c'][i-1] for i in categories ],
            'type': 'bar'
        }
    ]
    categories.reverse()
    series.append(
        {
            'name': 'z',
            'data': [ -fragments['z'][i-1] for i in categories ],
            'type': 'bar'
        }
    )
    return {'options': options,
            'series' : series  }



def make_highcharts(fasta, Q, raw_estimates, algos, make_plots = False):
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
                    precursors[ReactionsNo]['PTR']  += M['estimate']*PTRs
                    precursors[ReactionsNo]['ETnoD']+= M['estimate']*ETnoDs
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

    charts = {
        'precursor_intensities' : prepare_precursor_intensities(precursors_intensitites),
        'etnod_ptr_probability' : prepare_etnod_ptr_probabilities(algos),
        'branching_ratios'      : prepare_branching_ratios(branching_ratios, branching_ratio, pk),
        'fragmentation_probs'   : prepare_fragmentation_probs(fasta, algos),
        'fragment_pyramid'      : prepare_fragment_pyramid(fasta, fragments)
    }

    if make_plots:
        charts = map(make_highchart, charts)

    return charts
