from highcharts import Highchart


def make_precursor_intensity_plot(precursors):
    H = Highchart(width=750, height=600)
    charges = precursors.keys()
    Q_min = 1
    Q_max = max(charges)
    categories = range(Q_min,Q_max+1)
    values = [ int(precursors[q]) for q in categories ]
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


def make_etnod_ptr_probability_plot(res):
    categories = ('PTR', 'ETnoD')
    algos = ('basic_analysis', 'intermediate_analysis', 'advanced_analysis')
    H = Highchart(width=750, height=600)
    for cat in categories:
        H.add_data_set(
            data = [ res[a][0][cat]*100 for a in algos],
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
H = make_etnod_ptr_probability_plot(res)
H.save_file('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/etnod_ptr')
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
    H.add_data_set( data        = [(0,branching_ratio), (max(pk),branching_ratio) ],
                    type        = 'line',
                    dataLabels  = {'format': '{y:.2f}'},
                    name        = 'branching ratio aggregated over all reactions'     )
    return H

# import json
# H = make_branching_ratio_plot(branching_ratios, branching_ratio, pk)
# H.save_file('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/lines')
# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/lines_option.json','w') as h:
#     json.dump(H.option, h)
# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/lines_data.json','w') as h:
#     json.dump(H.data, h)


def make_fragmentation_prob_plot(fasta, res):
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
    algos = ('basic_analysis', 'intermediate_analysis', 'advanced_analysis')
    for algo in algos:
        probs, counts = res[algo]
        H.add_data_set(
            data = [probs[i]*100 for i in xrange(len(fasta))],
            name = algo.split('_')[0].title() ,
            type = 'column'
        )
    return H

# H = make_fragmentation_prob_plot(fasta, res)
# H.save_file('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/fragmentation_probs')
# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/fragmentation_probs_options.json','w') as h:
#     json.dump(H.option, h)
#
# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/fragmentation_probs_data.json','w') as h:
#     json.dump(H.data, h)
