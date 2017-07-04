from highcharts import Highchart


def make_precursor_intensity_plot(precursors):
    H = Highchart(width=750, height=600)
    charges = precursors.keys()
    Q_min = min(charges)
    Q_max = max(charges)
    categories = range(Q_min,Q_max+1)
    values = [ int(precursors[q]) for q in categories ]
    options = {
    	'title': {
            'text': 'Intensities Plot'
        },
        'xAxis': {
            'categories': categories,
            'title': {
                'text': 'Charge State'
            }
        },
        'yAxis': {
            'min': 0,
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
                'dataLabels': { 'enabled': True }
            }
        }
    }
    H.set_dict_options(options)
    H.add_data_set(values, 'bar', 'MassTodonPy: Precursors Charge')
    return H

# H = make_precursor_intensity_plot(precursors)
# H.save_file()
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
# H = make_etnod_ptr_probability_plot(res)
# H.save_file()
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


H = Highchart(width=750, height=600)
categories = ['0-4', '5-9', '10-14', '15-19',
  '20-24', '25-29', '30-34', '35-39', '40-44',
  '45-49', '50-54', '55-59', '60-64', '65-69',
  '70-74', '75-79', '80-84', '85-89', '90-94',
  '95-99', '100 + ']

options = {
    'chart': {
      'type': 'bar'
    },
    'title': {
      'text': 'Population pyramid for Germany, 2015'
    },
    'subtitle': {
      'text': 'Source: <a href="http://populationpyramid.net/germany/2015/">Population Pyramids of the World from 1950 to 2100</a>'
    },
    'xAxis': [{
      'categories': categories,
      'reversed': False,
      'labels': {
        'step': 1
      }
    }, {
      'opposite': True,
      'reversed': False,
      'categories': categories,
      'linkedTo': 0,
      'labels': {
        'step': 1
      }
    }],
    'yAxis': {
      'title': {'text': None},
      'labels': {
            'formatter': 'function () {\
                return Math.abs(this.value) + "%";\
            }'
        }
    },
    'plotOptions': {
      'series': { 'stacking': 'normal' }
    },
    'tooltip': {
        'formatter': 'function() {\
        return "<b>" + this.series.name + ", age " + this.point.category + "</b><br/>" +\
          "Population: " + Highcharts.numberFormat(Math.abs(this.point.y), 0);\
      }'
    }
}

H.set_dict_options(options)
H.add_data_set( data = [-2.2, -2.2, -2.3, -2.5, -2.7, -3.1, -3.2, -3.0, -3.2, -4.3, -4.4, -3.6, -3.1, -2.4, -2.5, -2.3, -1.2, -0.6, -0.2, -0.0, -0.0],
                type = 'bar',
                name = 'Male'     )

H.add_data_set( data = [2.1, 2.0, 2.2, 2.4, 2.6, 3.0, 3.1, 2.9, 3.1, 4.1, 4.3, 3.6, 3.4, 2.6, 2.9, 2.9,1.8, 1.2, 0.6, 0.1, 0.0 ],
                type = 'bar',
                name = 'Female'     )

H.save_file('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/age_pyramid')
