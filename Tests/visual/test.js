var option = {
        "loading": {},
        "subtitle": {
                "text": "as estimated by the MassTodon"
        },
        "xAxis": [{
                "reversed": false,
                "labels": {
                        "step": 1,
                        "format": "z{value}"
                },
                "categories": [76, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
        }, {
                "linkedTo": 0,
                "reversed": true,
                "labels": {
                        "step": 1,
                        "format": "c{value}"
                },
                "opposite": true
        }],
        "title": {
                "text": "Intensities of Fragments"
        },
        "series": {},
        "labels": {},
        "yAxis": {
                "labels": {
                        "formatter": function() {
                                return Math.abs(this.value);
                        }
                },
                "title": {
                        "text": null
                }
        },
        "chart": {
                "width": 750,
                "renderTo": "container",
                "type": "bar",
                "height": 1000
        },
        "tooltip": {
                "footerFormat": "</table>",
                "shared": true,
                "headerFormat": "<table>",
                "pointFormat": "<tr><td style='color:{series.color};padding:0'>{series.name}{point.x}: <b>{point.y:.0f}</b> </td></tr>",
                "useHTML": true,
                "valueDecimals": 0
        },
        "plotOptions": {
                "series": {
                        "stacking": "normal"
                }
        },
        "credits": {
                "enabled": false
        },
        "colors": {},
        "pane": {},
        "exporting": {},
        "drilldown": {},
        "navigation": {},
        "legend": {}
};
