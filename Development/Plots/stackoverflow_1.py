from bokeh.charts import Bar, show
import pandas as pd

data_dict = { 'numstudents' : [43, 22, 1,
                               2, 54, 9,
                               18, 10, 5,
                               14, 12, 3,
                               15, 11, 1,
                               14, 8, 2],
              'language' : ['Matlab','Matlab','Matlab',
                            'C/C++','C/C++','C/C++',
                            'Java','Java','Java',
                            'HTML/CSS','HTML/CSS','HTML/CSS',
                            'Python','Python','Python',
                            'Javascript','Javascript','Javascript'],
              'skill_level' : ['Beginner','Intermediate','Expert',
                               'Beginner','Intermediate','Expert',
                               'Beginner','Intermediate','Expert',
                               'Beginner','Intermediate','Expert',
                               'Beginner','Intermediate','Expert',
                               'Beginner','Intermediate','Expert']
            }

data_df = pd.DataFrame(data_dict)

p = Bar(data_df,
        values='numstudents',
        label='language',
        stack='skill_level',
        legend='top_right',
        title="ECEn 360 Student Self-Reported Programming Skills",
        tooltips=[('Students:', '@numstudents'), ('Language:', '@language')]
       )
show(p)
