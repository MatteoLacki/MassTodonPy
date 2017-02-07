from pyteomics import mzxml, auxiliary

path = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/'
file = path+'FRL_220715_ubi_952_ETD_40ms_04.mzXML'
file = '/Users/matteo/Documents/MassTodon/ReadingMzxml/data.mzXML'
file = path+'Ubiquitin_ETD_10 ms_1071.mzXML'

with mzxml.read(file) as reader:
    for spectrum in reader:
        print spectrum['intensity array']




mzxml.version_info(file)





with mzxml.read(file) as reader:
    for spectrum in reader:
        mz = spectrum['m/z array']
        intensity = spectrum['intensity array']

len(mz)
len(intensity)

len(intensity[intensity>0])
