def parse_path(path):
    '''Parsers path to the file.'''
    file_path, file_ext  = os.path.splitext(path)
    file_name = file_path.split('/')[-1]
    file_path = "/".join(file_path.split('/')[:-1]) + '/'
    return file_path, file_name, file_ext


def read_spectrum(  path    = None,
                    spectrum= None,
                    cut_off = None,
                    opt_P   = None,
                    digits  = 2   ):

    assert path or spectrum, "No path to mass spectrum or no mass spectrum."

    assert not (path and spectrum), "Please decide if you pass the spectrum as the argument (provide spectrum argument for method read_spectrum) or you want MassTodon to read the spectrum from file following the provided path."

    assert cutOff or opt_P, "Please decide if you want to apply an intensity based cut-off or to choose the optimal P-set of experimental peaks, e.g. a representative of the class of the smallest sets of peaks with a probability at least P."

    if path:
        file_path, file_name, file_ext = parse_path(path)
        file_ext = file_ext.lower()
        reader = {  '':         read_txt,
                    '.txt':     read_txt,
                    '.mzxml':   read_mzxml
        }[file_ext]
        spectrum = reader(path, digits)

    result = {}
    result['original spectrum'] = spectrum
    result['original total intensity'] = sum(spectrum[1])

    if cutOff:
        result['trimmed spectrum'] = trim_spectrum( *spectrum, cut_off=cut_off)
    else:
        result['trimmed spectrum'] = quantile_trim( *spectrum, perc = opt_P)

    result['total intensity of trimmed spectrum'] = result['trimmed spectrum'][1].sum()
    return result
