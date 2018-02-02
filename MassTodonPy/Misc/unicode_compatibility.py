def is_string(string):
    try:
        return isinstance(string, (str, unicode))
    except NameError:
        return isinstance(string, str)
