def listisize(func):
    def wrapper(args):
        args = set(args)
        results = {}
        for arg in args:
            results[arg] = func(arg)
        return results
    return wrapper
