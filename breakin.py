def breakin(fu):
    """ Decorator that loads the debugger in the beginning of the decorated function"""
    def wrap(*args, **kwargs):
        import pdb; pdb.set_trace()
        result = fu(*args, **kwargs)
        return result
    return wrap
