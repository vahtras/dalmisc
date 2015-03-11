import time

def time_me(func):
    def timed_function(*args, **kwargs):
        t1 = time.clock()
        result = func(*args, **kwargs)
        t2 = time.clock()
        print "Time used in %s: %f" % (func.__name__, t2 - t1)
        return result
    return timed_function




if __name__ == "__main__":

    @time_me
    def f(x):
        pass

    f(None)
