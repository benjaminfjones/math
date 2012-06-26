r"""
Module: tools
Author: Benjamin Jones <benjaminfjones@gmail.com>

DESCRIPTION:

Various useful decorators. In particular:

 * @decorator: applies update_wrapper to decorator
 * @n_ary: turns a binary function into an n_ary function
 * @trace: prints the trace of a recursive function
 * @disabled: set, e.g. trace = disabled to disable a decorator

"""
from functools import update_wrapper

def decorator(d):
    "Decorator to apply update_wrapper on decorators"
    def _d(fn):
        return update_wrapper(d(fn), fn)
    update_wrapper(_d, d)
    return _d

@decorator2
def n_ary(f):
    r"""
    Takes a binary function `f:X \times X \to X` and returns an n-ary version
    f : X^n \to X.
    """
    def _f(x, *args):
        if len(args) == 1:
            return f(x, args[0])
        if len(args) > 1:
            return f(x, _f(args[0], *args[1:]))
        else:
            raise ValueError("you must supply at least two input values")
    return _f

@n_ary
def f(x,y):
    "Simple binary addition -- used as an example for use with @n_ary decorator"
    return x+y

@decorator
def trace(f):
    """
    Trace decorator for use with recursive functions.
    """
    indent = '   '
    def _f(*args):
        signature = '%s(%s)' % (f.__name__, ', '.join(map(repr, args)))
        print '%s--> %s' % (trace.level*indent, signature)
        trace.level += 1
        try:
            result = f(*args)
            print '%s<-- %s == %s' % ((trace.level-1)*indent, 
                                      signature, result)
        finally:
            trace.level -= 1
        return result
    trace.level = 0
    return _f

@trace
def fib(n):
    "Fibonacci function -- Used as an example for use with @trace decorator"
    if n == 0 or n == 1:
        return 1
    else:
        return fib(n-1) + fib(n-2)

@decorator
def disabled(f): return f

# To disable trace, n_ary, etc.. uncomment the following lines
# trace = disabled
# n_ary = disabled
