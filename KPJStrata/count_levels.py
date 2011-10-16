import sys
from sage.all import *
from bucket_fill import *


def levels(n):
    """
    Returns a generator for all levels associated with quivers of total size `n`.
    """
    for (r,I) in quivers(n):
        for l in kpj_levels_gen(r,I):
            yield l


@parallel() # auto-detect # cpu's
def count_levels_byquiver(q):
    """
    INPUT:
        - q - a quiver (r,I)
    OUTPUT:
        - the number of levels associated to (r,I)

    EXAMPLES:
        sage: %timeit list(count_levels_byquiver([[[1, 1, 1, 1], [0, 1, 3]]]))
        5 loops, best of 3: 348 ms per loop

    """
    r,I = (q[0], q[1])
    s = Integer(0)
    for l in kpj_levels_gen(r,I):
        s += Integer(1)
    return s

def count_all_levels(n, pretty=False):
    count = dict()
    todo = []
    for (r,I) in quivers(n):
        str_r_I = "r = %s, I = %s" % (r,I)
        print '{0} ... '.format(str_r_I),
        if quiver_is_trivial(r,I):
            count[str_r_I] = 'TRIVIAL'
            print 'TRIVIAL'
        elif quiver_is_cotrivial(r,I):
            count[str_r_I] = 'CO-TRIVIAL'
            print 'CO-TRIVIAL'
        else:
            print 'todo'
            count[str_r_I] = Integer(0)
            todo.append([r,I])
    results = []
    # run kpj_levels in parallel for quivers in `todo`
    for X in count_levels_byquiver(todo):
        r,I = (X[0][0][0][0], X[0][0][0][1]) # resutls are unpacked
        str_r_I = "r = %s, I = %s" % (r,I) 
        count[str_r_I] = X[1]
    # print nice table of results
    if pretty:
        width=max(map(len,count.keys()))
        C = Integer(0)
        for s in count:
            print '{0:<{width}} #levels = {1}'.format(s,count[s],width=width)
            if type(count[s]) != str:
                C += count[s]
        print 'Count = {0}'.format(C)

##
## -- Parallel Level Enumeration --
##

# The next two functions count levels for each individual quiver
# in parallel. This is done by a special version of `kpj_levels_gen`
# which does the first (and largest) recursion level in parallel.
# The result is a much better distribution of multiple processes over
# the life of the computation.

@parallel()
def count_kpj_levels_initial(*args):
    """
    INPUT:
        - args should be (r, I, l):
        - r = list of integers
        - I = list of indicies in 0 ... len(r)-1
        - l = a partial `level` which is a list of len(r) bipartitions

    Counts levels for the quiver (r,I) whose initial segment is l.
    """
    (r,I,l) = args
    n = sum(r)
    t = len(r)
    count = Integer(0)
    for w in compatible_bipartitions(l[t-1], n, not (t-1 in I)):
        count += Integer(1)
    return count



def count_kpj_levels_parallel(r, I):
    """
    EXAMPLES:

        ! WAY slower than count_levels_byquiver

        sage: %timeit count_kpj_levels_parallel([1, 1, 1, 1], [0, 1, 3])
        5 loops, best of 3: 1.46 s per loop

    """
    t = len(r)
    if t == 0:
        return 0
    else:
        count = Integer(0)
        levels = kpj_levels_gen( r[:t-1], filter(lambda i:i < t-1, I))

        n = sum(r)
        def data_generator(r,I,levels):
            for l in levels:
                yield (r,I,l)

        data = data_generator(r,I,levels)

        for X in count_kpj_levels_initial(data):
            count += X[1]
        return count

# Notes:
# This counting function works well up to n=6
#

def count_all_levels2(n, pretty=False):
    count = dict()
    todo = []
    for (r,I) in quivers(n):
        str_r_I = "r = %s, I = %s" % (r,I)
        print '{0} ... '.format(str_r_I),
        if quiver_is_trivial(r,I):
            count[str_r_I] = 'TRIVIAL'
            print 'TRIVIAL'
        elif quiver_is_cotrivial(r,I):
            count[str_r_I] = 'CO-TRIVIAL'
            print 'CO-TRIVIAL'
        else:
            print 'todo'
            count[str_r_I] = Integer(0)
            todo.append([r,I])
    results = []
    # run special `kpj_levels_parallel` for quivers in `todo`
    print '\n*** Counting Non-trivial Levels ***\n'
    number_todo = len(todo)
    i = Integer(0)
    for (r,I) in todo:
        i += Integer(1)
        str_r_I = "r = %s, I = %s" % (r,I)
        print '{0} ({1}/{2}) ... '.format(str_r_I, i, number_todo),
        sys.stdout.flush()
        count[str_r_I] = count_kpj_levels_parallel(r,I)  # potentially a long time
        print '{0}'.format(count[str_r_I])
    # print nice table of results
    print '\n*** TABLE: ***\n'
    width=max(map(len,count.keys()))
    C = Integer(0)
    for s in count:
        print '{0:<{width}} #levels = {1}'.format(s,count[s],width=width)
        if type(count[s]) != str:
            C += count[s]
    print 'Count = {0}'.format(C)


#
# - main module code -

if __name__ == '__main__':
    n = int(raw_input('Count levels n = '))
    count_all_levels2(n, pretty=True)
