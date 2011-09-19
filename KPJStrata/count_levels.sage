import sys         # for sys.stdout.flush()


def levels(n):
    """
    Returns a generator for all levels associated with quivers of total size `n`.
    """
    for (r,I) in quivers(n):
        for l in kpj_levels_gen(r,I):
            yield l



# need generator version of kpj_levels for n >= 7 
# def count_levels(n, pretty=False):
#     count = dict()
#     for (r,I) in quivers(n):
#         str_r_I = "r = %s, I = %s" % (r,I) 
#         print '{0} ... '.format(str_r_I),
#         if quiver_is_trivial(r,I):
#             count[str_r_I] = 'TRIVIAL'
#             print 'TRIVIAL'
#         elif quiver_is_cotrivial(r,I):
#             count[str_r_I] = 'CO-TRIVIAL'
#             print 'CO-TRIVIAL'
#         else:
#             sys.stdout.flush()
#             s = 0
#             for l in kpj_levels_gen(r,I):
#                 s += 1
#             count[str_r_I] = s
#             print s
#     if pretty:
#         width=max(map(len,count.keys()))
#         C = 0
#         for s in count:
#             print '{0:<{width}} #levels = {1}'.format(s,count[s],width=width)
#             if type(count[s]) != str:
#                 C += count[s]
#         print 'Count = {0}'.format(C)

#
# next two functions count levels in parallel by assigning the count
# for each quiver to a separate process. This typically results in 
# lots of processes at first, dying down to a few long running 
# single processes that take the majority of the cpu time.
#

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
    s = 0
    for l in kpj_levels_gen(r,I):
        s += 1
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
            count[str_r_I] = 0
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
        C = 0
        for s in count:
            print '{0:<{width}} #levels = {1}'.format(s,count[s],width=width)
            if type(count[s]) != str:
                C += count[s]
        print 'Count = {0}'.format(C)

# The next two functions count levels for each individual quiver
# in parallel. This is done by a special version of `kpj_levels_gen`
# which does the first (and largest) recursion level in parallel.
# The result is a much better distribution of multiple processes over
# the life of the computation.

@parallel(2)
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
    count = 0
    for w in compatible_bipartitions(l[t-1], n, not (t-1 in I)):
        count += 1
    return count

@parallel()
def parallel_test(*args, **kwds):
    return 'args = {0}'.format(args)

def data_generator(r,I,levels):
    for l in levels:
        yield (r,I,l)

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
        count = 0
        levels = kpj_levels_gen( r[:t-1], filter(lambda i:i < t-1, I))

        n = sum(r)
        data = data_generator(r,I,levels)

        for X in count_kpj_levels_initial(data):
            count += X[1]
        return count

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
            count[str_r_I] = 0
            todo.append([r,I])
    results = []
    # run special `kpj_levels_parallel` for quivers in `todo`
    print '\n*** Counting Non-trivial Levels ***\n'
    number_todo = len(todo)
    i = 0
    for (r,I) in todo:
        i += 1
        str_r_I = "r = %s, I = %s" % (r,I)
        print '{0} ({1}/{2}) ... '.format(str_r_I, i, number_todo),
        sys.stdout.flush()
        count[str_r_I] = count_kpj_levels_parallel(r,I)
        print '{0}'.format(count[str_r_I])
    # print nice table of results
    print '\n*** TABLE: ***\n'
    width=max(map(len,count.keys()))
    C = 0
    for s in count:
        print '{0:<{width}} #levels = {1}'.format(s,count[s],width=width)
        if type(count[s]) != str:
            C += count[s]
    print 'Count = {0}'.format(C)

# output_N_7 = """
# sage: count_levels(7, pretty=True)
# r = [1, 1, 1, 1, 1, 1, 1], I = [] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [] ...  TRIVIAL
# r = [1, 1, 1, 2, 2], I = [] ...  TRIVIAL
# r = [1, 1, 1, 1, 3], I = [] ...  TRIVIAL
# r = [1, 2, 2, 2], I = [] ...  TRIVIAL
# r = [1, 1, 2, 3], I = [] ...  TRIVIAL
# r = [1, 1, 1, 4], I = [] ...  TRIVIAL
# r = [2, 2, 3], I = [] ...  TRIVIAL
# r = [1, 3, 3], I = [] ...  TRIVIAL
# r = [1, 2, 4], I = [] ...  TRIVIAL
# r = [1, 1, 5], I = [] ...  TRIVIAL
# r = [3, 4], I = [] ...  TRIVIAL
# r = [2, 5], I = [] ...  TRIVIAL
# r = [1, 6], I = [] ...  TRIVIAL
# r = [7], I = [] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 1, 1], I = [0] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0] ...  TRIVIAL
# r = [1, 1, 1, 2, 2], I = [0] ...  TRIVIAL
# r = [1, 1, 1, 1, 3], I = [0] ...  TRIVIAL
# r = [1, 2, 2, 2], I = [0] ...  TRIVIAL
# r = [1, 1, 2, 3], I = [0] ...  TRIVIAL
# r = [1, 1, 1, 4], I = [0] ...  TRIVIAL
# r = [1, 3, 3], I = [0] ...  TRIVIAL
# r = [1, 2, 4], I = [0] ...  TRIVIAL
# r = [1, 1, 5], I = [0] ...  TRIVIAL
# r = [1, 6], I = [0] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 1, 1], I = [0, 1] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0, 1] ...  TRIVIAL
# r = [1, 1, 1, 2, 2], I = [0, 1] ...  TRIVIAL
# r = [1, 1, 1, 1, 3], I = [0, 1] ...  TRIVIAL
# r = [1, 1, 2, 3], I = [0, 1] ...  TRIVIAL
# r = [1, 1, 1, 4], I = [0, 1] ...  TRIVIAL
# r = [1, 1, 5], I = [0, 1] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [5] ...  CO-TRIVIAL
# r = [1, 1, 1, 2, 2], I = [3] ...  59930
# r = [1, 2, 2, 2], I = [1] ...  8108
# r = [1, 1, 2, 3], I = [2] ...  2264
# r = [2, 2, 3], I = [0] ...  TRIVIAL
# r = [1, 2, 4], I = [1] ...  152
# r = [2, 5], I = [0] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 1, 1], I = [0, 1, 2] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0, 1, 2] ...  TRIVIAL
# r = [1, 1, 1, 2, 2], I = [0, 1, 2] ...  TRIVIAL
# r = [1, 1, 1, 1, 3], I = [0, 1, 2] ...  TRIVIAL
# r = [1, 1, 1, 4], I = [0, 1, 2] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0, 5] ...  3053812
# r = [1, 1, 1, 2, 2], I = [0, 3] ...  117732
# r = [1, 2, 2, 2], I = [0, 1] ...  TRIVIAL
# r = [1, 1, 2, 3], I = [0, 2] ...  4324
# r = [1, 2, 4], I = [0, 1] ...  TRIVIAL
# r = [1, 1, 1, 1, 3], I = [4] ...  CO-TRIVIAL
# r = [1, 1, 2, 3], I = [3] ...  CO-TRIVIAL
# r = [2, 2, 3], I = [2] ...  CO-TRIVIAL
# r = [1, 3, 3], I = [1] ...  284
# r = [3, 4], I = [0] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 1, 1], I = [0, 1, 2, 3] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0, 1, 2, 3] ...  TRIVIAL
# r = [1, 1, 1, 1, 3], I = [0, 1, 2, 3] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0, 1, 5] ...  4208248
# r = [1, 1, 1, 2, 2], I = [0, 1, 3] ...  154424
# r = [1, 1, 2, 3], I = [0, 1, 2] ...  TRIVIAL
# r = [1, 1, 1, 2, 2], I = [3, 4] ...  CO-TRIVIAL
# r = [1, 2, 2, 2], I = [1, 2] ...  14088
# r = [2, 2, 3], I = [0, 1] ...  TRIVIAL
# r = [1, 1, 1, 1, 3], I = [0, 4] ...  64474
# r = [1, 1, 2, 3], I = [0, 3] ...  4474
# r = [1, 3, 3], I = [0, 1] ...  TRIVIAL
# r = [1, 1, 1, 4], I = [3] ...  CO-TRIVIAL
# r = [1, 2, 4], I = [2] ...  CO-TRIVIAL
# r = [3, 4], I = [1] ...  CO-TRIVIAL
# r = [1, 1, 1, 1, 1, 1, 1], I = [0, 1, 2, 3, 4] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0, 1, 2, 3, 4] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0, 1, 2, 5] ...  6130472
# r = [1, 1, 1, 2, 2], I = [0, 1, 2, 3] ...  TRIVIAL
# r = [1, 1, 1, 2, 2], I = [0, 3, 4] ...  207308
# r = [1, 2, 2, 2], I = [0, 1, 2] ...  TRIVIAL
# r = [1, 1, 1, 1, 3], I = [0, 1, 4] ...  87714
# r = [1, 1, 2, 3], I = [0, 1, 3] ...  6014
# r = [1, 1, 2, 3], I = [2, 3] ...  CO-TRIVIAL
# r = [2, 2, 3], I = [0, 2] ...  1282
# r = [1, 1, 1, 4], I = [0, 3] ...  2200
# r = [1, 2, 4], I = [0, 2] ...  296
# r = [1, 1, 5], I = [2] ...  CO-TRIVIAL
# r = [2, 5], I = [1] ...  CO-TRIVIAL
# r = [1, 1, 1, 1, 1, 1, 1], I = [0, 1, 2, 3, 4, 5] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0, 1, 2, 3, 5] ... 8054712
# r = [1, 1, 1, 2, 2], I = [0, 1, 3, 4] ...  253444
# r = [1, 1, 1, 1, 3], I = [0, 1, 2, 4] ...  122772
# r = [1, 2, 2, 2], I = [1, 2, 3] ...  CO-TRIVIAL
# r = [1, 1, 2, 3], I = [0, 2, 3] ...  7462
# r = [1, 1, 1, 4], I = [0, 1, 3] ...  2888
# r = [1, 3, 3], I = [1, 2] ...  CO-TRIVIAL
# r = [1, 2, 4], I = [1, 2] ...  CO-TRIVIAL
# r = [1, 1, 5], I = [0, 2] ...  130
# r = [1, 6], I = [1] ...  CO-TRIVIAL
# r = [1, 1, 1, 1, 1, 1, 1], I = [0, 1, 2, 3, 4, 5, 6] ...  TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0, 1, 2, 3, 4, 5] ...  TRIVIAL
# r = [1, 1, 1, 2, 2], I = [0, 1, 2, 3, 4] ...  TRIVIAL
# r = [1, 1, 1, 1, 3], I = [0, 1, 2, 3, 4] ...  TRIVIAL
# r = [1, 2, 2, 2], I = [0, 1, 2, 3] ...  TRIVIAL
# r = [1, 1, 2, 3], I = [0, 1, 2, 3] ...  TRIVIAL
# r = [1, 1, 1, 4], I = [0, 1, 2, 3] ...  TRIVIAL
# r = [2, 2, 3], I = [0, 1, 2] ...  TRIVIAL
# r = [1, 3, 3], I = [0, 1, 2] ...  TRIVIAL
# r = [1, 2, 4], I = [0, 1, 2] ...  TRIVIAL
# r = [1, 1, 5], I = [0, 1, 2] ...  TRIVIAL
# r = [3, 4], I = [0, 1] ...  TRIVIAL
# r = [2, 5], I = [0, 1] ...  TRIVIAL
# r = [1, 6], I = [0, 1] ...  TRIVIAL
# r = [7], I = [0] ...  TRIVIAL
# r = [1, 1, 1, 4], I = [0, 1, 2]                      #levels = TRIVIAL
# r = [1, 2, 2, 2], I = [0, 1]                         #levels = TRIVIAL
# r = [3, 4], I = [0]                                  #levels = TRIVIAL
# r = [1, 1, 1, 4], I = []                             #levels = TRIVIAL
# r = [1, 3, 3], I = []                                #levels = TRIVIAL
# r = [7], I = [0]                                     #levels = TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0, 1, 2, 3, 4]          #levels = TRIVIAL
# r = [1, 2, 4], I = [0, 2]                            #levels = 296
# r = [1, 1, 1, 1, 1, 2], I = [5]                      #levels = CO-TRIVIAL
# r = [1, 1, 2, 3], I = []                             #levels = TRIVIAL
# r = [1, 1, 1, 1, 1, 1, 1], I = [0, 1, 2]             #levels = TRIVIAL
# r = [1, 2, 2, 2], I = [1]                            #levels = 8108
# r = [1, 2, 4], I = []                                #levels = TRIVIAL
# r = [1, 1, 2, 3], I = [0]                            #levels = TRIVIAL
# r = [1, 2, 2, 2], I = [0, 1, 2, 3]                   #levels = TRIVIAL
# r = [1, 1, 1, 4], I = [0, 3]                         #levels = 2200
# r = [3, 4], I = [1]                                  #levels = CO-TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0, 1, 2, 3, 5]          #levels = 8054712
# r = [2, 2, 3], I = [2]                               #levels = CO-TRIVIAL
# r = [1, 1, 1, 1, 1, 1, 1], I = [0, 1, 2, 3, 4]       #levels = TRIVIAL
# r = [1, 1, 5], I = [2]                               #levels = CO-TRIVIAL
# r = [1, 3, 3], I = [0, 1]                            #levels = TRIVIAL
# r = [1, 1, 1, 1, 1, 1, 1], I = []                    #levels = TRIVIAL
# r = [2, 2, 3], I = []                                #levels = TRIVIAL
# r = [1, 1, 1, 2, 2], I = [0, 3, 4]                   #levels = 207308
# r = [1, 2, 2, 2], I = [1, 2]                         #levels = 14088
# r = [1, 2, 4], I = [0]                               #levels = TRIVIAL
# r = [1, 1, 1, 4], I = [0]                            #levels = TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0, 1]                   #levels = TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0, 1, 2]                #levels = TRIVIAL
# r = [2, 2, 3], I = [0, 2]                            #levels = 1282
# r = [2, 5], I = [0]                                  #levels = TRIVIAL
# r = [1, 1, 1, 1, 3], I = [4]                         #levels = CO-TRIVIAL
# r = [1, 1, 2, 3], I = [0, 1, 2, 3]                   #levels = TRIVIAL
# r = [3, 4], I = []                                   #levels = TRIVIAL
# r = [1, 2, 2, 2], I = [0, 1, 2]                      #levels = TRIVIAL
# r = [1, 1, 1, 2, 2], I = [0, 1, 2]                   #levels = TRIVIAL
# r = [1, 1, 1, 2, 2], I = []                          #levels = TRIVIAL
# r = [1, 1, 1, 4], I = [0, 1, 3]                      #levels = 2888
# r = [1, 1, 1, 1, 1, 1, 1], I = [0]                   #levels = TRIVIAL
# r = [2, 5], I = [1]                                  #levels = CO-TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0, 1, 5]                #levels = 4208248
# r = [1, 6], I = []                                   #levels = TRIVIAL
# r = [1, 2, 4], I = [0, 1]                            #levels = TRIVIAL
# r = [1, 1, 1, 2, 2], I = [0, 1, 2, 3, 4]             #levels = TRIVIAL
# r = [1, 6], I = [0, 1]                               #levels = TRIVIAL
# r = [1, 1, 1, 2, 2], I = [3, 4]                      #levels = CO-TRIVIAL
# r = [1, 1, 1, 1, 3], I = [0, 1, 4]                   #levels = 87714
# r = [1, 2, 4], I = [2]                               #levels = CO-TRIVIAL
# r = [3, 4], I = [0, 1]                               #levels = TRIVIAL
# r = [1, 1, 1, 2, 2], I = [0, 3]                      #levels = 117732
# r = [1, 1, 1, 1, 1, 1, 1], I = [0, 1, 2, 3, 4, 5, 6] #levels = TRIVIAL
# r = [1, 1, 2, 3], I = [0, 2, 3]                      #levels = 7462
# r = [1, 1, 1, 2, 2], I = [0, 1, 2, 3]                #levels = TRIVIAL
# r = [1, 1, 1, 4], I = [0, 1, 2, 3]                   #levels = TRIVIAL
# r = [1, 1, 1, 1, 3], I = [0, 1, 2, 4]                #levels = 122772
# r = [1, 2, 4], I = [0, 1, 2]                         #levels = TRIVIAL
# r = [1, 1, 2, 3], I = [0, 3]                         #levels = 4474
# r = [1, 2, 2, 2], I = []                             #levels = TRIVIAL
# r = [2, 5], I = [0, 1]                               #levels = TRIVIAL
# r = [1, 1, 1, 4], I = [3]                            #levels = CO-TRIVIAL
# r = [1, 3, 3], I = [0]                               #levels = TRIVIAL
# r = [1, 1, 2, 3], I = [3]                            #levels = CO-TRIVIAL
# r = [1, 1, 2, 3], I = [0, 1, 3]                      #levels = 6014
# r = [2, 2, 3], I = [0, 1, 2]                         #levels = TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = []                       #levels = TRIVIAL
# r = [1, 1, 1, 4], I = [0, 1]                         #levels = TRIVIAL
# r = [2, 2, 3], I = [0, 1]                            #levels = TRIVIAL
# r = [2, 5], I = []                                   #levels = TRIVIAL
# r = [1, 1, 1, 2, 2], I = [0, 1]                      #levels = TRIVIAL
# r = [1, 3, 3], I = [0, 1, 2]                         #levels = TRIVIAL
# r = [1, 1, 5], I = []                                #levels = TRIVIAL
# r = [1, 1, 1, 1, 3], I = [0, 4]                      #levels = 64474
# r = [1, 1, 1, 1, 1, 2], I = [0, 5]                   #levels = 3053812
# r = [1, 1, 5], I = [0, 1, 2]                         #levels = TRIVIAL
# r = [1, 1, 1, 1, 3], I = [0, 1, 2, 3]                #levels = TRIVIAL
# r = [1, 3, 3], I = [1]                               #levels = 284
# r = [1, 1, 1, 2, 2], I = [0, 1, 3]                   #levels = 154424
# r = [1, 1, 2, 3], I = [0, 1, 2]                      #levels = TRIVIAL
# r = [1, 1, 1, 1, 3], I = [0, 1, 2]                   #levels = TRIVIAL
# r = [1, 6], I = [0]                                  #levels = TRIVIAL
# r = [1, 1, 2, 3], I = [2, 3]                         #levels = CO-TRIVIAL
# r = [1, 1, 1, 1, 1, 1, 1], I = [0, 1, 2, 3, 4, 5]    #levels = TRIVIAL
# r = [1, 1, 2, 3], I = [0, 1]                         #levels = TRIVIAL
# r = [1, 1, 1, 1, 3], I = [0, 1]                      #levels = TRIVIAL
# r = [1, 1, 1, 1, 3], I = [0, 1, 2, 3, 4]             #levels = TRIVIAL
# r = [1, 1, 1, 1, 1, 1, 1], I = [0, 1, 2, 3]          #levels = TRIVIAL
# r = [1, 1, 1, 1, 1, 1, 1], I = [0, 1]                #levels = TRIVIAL
# r = [1, 1, 1, 1, 3], I = [0]                         #levels = TRIVIAL
# r = [1, 2, 2, 2], I = [1, 2, 3]                      #levels = CO-TRIVIAL
# r = [1, 1, 1, 2, 2], I = [3]                         #levels = 59930
# r = [1, 1, 5], I = [0, 1]                            #levels = TRIVIAL
# r = [1, 2, 4], I = [1]                               #levels = 152
# r = [1, 1, 1, 1, 1, 2], I = [0, 1, 2, 3, 4, 5]       #levels = TRIVIAL
# r = [1, 1, 5], I = [0]                               #levels = TRIVIAL
# r = [1, 2, 2, 2], I = [0]                            #levels = TRIVIAL
# r = [1, 6], I = [1]                                  #levels = CO-TRIVIAL
# r = [1, 2, 4], I = [1, 2]                            #levels = CO-TRIVIAL
# r = [1, 1, 2, 3], I = [0, 2]                         #levels = 4324
# r = [1, 1, 1, 1, 1, 2], I = [0, 1, 2, 3]             #levels = TRIVIAL
# r = [7], I = []                                      #levels = TRIVIAL
# r = [1, 1, 1, 2, 2], I = [0, 1, 3, 4]                #levels = 253444
# r = [2, 2, 3], I = [0]                               #levels = TRIVIAL
# r = [1, 1, 5], I = [0, 2]                            #levels = 130
# r = [1, 1, 1, 2, 2], I = [0]                         #levels = TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0]                      #levels = TRIVIAL
# r = [1, 3, 3], I = [1, 2]                            #levels = CO-TRIVIAL
# r = [1, 1, 1, 1, 1, 2], I = [0, 1, 2, 5]             #levels = 6130472
# r = [1, 1, 1, 1, 3], I = []                          #levels = TRIVIAL
# r = [1, 1, 2, 3], I = [2]                            #levels = 2264
# Count = 22569008
# """
# output_N_6 = """
# r = [1, 1, 1, 3], I = [0, 1, 3]                #levels = 2888
# r = [3, 3], I = [0, 1]                         #levels = TRIVIAL
# r = [1, 5], I = [1]                            #levels = CO-TRIVIAL
# r = [1, 1, 1, 1, 1, 1], I = [0, 1, 2, 3, 4]    #levels = TRIVIAL
# r = [1, 5], I = [0, 1]                         #levels = TRIVIAL
# r = [1, 1, 2, 2], I = [0, 2, 3]                #levels = 7094
# r = [1, 1, 4], I = [0, 1, 2]                   #levels = TRIVIAL
# r = [1, 1, 1, 3], I = [3]                      #levels = CO-TRIVIAL
# r = [2, 4], I = [1]                            #levels = CO-TRIVIAL
# r = [1, 1, 1, 3], I = []                       #levels = TRIVIAL
# r = [1, 1, 2, 2], I = [0, 2]                   #levels = 4140
# r = [1, 1, 1, 3], I = [0, 3]                   #levels = 2200
# r = [1, 2, 3], I = [1, 2]                      #levels = CO-TRIVIAL
# r = [6], I = [0]                               #levels = TRIVIAL
# r = [1, 1, 1, 3], I = [0, 1]                   #levels = TRIVIAL
# r = [2, 4], I = []                             #levels = TRIVIAL
# r = [1, 1, 2, 2], I = [0]                      #levels = TRIVIAL
# r = [1, 1, 1, 1, 1, 1], I = [0, 1, 2, 3]       #levels = TRIVIAL
# r = [6], I = []                                #levels = TRIVIAL
# r = [1, 1, 2, 2], I = [2]                      #levels = 2166
# r = [1, 1, 2, 2], I = [0, 1, 2, 3]             #levels = TRIVIAL
# r = [1, 1, 1, 1, 1, 1], I = [0, 1, 2, 3, 4, 5] #levels = TRIVIAL
# r = [2, 4], I = [0, 1]                         #levels = TRIVIAL
# r = [1, 5], I = [0]                            #levels = TRIVIAL
# r = [1, 2, 3], I = [0]                         #levels = TRIVIAL
# r = [1, 2, 3], I = [2]                         #levels = CO-TRIVIAL
# r = [1, 2, 3], I = [0, 2]                      #levels = 296
# r = [1, 1, 1, 3], I = [0, 1, 2]                #levels = TRIVIAL
# r = [1, 1, 4], I = [0, 2]                      #levels = 130
# r = [1, 1, 4], I = [0, 1]                      #levels = TRIVIAL
# r = [1, 1, 1, 1, 2], I = [0, 1, 2, 3]          #levels = TRIVIAL
# r = [1, 1, 1, 1, 2], I = [0, 1, 2, 3, 4]       #levels = TRIVIAL
# r = [1, 1, 1, 1, 2], I = [0, 1]                #levels = TRIVIAL
# r = [1, 1, 1, 1, 2], I = []                    #levels = TRIVIAL
# r = [3, 3], I = []                             #levels = TRIVIAL
# r = [1, 1, 1, 1, 1, 1], I = [0, 1, 2]          #levels = TRIVIAL
# r = [1, 2, 3], I = []                          #levels = TRIVIAL
# r = [1, 1, 1, 1, 1, 1], I = []                 #levels = TRIVIAL
# r = [1, 1, 4], I = [2]                         #levels = CO-TRIVIAL
# r = [2, 4], I = [0]                            #levels = TRIVIAL
# r = [1, 1, 4], I = [0]                         #levels = TRIVIAL
# r = [1, 1, 1, 1, 2], I = [0, 1, 2]             #levels = TRIVIAL
# r = [1, 1, 2, 2], I = []                       #levels = TRIVIAL
# r = [2, 2, 2], I = [0, 1, 2]                   #levels = TRIVIAL
# r = [1, 1, 1, 1, 1, 1], I = [0]                #levels = TRIVIAL
# r = [1, 2, 3], I = [0, 1, 2]                   #levels = TRIVIAL
# r = [1, 1, 1, 1, 2], I = [0, 1, 4]             #levels = 84170
# r = [1, 1, 1, 1, 1, 1], I = [0, 1]             #levels = TRIVIAL
# r = [1, 1, 2, 2], I = [0, 1]                   #levels = TRIVIAL
# r = [2, 2, 2], I = [0]                         #levels = TRIVIAL
# r = [1, 1, 1, 3], I = [0]                      #levels = TRIVIAL
# r = [1, 1, 1, 3], I = [0, 1, 2, 3]             #levels = TRIVIAL
# r = [1, 2, 3], I = [0, 1]                      #levels = TRIVIAL
# r = [3, 3], I = [0]                            #levels = TRIVIAL
# r = [1, 1, 1, 1, 2], I = [4]                   #levels = CO-TRIVIAL
# r = [1, 2, 3], I = [1]                         #levels = 152
# r = [2, 2, 2], I = [0, 1]                      #levels = TRIVIAL
# r = [1, 1, 4], I = []                          #levels = TRIVIAL
# r = [1, 1, 1, 1, 2], I = [0, 1, 2, 4]          #levels = 117556
# r = [1, 1, 1, 1, 2], I = [0, 4]                #levels = 61922
# r = [2, 2, 2], I = []                          #levels = TRIVIAL
# r = [1, 1, 2, 2], I = [0, 1, 2]                #levels = TRIVIAL
# r = [1, 1, 1, 1, 2], I = [0]                   #levels = TRIVIAL
# r = [1, 1, 2, 2], I = [2, 3]                   #levels = CO-TRIVIAL
# r = [1, 5], I = []                             #levels = TRIVIAL
# 
# Count = 282714
# """
# 
# output_N_5 = """
# r = [1, 2, 2], I = []                    #levels = TRIVIAL
# r = [1, 1, 1, 1, 1], I = [0, 1]          #levels = TRIVIAL
# r = [1, 2, 2], I = [0, 1, 2]             #levels = TRIVIAL
# r = [2, 3], I = [0]                      #levels = TRIVIAL
# r = [1, 4], I = [0]                      #levels = TRIVIAL
# r = [2, 3], I = []                       #levels = TRIVIAL
# r = [1, 1, 1, 1, 1], I = []              #levels = TRIVIAL
# r = [1, 2, 2], I = [1]                   #levels = 146
# r = [2, 3], I = [0, 1]                   #levels = TRIVIAL
# r = [5], I = [0]                         #levels = TRIVIAL
# r = [1, 1, 1, 2], I = [0, 1, 2, 3]       #levels = TRIVIAL
# r = [1, 1, 1, 2], I = [0, 1, 2]          #levels = TRIVIAL
# r = [1, 1, 1, 2], I = [0]                #levels = TRIVIAL
# r = [1, 2, 2], I = [1, 2]                #levels = CO-TRIVIAL
# r = [1, 1, 3], I = [0, 2]                #levels = 130
# r = [1, 1, 1, 1, 1], I = [0]             #levels = TRIVIAL
# r = [1, 1, 1, 2], I = [0, 3]             #levels = 2124
# r = [1, 4], I = []                       #levels = TRIVIAL
# r = [1, 1, 1, 2], I = [0, 1]             #levels = TRIVIAL
# r = [2, 3], I = [1]                      #levels = CO-TRIVIAL
# r = [1, 4], I = [0, 1]                   #levels = TRIVIAL
# r = [1, 1, 3], I = []                    #levels = TRIVIAL
# r = [1, 2, 2], I = [0, 1]                #levels = TRIVIAL
# r = [1, 4], I = [1]                      #levels = CO-TRIVIAL
# r = [1, 1, 3], I = [0, 1, 2]             #levels = TRIVIAL
# r = [1, 1, 3], I = [0]                   #levels = TRIVIAL
# r = [1, 1, 1, 2], I = []                 #levels = TRIVIAL
# r = [1, 1, 3], I = [2]                   #levels = CO-TRIVIAL
# r = [1, 1, 1, 1, 1], I = [0, 1, 2]       #levels = TRIVIAL
# r = [1, 1, 3], I = [0, 1]                #levels = TRIVIAL
# r = [5], I = []                          #levels = TRIVIAL
# r = [1, 1, 1, 2], I = [3]                #levels = CO-TRIVIAL
# r = [1, 1, 1, 1, 1], I = [0, 1, 2, 3]    #levels = TRIVIAL
# r = [1, 1, 1, 2], I = [0, 1, 3]          #levels = 2784
# r = [1, 2, 2], I = [0]                   #levels = TRIVIAL
# r = [1, 1, 1, 1, 1], I = [0, 1, 2, 3, 4] #levels = TRIVIAL
# 
# Count = 5184
# """
#
# output_N_4 = """
# r = [1, 1, 2], I = [2]             #levels = CO-TRIVIAL
# r = [1, 3], I = [0]                #levels = TRIVIAL
# r = [1, 1, 1, 1], I = [0, 1, 2, 3] #levels = TRIVIAL
# r = [2, 2], I = [0, 1]             #levels = TRIVIAL
# r = [1, 1, 2], I = [0, 1, 2]       #levels = TRIVIAL
# r = [1, 3], I = [1]                #levels = CO-TRIVIAL
# r = [1, 1, 1, 1], I = [0, 1, 2]    #levels = TRIVIAL
# r = [1, 1, 2], I = [0, 2]          #levels = 130
# r = [1, 1, 1, 1], I = [0]          #levels = TRIVIAL
# r = [1, 1, 1, 1], I = []           #levels = TRIVIAL
# r = [1, 1, 2], I = [0]             #levels = TRIVIAL
# r = [1, 1, 2], I = [0, 1]          #levels = TRIVIAL
# r = [1, 1, 2], I = []              #levels = TRIVIAL
# r = [1, 1, 1, 1], I = [0, 1]       #levels = TRIVIAL
# r = [2, 2], I = [0]                #levels = TRIVIAL
# r = [2, 2], I = []                 #levels = TRIVIAL
# r = [1, 3], I = [0, 1]             #levels = TRIVIAL
# r = [4], I = [0]                   #levels = TRIVIAL
# r = [1, 3], I = []                 #levels = TRIVIAL
# r = [4], I = []                    #levels = TRIVIAL
# Count = 130
# """"
# 
#
# 
# output_N5_full_counts = """
# r = [1, 2, 2], I = []                    #levels = 73        
# r = [1, 1, 1, 1, 1], I = [0, 1]          #levels = 38768
# r = [1, 2, 2], I = [0, 1, 2]             #levels = 344
# r = [2, 3], I = [0]                      #levels = 20
# r = [1, 4], I = [0]                      #levels = 8
# r = [2, 3], I = []                       #levels = 10
# r = [1, 1, 1, 1, 1], I = []              #levels = 13862
# r = [1, 2, 2], I = [1]                   #levels = 146
# r = [2, 3], I = [0, 1]                   #levels = 34
# r = [5], I = [0]                         #levels = 2
# r = [1, 1, 1, 2], I = [0, 1, 2, 3]       #levels = 3524
# r = [1, 1, 1, 2], I = [0, 1, 2]          #levels = 2322
# r = [1, 1, 1, 2], I = [0]                #levels = 1082
# r = [1, 2, 2], I = [1, 2]                #levels = 252
# r = [1, 1, 3], I = [0, 2]                #levels = 130
# r = [1, 1, 1, 1, 1], I = [0]             #levels = 27724
# r = [1, 1, 1, 2], I = [0, 3]             #levels = 2124
# r = [1, 4], I = []                       #levels = 4
# r = [1, 1, 1, 2], I = [0, 1]             #levels = 1516
# r = [2, 3], I = [1]                      #levels = 20
# r = [1, 4], I = [0, 1]                   #levels = 12
# r = [1, 1, 3], I = []                    #levels = 34
# r = [1, 2, 2], I = [0, 1]                #levels = 216
# r = [1, 4], I = [1]                      #levels = 8
# r = [1, 1, 3], I = [0, 1, 2]             #levels = 154
# r = [1, 1, 3], I = [0]                   #levels = 68
# r = [1, 1, 1, 2], I = []                 #levels = 541
# r = [1, 1, 3], I = [2]                   #levels = 68
# r = [1, 1, 1, 1, 1], I = [0, 1, 2]       #levels = 59002
# r = [1, 1, 3], I = [0, 1]                #levels = 96
# r = [5], I = []                          #levels = 1
# r = [1, 1, 1, 2], I = [3]                #levels = 1082
# r = [1, 1, 1, 1, 1], I = [0, 1, 2, 3]    #levels = 86040
# r = [1, 1, 1, 2], I = [0, 1, 3]          #levels = 2784
# r = [1, 2, 2], I = [0]                   #levels = 146
# r = [1, 1, 1, 1, 1], I = [0, 1, 2, 3, 4] #levels = 120032
# """
