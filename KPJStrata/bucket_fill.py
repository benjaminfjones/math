# bucket_fill.sage
# Author: Benjamin Jones <benjaminfjones@gmail.com>
#
# load "/Users/jonesbe/Desktop/svn/sage_code/KPJStrata/bucket_fill.spyx"
# load "/Users/jonesbe/Desktop/svn/sage_code/KPJStrata/bucket_fill.sage"
#

from copy import copy
from sage.all import *
import bucket_fill_cython as bfc

# Fastest time so far: (for a Python generator version) 320 ms
# Fastest time so far: (for a Cython stored list version) 59 ms
#
# sage: timeit('C = list(bucket_vec_fill2_gen([5,4,3,2,1], 8))')
# 5 loops, best of 3: 320 ms per loop
# ==> about a 10% speedup compared to splitting out block_gen as a subroutine
#
def bucket_vec_fill2_gen(b, n):
    """
    Distribute n objects among buckets with capacity b[i].
    The result is a vector of length sum_i b[i] where each
    block of size b[i] is filled with entries 2, 1, -1, or 0
    --sorted in that order--
    having weights 2, 1, 1, and 0, respectively.

    The routine returns a generator object for such
    distributions.
    """
    # define local names to avoid lookup
    my_bucket_vec_fill = bucket_vec_fill2_gen
    my_block_fill = block_fill_c # implimented in cython
    my_range = range
    my_min = min
    my_max = max
    my_floor = floor
    l = len(b)
    # base case
    if l == 0:
        if n == 0:
            yield []
    # l > 0 case
    elif l > 1 or (l == 1 and 2*b[0] >= n):
        bz = b[l-1]
        bp = b[0:l-1]
        j = my_min(2*bz, n)
        for i in my_range(j+1):
            for v in my_bucket_vec_fill(bp, n-i):
                two_min = my_max(i-bz, 0)
                two_max = my_min(bz, my_floor(i/2))
                for h2 in my_range(two_min, two_max+1):
                    h1 = i-2*h2
                    r = bz-h2-h1
                    for k in my_range(h1+1):
                        yield v + my_block_fill(h2, k, h1-k, r)

def equal_part(p, q):
    """
    Assumes p and q are both partitions.

    Returns 'true' if they are equal partitions; they are allowed to have
    different numbers of trailing 0s.

    EXAMPLES::

        sage: equal_part([1],[1])
        True
        sage: equal_part([1],[1,0])
        True
        sage: equal_part([1,0,0],[1,0])
        True
        sage: equal_part([1,0,0],[1,0,1])
        False
    """
    lp = len(p)
    lq = len(q)
    if lp == lq:
        return p == q
    elif lp > lq:
        return p[0:lq] == q and p[lq] == 0
    else:
        return p == q[0:lp] and q[lp] == 0

def equal_bipart(p, q):
    """
    #
    # Assumes p and q are both bipartitions.  Returns 'true' if they are equal.
    #

    EXAMPLES::

        sage: equal_bipart([[1],[0,0,0]], [[1,0,0],[0]])
        True
        sage: equal_bipart([[3],[0]], [[3,1],[0]])
        False
    """
    return equal_part(p[0], q[0]) and equal_part(p[1], q[1])

def equal_level(L, M):
    """
    Determines if levels `L` and `M` are equal as lists of bipartitions.
    """
    return all(equal_bipart(p,q) for (p,q) in zip(L,M))

def uniqify_levels(L):
    """
    Returns a list of unique levels from the given list `L`.
    """
    R = []
    for l in L:
        inR = False
        for r in R:
            if equal_level(l,r):
                inR = True
                break
        if not inR:
            R.append(l)
    return R

def sum_partitions(mu, nu):
    """
    Returns the partition lam which is the sum of the partitions `mu` and `nu`
    which need not have the same number of parts.
    """
    lm = len(mu)
    ln = len(nu)
    if lm > ln:
        lam = [ mu[i] + nu[i] for i in range(ln) ]
        lam.extend(mu[ln:])
    else:
        lam = [ mu[i] + nu[i] for i in range(lm) ]
        lam.extend(nu[lm:])
    return lam

def bipartitionize( bp ):
    """
    Assumes bp is a pair of quasipartitions. The sum of its members need not be
    a partition.

    Returns the unique bipartition obtained by (1) sorting the rows so that
    the sum is indeed a partition, and then (2) moving individual rows left or
    right as necessary so that the members becomes honest partitions.

    This function is used to pick off the subordinate bipartitions of a signed
    quasibipartition.

    EXAMPLES::

        sage: bipartitionize([[1], [2,1]])
        [[1, 0], [2, 1]]
        sage: bipartitionize([[2,1], [1]])
        [[2, 1], [1, 0]]

        sage: bipartitionize([[1, 2], [1]])
        [[2, 2], [0, 0]]
        sage: bipartitionize([[3,2,3,2,1], [3,2,1,2]])
        [[3, 3, 3, 3, 1], [3, 1, 1, 1, 0]]
        sage: bipartitionize([[0], [3,2,1,2]])
        [[0, 0, 0, 0], [3, 2, 2, 1]]

    """
    # local names
    my_descent = bfc.descent_position2_c # Cython
    my_zero = bfc.zero_position2_c # Cython
    item_getter_0 = operator.itemgetter(0)
    item_getter_1 = operator.itemgetter(1)
    item_getter_2 = operator.itemgetter(2)
    # local copies to work on
    mu_t = copy(bp[0])
    nu_t = copy(bp[1])
    r = len(mu_t) - len(nu_t)
    # pad mu and nu
    if r >= 0:
        nu_t.extend([0]*r)
    else:
        mu_t.extend([0]*(-r))
    # sort mu_t and nu_t simultaneously by the sum of their rows
    l = [ (x,y,x+y) for (x,y) in zip(mu_t, nu_t) ]
    l = sorted(l, key=item_getter_2, reverse=True)
    mu_t = map(item_getter_0, l)
    nu_t = map(item_getter_1, l)
    # move rows right and left as neccesary
    j = my_descent(mu_t, nu_t)
    while j != -1:
        if mu_t[j] < mu_t[j+1]:
            nu_t[j] = nu_t[j] + mu_t[j] - mu_t[j+1]
            mu_t[j] = mu_t[j+1]
        else:
            mu_t[j+1] = mu_t[j+1] + nu_t[j+1] - nu_t[j]
            nu_t[j+1] = nu_t[j]
        j = my_descent(mu_t, nu_t)
    j = my_zero(mu_t, nu_t)
    # remove all `zero` rows
    if j >= 0:
        mu_t = mu_t[0:j]
        nu_t = nu_t[0:j]
    return [ mu_t, nu_t ]

def compatible_bipartitions(p, n, e, debug=False):
    """
    Generator object for bipartitions `r` such that there exists a signed
    quasi-bipartition `xi` having `p` and `r` as subordinate bipartitions.
    If e == True then `p` is the + bipartition and `r` is the -, otherwise the
    opposite.

    Note: repetitions can occur and are not checked for

    INPUT:
        - p - a bipartition [ part, part ]
        - n - a positive integer of sufficient size. n >= size(p) is garunteed to work
        - e - a boolean flag

    OUTPUT:
        - Generator object for bipartitions `r` = [ part, part ] (repetitions
          do occur)

    EXAMPLES:

        The following lists check with the GAP algorithm after removing
        repetitions::

        sage: P = list(compatible_bipartitions( [[1], [0]], 4, True ))
        sage: for p in P:
        ....:     print p
        ....:
        [[1, 0, 0], [1, 1, 1]]
        [[0, 0, 0, 0], [1, 1, 1, 1]]
        [[1, 1, 1, 1], [0, 0, 0, 0]]
        [[0, 0, 0, 0], [1, 1, 1, 1]]

        sage: P = list(compatible_bipartitions( [[1], [0]], 4, False ))
        sage: for p in P:
            print p
        ....:
        [[2, 1, 1], [0, 0, 0]]
        [[2, 1, 1], [0, 0, 0]]
        [[1, 1, 1, 1], [0, 0, 0, 0]]
        [[1, 1, 1, 1], [0, 0, 0, 0]]

        sage: P = list(compatible_bipartitions( [[2,1],[1,0]], 8, True ))
        sage: for p in P: print p
        ....:
        [[2, 1, 0, 0], [2, 1, 1, 1]]
        [[2, 0, 0, 0, 0], [2, 1, 1, 1, 1]]
        [[2, 1, 1, 1, 1], [2, 0, 0, 0, 0]]
        [[1, 1, 0, 0, 0], [2, 1, 1, 1, 1]]
        [[2, 1, 0, 0, 0], [1, 1, 1, 1, 1]]
        [[2, 0, 0, 0, 0], [2, 1, 1, 1, 1]]
        [[1, 0, 0, 0, 0, 0], [2, 1, 1, 1, 1, 1]]
        [[1, 1, 1, 1, 1, 1], [2, 0, 0, 0, 0, 0]]
        [[2, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1]]
        [[2, 1, 1, 1, 1, 1], [1, 0, 0, 0, 0, 0]]
        [[1, 1, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1]]
        [[1, 0, 0, 0, 0, 0], [2, 1, 1, 1, 1, 1]]
        [[2, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1]]
        [[1, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1]]
        [[1, 1, 1, 1, 1, 1, 1], [1, 0, 0, 0, 0, 0, 0]]
        [[1, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1]]

        sage: %timeit list(compatible_bipartitions( [[5,4,3,2,1],[1,1,1,1,1]], 30, True ))
        5 loops, best of 3: 2.69 s per loop
        sage: time P = list(compatible_bipartitions( [[5,4,3,2,1],[1,1,1,1,1]], 30, True ))
        CPU times: user 2.82 s, sys: 0.00 s, total: 2.66 s

        After optimizing subord_bipartitions_c (better than GAP!)::

        sage: %timeit list(compatible_bipartitions( [[5,4,3,2,1],[1,1,1,1,1]], 30, True ))
        5 loops, best of 3: 1.02 s per loop

    """
    q = [ copy(p[0]), copy(p[1])]
    lam = sum_partitions(q[0], q[1])
    xi = [ 2*x-1 for x in lam ] # construct partition xi to add boxes to
    r = n + sum(lam) - sum(xi) # number of extra boxes to add
    q[0].extend([0 for i in range(r+len(lam)-len(q[0]))])
    xi.extend([0 for i in range(r)])
    # if debug: print "lam = ", lam
    # if debug: print "xi = ", xi

    # count occurences of each distinc part of xi -> mults
    dxi = dict() # create a dictionary of occurences
    for x in xi:
        dxi[x] = dxi.get(x, 0) + 1
    occur = sorted(dxi.iteritems(), key=operator.itemgetter(0), reverse=True)
    mults = [ x2 for (x1,x2) in occur ]
    lm = len(mults)
    # if debug: print "mults = ", mults, " r = ", r

    # construct the different fillings of buckets (buckets = blocks of distinct parts of xi)
    add = []; add_append = add.append
    for i in range(min(r,mults[lm-1])+1):
        for v in bfc.bucket_vec_fill1_c(mults[lm-1], i):
            for t in bfc.bucket_vec_fill2_c(mults[0:lm-1], r-i):
                add_append(t+v)
    # if debug:
    #        print "size(add) = %d" % len(add)

    # if debug: print "START v loop:\n"
    for v in add:
        # if debug: print "v = ", v
        theta = filter(bfc.non_zero_c, [xi[i]+abs(v[i]) for i in range(len(xi))])
        # if debug: print "theta = ", theta
        # if debug: print "q[0] = ", q[0]
        # Construct mu, nu, and epsilon partitions
        mu = []; nu = []; ep = []
        for k in range(len(theta)):
            if q[0][k] == 0:
                mu.append(0)
                if e:
                    if v[k] >= 1: ep.append(-1)
                    else: ep.append(1)
                else:
                    if v[k] >= 1: ep.append(1)
                    else: ep.append(-1)
            else:
                if e:
                    if v[k] >= 1:
                        mu.append(2*q[0][k])
                        ep.append(-1)
                    else:
                        mu.append(2*q[0][k]-1)
                        ep.append(1)
                else:
                    if v[k] >= 1:
                        mu.append(2*q[0][k]+1)
                        ep.append(1)
                    else:
                        mu.append(2*q[0][k])
                        ep.append(-1)
            if mu[k] > theta[k]: mu[k] = theta[k]
            nu.append(theta[k] - mu[k])
        # END k loop
        # if debug:
            # print "mu = ", mu, " nu = ", nu, " ep = ", ep

        # form the subordinate bipartitions from the SQBP (mu, nu, ep)
        sb0 = bfc.subord_bipartitions_c(mu, nu, ep)
        sb = [ bipartitionize([sb0[0], sb0[1]]), bipartitionize([sb0[2], sb0[3]]) ]
        # if debug: print "sb = ", sb
        # if debug: print "\n"
        if e and equal_bipart(q, sb[0]):
            yield sb[1]
        elif (not e) and equal_bipart(q, sb[1]):
            yield sb[0]
        # (still in "v" loop)
        # include choices where the last block of +- 1's is on the mu side
        if not e:
            for k in [ x for x in range(len(theta)) if (q[0][x] == 0 and v[x] == 1) ]:
                mu[k] = min(1, theta[k])
                nu[k] = theta[k] - mu[k]
            sb0 = bfc.subord_bipartitions_c(mu, nu, ep)
            sb = [ bipartitionize([sb0[0], sb0[1]]), bipartitionize([sb0[2], sb0[3]]) ]
            if equal_bipart(q, sb[1]):
                yield sb[0]
    # END "v" loop

def bipartition_to_quiver(mu, nu):
    """
    Takes a bipartition `(mu, nu)` and returns the associated quiver data.
    Quiver data consists of a list $R = (R_i | i = 0, \ldots, t-1 )$ and an index
    set $I \subset \{0, \ldots, t-1\}$ which are determined by a bipartition
    $(\mu, \nu)$ as described below.

    ALGORITHM:

        - lam := mu + vu (a partition)
        - t := lam[0] (largest part of lam)
        - Construct R[0], R[1], ... the list of column sizes of
          lam in increasing order.
        - I := subset of {0, ... , t-1} where i in I iff.
          R[i] is a column length of `mu` and R[i] = R[i+1] with i+1 in I
          forces i in I (so columns lengths of `mu` appear first in the
          sequence of R[j]'s)
        - return [ R, I ]

    EXAMPLES::

    In this example, `lam' is [4,2], `lam` conjugate is [2,2,1,1]::

        sage: bipartition_to_quiver([3,2],[1])
        [[1, 1, 2, 2], [0, 2, 3]]

    In this example, `lam` is [6,4,2], `lam` conjugate is
    [3,3,2,2,1,1]::

        sage: bipartition_to_quiver([3,2,1],[3,2,1])
        [[1, 1, 2, 2, 3, 3], [0, 2, 4]]

    """
    if type(mu) != sage.combinat.partition.Partition_class:
        muP = Partition(mu)
        nuP = Partition(nu)
    else:
        muP = mu
        nuP = nu
    lam = Partition(sum_partitions(muP, nuP))
    t = lam[0]
    R = sorted(lam.conjugate().to_list())
    mu_cols = sorted(muP.conjugate().to_list())
    I = []
    for i in range(t):
        if R[i] in mu_cols:
            I.append(i)
            j = mu_cols.index(R[i])
            del(mu_cols[j])
    return [R, I]

def quivers(n):
    """
    Return a list of all "quiver data: of total size `n`.
    See the function `bipartition_to_quiver`.

    INPUT:

        - n - an integer > 0

    OUTPUT:

        - A generator object yielding quivers of size `n`

    EXAMPLES::

        sage: for s in quivers(3):
        ....:     print s
        ....:
        [[1, 1, 1], []]
        [[1, 2], []]
        [[3], []]
        [[1, 1, 1], [0]]
        [[1, 2], [0]]
        [[1, 1, 1], [0, 1]]
        [[1, 2], [1]]
        [[1, 1, 1], [0, 1, 2]]
        [[1, 2], [0, 1]]
        [[3], [0]]

        sage: time L = list(quivers(20))
        CPU times: user 1.76 s, sys: 0.00 s, total: 1.76 s
        Wall time: 1.76 s

    """
    for i in range(n+1):
        for (mu,nu) in CartesianProduct(Partitions(i), Partitions(n-i)):
            yield bipartition_to_quiver(mu, nu)

def quiver_is_trivial(r,I):
    """Returns True if I is [ 0, 1, ..., len(I)-1 ]"""
    return I == range(len(I))

def quiver_is_cotrivial(r,I):
    """Returns True if I is a consequtive list: [ len(r)-len(I), ..., len(r) ]"""
    return I == range(len(r)-len(I), len(r))

def kpj_levels( r, I, unique=False ):
    """
    Returns a list of all KPJ levels associated to the quiver data (r, I).
    A "level" is a list of sequentially compatible bipartitions (see
    `compatible_bipartitions`)

    INPUT:
        - r - a list of `n` positive integers
        - I - a list of

    EXAMPLES::

        sage: for l in kpj_levels([1,1],[0,1], unique=True):
            print l
        ....:
        [[[], []], [[0], [1]], [[0], [2]]]
        [[[], []], [[0], [1]], [[0, 0], [1, 1]]]
        [[[], []], [[0], [1]], [[1, 1], [0, 0]]]
        [[[], []], [[1], [0]], [[2], [0]]]
        [[[], []], [[1], [0]], [[1, 1], [0, 0]]]

        sage: for l in kpj_levels([1,1],[0], unique=True):
            print l
        ....:
        [[[], []], [[0], [1]], [[0], [2]]]
        [[[], []], [[0], [1]], [[0, 0], [1, 1]]]
        [[[], []], [[1], [0]], [[1], [1]]]
        [[[], []], [[1], [0]], [[0, 0], [1, 1]]]
        [[[], []], [[1], [0]], [[1, 1], [0, 0]]]

        sage: for l in kpj_levels([1,1],[], unique=True):
            print l
        ....:
        [[[], []], [[0], [1]], [[0], [2]]]
        [[[], []], [[0], [1]], [[0, 0], [1, 1]]]

        sage: for l in kpj_levels([1,1],[1], unique=True):
            print l
        ....:
        [[[], []], [[0], [1]], [[0], [2]]]
        [[[], []], [[0], [1]], [[0, 0], [1, 1]]]
        [[[], []], [[0], [1]], [[1, 1], [0, 0]]]

    """
    t = len(r)
    if t == 0:
        return [ [ [ [ ], [ ] ] ] ]
    levels0 = kpj_levels( r[:t-1], filter(lambda i:i < t-1, I))
    n = sum(r)
    result = []; result_append = result.append
    compats = []
    for L in levels0:
        for w in compatible_bipartitions(L[t-1], n, not (t-1 in I)):
            result_append(L + [w])
    if unique:
        return uniqify_levels(result)
    return result

def kpj_levels_gen(r, I):
    """
    Same as `kpj_levels`, but returns a generator object.
    """
    t = len(r)
    if t == 0:
        yield [ [ [ ], [ ] ] ]
    else:
        levels0 = kpj_levels_gen( r[:t-1], filter(lambda i:i < t-1, I))
        n = sum(r)
        for L in levels0:
            for w in compatible_bipartitions(L[t-1], n, not (t-1 in I)):
                yield L + [w]
