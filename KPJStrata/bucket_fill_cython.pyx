# LOG
# TODO 3/17/2011 optimize subord_bipartitions_c more (see profile)

def block_fill_c(int n1, int n2, int n3, int n4):
    cdef int x
    L = [2 for x in range(n1)]
    L.extend([1 for x in range(n2)])
    L.extend([-1 for x in range(n3)])
    L.extend([0 for x in range(n4)])
    return L

# Fastest time so far:
# sage: timeit('C = bucket_vec_fill2_c([5,4,3,2,1],8)')
# 5 loops, best of 3: 59 ms per loop
#
# After eliminating function call to block_fill_c:
#
# sage: %timeit C = bucket_vec_fill2_c([5,4,3,2,1],8)
# 25 loops, best of 3: 27.7 ms per loop
#
def bucket_vec_fill2_c(b, int n):
    """
    Distribute n objects among buckets with capacity b[i].
    The result is a vector of length sum_i b[i] where each
    block of size b[i] is filled with entries 2, 1, -1, or 0
    --sorted in that order--
    having weights 2, 1, 1, and 0, respectively.

    EXAMPLES::

    sage: from bucket_fill_cython import bucket_vec_fill2_c
    sage: bucket_vec_fill2_c([3,2,1], 10)
    [[2, 2, 2, 2, 2, 0],
    [2, 2, 2, 2, -1, -1],
    [2, 2, 2, 2, -1, 1],
    [2, 2, 2, 2, 1, -1],
    [2, 2, 2, 2, 1, 1],
    [2, 2, -1, 2, 2, -1],
    [2, 2, -1, 2, 2, 1],
    [2, 2, 1, 2, 2, -1],
    [2, 2, 1, 2, 2, 1],
    [2, 2, 2, -1, -1, 2],
    [2, 2, 2, 1, -1, 2],
    [2, 2, 2, 1, 1, 2],
    [2, 2, 2, 2, 0, 2],
    [2, 2, -1, 2, -1, 2],
    [2, 2, -1, 2, 1, 2],
    [2, 2, 1, 2, -1, 2],
    [2, 2, 1, 2, 1, 2],
    [2, -1, -1, 2, 2, 2],
    [2, 1, -1, 2, 2, 2],
    [2, 1, 1, 2, 2, 2],
    [2, 2, 0, 2, 2, 2]]
    """
    cdef int l = len(b)
    cdef int bz, i, two_min, two_max, h1, h2, x
    # base case
    if l == 0:
        if n == 0:
            return [[]]
        else:
            return []
    # l > 0 case
    elif l > 1 or (l == 1 and 2*b[0] >= n):
        R = []; R_append = R.append
        bz = b[l-1]
        bp = b[0:l-1]
        t = [0 for x in range(bz)]
        j = min(2*bz, n)
        for i in range(j+1):
            for v in bucket_vec_fill2_c(bp, n-i):
                two_min = max(i-bz, 0)
                if i % 2 == 0:
                    two_max = min(bz, i/2)
                else:
                    two_max = min(bz, (i-1)/2)
                for h2 in range(two_min, two_max+1):
                    h1 = i-2*h2
                    for k in range(h1+1):
                        for x in range(0,h2):
                            t[x] = 2
                        for x in range(h2, h2+k):
                            t[x] = 1
                        for x in range(h2+k, h2+h1):
                            t[x] = -1
                        for x in range(h2+h1, bz):
                            t[x] = 0
                        R_append(v + t)
        return R
    else:
        return []

def bucket_vec_fill1_c(int b, int n):
    """
    Distribute n +/- 1's `b` spaces. The result is a vector of length `b`
    with entries 1, -1, or 0 --sorted in that order--

    The routine returns a list of such vectors.

    EXAMPLES::

        sage: bucket_vec_fill1_c(3,3)
        [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [1, 1, 1]]
        sage: bucket_vec_fill1_c(3,4)
        []
        sage: bucket_vec_fill1_c(5,3)
        [[-1, -1, -1, 0, 0], [1, -1, -1, 0, 0], [1, 1, -1, 0, 0], [1, 1, 1, 0, 0]]

    """
    cdef int i, j
    if b == 0:
        if n == 0:
            return [[]]
        else:
            return []
    elif n <= b:
        R = []; R_append = R.append
        for i in range(n+1):
            t = []; t_append = t.append
            for x in range(i):
                t_append(1)
            for x in range(n-i):
                t_append(-1)
            for x in range(b-n):
                t_append(0)
            R_append(t)
        return R
    else:
        return []

# sage: timeit subord_bipartitions_c( [2,3,2,1,2,0,1], [4,2,3,3,1,2,0], [-1,1,-1,1,-1,-1,1] )
# 625 loops, best of 3: 71.8 µs per loop
#
# After eliminating function call to floor_over_two_c:
# sage: %timeit subord_bipartitions_c( [2,3,2,1,2,0,1], [4,2,3,3,1,2,0], [-1,1,-1,1,-1,-1,1] )
# 625 loops, best of 3: 9.67 µs per loop

def subord_bipartitions_c( mu, nu, ep):
    """
    Returns the + and - subordinate bipartitions associated to the signed
    quasi-bipartition (mu, nu, ep).

    EXAMPLES::

        sage: subord_bipartitions_c([3],[2],[1])
        [[[2], [1]], [[1], [1]]]
        sage: subord_bipartitions_c( [2,3], [4,2], [-1, 1] )
        [[[2, 2], [1, 1]], [[1, 1], [2, 1]]]
        sage: subord_bipartitions_c( [2,3,2,1,2,0,1], [4,2,3,3,1,2,0], [-1,1,-1,1,-1,-1,1] )
        [[[2, 2, 1, 1, 1, 1, 1], [1, 1, 1, 1, 0, 0, 0]], [[1, 1, 1, 1, 1, 0], [2, 2, 1, 1, 1, 1]]]

    """
    cdef int i
    cdef int n = len(mu)
    mup = [0 for i in range(n)]
    nup = [0 for i in range(n)]
    mum = [0 for i in range(n)]
    num = [0 for i in range(n)]
    mui = [ <int>mu[i] for i in range(n) ]
    nui = [ <int>nu[i] for i in range(n) ]
    epi = [ <int>ep[i] for i in range(n) ]
    for i in range(n):
        if mui[i] > 0 and ((epi[i] == 1 and mui[i] % 2 == 0) or
                           (epi[i] == -1 and mui[i] % 2 == 1)):
            mup[i] = mui[i]/2
            nup[i] = (mui[i] + nui[i] + (1+epi[i])/2)/2 - mup[i]
            mum[i] = (mui[i]-1)/2
            num[i] = (mui[i] + nui[i] + (1-epi[i])/2)/2 - mum[i]
        else:
            mup[i] = (mui[i]+1)/2
            nup[i] = (mui[i] + nui[i] + (1+epi[i])/2)/2 - mup[i]
            mum[i] = mui[i]/2
            num[i] = (mui[i] + nui[i] + (1-epi[i])/2)/2 - mum[i]
    return [ mup, nup, mum, num ]

## Beware: undocumented helper functions below:

cpdef int descent_position_c(L):
    cdef int i
    for i in range(len(L)-1):
        if L[i] < L[i+1]:
            return i
    return -1

cpdef int descent_position2_c(L,M):
    cdef int i
    for i in range(len(L)-1):
        if L[i] < L[i+1] or M[i] < M[i+1]:
            return i
    return -1

cpdef int zero_position2_c(L,M):
    cdef int i
    for i in range(len(L)):
        if L[i] == 0 and M[i] == 0:
            return i
    return -1

cpdef non_zero_c(int n):
    if n == 0:
        return False
    else:
        return True
