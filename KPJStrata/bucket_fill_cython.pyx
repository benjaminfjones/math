# bucket_fill_cython.pyx
# Author: Benjamin Jones <benjaminfjones@gmail.com>
#
# This file contains Cython optimized routines which are used in block_fill.py.
#

from libc.stdlib cimport calloc, malloc, free

cdef void block_fill_c(int *s, int n1, int n2, int n3, int n4):
    """
    Fill array `s` with the indicated number of 2's, 1's, -1's, and 0's.
    """
    cdef int x, i
    x = n1+n2+n3+n4
    for i from 0 <= i < n1:
        s[i] = 2
    for i from n1 <= i < n1+n2:
        s[i] = 1
    for i from n1+n2 <= i < x-n4:
        s[i] = -1
    for i from x-n4 <= i < x:
        s[i] = 0


def bucket_vec_fill2_c(b, int n):
    """
    Distribute n objects among buckets with capacity b[i].
    The result is a vector of length sum_i b[i] where each
    block of size b[i] is filled with entries 2, 1, -1, or 0
    --sorted in that order--
    having weights 2, 1, 1, and 0, respectively.

    The routine returns a list of all such distributions.

    TIMING: timeit('C = bucket_vec_fill2_c([5,4,3,2,1],8)')

    Original code:
    5 loops, best of 3: 59 ms per loop

    After eliminating function call to block_fill_c:
    25 loops, best of 3: 27.7 ms per loop

    After factoring out block_fill_c to use C int array
    25 loops, best of 3: 25.6 ms per loop
    """
    cdef int l  # length of b
    cdef int bz # last entry of b
    cdef int *s # temporary storage for block_fill_c
    cdef i, j, two_min, two_max, h1, h2, x, k
    cdef list R, v # python lists for return results

    l = len(b)
    # base case
    if l == 0:
        if n == 0:
            return [[]]
        else:
            return []
    # l > 0 case
    elif l > 1 or (l == 1 and 2*b[0] >= n):
        R = []
        R_append = R.append # avoid lookup
        bz = b[l-1]
        bp = b[0:l-1]
        j = min(2*bz, n)
        s = <int*> malloc(bz * sizeof(int)) # allocate temp storage
        for i from 0 <= i < j+1:
            for v in bucket_vec_fill2_c(bp, n-i):
                two_min = max(i-bz, 0)
                if i % 2 == 0:
                    two_max = min(bz, i/2)
                else:
                    two_max = min(bz, (i-1)/2)
                for h2 in range(two_min, two_max+1):
                    h1 = i-2*h2
                    for k in range(h1+1):
                        block_fill_c(s, h2, k, h1-k, bz-h2-h1) # factored out
                        R_append(v + [ s[x] for x in range(bz) ])
        free(s)
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

    TIMING:

        timeit('C = bucket_vec_fill1_c(200,100)')
        625 loops, best of 3: 537 µs per loop
    """
    cdef int i, x
    cdef list R
    if b == 0:
        if n == 0:
            return [[]]
        else:
            return []
    elif n <= b:
        R = []
        R_append = R.append
        s = <int*> malloc(b * sizeof(int))
        for i from 0 <= i < n+1:
            for x from 0 <= x < i:
                s[x] = 1
            for x from i <= x < n:
                s[x] = -1
            for x from n <= x < b:
                s[x] = 0
            R_append([ s[x] for x in range(b) ])
        return R
    else:
        return []


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

    TIMING: %timeit subord_bipartitions_c( [2,3,2,1,2,0,1], [4,2,3,3,1,2,0], [-1,1,-1,1,-1,-1,1] )

        Original code:
        625 loops, best of 3: 71.8 µs per loop

        After eliminating function call to floor_over_two_c:
        625 loops, best of 3: 9.67 µs per loop

        After allocating int* arrays manually:
        625 loops, best of 3: 6.78 µs per loop
    """
    cdef int i
    cdef int n = len(mu)
    cdef int *mup, *nup, *mum, *num, *mui, *nui, *epi
    cdef list r_mup, r_nup, r_mum, r_num
    mup = <int*> malloc(n * sizeof(int))
    nup = <int*> malloc(n * sizeof(int))
    mum = <int*> malloc(n * sizeof(int))
    num = <int*> malloc(n * sizeof(int))
    mui = <int*> malloc(n * sizeof(int))
    nui = <int*> malloc(n * sizeof(int))
    epi = <int*> malloc(n * sizeof(int))
    for i from 0 <= i < n:
        mui[i] = <int>mu[i]
        nui[i] = <int>nu[i]
        epi[i] = <int>ep[i]
    for i from 0 <= i < n:
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
    r_mup = [ mup[i] for i in range(n) ]
    r_nup = [ nup[i] for i in range(n) ]
    r_mum = [ mum[i] for i in range(n) ]
    r_num = [ num[i] for i in range(n) ]
    free(mup)
    free(nup)
    free(mum)
    free(num)
    free(mui)
    free(nui)
    free(epi)
    return [ r_mup, r_nup, r_mum, r_num ]

# Below are various small helper functions, these could be inline.

cpdef int descent_position_c(L):
    cdef int i
    for i from 0 <= i < len(L)-1:
        if L[i] < L[i+1]:
            return i
    return -1

cpdef int descent_position2_c(L,M):
    cdef int i
    for i from 0 <= i < len(L)-1:
        if L[i] < L[i+1] or M[i] < M[i+1]:
            return i
    return -1

cpdef int zero_position2_c(L,M):
    cdef int i
    for i from 0 <= i < len(L):
        if L[i] == 0 and M[i] == 0:
            return i
    return -1

cpdef non_zero_c(int n):
    if n == 0:
        return False
    else:
        return True

