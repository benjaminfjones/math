def block_gen(b, n):
    """
    Fill in `b` spaces with elements 2, 1, -1, 0 so that 
    their absolute values add to `n`. The resulting blocks
    are sorted by descending absolute value. 
    
    The routine returns a generator object for such 
    blocks.
    """
    if b > 0:
        # define local names to avoid lookup
        myrange = range
        my_block_fill = block_fill_c # implimented in cython
        # when n > b, min number of 2's `m` satisfies: m >= ceil((n-b)/2)
        # and satisfies: n-2m <= b-m  i.e. m >= n-b
        # On the othe hand, the max number of 2's `j` satisfies:
        # m <= b and m <= floor(n/2)
        two_min = max(n-b, 0)
        two_max = min(b, floor(n/2))
        for h2 in myrange(two_min, two_max+1):
            h1 = n-2*h2
            r = b-h2-h1
            for k in myrange(h1+1):
                yield my_block_fill(h2, k, h1-k, r)
                # yield [2]*h2 + [1]*k + [-1]*(h1-k) + [0]*(b-h2-h1)
                # note: [2 for x in myrange(h2)] is a factor of ~5x faster than [2]*h2
                # yield [2 for x in myrange(h2)] + [1 for x in myrange(k)] + [-1 for x in myrange(h1-k)] + [0 for x in myrange(b-h2-h1)]

#
# fastest version so far:
#
# sage: timeit('C = list(bucket_vec_fill_gen([5,4,3,2,1], 8))')
# 5 loops, best of 3: 359 ms per loop
#
def bucket_vec_fill_gen(b, n):
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
    my_block_gen = block_gen
    my_bucket_vec_fill = bucket_vec_fill_gen
    my_range = range
    l = len(b)
    # base case
    if l == 0:
        if n == 0:
            yield []
    # l > 0 case
    elif l > 1 or (l == 1 and 2*b[0] >= n):
        bz = b[l-1]
        bp = b[0:l-1]
        j = min(2*bz, n)
        for i in my_range(j+1):
            for v in my_bucket_vec_fill(bp, n-i):
                for t in my_block_gen(bz, i):
                    yield v+t

def subord_bipartitions( mu, nu, ep):
    """
    Returns the + and - subordinate bipartitions associated to the signed
    quasi-bipartition (mu, nu, ep).
    
    EXAMPLES::
    
        sage: subord_bipartitions([3],[2],[1])
        [[[2], [1]], [[1], [1]]]
        sage: subord_bipartitions( [2,3], [4,2], [-1, 1] )
        [[[2, 2], [1, 1]], [[1, 1], [2, 1]]]
        sage: subord_bipartitions( [2,3,2,1,2,0,1], [4,2,3,3,1,2,0], [-1,1,-1,1,-1,-1,1] )
        [[[2, 2, 1, 1, 1, 1, 1], [1, 1, 1, 1, 0, 0, 0]], [[1, 1, 1, 1, 1, 0], [2, 2, 1, 1, 1, 1]]]
        
    """
    # local names
    my_floor_over_two = floor_over_two_c # Cython: returns floor(x/2)
    my_floor = floor
    my_sub = operator.sub
    mu_t = copy(mu)
    nu_t = copy(nu)
    mr = range(len(mu))
    for i in mr:
        if mu[i] > 0 and ((ep[i] == 1 and mu[i] % 2 == 0) or
                          (ep[i] == -1 and mu[i] % 2 == 1)):
            mu_t[i] -= 1
            nu_t[i] += 1
    lp = [ my_floor_over_two(mu_t[i] + nu_t[i] + 1/2 + ep[i]/2) for i in mr ]
    lm = [ my_floor_over_two(mu_t[i] + nu_t[i] + 1/2 - ep[i]/2) for i in mr ]
    mup = [ my_floor_over_two(mu_t[i] + 1) for i in mr ]
    nup = map(my_sub, lp, mup)
    mum = [ my_floor_over_two(mu_t[i]) for i in mr ]
    num = map(my_sub, lm, mum)
    return [ bipartitionize([mup, nup]), bipartitionize([mum, num]) ]