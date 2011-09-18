import networkx
from copy import copy

#################################
# Partition funcitons
#################################

def dcoeff(l,m,N):	
    """Returns the 'D' coefficient for partitions l, m having exactly N parts"""
    # make padded copies of l, m
    l1 = copy(l)
    m1 = copy(m)
    l1.extend([0]*(N-len(l)))
    m1.extend([0]*(N-len(m)))
    # for the defining matrix
    M = Matrix( [ [ binomial(l1[i]+N-(i+1), m1[j]+N-(j+1)) for j in range(N) ] for i in range(N) ])
    return det(M)

def dual(lam, l, r):
	"""Returns dual of partition lam in l x r rectangle"""
	l1=copy(lam)
	l1.extend([0]*(r-len(lam)))
	dl=[0]*l
	for i in range(l):
		dl[i]=r-l1[i]
	dl.reverse()
	return trunc_part(dl)

def conj(lam):
	"""Returns the conjugate of partition lam"""
	(za,zb) = to_zel_part(lam)
	za.reverse()
	zb.reverse()
	return from_zel_part(zb,za)

def is_sub_part(lam,mu):
	"""Returns True if lam >= mu"""
	s=len(lam)
	if s<len(mu):
		return False
	mu1=copy(mu)
	mu1.extend([0]*(s-len(mu)))
	for i in range(len(lam)):
		if lam[i] < mu1[i]:
			return False
	return True

def trunc_part(p):
    """Truncates trailing zeros from a partition"""
    if 0 in p:
        sp = p[:p.index(0)]
    else:
        sp = p[:]
    return sp


def flagged_partition_iterator(l,r):
	"""Returns a 2 part flag of partitions in a l x r rectangle"""
	PB = PartitionsInBox(l,r)
	PB2 = PartitionsInBox(l,r)
	for lam in PB:
		for mu in PB2:
			llam = list(lam)
			lmu = list(mu)
			if is_sub_part(llam,lmu):
				yield (llam,lmu)

#################################
# Functions specific to Zelevinsky 
# resolution combinatorics
#################################

def to_zel_part(part):
    """Returns the Zelevinsky data for the partition 'part'"""
    sp = trunc_part(part)
    if sp == []:
        return ([],[])
    
    sp.reverse()
    b = [ sp[0] ]
    a = []
    acount = 0
    for i in range(len(sp)):
        if sp[i] > sum(b):
            b.append(sp[i]-sum(b))
            a.append(acount)
            acount=1
        else:
            acount += 1
    a.append(acount)
    return (a,b)

def from_zel_part(za,zb):
    """Return a usual partition from Zelevinsky data"""
    p = []
    for i in range(len(za)):
        b = sum(zb[:i+1])	
        p = p + [ b for j in range(za[i]) ]
    p.reverse()	
    return trunc_part(p)


def remove_peak(part, j):
	"""Returns the partition P(j) obtained by removing peak 'j' from
	the partition P"""
	za, zb = to_zel_part(part)
	l=len(za)
	if l == 1:
		return [0]
	elif j == 1:
		za = za[1:]
		zb[j-1:j+1] = [ zb[j-1]+zb[j] ]
	elif j == l:
		za[j-2:j]   = [ za[j-2]+za[j-1] ]
		zb = zb[:-1]
	else:
		za[j-2:j]   = [ za[j-2]+za[j-1] ]
		zb[j-1:j+1] = [ zb[j-1]+zb[j] ]
	return from_zel_part(za,zb)
			
def incidence_graph(plist, n=0):
    """Creates an incidence graph from: n = 'initial chain size'
    and plist = 'a permutation'"""
    if n == 0:
        n = len(plist)+1
    # First, create the inital chain
    d = {}
    for i in range(n-1):
        d.update( { i:[i+1] } )
    g = networkx.DiGraph(d)
    # node list
    nlist = range(n)
    # Next, add extra nodes
    for i in range(len(plist)):
        ni = i + n
        p = plist[i]
        # add new nodes and edges
        g.add_node(ni)
        g.add_path([ nlist[p-1], ni, nlist[p] ])
        # modify the node list
        nlist[p-1:p+1] = [ni]
    return g


def ig_weights(part, plist):
	"""Returns a list of weights, one for each node (in order)
	of the graph incidence_graph( n_peaks, plist)."""
	part_list = [ part ]
	z = [ to_zel_part(part) ]
        # initial weights
        w = [ [ sum(z[-1][0][:i]) + sum(z[-1][1][:i]) for i in range(len(z[-1][0])+1) ] ]
        # weights for the nodes off the initial chain
	for j in range(len(plist)):
            part_list.append( remove_peak( part_list[-1], plist[j] ) )
            z.append( to_zel_part( part_list[-1] ) )
            # shift weights more if the 1st peak is removed
            if plist[j] == 1:
                shift = w[-1][0] + z[-2][0][0]
            # standard shift
            else:
                shift = w[-1][0]
            w.append([ sum(z[-1][0][:i]) + sum(z[-1][1][:i]) + shift \
                       for i in range(len(z[-1][0])+1) ])
        # weights for the nodes in initial chain
        d = copy(w[0])
        # weights for the nodes below the chain
        for i in range(len(plist)):
            d.append(w[i+1][plist[i]-1])
	return d


def get_bounds(g,i):
    """Return the first predecessor and successor of node i in g among
       all those nodes connected to i with smaller node number.
    """
    pred_list = g.predecessors(i)
    succ_list = g.successors(i)
    pred_list = [ j for j in pred_list if j < i ]
    succ_list = [ j for j in succ_list if j < i ]
    return (pred_list[0], succ_list[0])



def tensor_term_iterator(g,w):
	"""Returns a coefficient, a list of S partitions, and a list of Q
	partitions which parametrize one term in the chern class of the
	tensor product of the bundles specified by the input data.
	INPUT: g = incidence_graph( ... )
           w = ig_weights( ... , ... )
    OUTPUT: a tuple ( coeff, S_list, Q_list )
	"""
	n = (len(w)+1)/2
	l = []
	r = []
	FL_part = []
	for i in range(n-1):
		(pred, succ) = get_bounds(g,n+i)
		node_weight = w[n+i]
		pred_weight = w[pred]
		succ_weight = w[succ]
		l.append(node_weight - pred_weight)
		r.append(succ_weight - node_weight)
		# generate 2 part flags in the rectangle l_i x r_i
		FL_part.append(flagged_partition_iterator(l[-1],r[-1]))
	# now for each node, we choose a flag of partitions and form an S_list, Q_list with coeff.

	CPI = cartesian_product_iterator( [ list(a) for a in FL_part ] )
	for c in CPI:
		S_list=[]
		Q_list=[]
		coeff=1
		for i in range(n-1):
			(lam,mu) = c[i]
			coeff = coeff * dcoeff(lam,mu,r[i])
			S_list.append(mu)
			Q_list.append(dual(conj(lam),l[i],r[i]))
		yield (coeff,S_list,Q_list)


#################################
# Functions specific to localization
#################################

def ic_fixed_points(w):
    """Returns the unique fixed point data assoc. to the nodes in an abstract
    initial chain.

    INPUT: n = init. chain size, w = list of weights for nodes in init. chain

    OUTPUT: one fixed point, i.e. a list of index sets, one for each ndoe in the
    init. chain.
    """
    return [ map(lambda x:x+1, range(i)) for i in w ]

def fixed_points(g, w, n, i):
    """Returns a list of fixed points of the (projection of the) incidence
    variety assoc. to g,w (where we project to include only nodes 1..i). This
    allows us to build a list of fixed points for the whole incidence variety
    inductively.

    INPUT:
    g = incidence graph (result of calling incidence_graph(n,plist)
    w = weights for the nodes (result of calling ig_weights(part,plist)
        for the same plist)
    n = length of init. chain (optional). defaults to (len(w)+1)/2
    i = How many nodes to include

    OUTPUT: [ fp_1, fp_2, ... ] where each fp_i is itself a list of
    index sets, one for each node from 1..i. An index set is a list of
    indices which determines the character of the $T$ action on the
    fiber of the corresponding universal bundle at the fixed point
    fp_i.

    NOTES:
    Tends to bog down around 8 levels of recursion, i.e.
    p=[8,7,6,5,4,3,2,1] # part and plist will be the same here
    g=incidence_graph(9,p)
    w=ig_weights(p,p)
    FP = fixed_points(g,w,9,7)   # runs for 10 sec
    FP = fixed_points(g,w,9,8)   # runs for a LONG time

    or...

    sage: plist=[8,7,6,5,4,3,2,1]
    sage: g = incidence_graph(plist)
    sage: w = ig_weights(plist,plist)
    sage: time FP2 = fixed_points(g,w,len(plist)+1,len(plist)-2)
    sage: time FP1 = fixed_points(g,w,len(plist)+1,len(plist)-1)
    sage: time FP0 = fixed_points(g,w,len(plist)+1,len(plist))
    Time: CPU 0.54 s, Wall: 0.54 s
    Time: CPU 9.14 s, Wall: 9.15 s
    Time: CPU 557.46 s, Wall: 558.66 s
    sage: len(FP2)
    5040
    sage: len(FP1)
    40320
    sage: len(FP0)
    362880
    """
    ic_size = len(w[:n])
    ic_fix = ic_fixed_points(w[:n])
    
    (pred,succ) = get_bounds(g,n+i-1)
    node_weight = w[n+i-1]
    pred_weight = w[pred]
    weight_diff = node_weight - pred_weight

    # base of induction
    if i == 1:
        pred_index_set = set(ic_fix[pred])
        succ_index_set = set(ic_fix[succ])
        diff_index_set = succ_index_set.difference(pred_index_set)
        C = combinations_iterator( list(diff_index_set), weight_diff )
        poss = [ [ic_fix[pred]+c]  for c in C ]
        return poss
    elif i > 1:
        # get the result for i-1
        induct = fixed_points(g,w,n,i-1)
        poss = []
        for fp in induct:
            fpp = ic_fix + fp  # tack on the initial data
            pred_index_set = set(fpp[pred])
            succ_index_set = set(fpp[succ])
            diff_index_set = succ_index_set.difference(pred_index_set)
            C = combinations_iterator( list(diff_index_set), weight_diff )
            newfps = [ fp+[fpp[pred]+c] for c in C ]
            poss = poss + newfps
        return poss
    else:
        return None

def fixed_points_iterator(g, w, n, i):
    """Same as fixed_points() but acts as an interator.
    See doc string for fixed_points().

    NOTES: The iterator takes dramatically less time than
    fixed_points() when i gets large. For example:

    sage: ## TEST THE ITERATOR
    sage: plist=[8,7,6,5,4,3,2,1]
    sage: g = incidence_graph(plist)
    sage: w = ig_weights(plist,plist)
    sage: FP = fixed_points_iterator(g,w,len(plist)+1,len(plist))
    sage: t0 = cputime()
    sage: count=0
    sage: for i,fp in enumerate(FP):
    ...       if i%50000 == 0:
    ...           print cputime()-t0,",",
    ...
    sage: print "Total Time (iterator):",cputime()-t0
    0.000598000000011 , 6.469392 , 12.936641 , 19.405214 ,
    25.886539 , 32.357304 , 38.827102 , 45.297906 , Total Time
    (iterator): 46.964347

    Compare to 557.46 s for the same call to fixed_points()
    """
    # initial chain length: node labels in g are shifted by this amount
    ic_size = len(w[:n])
    ic_fix = ic_fixed_points(w[:n])
    
    (pred,succ) = get_bounds(g,n+i-1)
    node_weight = w[n+i-1]
    pred_weight = w[pred]
    weight_diff = node_weight - pred_weight

    # base of induction
    if i == 1:
        pred_index_set = set(ic_fix[pred])
        succ_index_set = set(ic_fix[succ])
        diff_index_set = succ_index_set.difference(pred_index_set)
        C = combinations_iterator( list(diff_index_set), weight_diff )
        for c in C:
            yield [ic_fix[pred] + c]

    elif i > 1:
        # get a fixed point associated to i-1
        FPI = fixed_points_iterator(g,w,n,i-1)
        for fp in FPI:
            fpp = ic_fix + fp  # tack on the initial data
            pred_index_set = set(fpp[pred])
            succ_index_set = set(fpp[succ])
            diff_index_set = succ_index_set.difference(pred_index_set)
            C = combinations_iterator( list(diff_index_set), weight_diff )
            for c in C:
                yield fp+[fpp[pred]+c]
    else:
        yield None

def evaluate_schur(lamb,X):
    """Evaluates s_\lambda(B) where B has C^* character X"""
    # print "evaluate_schur(",lamb,",",X,")"
    s=SFASchur(QQ)
    if trunc_part(lamb) == []:
        return 1
    else:
        f=s(lamb)
        fp=f.expand(len(X))
        return fp(X)

def bott_denominator(g,w,fp,X):
    """Calculate the top Chern class of the normal bundle to an isolated 
       fixed point

    INPUT: g = incidence graph
           w = node 'weights' (result of ig_weights() call)
           fp = fixed point data for a single fixed point
           X = character of C^* on C^n (a list of 'n' integers)

    OUTPUT: an integer
    """
    n = (len(w)+1)/2   # size of the initial chain in g
	ic_fix = ic_fixed_points(w[:n])
    denominator = 1

    # loop over the nodes
    for i in range(n-1):
         # get succ and pred
 		(pred, succ) = get_bounds(g,n+i)
        # tack on initial fixed point data
        fpp = ic_fix + fp
        Q_char = [ X[j-1] for j in fpp[succ] if j not in fpp[n+i] ]
        S_char = [ -X[j-1] for j in fpp[n+i] if j not in fpp[pred] ]
        # First, calculate the denominator
        StQ_char = [ q+s for s in S_char for q in Q_char ]
        denominator = denominator * prod(StQ_char)
    return denominator

def bott_term(g,w,fp,X,lambda_S,lambda_Q):
    """Calculate the localization term at an isolated 
       fixed point
    INPUT: g = incidence graph
           w = node 'weights' (result of ig_weights() call)
           fp = fixed point data for a single fixed point
           X = character of C^* on C^n (a list of 'n' integers)
           lambda_S = list of partitions for s_\lambda(S),one for each node
           lambda_Q = list of partitions for s_\lambda(Q),one for each node
    OUTPUT: an integer
    """
    n = (len(w)+1)/2   # size of the initial chain in g
    ic_fix = ic_fixed_points(w[:n])
    numerator = 1
    denominator = 1
    for i in range(n-1):
        (pred, succ) = get_bounds(g,n+i)
        # tack on initial fixed point data
        fpp = ic_fix + fp
        Q_char = [ X[j-1] for j in fpp[succ] if j not in fpp[n+i] ]
        S_char = [ -X[j-1] for j in fpp[n+i] if j not in fpp[pred] ]
        # First, calculate the denominator
        StQ_char = [ q+s for s in S_char for q in Q_char ]
        # print "S_char=",S_char
        # print "Q_char=",Q_char
        # print "StQ_char=",StQ_char
        denominator = denominator * prod(StQ_char)
        # Second, calculate the numerator
        f1 = evaluate_schur(lambda_S[i],S_char)
        f2 = evaluate_schur(lambda_Q[i],Q_char)
        numerator = numerator * f1 * f2
	return numerator / denominator
	
def bott_term_general(g,w,fp,X,lambda_S,lambda_Q):
    """Calculate the general localization term at an isolated 
       fixed point. This procedure counts only terms of the correct (top) degree
	   and allows for a list of s_\lambda(B) factors in the numerator where
	   the imput list "L" is a list of tuples ( bundle_type , bundle_id, lambda )
    INPUT: g = incidence graph
           w = node 'weights' (result of ig_weights() call)
           fp = fixed point data for a single fixed point
           X = character of C^* on C^n (a list of 'n' integers)
		   L = list of tuples ( bundle_type , bundle_id, lambda ) where:
			   bundle_type = ( 's' | 'q' )
			   bundle_id   = i in { 1, ... , n}; n = number of nodes
			   lambda      = partition
    OUTPUT: a rational number
    """
    n = (len(w)+1)/2   # size of the initial chain in g
    ic_fix = ic_fixed_points(w[:n])
    numerator = 1
    for (bundle_type, i, lam) in L:
        (pred, succ) = get_bounds(g,n+i)
        # tack on initial fixed point data
        fpp = ic_fix + fp
		if bundle_type == 's':
        	B_char = [ -X[j-1] for j in fpp[n+i] if j not in fpp[pred] ]
		else if bundle_type == 'q':
        	B_char = [ X[j-1] for j in fpp[succ] if j not in fpp[n+i] ]
		else:
			raise TypeError, 'bundle_type'
		ev = evaluate_schur(lam,B_char)
        numerator = numerator * ev
	denominator = bott_denominator(g,w,fp,X)
	return numerator / denominator

