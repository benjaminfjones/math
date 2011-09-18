#
# bott_residue.py -- Funcitons for computing with the Bott Residue 
#                    fomrula and specifically for Z resolutions
# Author: benjamin F Jones
# Version: 7-15-08
#
from csm_partition import *
from Z_resolution import *
from sage.combinat.sf.sfa import SFASchur

class ChernClass:
    """
    Used as a C-like struct for holding 4 pieces of data which
    determine a certain type of Chern class.
    DATA:
        coeff       = integer coefficient
        bundle_type = either 's' = sub, 'q' = quotient, or 'e' = trivial
        bundle_id   = between 1 and Z.num_nodes (Z is a Zresolution object)
        part        = partition
    """
    def __init__(self, c=1,t='e',i=0,p=[]):
        self.coeff = c      
        self.bundle_type = t 
        self.bundle_id = i   
        self.part = p
        
    def __str__(self):
        return '%d * s_%s(%s[%d])' % (self.coeff,str(self.part),str(self.bundle_type.upper()),\
                                  self.bundle_id)
    def degree(self):
        """Returns the degree of the Chern class."""
        return sum( self.part )
        
class ChernClassMonomial(list):
    """
    Keeps track of a product of Chern classes
    Methods: __init__, __str__, and append are overridden
    INIT: a list of ChernClass objects. The coeff's get multiplied into the 
          'overall' coeff which is an attribute of this class
    """
    def __init__(self, data=[]):
        list.__init__(self)
        self.extend(data)
        self.coeff = 1
        for d in self:
            self.coeff *= d.coeff
            
    def __str__(self):
        s = "%d * ( " % self.coeff
        for d in self:
            if trunc_part(d.part) != []:
                s += "s_%s(%s[%d]) * " % (str(d.part),str(d.bundle_type.upper()), d.bundle_id)
        s += ")"
        return s
    def _latex_(self):
        """
        Return \LaTeX representation of X.

        EXAMPLES:
            sage: a = X(1,2)
            sage: latex(a)
            '\\frac{1}{2}'
        """
        if self.coeff != 1:
            s = "%s " % self.coeff
        else:
            s = ""
        for d in self:
            if trunc_part(d.part) != []:
                s += "s_{%s}(\underline{%s}_{%s}) " % (d.part,d.bundle_type.upper(), d.bundle_id)
        return s
    def append(self, item):
        """Appends an item and multiplies the overall coefficient by item.coeff"""
        if trunc_part(item.part) != []:  # ignore factors of 1
            list.append(self, item)
            self.coeff *= item.coeff
    def degree(self):
        """Rerturns the total degree of the monomial, i.e. the sum of the 
        sizes of all partitions"""
        return sum( [ d.degree() for d in self ] )

#
# Now in class Zresolution
#            
# def tensor_term_iterator(Z):
#     """Returns a coefficient, a list of S partitions, and a list of Q
#     partitions which parametrize one term in the chern class of the
#     tensor product of the bundles specified by the input data.
#     INPUT: Z = instance of class Zresolution
#     OUTPUT: a list of tuples ( coeff, bundle_type, bundle_id, lam ) where:
#             coeff = coefficient of the term (an integer)
#             bundle_type = ( 's' | 'q' ), s = sub, q = quotient
#             bundle_id   = the 'i' in S_i or Q_i
#             lam         = the lambda in s_lambda(B)
#     """
#     l = []
#     r = []
#     FL_part = []
#     n = Z.num_nodes
#     for i in range(n):
#         (pred, succ) = Z.get_bounds(i)
#         node_weight = Z.weights[Z.ic_size+i]
#         pred_weight = Z.weights[pred]
#         succ_weight = Z.weights[succ]
#         l.append(node_weight - pred_weight)
#         r.append(succ_weight - node_weight)
#         # generate 2 part flags in the rectangle l_i x r_i
#         FL_part.append(flagged_partition_iterator(l[-1],r[-1]))
# 
#     CPI = cartesian_product_iterator( [ list(a) for a in FL_part ] )
#     for c in CPI:
#         M = ChernClassMonomial()
#         for i in range(n):
#             (lam,mu) = c[i]
#             M.append( ChernClass( dcoeff(lam,mu,l[i]), 's', i+1, mu ) )
#             M.append( ChernClass( 1, 'q', i+1, conj( dual(lam,l[i],r[i]) )))
#         yield M

# def evaluate_schur(lamb,X):
#     """Evaluates s_\lambda(B) where B has C^* character X"""
#     # print "evaluate_schur(",lamb,",",X,")"
#     s=SFASchur(QQ)
#     if trunc_part(lamb) == []:
#         return 1
#     else:
#         f=s(lamb)
#         fp=f.expand(len(X))
#         return fp(X)

def evaluate_schur(Z,lamb,X):
    """Evaluates s_\lambda(B) where B has C^* character X

	TODO: This should probably be a class method of Z
	"""
    t=(tuple(trunc_part(lamb)),len(X))
    if t in Z.schur_expand:
        return Z.schur_expand[ t ](X)
    else:
        raise KeyError, "invalid schur function"
    
def bott_denominator(Z,fp,X):
    """Calculate the top Chern class of the normal bundle to an isolated 
       fixed point

    INPUT: Z = Zresolution object
           fp = fixed point data for a single fixed point
           X = integral character of C^* on C^n (a list of 'n' integers)

    OUTPUT: an integer

	TODO: This should probably be a class method of Z	
    """
    n = Z.num_nodes
    ic_fix = Z.ic_fix
    denominator = 1

    for i in range(n):
         # get succ and pred
        (pred, succ) = Z.get_bounds(i)
        # tack on initial fixed point data
        fpp = ic_fix + fp
        Q_char = [ X[j-1] for j in fpp[succ] if j not in fpp[Z.ic_size + i] ]
        S_char = [ -X[j-1] for j in fpp[Z.ic_size + i] if j not in fpp[pred] ]
        # First, calculate the denominator
        StQ_char = [ q+s for s in S_char for q in Q_char ]
        denominator *= prod(StQ_char)
        
    return denominator

# !!! needs to be updated
# def bott_term(g,w,fp,X,lambda_S,lambda_Q):
#     """Calculate the localization term at an isolated 
#        fixed point
#     INPUT: g = incidence graph
#            w = node 'weights' (result of ig_weights() call)
#            fp = fixed point data for a single fixed point
#            X = character of C^* on C^n (a list of 'n' integers)
#            lambda_S = list of partitions for s_\lambda(S),one for each node
#            lambda_Q = list of partitions for s_\lambda(Q),one for each node
#     OUTPUT: an integer
#     """
#     n = (len(w)+1)/2   # size of the initial chain in g
#     ic_fix = ic_fixed_points(w[:n])
#     numerator = 1
#     denominator = 1
#     
#     for i in range(n-1):
#         (pred, succ) = get_bounds(g,n+i)
#         # tack on initial fixed point data
#         fpp = ic_fix + fp
#         Q_char = [ X[j-1] for j in fpp[succ] if j not in fpp[n+i] ]
#         S_char = [ -X[j-1] for j in fpp[n+i] if j not in fpp[pred] ]
#         # First, calculate the denominator
#         StQ_char = [ q+s for s in S_char for q in Q_char ]
#         # print "S_char=",S_char
#         # print "Q_char=",Q_char
#         # print "StQ_char=",StQ_char
#         denominator = denominator * prod(StQ_char)
#         # Second, calculate the numerator
#         f1 = evaluate_schur(lambda_S[i],S_char)
#         f2 = evaluate_schur(lambda_Q[i],Q_char)
#         numerator = numerator * f1 * f2
#         
#     return numerator / denominator
    
def bott_term_general(Z,fp,X,M):
    """Calculate the general localization term at an isolated 
       fixed point. This procedure counts only terms of the top degree (other
       terms integrate to zero).
    INPUT: Z = Zresolution object
           fp = a fixed point datum (e.g. from Z.fixed_points() )
           X = character of C^* on C^n (a list of 'n' integers)
           M = Chern class monomial object
    OUTPUT: a rational number

	TODO: This should probably be a class method of Z
    """
    # check that the dimension == degree
    if Z.dimension != M.degree():
        return 0
        
    n = Z.num_nodes
    ic_fix = Z.ic_fix
    numerator = 1

    for c in M:
        fpp = ic_fix + fp
        if c.bundle_type == 'e':
            # 'e' means the quotient bundle on the ambient grassmannian of the image
            # of the resolution
            (pred,succ) = (0,Z.ic_size-1)
            B_char = [ X[j-1] for j in fpp[succ] if j not in fpp[-1] ]
            #print "case 'e': pred=",pred,"succ=",succ,"B_char=",B_char
        elif c.bundle_type == 's':
            (pred, succ) = Z.get_bounds(c.bundle_id - 1)
            B_char = [ -X[j-1] for j in fpp[Z.ic_size + c.bundle_id - 1] if j not in fpp[pred] ]
        elif c.bundle_type == 'q':
            (pred, succ) = Z.get_bounds(c.bundle_id - 1)
            B_char = [ X[j-1] for j in fpp[succ] if j not in fpp[Z.ic_size + c.bundle_id - 1] ]
        else:
            raise TypeError, 'c.bundle_type'
        # ev = evaluate_schur( c.part, B_char )
        ev = evaluate_schur( Z, c.part, B_char )
        numerator *= ev
        
    denominator = bott_denominator(Z,fp,X)
    return M.coeff * numerator / denominator

# def csm_coeff(Z,part):
# 	"""
# 	INPUT:
# 	OUTPUT:
# 	TODO: This should probably be a class method of Z
# 	"""
# 	if trunc_part(Z.part) == []:
# 		if trunc_part(part) == []:
#             return 1
#         else:
#             return 0
#     N = len(Z.part) + Z.part[0]
#     X = list(range(N))
#     c = ChernClass(1,'e',0,conj(part))
#     s = 0
#     for t in Z.tensor_term_iterator():
#         t.append(c)
#         for fp in Z.fixed_points_iterator():
#             s += bott_term_general(Z,fp,X,t)
#     return s