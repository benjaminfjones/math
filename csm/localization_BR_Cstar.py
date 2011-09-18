from copy import *

def loc_factor(sp, weights):
	"""Returns the result of evaluating a given symmetric poly at the indicated values
	
	INPUT:
		sp = element of a symmetric function algebra SFA(...)
		weights = list of values to evaluate sp on.
		   
	OUTPUT: 
		A rational number
	"""
	f = sp.expand(len(weights))
	# to deal with SFAElemantary bug 2
	#if len(weights) < len(f.parent().gens()):
	#	return f(weights+[0 for i in range(len(f.parent().gens())-len(weights))])
	return f(weights)
	
def bott_residue_term( splist, wlist, sp_denom, weights_denom):
	"""Returns a term from the Bott residue formula. Each term has the form:
	$$ \frac{ \prod_i c_{\lambda_i}(B_i) }{ c_d( N ) } $$ 
	
	INPUT:
		splist = list of symmetric polys
		wlist  = list of weights lists, one weight list per poly
		sp_denom = symmetric poly for the denominator
		weights_denom = weights for sp_denom
		
	OUTPUT:
		A rational number
	"""
	res = 1
	for i in range(len(splist)):
		res = res * loc_factor( splist[i], wlist[i] )
	return res / loc_factor( sp_denom, weights_denom )

def loc_term(e1,e2,f1,f2):
    """Returns a term from the Bott-Residue formula

    INPUT: e1 = integer weight
           e2 = integer weight
           f1 = list of weights, should contain e1
           f2 = list of weights, should contain e1,e2

    OUTPUT: integer (sum of rational numbers?)
    """
    flag2 = copy(f2)
    flag2.remove(e1)
    flag2.remove(e2)

    flag1 = copy(f1)
    flag1.remove(e1)
    
    denom = prod( [ (f-e1) for f in flag1] ) * prod( [ (f-e2) for f in flag2 ] )
    return (-e1*e2*(e1+e2))/( denom )

def fixed_points_24():
    """Finds the fixed points of T on ..."""
    c4=[1,2,3,4]
    c2=[1,2]
    G24T = combinations(c4,2)
    G12T = combinations(c2,1)
 

#c4weights = [2,1,-1,-2]
#c4weights = [1,2,3,4]
#c4weights = [-2,-1,1,2]

def intZ(w):
	s=SFASchur(QQ)
	e=SFAElementary(QQ)
	c4weights = w
	c2weights = c4weights[0:2]

	slist = []
	for i in c2weights:
		temp_c2=copy(c2weights)
		temp_c2.remove(i)
    		temp_c4 = copy(c4weights)
		temp_c4.remove(i)
    		for j in temp_c4:
			nw=[i,j]
			temp2_c4=copy(temp_c4)
			temp2_c4.remove(j)
			dw=[f-i for f in temp_c2] + \
			   [f-j for f in temp2_c4]
			t=bott_residue_term( [-s([1]),s([1,1])], [ [i], temp2_c4 ], \
					     e([3]), dw )
        		slist.append(t)
			print "i=",i+1,"j=",j+1,"term=",t

	return sum( slist )

intZ([2,1,-1,-2])
