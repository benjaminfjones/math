R.<x,y,z,w> = PolynomialRing(QQ)
Q = R.fraction_field()

def loc_term(e1,e2,glist):
    """Returns a term from the Bott-Residue formula

    INPUT: e1 = generator of R
           e2 = generator of R
           glist = R.gens

    OUTPUT: element of Q
    """
    nlist = glist
    nlist.remove(e1)
    nlist.remove(e2)
    denom = (e2 - e1) * prod( [ (e2 - f) for f in nlist ] )
    return (e1*e2*(e1+e2))/( denom )

s1 = sum([ loc_term(x,t,[x,y,z,w]) for t in [y,z,w] ])
s2 = sum([ loc_term(y,t,[x,y,z,w]) for t in [x,z,w] ])
s = s1 + s2; s
