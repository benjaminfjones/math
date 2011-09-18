# demazure.sage
# 1-12-2011
# Benjamin Jones (jonesbe@uwstout.edu)
# file = "/Users/jonesbe/Desktop/svn/sage_code/demazure/demazure.sage"
#

#
# enhancements to the WeightRingElement class
#
def wr_degree(self):
    """
    The degree of the character, that is, the dimension of module.
    
    EXAMPLES::
        
        sage: b = WeightRing(WeylCharacterRing(['B',2]))
        sage: ch = 3*b((3,1)) + 2*b((2,1)) + b((1,1))
        sage: ch.degree()
        6
    """
    return sum(self._mdict[k] for k in self._mdict)
    
def sl3_2d_coords(v):
    """
    INPUT:
    v = 3d vector lying in the plane x+y+z=0

    OUTPUT:
    2d vector of coordinates with respect to basis vectors: 
    u1 = vector((1/2, -1/2, 0))
    u2 = vector((1/2*sqrt(1/2)*sqrt(2/3), 1/2*sqrt(1/2)*sqrt(2/3), 
          -sqrt(1/2)*sqrt(2/3)))
    """
    u1 = vector((1/2, -1/2, 0))
    u2 = vector((1/2*sqrt(1/2)*sqrt(2/3), 1/2*sqrt(1/2)*sqrt(2/3), 
          -sqrt(1/2)*sqrt(2/3)))
    return vector((v*u1, v*u2))

def g2_2d_coords(v):
    """
    INPUT:
    v = 3d vector lying in the plane x+y+z=0

    OUTPUT:
    2d vector of coordinates with respect to basis vectors: 
    u1 = vector((0, 1, -1))
    u2 = vector((1, -1/2, -1/2))
    """
    u1 = vector((0, 1, -1))
    u2 = vector((1, -1/2, -1/2))
    return vector((v*u1, v*u2))


def wr_plot(self, plot_roots=True, arrow_style=None, point_style=None,
            mult_scale=None):
    """
    Plot the character if cartan type is rank 2.
    
    EXAMPLES::
    
        sage: A = WeylCharacterRing(['B',2])
        sage: a = WeightRing(A)
        sage: lam = a.space().fundamental_weights()
        sage: M = demazure(a(4*lam[1]+4*lam[2]), [1,2])
        sage: M.plot()
    """
    # plot parameters
    if arrow_style == None:
        arrow_style={'rgbcolor':(0,0,1),'width':1,'linestyle':'dotted'}
    if point_style == None:
        point_style={'rgbcolor':(1,0,0)}
    if mult_scale == None:
        mult_scale = lambda m: m^2*10

    # Type A2
    if tuple(self.cartan_type()) == ('A',2):
        A2 = WeylCharacterRing(['A',2])
        roots = A2.space().roots()
        
        # compute coordinates of the roots
        root_coords = [ sl3_2d_coords(vector(A2.coerce_to_sl(r)))
                      for r in roots ]
        # compute coordinates of the weights and store with multiplicity
        wt_coords = [ (tuple(sl3_2d_coords(vector(A2.coerce_to_sl(wt)))),m) 
                      for (wt,m) in self.mlist() ]
    # Type B2: alpha[2] is short
    elif tuple(self.cartan_type()) == ('B',2):
        B2 = WeylCharacterRing(['B',2])
        roots = B2.space().roots()
        # Note: we reflect about the line y=x here to be consistent with
        #       the std. picture of B2 where alpha[2] points east and
        #       alpha[1] points northwest
        root_coords = [ vector((r.coefficient(1), r.coefficient(0))) 
                        for r in roots ]
        wt_coords = [ ((wt.coefficient(1), wt.coefficient(0)), m) 
                      for (wt,m) in self.mlist() ]
    # Type G2: alpha[1] is short
    elif tuple(self.cartan_type()) == ('G',2):
        G2 = WeylCharacterRing(['G',2])
        roots = G2.space().roots()
        root_coords = [ g2_2d_coords(vector(map(r.coefficient,range(3))))
                      for r in roots ]
        wt_coords = [ (tuple(g2_2d_coords(vector(map(wt.coefficient,range(3))))), m) 
                      for (wt,m) in self.mlist() ]
    else:
        raise NotImplementedError('Plotting is not implemented for this type')
        
    root_plot = sum(plot(v,**arrow_style) for v in root_coords)
    wt_plot = sum([point(wt, size=mult_scale(m), **point_style) 
                   for (wt,m) in wt_coords])
    if plot_roots:
        return wt_plot+root_plot
    else:
        return wt_plot

# add methods to the existing class WeightRingElement
WeightRingElement.degree = wr_degree
WeightRingElement.plot = wr_plot

# 
# Demazure operator
#
def demazure(ch, i):
    """
    EXAMPLES:
        # Check the Demazure character formula in A2
        sage: A = WeylCharacterRing(['A',2])
        sage: a = WeightRing(A)
        sage: space = a.space()
        sage: lam = space.fundamental_weights()
        sage: rho = sum(list(lam))
        sage: wch = A(2*rho) # character of the irred. representation
        sage: ch = a(2*rho)  # element in the weight ring
        sage: dch = demazure(ch,[2,1,2]) # apply 3 Demazure operators
        sage: sorted(wch.mlist()) == sorted(dch.mlist())
        True
        
        # Check the Demazure character formula in B2
        sage: A = WeylCharacterRing(['B',2])
        sage: a = WeightRing(A)
        sage: space = a.space()
        sage: lam = space.fundamental_weights()
        sage: rho = sum(list(lam))
        sage: wch = A(2*rho) # character of the irred. representation
        sage: ch = a(2*rho)  # element in the weight ring
        sage: dch = demazure(ch,[1,2,1,2]) # apply 4 Demazure operators
        sage: sum([m for (wt,m) in dch.mlist()]) # degree should be 81
        81
        sage: sorted(wch.mlist()) == sorted(dch.mlist())
        True
        
    """
    a = ch.parent()
    # passing a list applies the operators in the indicated order
    # elements of i index simple roots
    if type(i) == list:
        ch_temp = copy(ch)
        for j in i:
            ch_temp = demazure(ch_temp, j)
        return ch_temp
    # otherwise i is an index for a simple root
    space = a.space()
    if i not in space.index_set(): 
        raise IndexError('%s does not index a simple root' % i.__str__())
    alpha = space.simple_root(i) # if i is in space.index_set()
    ch_temp=a(0)
    for (wt,wtm) in ch.mlist():
        n = wt.scalar(2*alpha/alpha.scalar(alpha))
        if n >= 0:
            ch_temp += wtm * sum([a(wt-i*alpha) for i in range(n+1)])
        if n == -1: # for clarity
            ch_temp += 0
        if n <= -2:
            ch_temp += wtm * sum([-1*a(wt+i*alpha) for i in range(1,-n)])
    return ch_temp

def weight_string(M,wt,alpha):
    """
    INPUT:
    M = element of WeightRing(...)
    wt = a weight from M of multiplicity at least 1
    alpha = positive root of the corresponding RootSystem, 
    i.e. M.parent().positive_roots()
    
    OUTPUT:
    Returns the alpha-string of weights in M through wt. This is
    a list of the form [ (wt1, m1), (wt2, m2), ... ] where wt's are weights
    and m's are multiplcities. 
    """
    # sanity check
    D = M._mdict
    if wt not in D.keys():
        raise ValueError('Specified weight is not in the given ' +
                         'WeightRingElement')
    # find the upper (n) and lower (m) bounds of the string
    n=-1
    wt_temp = wt
    while wt_temp in D.keys():
        wt_temp += alpha
        n += 1
    m=-1
    wt_temp = wt
    while wt_temp in D.keys():
        m = m+1
        wt_temp -= alpha
    return [ (wt + i*alpha, D[wt + i*alpha]) for i in range(-m,n+1) ]
    
def c_kernel(M,alpha):
    """
    INPUT:
    M = element of WeightRing(...)
    alpha = positive root of the corresponding RootSystem, 
    i.e. M.parent().positive_roots()
    
    OUTPUT:
    The number of weights in M (with multiplicity) that are anhilated by alpha
    """
    return sum([m for (wt,m) in M.mlist() if wt+alpha not in M._mdict.keys()])
    
 
#
# __main__
#
CT = ['B',2]
A = WeylCharacterRing(CT); A2=A
a = WeightRing(A)
space = a.space()
lam = space.fundamental_weights()
M = demazure(a(4*lam[1]+4*lam[2]),[1,2])
M.plot()



