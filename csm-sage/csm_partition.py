#
# csm_partition.py
# Author: Benjamin F Jones
# Version: 7-13-08
#
#################################
# Partition funcitons
#################################

from copy import copy
from sage.misc.functional import det
from sage.rings.arith import binomial
from sage.combinat.partition import *
from sage.matrix.constructor import *

def dcoeff(l,m,N):	
    """Returns the 'D' coefficient for partitions l, m having exactly N parts. See 
	   \cite{Fulton:Intersection_Theory}"""
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
    if len(lam) > l: raise IndexError, "lam is longer than l"
    dl = [ r-i for i in lam ]
    dl.extend([r]*(l-len(lam)))
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
				
#
# Partition functions specific to Zelevinsky notation and resolutions
#
				
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