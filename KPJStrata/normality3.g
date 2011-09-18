###########################################################################
#
# GAP3/4 code for studying the normality of enhanced quiver varieties.
#
# The goal is to enumerate the Kraft-Procesi-Johnson strata of enhanced 
# quiver varieties (they have a combinatorial parametrization) and then
# compute a combinatorial bound on their dimensions.
#
# Conventions:
#
# A "partition" is always a list of integers in nonincreasing order, but
#   trailing 0s are allowed.
# A "quasipartition" is a list of integers in which a later entry may exceed
#   an earlier entry by at most 1.
# A "quasibipartition" is a pair of quasipartitions with the same length,
#   whose termwise sum is a partition.  Note that the same length requirement
#   often forces one member of the pair to have trailing 0s.  It is allowed
#   for both members to have trailing 0s.
# A "signed quasibipartition" is a list of length 3, say [a,b,c], such that
#   [a,b] is a quasibipartition, and c is list of the same length as a and b 
#   whose members are all +1s and -1s.  Moreover, the entries of c obey the
#   following rules:
#   - if a[i] > 0, then c[i] = -(-1)^(a[i])
#   - if a[i] = a[j] = 0 and b[i] = b[j] and
# A "stratum" is list of signed quasibipartitions.  They satisfy conditions 
#   determined by a pair (r_i),I.  The r_i determine the total sizes of all the
#   signed quasibipartitions in the list.  Each signed quasibipartition in the
#   list is required to have a subordinate bipartition (either + or -) that
#   agrees with one of the subordinate bipartitions of the next term in the
#   list; the specific choices of + and - are determined by I.  HOWEVER, the
#   (r_i),I are not carried along as part of the data structure of a stratum.
#
# The first part of the file is fairly complete program for testing the 
# dimension conjecture by enumerating all strata of all enhanced quiver
# varieties associated to bipartitions of a given total size.  The main
# function is:
#
#   TestAllQuivers( n );
#
# However, enumerating all strata quickly becomes far too slow to do.  The
# dimension estimate depends only on the sequence of subordinate bipartitions,
# so a more efficient approach is to enumerate only all possible sequences
# of subordinate bipartitions.  There are far fewer of those than there are
# strata.  See the comments in the farther down in the file.
#
# Read("/Users/jonesbe/Desktop/svn/sage_code/KPJStrata/normality3.g");
#
###########################################################################

###########################################################################
#
# For use in GAP3, uncomment the following lines:
#StructuralCopy := Copy;
#fail := false;
#
###########################################################################

KPJ_SKIPTRIVIAL := true;
#
# If set to true, the various "TestAll..." functions will skip quivers in
# which I = [0 .. k] for some k, since those quivers correspond to enhanced
# nilpotent orbits already known to have normal closures.
#
KPJ_SKIPCOTRIVIAL := true;
#
# If set to true, skips those with I = [k .. t-1].
#

PartN := function( p )
#
# Returns n(p), where p is a partition.
#
  return Sum(List([1..Length(p)], i -> (i-1)*p[i]));
end;

TransposedPart := function( p )
#
# Returns the transpose of the partition p.
#
  return List([1..Maximum(p)], i -> Number(p, j -> j >= i));
end;

QuasiIncreasing := function( p )
#
# Returns 'true' if p is a quasipartition.
#
  local l;
  l := Length(p);
  return ForAll([1..l-1], i -> ForAll([i+1..l], j -> p[i] >= p[j]-1));
end;

AllQuasibipartitions := function( p )
#
# Assumes p is a partition.
# Returns the list of all quasibipartitions whose sum is the partition p.
#
  local g;
  g := List([0..p[1]], a -> [a, p[1]-a]);
  if Length(p) = 1 then
    return  List([0..p[1]], a -> [[a], [p[1]-a]]);
  else
    return Filtered( List( 
      Cartesian(g, AllQuasibipartitions( p{[2..Length(p)]} )),
        w -> [Concatenation([w[1][1]], w[2][1]),
          Concatenation([w[1][2]], w[2][2])] ),
      b -> QuasiIncreasing(b[1]) and QuasiIncreasing(b[2]) );
  fi;
end;

AllSignedQuasibipartitions := function( p )
#
# Assumes p is a partition.
# Returns the list of all signed quasibipartitions whose sum is p.
#
  local qb, a, x, y, c, d, e, r;
  qb := AllQuasibipartitions(p);
  a := [];
  for x in qb do
# for an element x, think: x[1] is \mu, and x[2] is \nu
    c := Filtered([1..Length(x[1])], i -> x[1][i] = 0);
# c is the list of row numbers where \mu is empty.  These are the only rows
# where there is a choice of sign.
    d := Collected(x[2]{c});
# d is the list of rows lengths and multiplicities in \nu corresponding to
# empty rows in \mu.  For each entry [ length, mult ] in d, we have to put
# anywhere from 0 to mult +1s in the corresponding entries in the list of signs,
# and fill in the rest with -1s.
    e := List(d, t -> Filtered(c, i -> x[2][i] = t[1]));
# e is a list of lists of row numbers of \nu, grouped by row length.
    r := List(Cartesian(List(e, l -> List([0..Length(l)], j -> l{[1..j]}))),
           t -> Concatenation(t));
# r is the list whose members are all possible lists formed in this way: take
# an initial segment (possibly empty) of each member of e, and concatenate them.
# Each member of r corresponds to one possible choice of sign assignment to the
# (unsigned) quasibipartition x.
    for y in r do
      Add(a, [ x[1], x[2], List([1..Length(x[1])], function(i)
        if x[1][i] > 0 then return -(-1)^x[1][i];
        elif i in y then return +1;  else return -1;  fi;  end) ]);
    od;
  od;
  return a;
end;

Bipartitionize := function( bp )
#
# Assumes bp is a pair of quasipartitions.  The sum of its members need not be
#  a partition.
# Returns the unique bipartition obtained by (1) sorting the rows so that the
#  sum is indeed a partition, and then (2) moving individual rows left or right
#  as necessary so that the members becomes honest partitions.
# This function is used to pick off the subordinate bipartitions of a signed
#  quasibipartition.
#
  local mu, nu, j, a;
  mu := StructuralCopy(bp[1]);
  nu := StructuralCopy(bp[2]);

# First sort the rows
  j := PositionProperty([1..Length(mu)-1], i -> mu[i]+nu[i] < mu[i+1]+nu[i+1]);
  while j <> fail do
    a := mu[j];  mu[j] := mu[j+1];  mu[j+1] := a;
    a := nu[j];  nu[j] := nu[j+1];  nu[j+1] := a;
    j := PositionProperty([1..Length(mu)-1],
      i -> mu[i]+nu[i] < mu[i+1]+nu[i+1]);
  od;
#1/0;
  j := PositionProperty([1..Length(mu)-1],
    i -> mu[i] < mu[i+1] or nu[i] < nu[i+1]);
  while j <> fail do
    if mu[j] < mu[j+1] then
      nu[j] := mu[j] + nu[j] - mu[j+1];  mu[j] := mu[j+1];
    else
      mu[j+1] := mu[j+1] + nu[j+1] - nu[j];  nu[j+1] := nu[j];
    fi;
    j := PositionProperty([1..Length(mu)-1],
      i -> mu[i] < mu[i+1] or nu[i] < nu[i+1]);
  od;
# NEW Feb 2010: remove all-0 rows from mu and nu
  j := PositionProperty([1..Length(mu)], i -> mu[i] = 0 and nu[i] = 0);
  if j <> fail then
    mu := mu{[1..j-1]};  nu := nu{[1..j-1]};
  fi;
  return [ mu, nu ];
end;

SignedQbpSubord := function( sqp )
#
# Assumes sqp is a signed quasibipartition.
# Returns a pair consisting of the "+" and "-" subordinate bipartitions.
#
  local sq, j, m, lp, lm, mup, nup, mum, num;

  sq := StructuralCopy(sqp);
  m := Length(sq[1]);
  for j in [1..m] do 
      if sq[1][j] > 0 and sq[3][j]<>-(-1)^sq[1][j] then
        sq[1][j] := sq[1][j] - 1; sq[2][j] := sq[2][j] + 1;
      fi; 
  od;
  lp := List([1..m], i -> Int((sq[1][i] + sq[2][i] + 1/2 + sq[3][i]/2)/2) );
  lm := List([1..m], i -> Int((sq[1][i] + sq[2][i] + 1/2 - sq[3][i]/2)/2) );
  mup := List([1..m], i -> Int((sq[1][i]+1)/2) );
  nup := lp - mup;
  mum := List([1..m], i -> Int( sq[1][i] / 2 ) );
  num := lm - mum;
  return [ Bipartitionize([mup, nup]), Bipartitionize([mum, num]) ];
end;

# Quasibipartitions := function( p, n )
# # THIS FUNCTION IS NOT USED ELSEWHERE AND CAN BE DELETED.
# # p is a partition
# # n is an integer such that the "mu"s (left-hand partitions) have size n
#   local g, m;
#   g := List([0..p[1]], a -> [a, p[1]-a]);
#   m := Sum(p) - n;
#   if m < 0 then  return [ ];
#   elif m = 0 then  return [ [ p, List(p, i -> 0) ] ];
#   elif n = 0 then  return [ [ List(p, i -> 0), p ] ];
#   elif Length(p) = 1 then  return [ [ n ], [ m ] ];
#   else
#     return Concatenation(List([0..n], a ->
#       Filtered(List( Quasibipartitions( p{[2..Length(p)]}, n - a ), bp ->
#       [ Concatenation([a],bp[1]), Concatenation([p[1]-a],bp[2]) ] ),
#       b -> QuasiIncreasing(b[1]) and QuasiIncreasing(b[2])) ));
#   fi;
# end;

SignedQuasibipartitions := function( p, n )
#
# Returns the list of signed quasibipartitions whose sum is n and whose "+"
# subordinate bipartition has total size n.
# 
  return Filtered(AllSignedQuasibipartitions(p), function(bp)
    local t;  t := SignedQbpSubord(bp);  return Sum(t[1][1])+Sum(t[1][2])=n;
    end);
end;

# CheckSubord := function( sbp )
# #
# # I DONT KNOW WHAT THIS DOES.  IT DOES NOT MATTER, SINCE IT IS NOT USED.
# #
#   local lp, lm, r;
#   lp := sbp[1][1] + sbp[1][2];
#   lm := sbp[2][1] + sbp[2][2];
#   r := Sum(lp) - Sum(lm);
#   if r <=0 or Number(lp, i -> i > 0) < r then  return 0;  fi;
#   return Sum(sbp[1][1]) - Sum(sbp[2][1])
#     <= r + PartN(lp) - PartN(lm) - r*(r-1)/2;
# end;
# 
# CheckAll := function( N, k )
# #
# # I DONT KNOW WHAT THIS DOES.  IT DOES NOT MATTER, SINCE IT IS NOT USED.
# #
#   local P, S, SB;
#   P := Partitions(N);
#   S := Concatenation(List(P, SignedQuasibipartitions));
#   SB := List(S, SignedQbpSubord);
#   if k = 0 then return List(SB, CheckSubord); else return S[k]; fi;
# end;
# 
# Choose := function( s, k )
# # THIS FUNCTION IS NOT USED ELSEWHERE AND CAN BE DELETED.
# # returns the set of sublists of length k from the list s
# #
#   local l;
#   l := Length(s);
#   if k = 0 then  return [ [ ] ];
#   elif k = 1 then  return List(s, a -> [a]);
#   elif k = s then  return [ s ];
#   elif k < 0 or k > s then  return [ ];
#   else
#     return Concatenation(List([1..l-k], i -> List(Choose(s{[i+1..l]}, k-1),
#       u -> Concatenation([s[i]], u) )) );
#   fi;
# end;

EqPart := function( p, q )
#
# Assumes p and q are both partitions.
# Returns 'true' if they are equal partitions.  (They are allowed to have
#   different numbers of trailing 0s.)
#
  if Length(p) = Length(q) then
    return p = q;
  elif Length(p) > Length(q) then
    return p{[1..Length(q)]} = q and
      ForAll(p{[Length(q)+1..Length(p)]}, x -> x = 0);
  else
    return q{[1..Length(p)]} = p and
      ForAll(q{[Length(p)+1..Length(q)]}, x -> x = 0);
  fi;
end;

EqBipart := function( p, q )
#
# Assumes p and q are both bipartitions.  Returns 'true' if they are equal.
#
  return EqPart(p[1],q[1]) and EqPart(p[2],q[2]);
end;

KPJStrata := function( rl, I )
#
# Returns the list of all strata associated to the data (r_i),I
#
# Warning: SignedQuasibipartions and SignedQbpSubord work with +s and -s.
# Must be careful about which goes with which, using I.
#
  local t, kpj, bp, left, tail, b, ind, N, j, k, a;

  t := Length(rl);
  if t = 0 then  return [ [ ] ];  fi;
  kpj := KPJStrata( rl{[1..t-1]}, Filtered(I, i -> i < t-1) );
# kpj is the list of strata for a quiver with one less vertex
  if t - 1 in I then
    bp := Concatenation( List(Partitions( 2*Sum(rl) - rl[t] ),
      p -> SignedQuasibipartitions(p, Sum(rl)) ) );
    left := List(bp, s -> SignedQbpSubord(s)[2]);
  else
    bp := Concatenation( List(Partitions( 2*Sum(rl) - rl[t] ),
      p -> SignedQuasibipartitions(p, Sum(rl) - rl[t]) ) );
    left := List(bp, s -> SignedQbpSubord(s)[1]);
  fi;
# bp is the list of all possible signed quasibipartitions associated to the
# last edge in the quiver, and left is the list of subordinate bipartitions
# attached to the penultimate vertex.
  if t = 1 then  return List(bp, p -> [p]);  fi;
  N := Length(left);
  a := [];
  if t - 2 in I then  tail := 1;  else  tail := 2;  fi;
  for j in kpj do
# for each stratum in kpj, tack on all possible signed quasibipartitions
# with matching subordinate
    b := SignedQbpSubord( j[t-1] )[tail];
    ind := Filtered([1..N], i -> EqBipart(b, left[i]));
    for k in ind do
      Add(a, Concatenation(j, [bp[k]]));
    od;
  od;
  return a;
end;    

KPJSubord := function( kp, I, n )
#
# Assumes that kp is a stratum, I is part of the (rl,I) data, and n is the
#   index of some vertex in the quiver, in the range 0 .. Length(kp).
# Returns the bipartition attached to the nth vertex.  This bipartition is 
#   subordinate to both the nth and (n+1)th terms signed quasibipartitions in
#   kp.
#
  local t, a, b;
  t := Length(kp);
  if n > 0 then
    if n-1 in I then  a := SignedQbpSubord(kp[n])[1];
    else  a:= SignedQbpSubord(kp[n])[2];
    fi;
  fi;
  if n < t then
    if n in I then  b := SignedQbpSubord(kp[n+1])[2];
    else  b := SignedQbpSubord(kp[n+1])[1];
    fi;
  fi;
  if n = 0 then  return b;
  elif n = t then  return a;
  else
    if not EqBipart(a,b) then Error("Error!"); else return a; fi;
  fi;
end;

KPDelta := function( sbp )
#
# Assumes sbp is a signed quasibipartition.
# Forgets to a signed bipartition, then computes the number "Delta"
#   from [KP] Proposition 5.3.
#
  local p, n, a, b, i, k;
  p := Length(sbp[1]);
  n := sbp[1][1] + sbp[2][1];  # max row length
  if n mod 2 = 0 then  n := n + 1;  fi;
  a := List([1..n], i -> 0);
  b := StructuralCopy(a);
  for i in [1..p] do
    k := sbp[1][i] + sbp[2][i];
    if sbp[3][i] = 1 then a[k] := a[k] + 1;
      else  b[k] := b[k] + 1;  fi;
  od;
  return Sum(List([1,3..n], i -> a[i]*b[i]));
end;

KPJEstdim := function( rl, I, kp )
#
# Assumes kp is a stratum corresponding to the data (rl,I).
# Returns the combinatorial estimate of the dimension of that stratum, minus
#  the expected dimension of the entire variety.  Thus, if our conjectures
#  are true, this function should always return nonpositive numbers.
#
  local t, A, B, fp;
  t := Length(kp);
  A := Filtered([0..t-1], i -> (i in I) and not (i-1 in I));
  B := Filtered([1..t], i -> (i-1 in I) and not (i in I));
  fp := KPJSubord(kp,I,t);
  return PartN(TransposedPart(rl)) - PartN(fp[1]+fp[2])
    + Sum(List(B, i -> Sum(KPJSubord(kp,I,i)[1])))
    - Sum(List(A, i -> Sum(KPJSubord(kp,I,i)[1])))
    - Sum(List(I, i -> rl[i+1]))
    - Sum(List(kp, x -> KPDelta(x)));
end;

ENOrbitdim := function( bp )
#
# Assumes bp is a bipartition.
# Returns the dimension of the corresponding enhanced nilpotent orbit.
#
  local e, f;
  e := Sum(bp[1]);
  f := Sum(bp[2]);
  return (e+f)^2 - (e+f) - 2* PartN(bp[1]+bp[2]) + e;
end;
  
ENPOrbitdim := function( sbp )
#
# Assumes sbp is a signed quasibipartition.
# Returns the dimension of the corresponding enhanced nilpotent orbit,
#  using Johnsons formula.
#
  local n, bp, m;
  n := Sum(sbp[1]) + Sum(sbp[2]);
  bp := Bipartitionize([sbp[1], sbp[2]]);
  m := Sum(List([1..Length(sbp[1])], i ->
     Int(sbp[1][i]/2 + 1/4 + sbp[3][i]/4) ));
  return n*(n-1)/2 - PartN(bp[1]+bp[2]) + m;
end;

KPJActualdim := function( rl, I, kp )
#
# Assumes kp is a stratum corresponding to the data (rl,I).
# Returns the actual dimension of that stratum.
#
  return Sum(List(kp, s -> ENPOrbitdim(s)))
    - Sum(List([1..Length(kp)-1], i -> ENOrbitdim(KPJSubord(kp,I,i))));
end;

KPJSmooth := function( rl, I, kp )
#
# Assumes kp is a stratum corresponding to the data (rl,I).
# Returns 'true' if all 'A_i' are injective, all 'B_i' are surjective,
#  and if i is in I then u_i is not in the image of A_i
#
  local t, sq, p, a;
  t := Length(kp);
  sq := kp[t];
  p := Length(sq[1]);
  if t-1 in I then
    a := (ForAll(sq[3], x -> x = 1) # A_{t-1} is injective
      and ForAny([1..p], i-> sq[2][i] = 0 and sq[1][i] > 0)) # u_i not in image
      or ForAll([1..p], i -> sq[3][i] = 
      -(-1)^(sq[1][i] + sq[2][i])); # B_{t-1} is surjective
  else
    a := ForAll(sq[3], x -> x = -1) # A_{t-1} is injective
      or ForAll([1..p], i -> sq[3][i] = 
      (-1)^(sq[1][i] + sq[2][i])); # B_{t-1} is surjective
  fi;
  if t = 1 then
    return a;
  else
    return a and KPJSmooth( rl{[1..t-1]}, I, kp{[1..t-1]});
  fi;
end;

KPJEstdimlist := function( rl, I )
#
# Returns the list of all dimension estimates (minus expected dimension of the
# variety) for all strata associated to the data (rl,I).
#
  return List(KPJStrata(rl, I), k -> KPJEstdim(rl, I, k));
end;

AllQuivers := function( n )
#
# Returns the list of all possible pairs (rl,I) where the sum of the entries 
# of rl is n.
#
  local ptn, bip, i, p, a, rl, I, q;
  ptn := List([0..n], i -> Partitions(i));
  bip := [ ];
  for i in [0..n] do
    for p in ptn[i+1] do
      Append(bip, List(ptn[n-i+1], q -> [p,q]));
    od;
  od;
  a := [ ];
  for p in bip do
    rl := Concatenation(p[1], p[2]);
    Sort(rl);
    I := [ ];
    for q in p[1] do
      Add(I, First([0..Length(rl)-1], i -> not(i in I) and rl[i+1] = q));
    od;
    Add(a, [rl, I]);
  od;
  return a;
end;

TestAllQuivers := function( n )
#
# Main testing function.  Generates all quivers of total size n, and for each,
# enumerates all strata, and estimates their dimension.  If any such dimension
# estimate comes out positive, it reports the error.
#
  local aq, i, kp, m, j, dim;
  aq := AllQuivers(n);
  Print(Length(aq)," quivers\n");
  for i in [1..Length(aq)] do
    Sort(aq[i][2]);
    if KPJ_SKIPTRIVIAL and (Length(aq[i][2]) = 0 or 
      aq[i][2] = [0.. Length(aq[i][2])-1] or (KPJ_SKIPCOTRIVIAL and
      aq[i][2] = [Length(aq[i][1])-Length(aq[i][2])..Length(aq[i][1])-1])) then
      Print("Q", i, ": trivial\n");
    else
    kp := KPJStrata(aq[i][1],aq[i][2]);
    m := Length(kp);
    for j in [1..m] do
      if j mod 500 = 0 then Print("Q", i, ": ", j, " of ", m, " strata\r"); fi;
      dim := KPJEstdim(aq[i][1], aq[i][2], kp[j]);
      if dim > 0 then
        Error("\nFalse: Quiver ", aq[i], "\n Stratum ", kp[j], "\n");
      elif dim >= -1 and not KPJSmooth(aq[i][1],aq[i][2],kp[j]) then
        Error("\nNot smooth: Quiver ", aq[i], "\n Stratum ", kp[j], "\n");
      fi;
    od;
    Print("Q", i, ": ", m, " strata\n");
  fi; od;
#  return
#    ForAll( AllQuivers(n), q -> ForAll(KPJEstdimlist(q[1],q[2]), i -> i<=0) );
end;

TestAllQuiversActual := function( n )
#
# Like TestAllQuivers, but uses actual dimensions of strata instead of the
# estimate.
#
  local aq, i, rl, I, t, dri, kp, m, j, dim;
  aq := AllQuivers(n);
  Print(Length(aq)," quivers\n");
  for i in [1..Length(aq)] do
    Sort(aq[i][2]);
    if KPJ_SKIPTRIVIAL and (Length(aq[i][2]) = 0 or 
      aq[i][2] = [0.. Length(aq[i][2])-1] or (KPJ_SKIPCOTRIVIAL and
      aq[i][2] = [Length(aq[i][1])-Length(aq[i][2])..Length(aq[i][1])-1])) then
      Print("Q", i, ": trivial\n");
    else
    rl := aq[i][1];
    I := aq[i][2];
    t := Length(rl);
    kp := KPJStrata(rl, I);
    m := Length(kp);
    dri := Sum(List([1..t-1], i -> Sum(rl{[1..i]})^2))
      + 2 * Sum(List([1..t-1], i -> rl[i]*Sum(rl{[i+1..t]})))
      + Sum(rl{List(I, i -> i+1)});
    for j in [1..m] do
      if j mod 500 = 0 then Print("Q", i, ": ", j, " of ", m, " strata\r"); fi;
      dim := KPJActualdim(rl, I, kp[j]) - dri;
      if dim > 0 then
        Error("\nFalse: Quiver ", aq[i], "\n Stratum ", kp[j], "\n");
      elif dim >= -1 and not KPJSmooth(aq[i][1],aq[i][2],kp[j]) then
        Error("\nNot smooth: Quiver ", aq[i], "\n Stratum ", kp[j], "\n");
      fi;
    od;
    Print("Q", i, ": ", m, " strata\n");
  fi; od;
end;

CountStrata := function( n )
#
# Returns the total number of strata for all quivers of total size n.
#
  local aq, i, a;
  aq := AllQuivers(n);
  a := 0;
  Print(Length(aq)," quivers\n");
  for i in [1..Length(aq)] do
    Print("Q", i, " of ", Length(aq), "\r");
    a := a + Length(KPJStrata(aq[i][1],aq[i][2]));
  od;
  return a;
end;

###########################################################################
#
# SECOND PART.  A pair of bipartitions is called "compatible" if there exists
# a signed quasibipartition whose two subordinate bipartitions are the given
# ones.  (To be more precise, compatibility is not a symmetric relation, since
# the + and - sides are different.)  
#
# A "level" is a sequence of bipartitions in which consecutive members are
#   compatible.  A level is associated to data (r_i),I, just like the KPJ
#   strata.
#
# This section is a program to enumerate all levels and estimate their
# dimensions.  The main function is  
#
#    TestAllQuiversLev( n )
#
# whose use is just like TestAllQuivers above.
#
###########################################################################

#
# The various BucketFill functions enumerate ways to distribute n objects
# into a bunch of "buckets" of various capacities.  Only BucketVecFill3
# is actually used; the rest were written as "practice" (or because I thought
# I would use them and then changed the way I wanted to invoke the function).
# See CompatibleBipartitions for how it is used.
#

BucketFill := function( buc, n )
#
# buc is a vector of integers (bucket sizes).  n is an integer.
# Returns a list of ways to distribute n into the buckets.
#
  local b, a, bp, i, ap, l;
  b := Length(buc);
  if b = 0 then  return [ ];
  elif b = 1 then  if buc[1] >= n then return [ [n] ]; else return [ ]; fi;
  else
    a := [ ];
    bp := buc{[1..b-1]};
    for i in [0..Minimum(buc[b],n)] do
      ap := BucketFill(buc{[1..b-1]}, n-i);
      for l in ap do Add(l,i); od;
      Append(a, ap);
    od;
    return a;
  fi;
end;

BucketVecFill := function( buc, n )
#
# Same as above but returns answers as blocks of 1s and 0s
#
  local b, a, bp, i, ap, l, j;
  b := Length(buc);
  if b = 0 then  return [ ];
  elif b = 1 then
    if buc[1] >= n then return [ Concatenation(List([1..n], k->1),
      List([n+1..buc[1]], k->0)) ]; else return [ ]; fi;
  else
    a := [ ];
    bp := buc{[1..b-1]};
    for i in [0..Minimum(buc[b],n)] do
      ap := BucketVecFill(buc{[1..b-1]}, n-i);
      j := Concatenation(List([1..i], k->1), List([i+1..buc[b]], k->0));
      for l in ap do Append(l,j); od;
      Append(a, ap);
    od;
    return a;
  fi;
end;
  
BucketVecFill2 := function( buc, n )
#
# Same as above but returns answers as blocks of 1s, -1s, and 0s
#
  local b, a, bp, i, ap, l, jp;
  b := Length(buc);
  if b = 0 then  return [ ];
  elif b = 1 then
    if buc[1] >= n then return [ Concatenation(List([1..n], k->1),
      List([n+1..buc[1]], k->0)) ]; else return [ ]; fi;
  else
    a := [ ];
    bp := buc{[1..b-1]};
    for i in [0..Minimum(buc[b],n)] do
      ap := BucketVecFill2(buc{[1..b-1]}, n-i);
      jp := List([0..i], k -> Concatenation(List([1..k], h->1),
        List([k+1..i], h->-1), List([i+1..buc[b]], h->0)) );
      for l in Cartesian(ap, jp) do Add(a, Concatenation(l)); od;
    od;
    return a;
  fi;
end;
  
BucketVecFill3 := function( buc, n )
#
# Same as above but returns answers as blocks of 2s, 1s, -1s, and 0s
#
  local b, a, bp, i, ap, l, jp, two_pars;
  b := Length(buc);
  if b = 0 then  if n = 0 then return [ [ ] ]; else return [ ]; fi;
  elif b = 1 and 2*buc[1] < n then  return [ ];
  else
    a := [ ];
    bp := buc{[1..b-1]};
    for i in [0..Minimum(2*buc[b],n)] do
      ap := BucketVecFill3(buc{[1..b-1]}, n-i); 
      two_pars := Filtered(Partitions(i+2,2), h->h[1]<=buc[b]+1);
      # for h in two_pars:
      # 	for k from 0 to h[1]-h[2] (which is less eq. to buc[b])
      # 		make a list with h[2]-1 2's and k 1's and h[1]-h[2]-k -1's 
      # 		and the rest 0's. For all `k`, jp has length: 
      # 		h[2]-1+k+h[1]-h[2]-k+buc[b]-h[1]+1 = buc[b]
      jp := Concatenation(List(two_pars, h -> List([0..h[1]-h[2]],
                                         k -> Concatenation(List([1..h[2]-1], c->2), 
                                                            List([1..k], c->1), 
															List([k+1..h[1]-h[2]], c->-1),
															List([h[1]..buc[b]], c->0)))));
      for l in Cartesian(ap, jp) do Add(a, Concatenation(l)); od;
    od;
    return a;
  fi;
end;
  

CompatibleBipartitions := function( p, n, e )
#
# Assumes: p is a bipartition; n is an integer bigger than the size of p;
#  and e is "true" or "false."
# WARNING: p may not have trailing 0s in both partitions.  Should be ok, since
#  SignedQbpSubord (via Bipartitionize) does not return bipartitions in which
#  both members have trailing 0s.
#
# Returns the list of all bipartitions q of total size n such that there
# exists a signed quasibipartition whose subordinate bipartitions are p
# and q.  If e is true, p should be the + subordinate, and q the - subordinate;
# the other way around if e is false.
#
# PROFILE:
# gap> C := CompatibleBipartitions( [[3,3,3,2,2,2,1,1,1],[3,2,1,1,0,0,0,0,0]], 25, true);;
# gap> DisplayProfile();
#   count  self/ms  chld/ms  function                                           
#  180087      195       15  ADD_LIST                                           
#  411574      288      -10  APPEND_LIST                                        
#  353957      357      -13  ListOp: for a collection that is a list            
#  382391      363        1  SHALLOW_COPY_OBJ: for a presentation in default *  
# 1740952      702      -39  Iterator: for a mutable iterator                   
# 1376604      772       18  AbsoluteValue: rationals                           
# 1740952     1096        6  Int: for a rational                                
# 1030935     3811      409  PositionProperty: for dense list and function      
# 3481904     3005     1829  Int                                                
# 1030935      856     4211  PositionProperty                                   
# 124360*    50005     -360  Size: for a list that is a collection              
#            61451           TOTAL                                              
#
  local lam, q, xi, x, r, parts, add, add0, pmult, v, theta, mu, nu, ep, k, sb, a;

#  if e then q := StructuralCopy(p); else q := [p[2],p[1]]; fi;
  q := StructuralCopy(p);
  if Length(q[1]) = 0 then lam := [ ]; else lam := q[1] + q[2]; fi;
  xi := List(lam, i -> 2*i-1);
# xi should be thought of like this: draw p with '+'s in the boxes, and insert
# boxes containing '-'s between every pair of consecutive boxes of p.
  r := n + Sum(lam) - Sum(xi);
# Any putative signed quasibipartition having p as a subordinate will look like
# the diagram xi with r additional '-' boxes added on at the ends of rows.  The
# main task is to determine all possible ways to attach those extra boxes.  This
# is how BucketVecFill3 will be used.
  Append(xi, List([1..r], i -> 0)); # append 'r' zeros
  Append(q[1], List([1..r], i -> 0));
  parts := Reversed(Collected(xi)); # count freq's and reverse
  x := Length(xi);
  #Print("xi = ", xi, "\n");
  # (BJ) following can be combined into 1 highly optimized helper function
  # e.g. if `xi` has row 1 w/ mult m1, then we're looking
  # distribute some of the `r` boxes to the ends of the first m1 rows, but no
  # more than m1 boxes can go here. 
  # Then filter so that ?? 
  pmult := List(parts, p->p[2]);
  add0 := BucketVecFill3( pmult, r );
  add := Filtered(add0, v -> v[x-r+1]<>2);
  #Print("mult = ", pmult, "\n");
      
  #
  #Print("START v loop:\n\n");
  a := [];
  ### above lines def'ing add and this one can be combined in a generator loop
  for v in add do    
    theta := Filtered( List([1..x], i -> xi[i] + AbsInt(v[i])),  # AbsInt = abs(x)
               i -> i<>0 ); # filter to remove rows w/ length 0
    ### construct mu and nu from theta
    mu := []; nu := []; ep := [];
    for k in [1..Length(theta)] do
      if q[1][k] = 0 then  mu[k] := 0;
      else
        if e then  
          if v[k]>=1 then mu[k] := 2*q[1][k]; else mu[k] := 2*q[1][k]-1; fi;
        else
          if v[k]>=1 then mu[k] := 2*q[1][k]+1; else mu[k] := 2*q[1][k]; fi;
        fi;
      fi;
      if mu[k] > theta[k]  then  mu[k] := theta[k];  fi;
      nu[k] := theta[k] - mu[k];
      if e then if v[k]>=1 then ep[k] := -1; else ep[k] := 1; fi;
      else if v[k]>=1 then ep[k] := 1; else ep[k] := -1; fi;
      fi;
    od; 
    #
    #Print("v = ", v, "\n");
    #Print("mu = ", mu, "\n");
    #Print("nu = ", nu, "\n");
    #Print("ep = ", ep, "\n");
    sb := SignedQbpSubord([ mu, nu, ep ]); ### form the SQBP from mu, nu, ep
    #Print("sb = ", sb, "\n\n");
    if not e then  sb := [sb[2], sb[1]];  fi;
    # check if we should add it to `a` and avoid repetitions
    if EqBipart(q,sb[1]) and not ForAny(a,b->EqBipart(b,sb[2])) # ? can be skipped
      then Add(a, sb[2]); fi;
    # In the e false case, also add the one where the length 1 +s are in mu
    if not e then
      for k in Filtered([1..Length(theta)], i->q[1][k]=0 and v[k]=1) do
        mu[k] := Minimum(1,theta[k]);  nu[k] := theta[k] - mu[k];
      od;
      sb := SignedQbpSubord([ mu, nu, ep ]);
      if EqBipart(q,sb[2]) and not ForAny(a,b->EqBipart(b,sb[1])) # ? can be skipped
        then Add(a, sb[1]); fi;
    fi;
  od;
  return a;
end;


KPJLevels := function( rl, I )
#
# Returns the list of all levels associated to the data (rl,I).
#
  local t, kpj, a, penum, compats, n, v, k, w;

  t := Length(rl);
  if t = 0 then  return [ [ [ [ ], [ ] ] ] ];  fi;
  kpj := KPJLevels( rl{[1..t-1]}, Filtered(I, i -> i < t-1) );
  a := [];
  penum := [];
  compats := [];
  n := Sum(rl);
  for v in kpj do
    k := Position(penum, v[t]);
    if k = fail then
      Add(penum, v[t]);
      Add(compats, CompatibleBipartitions(v[t], n, not (t-1 in I)));
      k := Length(penum);
    fi;
    for w in compats[k] do  Add(a, Concatenation(v, [w]));  od;
  od;
  return a;
end;

KPJLevdim := function( rl, I, kp )
#
# Assumes that kp is a level associated to the data (rl,I).
# Returns the dimension estimate for that level, minus the expected dimension
# of the whole variety.
#
  local t, A, B;
  t := Length(kp)-1;
  A := Filtered([0..t-1], i -> (i in I) and not (i-1 in I));
  B := Filtered([1..t], i -> (i-1 in I) and not (i in I));
  return PartN(TransposedPart(rl)) - PartN(kp[t+1][1]+kp[t+1][2])
    + Sum(List(B, i -> Sum(kp[i+1][1])))
    - Sum(List(A, i -> Sum(kp[i+1][1])))
    - Sum(List(I, i -> rl[i+1]));
end;

KPJLevdimlist := function( rl, I )
#
# Analogue for levels of KPJEstdimlist.
#
  return List(KPJLevels(rl, I), k -> KPJLevdim(rl, I, k));
end;

TestAllQuiversLev := function( n )
#
# Main testing function.  Generates all quivers of total size n, and for each,
# enumerates all levels, and estimates their dimension.  If any such dimension
# estimate comes out positive, it reports the error.
#
  local aq, i, kp, m, j;
  aq := AllQuivers(n);
  Print(Length(aq)," quivers\n");
  for i in [1..Length(aq)] do
    Sort(aq[i][2]);
    if KPJ_SKIPTRIVIAL and (Length(aq[i][2]) = 0 or 
      aq[i][2] = [0.. Length(aq[i][2])-1] or (KPJ_SKIPCOTRIVIAL and
      aq[i][2] = [Length(aq[i][1])-Length(aq[i][2])..Length(aq[i][1])-1])) then
      Print("Q", i, ": trivial\n");
    else
    kp := KPJLevels(aq[i][1],aq[i][2]);
    m := Length(kp);
    for j in [1..m] do
      if j mod 500 = 0 then Print("Q", i, ": ", j, " of ", m, " levels\r"); fi;
      if KPJLevdim(aq[i][1], aq[i][2], kp[j]) > 0 then
        Error("\nFalse: Quiver ", aq[i], "\n Level ", kp[j], "\n");
      fi;
    od;
    Print("Q", i, ": ", m, " levels\n");
  fi; od;
#  return
#    ForAll( AllQuivers(n), q -> ForAll(KPJEstdimlist(q[1],q[2]), i -> i<=0) );
end;

CountLevels := function( n )
#
# Returns the total number of levels for all quivers of total size n.
#
  local aq, i, a;
  aq := AllQuivers(n);
  a := 0;
  Print(Length(aq)," quivers\n");
  for i in [1..Length(aq)] do
    Print("Q", i, " of ", Length(aq), "\r");
    a := a + Length(KPJLevels(aq[i][1],aq[i][2]));
  od;
  return a;
end;
