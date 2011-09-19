#
# File: quadratic_irrational_period.sage
# Author: Benjamin Jones <benjaminfjones@gmail.com>
#
#
# Tests: From Rosen pg. 502: period of fractions sqrt(D) for various D:
# for d in [23, 31, 46, 76, 97]:
#     p = quadratic_irrational_period(0,d,1)
#     print "sqrt(%d), period = %d" % (d, p)



## TODO: Not sure this is correct if c < 0 ???
def is_periodic(a, b, c):
    return bool((b > c*c - 2*a*c + a*a) and (b < c*c + 2*a*c + a*a) and (b > a*a))

## TODO: doesn't work in many cases. Do we need to check test(i) at all test points?
def qs_floor(a, b, c):
    """
    Returns the floor of (a + sqrt(b))/c without calling sqrt() to any
    great precision.
    """
    approx = floor((float(a)+sqrt(float(b)))/float(c))
    if c > 0:
        test = lambda i: ((i + approx)*c - a)**2 > b
    else:
        test = lambda i: ((i + approx)*c - a)**2 < b
    for i in range(-1, 2, 1):
        if test(i):
            break
    return i-1+approx

def quadratic_irrational_period(a, b, c, debug=False):
    """
    Returns the period of the continued fraction representation of a quadratic irrational number
    of the form $\frac{a + \sqrt{b}}{c}$.

    INPUT:

        - a - an integer
        - b - a natural number
        - c - an integer

    OUTPUT:

        A finite list [p_1, p_2, ... ] equaling the period of the continued fraction representation of `\frac{a + \sqrt{b}}{c}`.

    ALGORITHM:

        1. Input: $(a,b,c)\in Z \times N\times Z$  and $b$  not a perfect square. This triple represents $\frac{a+\sqrt{b}}{c}$.

        2. If $c$  does not divide $b-a^{2}$ then replace $(a, b, c)$ by the new triple
           $(a',b',c') = (|c| a, c^2 b, |c| c)$

        3. Set $P_{0}=a'$ , $Q_{0}= c'$, $d=b'$.

        4. For $k\geq0$  set

        (a) $\alpha_{k}=\frac{P_{k}+\sqrt{d}}{Q_{k}}$

        (b) $a_{k}=\lfloor\alpha_{k}\rfloor$

        (c) $P_{k+1}=a_{k}Q_{k}-P_{k}$

        (d) $Q_{k+1}=\frac{b-P_{k+1}^{2}}{Q_{k}}$

        5.  Then $\alpha_{0}=\frac{a+\sqrt{b}}{c}=[a_{0},a_{1},a_{2},...]$

        6. Let $k< l$  be the first pair of integers which satisfy that $$(P_{k},Q_{k})=(P_{l},Q_{l})$$  then the period of $\alpha_{0}$  is $(a_{k},\dots,a_{l-1})$ .

    """
    if b.is_square():
        raise ValueError("INPUT `b` should not be a perfect square")

    # ensure that c divides b-a**2
    if not c.divides(b-a**2):
        if debug: print "replacing ", (a,b,c), "by ",
        a,b,c = (abs(c)*a, c**2*b, c**2)
        if debug: print (a,b,c)

    # Main algorithm:
    #     1. Unwind the continued fraction until we find a purely period surd.
    #     2. At that point we record P_k and Q_k
    #     3. Continue unwinding until we find a repeat (P_l, Q_l) = (P_k, Q_k)
    cont = True
    found_periodic = False
    pre_period_counter = 1
    period_counter = 0
    # initial values of P_i, Q_i
    P_i = a
    Q_i = c
    while (cont):
        alpha_i = (P_i + sqrt(b))/Q_i  # exact expression
        ## Below, when floor() is called on an exact expression like alpha_i
        ## real interval field arithmetic with automatically adjusting precision up to 20000 bits is used if needed.
        ## Note: 20000 can be replaced with an adjustable parameter if wanted.
        #a_i = floor(alpha_i)
        a_i = qs_floor(P_i, b, Q_i)

        if debug:
            print "PQaa: ", (P_i, Q_i, alpha_i, a_i)
        if is_periodic(P_i,b,Q_i) and (not found_periodic):
            P_k = P_i
            Q_k = Q_i
            found_periodic = True
            if debug: print "Found periodic at ", pre_period_counter

        pre_period_counter += 1
        P_i=a_i*Q_i-P_i
        Q_i=(b-P_i**2)/Q_i  # note that this line uses the `new` P_i just computed, not the old one

        if found_periodic:
            period_counter += 1
            if (P_k, Q_k) == (P_i, Q_i):
                cont = False

    return period_counter



#
#  - OLD CODE -
#

# def conjugate(a,b,c):
#     """ Returns the conjugate of $(a + \sqrt(b))/c$ as a list (p,q,r)."""
#     return (-a, b, -c)

# def inverse(a,b,c):
#     """
#         Returns the inverse of $(a + \sqrt(b))/c$ as a list (p,q,r).
#         Note: $c/(a+\sqrt(b)) = c*(a-sqrt(b))/(a^2-b) = (-c*a + sqrt(c^2*b))/(b-a^2)
#     """
#     return (-c*a, c^2*b, b-a^2)

# def value(a,b,c):
#     """ Returns the expression (a+sqrt(b))/c """
#     return (a+sqrt(b))/c

# def is_periodic(a,b,c):
#     """ Returns true iff. value(a,b,c) > 1 and value(*conjugate(a,b,c)) is in (-1,0). """
#     vc = value(*conjugate(a,b,c))
#     return bool((value(a,b,c) > 1) and (vc > -1) and (vc < 0))
