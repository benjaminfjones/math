#d=3 example
(t1,t2,t3)=var('t1,t2,t3')
(u1,u2,u3)=var('u1,u2,u3')

s=1
T = [t1, t2,t3]
U = [u1, u2,u3]
#T.reverse()
#U.reverse()
for i in range(3):
    for j in range(3):
        if i < j:
            s *= (T[i]-T[j])*(U[i]-U[j])/(1-2*T[j]+T[i]*T[j])
        s *= (1-T[i])/(1-T[i]*(1+U[j]))

UPARTS = [\
[3,3,3],\
[3,3,2],\
[3,2,2],\
[3,3,1],\
[2,2,2],\
[3,2,1],\
[2,2,1],\
[3,1,1],\
[2,1,1],\
[1,1,1],
[3,3,0],\
[3,2,0],\
[2,2,0],\
[3,1,0],\
[2,1,0],\
[1,1,0],
[3,0,0],
[2,0,0],
[1,0,0],
[0,0,0] ]

for tpart in UPARTS:
    for upart in UPARTS:
        print "tpart =", tpart, "upart =", upart
        f=s
        vdict = { t1:3+tpart[0], t2:2+tpart[1], t3:1+tpart[2], \
                  u1:3+upart[0], u2:2+upart[1], u3:1+upart[2] }
        for v in U+T:
            nf=f.taylor(v,0,vdict[v]+1)
            f=nf.coeff(v,vdict[v])
        print "coeff =", latex(f)
