# The full exceptional fibre in a quadratic leading coefficient

**Theorem KellerQuadraticMateFullExceptionalFibre.**
Let \(r=g(x)+x^m v\) with \(m\ge2\) over a characteristic-zero
field. If \((P,Q,r)\) is a polynomial Keller triple and \(Q\)
is quadratic in \(y\), then its leading \(y\)-coefficient
is divisible by \(r-g(0)\) in the polynomial ring.
There is no ordinary degree bound on the triple.

Work over \(\mathbb C\) by embedding the finitely generated
coefficient field. Since \(r\) is a submersion, \(g'(0)\ne0\).
Subtract \(g(0)\) and write
\[
                    r=xH,\qquad H=g(x)/x+x^{m-1}v .
\]
Put \(B_0=\mathbb C[x,v]\), and let \(q,b,c\) be the coefficients
of \(Q=qy^2+by+c\). Set \(K=c-b^2/(4q)\) and \(Y=y+b/(2q)\).
The bracket along \(r\) is \([P,Q]_D=\kappa\ne0\), where
\(D=-x^m\partial_x\) at fixed \(r\).

By the
[quadratic-mate theorem](KELLER_REPEATED_PARALLEL_SCALAR_QUADRATIC_MATE.md),
all prime factors of \(q\) occur in scalar fibres, and \(q\)
is a square modulo \(\mathbb C(r)^\times\) in the generic field.
The only reducible scalar fibre is \(r=0\), with distinct
factors \(x,H\). Hence \(q=e(r)x^i H^j\) with \(e(0)\ne0\),
\(i,j\ge0\) and \(i-j\) even. If \(r\) does not divide \(q\),
at least one of \(i,j\) is zero. The two possible forms are
\(q=e(r)x^{2l}\) and \(q=e(r)H^{2l}\), \(l\ge0\).

The even part of \(P\) in \(Y\) is \(\Psi(Q,r)\), with
coefficients in \(B_0[q^{-1}]\cap\mathbb C(r)\). These
coefficients have no pole at \(r=0\). Indeed, if a reduced
fraction in \(r\) has such a pole and belongs to \(B_0[q^{-1}]\),
the substituted Bézout identity forces \(r\) to divide a
power of \(q\). Since \(r=xH\) is squarefree, this would give
\(r\mid q\). It follows that the odd part \(O=P-\Psi(Q,r)\)
is regular along both components of \(r=0\).

Write \(f\) for the coefficient of \(Y\) in \(O\).
The constant bracket coefficient gives \(-fDK=\kappa\);
thus \(f\) is a unit of \(B_0[q^{-1}]\). Write \(q=e(r)w^2\)
using \(w=x^l\) in the first case and \(w=H^l\) in the second.
For \(Z=wY\) the odd-part formula of the quadratic-mate
theorem is
\[
               O=\int_0^Z F(Q-eu^2)\,du,\qquad F(K)=f/w ,
                                                               \tag{1}
\]
where \(F\) is a nonzero polynomial over \(\mathbb C(r)\).
Use \(d=1+\deg F\).

## Coefficient supported at the line \(x=0\)

Here \(f=f_\star(r)x^k\), where \(k\) is an integer and
\(f_\star\) is a rational unit at \(r=0\). Integrating
\(K_x=\kappa/(x^m f)\) gives
\[
        K=\alpha(r)x^t+\beta(r),\qquad t=1-m-k ,
\]
where \(\alpha\) is a unit at \(r=0\). The exponent \(t=0\)
is impossible because it would require a rational primitive
of a nonzero multiple of \(x^{-1}\).

Since \(F(K)=f_\star x^{k-l}\) is a monomial in \(x\)
over \(\mathbb C(r)\), comparison of Laurent powers gives
\(F(T)=\eta(r)(T-\beta)^{d-1}\), with
\(t(d-1)=k-l\) and \(\eta\) a unit at \(r=0\).
In particular \(td=1-m-l<0\).
On \(H=0\) the functions \(x\) and \(q\) are units and \(K\)
is regular. Thus \(\beta\) has no pole at \(r=0\). On \(x=0\),
in contrast, \(K\) has a pole since \(t<0\). Set \(y=0\),
and write \(Z_0=b/(2ew)\); then \(eZ_0^2=c-K\), so \(Z_0\)
has a pole. In (1) at \(y=0\), the coefficient of the
highest power \(Z_0^{2d-1}\) is a unit, and all other
coefficients are regular. The resulting pole contradicts
regularity of \(O\) along \(x=0\).

## Coefficient supported at the component \(H=0\)

Here \(f=f_\star(r)H^k\), with \(f_\star\) again a unit at
\(r=0\). Since \(H=r/x\), integration of \(K_x\) gives
\[
          K=\alpha(r)x^t+\beta(r),\qquad
          t=k-m+1\ne0,\qquad \operatorname{ord}_0\alpha=-k .
\]
Now \(F(K)=f_\star r^{k-l}x^{l-k}\), so the same comparison gives
\[
 \begin{aligned}
 F(T)&=\eta(r)(T-\beta)^{d-1},&
 t(d-1)&=l-k,&
 \operatorname{ord}_0\eta&=dk-l=(d-1)(m-1).
 \end{aligned}                                               \tag{2}
\]
At \(x=0\) the coefficient \(q=e(r)H^{2l}\) is a unit
because \(H|_{x=0}=g'(0)\ne0\), so \(K\) is regular there.
The term \(\alpha x^t\) has order \(t-k=1-m\) along
\(x=0\). Therefore \(\beta\) has a pole of order exactly
\(m-1\) at \(r=0\).

At \(H=0\), \(x\) and \(e\) are units, and \(r\) has order
one. Evaluate (1) at \(y=0\) and put \(\sigma=c-\beta\),
so \(eZ_0^2=c-K\).
If \(t>0\), then \(k>m-1\) and \(\alpha x^t\) has strictly
greater pole order than \(\beta\). Thus \(eZ_0^2\) has pole
order \(k\), strictly larger than that of \(\sigma\).
The highest term in \(Z_0\) of the integral dominates
and has order
\[
              \operatorname{ord}_0\eta-(d-\tfrac12)k<0 .
\]
If \(t<0\), then \(\beta\) dominates \(\alpha x^t\), so
\(eZ_0^2\) and \(\sigma\) have equal leading term and pole
order \(m-1\). In this case write (1) as
\[
 O|_{y=0}=\eta Z_0\,\sigma^{d-1}
       \sum_{j=0}^{d-1}\frac{(-1)^j}{2j+1}
                \binom{d-1}{j}\left(\frac{eZ_0^2}{\sigma}\right)^j .
\]
The sum is a unit because its residue is
\(\int_0^1(1-u^2)^{d-1}\,du\ne0\).
By (2), the order of the displayed expression is then
\(-(m-1)/2<0\). Both cases give a pole of \(O|_{y=0}\)
at \(H=0\), contrary to its regularity.
