# Angular separation near two prime-power private exponents

Suppose an ordinary two-monomial profile has exactly one private
variable on each side, with exponents \(A=p^j\) and \(B=p^\ell\)
for one prime \(p\) and integers \(j,\ell\ge0\). Either exponent
can be one. Equal rows can be omitted. Let
\[
                V=\sum_{\rm shared}
                     \frac{(b_i-a_i)^2}{\min(a_i,b_i)}.
\]
There is an explicit positive bound \(v_{J,p,A,B}\), given below,
such that consecutive moments have disjoint root arguments in
every degree if \(V<v_{J,p,A,B}\). The shared rows may point in
either direction and there is no restriction on their signed total
gap. The proof uses the
[absolute ordinary eventual angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
\(J\), whose numerical value is unevaluated. The pure private
profile \(V=0\) has all-degree angular separation directly,
independently of that cutoff.

The small-tail condition is an inequality in the exact
[matched tanh data](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid).
At any reference consecutive pair \(n,n+1\), put
\[
 u=\left(\frac{n}{T'}-\frac{n-1}{T}\right)^{-1}
\]
on every matched Gauss slot belonging to a shared row. That row
has \(d=|b_i-a_i|\) slots, each with \(u=d/(2\min(a_i,b_i)+d)\).
Thus \(V=\sum_{\rm shared\ slots}2u/(1-u)\).
The criterion uses the exact private factors and bounds the shared tail.

## Modular angular separation of the exact private profiles

For \(D=p^h\) define the positive integer
\[
                  a_{D,k}=\frac{(Dk)!}{k!\,p^{k(D-1)/(p-1)}}.
\]
Repeatedly stripping the multiples of \(p\) proves this is a
\(p\)-adic unit. More precisely, for \(L\ge0\),
\[
 \frac{(pL)!}{L!p^L}
     =\prod_{r=0}^{L-1}\prod_{i=1}^{p-1}(pr+i)
                \equiv(-1)^L\pmod p.
\]
Wilson's theorem therefore gives
\(a_{D,k}\equiv\sigma_D^k\pmod p\), with
\(\sigma_D=(-1)^{(D-1)/(p-1)}\), including \(D=1\).
Set
\[
     G_m(z)=\sum_{k=0}^m g_{m,k} z^k,\qquad
             g_{m,k}=a_{A,m-k}a_{B,k}.
\]
The two private factorials cancel the binomial coefficient up to
\(m!\), after extracting \(k!\) and \((m-k)!\). Thus \(G_m\) is
the pure private factorial polynomial up to a positive scalar and
positive rescaling of \(z\).

Every root of \(G_m\) is a \(p\)-adic unit, since all coefficients
are \(p\)-integral and its endpoints are units. Its reduction is the
geometric polynomial
\[
      G_m(z)\equiv\sigma_A^m\sum_{k=0}^m
                  (\sigma_B\sigma_A^{-1} z)^k\pmod p.
\]
The quotient of the reductions of any two roots consequently has
multiplicative order dividing the prime-to-\(p\) part of \(m+1\).

Suppose two pairs of roots of \(G_n,G_N\), \(N=n+1\), have equal
quotients, and the two roots in each pair are distinct. Work in a
single algebraic closure with one extension of the \(p\)-adic
valuation. Both quotients reduce to one, since their reduced orders
divide the coprime numbers \(n+1,n+2\). One of those two numbers is
prime to \(p\). At that degree the reduced geometric polynomial is
squarefree, so two distinct roots cannot have identical reductions.
This is a contradiction.

For a common nonreal argument, take the two pairs to be each root
and its complex conjugate. Positive coefficients exclude the positive
real ray. On the negative ray the even pure private moment is
strictly positive, from its expectation as an even power of the
difference of two independent positive monomial values. This proves
the pure angular separation, and also shows that the two sets of
pure values \(z/\bar z\) are disjoint.

## A uniform margin for the bounded degrees

Take an integer \(J\ge2\) at least the absolute cutoff. Put
\[
  D_0=\max(A,B),\quad M=(D_0J)^{D_0},\quad
  \gamma=\tfrac12(2M^{2J+2})^{-J^4},\quad r=\frac{\gamma}{8M},\quad H=(J+1)^2(2M^2)^J.
\]
We use the sufficient bound
\[
    v_{J,p,A,B}=
        \min\left\{\frac1{2(J+1)},
           \frac{r}{H\,[M(J+1)]^{J-1}}\right\}.    \tag{1}
\]

For \(1\le m\le J\), the ratio \(a_{D,k+1}/a_{D,k}\) is the integer
\(D\prod_{i=1}^{D-1}(Dk+i)/p^{(D-1)/(p-1)}\) for \(k<m\).
Its numerator has \(p\)-adic valuation \((D-1)/(p-1)\), because
\(v_p(Dk+i)=v_p(i)\) for \(1\le i<D\). Hence this ratio lies in
\([1,M]\), and all \(g_{m,k}\le M^J\). Adjacent ratios in \(G_m\)
lie in \([M^{-1},M]\), so its roots have moduli in that interval
by Eneström–Kakeya.

Take pure roots \(\alpha,\beta\) in consecutive degrees \(n,N\le J\).
With \(X=g_{n,n}\alpha\) and \(Y=g_{N,N}\beta\), both \(X,Y\) and
their ordinary complex conjugates are algebraic integers. All their
conjugate moduli are at most \(M^{J+1}\), since their monic
polynomials are obtained by substituting the variable divided
by the leading coefficient. Put \(Z=X\bar Y-\bar X Y\). The
modular argument proves \(Z\ne0\) if both roots are nonreal.
If one is real, \(Z\ne0\) by strict positivity of the even
moment on the real line. Each of \(X,\bar X,Y,\bar Y\) has degree
at most \(J\), so their joint field has degree at most \(J^4\).
The norm of the nonzero algebraic integer \(Z\) has absolute
value at least one, and each conjugate of \(Z\) has modulus
at most \(2M^{2J+2}\). It follows that
\[
    \left|\frac{\alpha}{\bar\alpha}-\frac{\beta}{\bar\beta}\right|
              \ge (2M^{2J+2})^{-J^4}.
\]
The circular distance of their arguments is therefore at least
\(\gamma\), since the double-angle exponential is \(2\)-Lipschitz.

The pure roots within a degree also have a uniform mesh. After
cancelling the private units, the coefficients before any positive
rescaling are moments under the left endpoint tilt of a positive
ratio. The
[geometric moment angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry)
apply to this ratio. The positive boxes have gap greater than
\(\pi/m\) between successive labels. At the positive cut their
reflection gap exceeds \(2\pi/m\). The last box has its upper endpoint
at most \(\pi-\pi/(m+1)\) in even degree and
\(\pi-2\pi/(m+1)\) in odd degree, with the negative real root
present only in odd degree. Thus the circular mesh is at least
\(\pi/(m+1)\). The constant ratio case \(A=B=1\) is the geometric
polynomial and has this mesh as well. Degree one has only one root.

Given this mesh and both moduli at least \(1/M\), the distance
between any two pure roots in degree \(m\) is at least
\(2/[M(J+1)]\). In fact for moduli \(s,t\) and circular angular
distance \(\theta\) the squared distance is
\((s-t)^2+4st\sin^2(\theta/2)\), and
\(\sin(\theta/2)\ge\theta/\pi\) on \([0,\pi]\).
Thus the closed disks of radius \(r\) about the pure roots
are disjoint. At every boundary point, the distance to any
other pure root is at least \(1/[M(J+1)]\), since
\(r\le1/[M(J+1)]\). Factoring the pure polynomial and using its
integer leading coefficient gives the lower bound
\(r[M(J+1)]^{-(J-1)}\) on all those boundaries.

## Adding the exact ordinary shared factors

Write \(d_i=b_i-a_i\) on a shared row. For a degree \(m\) and
index \(k\), put \(x=m a_i\), \(q=d_i k\). Strip the positive
scalar \(x!\) and the power \(x^q\) from \((x+q)!\).
The logarithm of the remaining factor is
\[
 E_i(m,k)=
 \begin{cases}
    \sum_{r=1}^q \log(1+r/x), &q\ge0,\\
    -\sum_{r=0}^{-q-1}\log(1-r/x), &q<0.
 \end{cases}
\]
For each logarithm bound the absolute value by
\(r/[m\min(a_i,b_i)]\), using integration of \(1/t\) on the
interval between \(x-r\) and \(x+r\) as appropriate.
It follows that for \(E(m,k)=\sum_i E_i(m,k)\),
\[
 0\le E(m,k)\le \frac{V(k^2+k)}{2m}\le\frac{m+1}{2}V,\qquad
 0\le E(m,k+1)-E(m,k)\le V .
\]
In the second bound each new summand has \(r\le |d_i|(k+1)\)
and there are \(|d_i|\) new terms.

All the stripped index powers are absorbed in a positive
rescaling of \(z\), independently at each degree. These scalings
include \(m^{\sum d_i}\), which is why no signed total gap
condition is necessary for angular separation. The exact moment
polynomial in the new variable is
\[
                          \widetilde G_m(z)
                  =\sum_{k=0}^m g_{m,k} e^{E(m,k)} z^k .
\]
Its roots have modulus at most \(M\), because its consecutive
coefficient ratios increase by a factor \(e^{E(m,k+1)-E(m,k)}\ge1\).
For \(m\le J\), (1) makes \(E(m,k)\le1/4\) and
\(e^{E(m,k)}-1\le(J+1)V\). On \(|z|\le2M\), which contains
all the closed pure-root disks, we have
\(\lvert\widetilde G_m(z)-G_m(z)\rvert\le HV\).
By (1) this is strictly below the pure modulus on each boundary.
Rouché gives one actual root in each pure-root disk. Each pure
root has modulus at least \(1/M\), so its argument moves by
less than \(2Mr=\gamma/4\), using \(|\log(1+w)|\le2|w|\)
at \(|w|\le1/2\). The pure consecutive angular margin
therefore persists.
All \(N>J\) are covered by the absolute angular cutoff for
ordinary polynomials. The argument uses the exact factorial
arrays and has no finite root computations.
