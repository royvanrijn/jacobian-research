# Prime-power Frobenius reduction for a single shared mirrored pair

Fix positive integers \(\ell,d\) and take the unaugmented
ordinary monomials \(X^{\ell+d}Y^\ell\) and \(X^\ell Y^{\ell+d}\).
At an upper order \(Q=p^q\equiv1\pmod{2d}\), with prime
\(p>2\ell+d\), the exact factorial pair has a residue-field
coprimality test of degrees bounded in \(\ell,d\), independently
of \(q\). For every sufficiently large \(p\), an explicit
arithmetic progression of the exponents \(q\) gives coprimality.
When \(\ell=1\), the simpler condition \(Q\equiv1\pmod{4d(d+1)}\)
already suffices for every \(p>d+2\). Equal rows contribute
only scalars.

The bounded-degree reduction uses coefficient valuations and
Frobenius on a common residue root. It concerns the exact
ordinary factorial arrays, with no limiting coefficient
comparison or checks at individual orders.

## Two compressed residue polynomials

Normalize the ordinary moment polynomial at degree \(m\)
to be monic, with coefficient
\[
 c_{m,k}=\binom mk
    \frac{(m\ell+dk)!\,(m(\ell+d)-dk)!}
         {(m\ell)!\,(m(\ell+d))!}.
\]
Put \(M=(Q-1)/d\), which is even, and set
\(\epsilon=(-1)^{d+1}\). In this proof \(z\) is the polynomial
variable and \(X\) is an auxiliary variable to be evaluated
at \(z^M\). Define the rational constants
\[
                b_u=(-1)^u
                 \frac{d(\ell+u)!\,(\ell+d-u-1)!}
                      {\ell!\,(\ell+d)!},
                   \qquad 0\le u\le d-1 .
\]
The upper residue polynomial \(\overline P_Q\) satisfies
\[
 \begin{aligned}
 (1-\epsilon z)\overline P_Q(z)&={\cal U}(z,z^M),\\
 {\cal U}(z,X)&=1+(b_0-1)\epsilon z+
       \epsilon z\sum_{u=1}^{d-1}(b_u-b_{u-1})X^u\\
   &\hspace{22mm}+(1-\epsilon b_{d-1})zX^d-\epsilon z^2 X^d .
 \end{aligned}                                                   \tag{1}
\]
Here and below rational coefficients in residue equations are
reduced modulo \(p\).

For the lower residue, use
\[
 {\cal F}(z,X)=\sum_{u,r}
      (-1)^v\,\frac{(\ell+u-1)!\,(\ell+d-u-1)!}
                      {(\ell-r-1)!\,(\ell+d+r-1)!}\,z^v X^u ,
 \tag{2}
\]
where
\[
 0\le u\le d,\quad -\ell-d+1\le r\le\ell-1,\quad
 v=(u+r)/d\ \hbox{is an integer}.
\]
At \(u=0\) impose \(v\ge0\), and at \(u=d\) impose \(v\le0\).
Then \(\overline P_{Q-1}(z)={\cal F}(z,z^M)\). This
substitution has no negative powers even when (2) does.

To verify (1), factor the coefficient as
\[
 c_{Q,k}
  =\frac{\binom Qk}{\binom {dQ}{dk}}\,
       \frac{\binom{\ell Q+dk}{\ell Q}
             \binom{(\ell+d)Q-dk}{\ell Q}}
            {\binom{(\ell+d)Q}{\ell Q}} .
\]
For \(0<k<Q\) write \(dk=uQ+s\), \(0<s<Q\). The first
ratio equals \(\binom{Q-1}{k-1}/\binom{dQ-1}{dk-1}\),
so Lucas gives the unit ratio
\((-1)^{k-s}/\binom{d-1}{u}\). The other three binomials
reduce to \(\binom{\ell+u}{\ell}\),
\(\binom{\ell+d-u-1}{\ell}\), and \(\binom{\ell+d}{\ell}\).
The product reduces to \(b_u\epsilon^k\).
The band of indices with this \(u\) is
\(uM<k\le(u+1)M\). Using both normalized endpoints and
\(\epsilon^M=1\) proves (1).

For the lower order \(Q-1\), the binomial factor is always a
unit and reduces to \((-1)^k\) by Lucas. At a level \(P=p^h\)
with \(1\le h\le q\), put \(t=dk\bmod P\). Subtracting the two
endpoint factorial floors from the two coefficient floors
gives the contribution
\[
                   1-\mathbf1_{t<\ell}
                       -\mathbf1_{t>P-\ell-d}.
\]
Indeed the changing floor offsets from multiples of \(P\)
are \(t-\ell\) and \(-t-\ell-d\), and both endpoint offsets
are negative and greater than \(-P\). Since \(p>2\ell+d\),
these contributions are nonnegative. There are no higher
levels because all factorial arguments are less than \(pQ\).
Having contribution zero at \(P=Q\) is equivalent to
\(dk=uQ+r\) with \(-\ell-d+1\le r\le\ell-1\).
It then vanishes at every lower level as well, because \(r\)
is smaller than \(p\) in absolute value.
For such \(k\), the two factorial arguments are
\((\ell+u)Q-(\ell-r)\) and
\((\ell+d-u)Q-(\ell+d+r)\). The analogous endpoint
arguments have \(u=r=0\). After taking the unit parts
of the factorials, the high base-\(p\) digits give the
top factorials in (2) divided by the corresponding endpoint
factorials. The least digits give the opposite endpoint
factor and the bottom factorials shown in (2). The other
digits are \(p-1\) in all four factorials.

For clarity, a factorial unit at \(p\) equals
\((-1)^{v_p(n!)}\) times the product of the digit factorials
of \(n\) modulo \(p\). This follows by repeated division by
\(p\) and Wilson's theorem on each full block. The valuation
signs cancel in the ratio here. In the least digits use
\((p-1-h)!\equiv(-1)^{h+1}/h!\pmod p\); their signs also
cancel. These are the usual Lucas and unit-factorial
congruences. [Granville's account](https://www.cecm.sfu.ca/organics/papers/granville/paper/binomial/html/node2.html)
Now \(k=uM+v\), with \(v=(u+r)/d\). The two boundary
restrictions in (2) are exactly the requirements
\(0\le k\le Q-1\) at \(u=0,d\). For interior \(u\),
they hold automatically since \(Q>2\ell+d\).
The binomial sign is \((-1)^v\), completing the proof
of (2). Both residue polynomials have unit constant
and leading coefficients in the original variable \(z\).

## A nonzero bounded elimination polynomial

Clear negative powers in (2) by multiplying by \(z^{\ell+d}\),
and set
\[
 {\cal R}(z)=\operatorname{Res}_X
          \bigl(z^{\ell+d}{\cal F}(z,X),{\cal U}(z,X)\bigr).
\]
This is a nonzero polynomial over \(\mathbb F_p\), with
\(\deg{\cal R}\le D:=2d(\ell+d+1)\).
The degree estimate follows from the Sylvester determinant,
using \(\deg_X{\cal F},\deg_X{\cal U}\le d\),
\(\deg_z(z^{\ell+d}{\cal F})\le2(\ell+d)\) and
\(\deg_z{\cal U}\le2\).

For nonvanishing, work over the algebraic closure of
\(\mathbb F_p((z))\) with the \(z\)-adic valuation. Any root
\(X\) of \({\cal U}\) has valuation \(-1/d\). Indeed its
constant coefficient is a unit, its lower nonconstant
coefficients have valuation at least one, and its highest
coefficient has valuation one because
\(1-\epsilon b_{d-1}=\ell/(\ell+d)\). A valuation
higher than \(-1/d\) makes the constant term uniquely
smallest in valuation. A valuation lower than \(-1/d\)
does the same for the highest term.

For such an \(X\), the term of \({\cal F}\) with powers
\(z^vX^u\) has valuation \(v-u/d\). Write \(\ell-1=dL+r_0\),
\(0\le r_0<d\). At \(1\le u\le r_0\), the least admissible
\(v\) is \(-L-1\). At \(r_0<u\le d\) it is \(-L\).
At \(u=0\) all admissible \(v\) are nonnegative. Thus
the unique least valuation occurs at \(u=r_0,\ v=-L-1\)
if \(r_0>0\), and at \(u=d,\ v=-L\) otherwise.
Its coefficient is a \(p\)-unit and \({\cal F}(z,X)\ne0\).
This proves the asserted nonvanishing of \({\cal R}\).

Any common root of the original pair can be reduced to a
nonzero common root of the monic integral residue polynomials.
It is then a root of \({\cal R}\). To test for such a root
one may remove factors \(z\) from \({\cal R}\), compute
\(z^M\bmod{\cal R}\) by modular exponentiation, and take
the gcd of \({\cal R}\) with
\(z^{\ell+d}{\cal F}(z,z^M)\) and \({\cal U}(z,z^M)\).
A constant gcd proves coprimality of the original pair.
The multiplied factor \(1-\epsilon z\) in (1) can introduce
a spurious candidate in this sufficient test.

## Explicit infinite sequences of terminal orders

Form the two rational numbers
\[
 A_{\ell,d}={\cal F}(-1,1)>0,\qquad
 B_{\ell,d}=2-d^{-1}\sum_{u=0}^{d-1}b_u>1 .
\]
The first is positive because every sign \((-1)^v\) cancels
at \(z=-1\). For the second, the moduli of the \(b_u\)
are at most \(d/(\ell+d)<1\), since their factorial product
is log-convex in \(u\) and its maxima are at the endpoints.
Choose \(p>2\ell+d\) avoiding the two numerators, and let
\(q\) be a positive multiple of the multiplicative order
of \(p\) modulo
\[
                  2d\,\operatorname{lcm}_{1\le g\le D}(p^g-1).
 \tag{3}
\]
This modulus is prime to \(p\), and the resulting \(Q=p^q\)
is congruent to one modulo \(2d\). If there were a nonzero
common residue root \(z\), it would be a root of \({\cal R}\)
and hence belong to some \(\mathbb F_{p^g}\) with \(g\le D\).
By (3), \(M=(Q-1)/d\) is divisible by \(p^g-1\), so \(X=z^M=1\).
Equation (1) gives
\[
                        {\cal U}(z,1)=(1-\epsilon z)(1+z).
\]
The candidate \(-1\) is excluded by \(A_{\ell,d}\).
The only other candidate is \(+1\) when \(\epsilon=1\),
which can be an artifact of multiplying the upper polynomial
by \(1-z\). The actual upper residue at \(+1\) is
\(2+M\sum_u b_u\equiv B_{\ell,d}\ne0\), since each band
has length \(M\) and now \(\epsilon=1\).
Thus there is no common residue root.

For \(\ell=1\), (2) reduces to \({\cal F}=1+X+\cdots+X^d\).
In fact the only allowable term for \(u\) has \(r=-u,v=0\)
and the factorials cancel. Suppose \(Q\equiv1\pmod{4d(d+1)}\).
At a common residue root \(z\), \(X^{d+1}=1\) and
\(X\ne1\). Frobenius fixes \(X\) and sends \(z\) to \(z/X\),
because \(z^Q=zX^d=z/X\). Equation (1), viewed as a quadratic
in \(z\) at this \(X\), has both \(z\) and \(z/X\) as
distinct roots. Its constant is one and its leading term
is \(-\epsilon X^d\). Their product therefore gives
\(z^2=-X^2/\epsilon\). Frobenius fixes this square, so
\(z^{2Q}=z^2/X^2=z^2\) forces \(X^2=1\). Hence \(X=-1\),
\(d\) is odd, and \(\epsilon=1\). Now \(z^2=-1\) and
\(Q\equiv1\pmod4\) force \(z^Q=z\), contradicting \(z^Q=-z\).
This proves coprimality in the simpler progression,
with no additional prime exclusions.

These residue exclusions prove coprimality of the actual
ordinary polynomials over \(\mathbb Q\) and hence over every
characteristic-zero field. They apply to the indicated
terminal orders and leave the other orders open.
