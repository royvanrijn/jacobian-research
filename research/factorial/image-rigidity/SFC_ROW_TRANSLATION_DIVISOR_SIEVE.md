# A divisor sieve for factorial windows under row translation

Consider \(s\ge2\) distinct monomial exponent columns. In one
coordinate their exponents are \(L,L+e_1,\ldots,L+e_{s-1}\),
where \(L\) is a nonnegative integer and all \(e_i>0\).
The other coordinates are fixed arbitrary nonnegative integers.
Let \(I\) be any set of \(s\) distinct positive orders,
take \(j\in I\), and write \(I'=I\setminus\{j\}\) and
\(D=\prod_{q\in I'}q\).
The terminal \(j=\max I\) always has the property that it
does not divide any \(q\in I'\).

Write \(Z=(Z_0,\ldots,Z_{s-1})\) for the coefficient variables.
For \(|\alpha|=q\) let \(C_{q,\alpha}\) be the multinomial
coefficient times the factorials contributed by all unshifted
coordinates. Set \(e_0=0\). Divide the order-\(q\) moment by
\((qL)!\) to obtain the integer polynomial form
\[
 P_q(L,Z)=\sum_{|\alpha|=q}
 C_{q,\alpha}(qL+1)_{\sum_i e_i\alpha_i}Z^\alpha,
 \qquad
 (u)_k=\prod_{\nu=0}^{k-1}(u+\nu).
 \tag{1}
\]
For \(a=\min_{i>0}e_i\) and \(w=\max_{i>0}e_i-a\),
define \(U_q(A,Z')\), \(Z'=(Z_1,\ldots,Z_{s-1})\),
by deleting column zero and using exponents
\(A+e_i-a\) in the chosen coordinate, with the same factorial
normalization at \(A\). Thus
\[
 P_q(L,0,Z')=(qL+1)_{qa}\,U_q(L+a,Z').
 \tag{2}
\]
Let \(R(A)=\operatorname{Res}_{Z'}(U_q(A,Z'):q\in I')\)
be its homogeneous resultant. For one variable, the resultant of
the one form \(cZ_1^q\) means \(c\).

## 1. A necessary divisibility at every admissible mark

**Theorem.** For \(1\le h\le a\) with \(\gcd(h,j)=1\), define
\[
 B_{q,h}=\prod_{\nu=1}^{qa}(j\nu-qh),\qquad
 W_{q,h}(Z')=j^{qw}U_q(a-h/j,Z'),\qquad
 R_h=\operatorname{Res}_{Z'}(W_{q,h}:q\in I').
\]
All these quantities have integer coefficients. Put
\(c_j=C_{j,j\mathbf e_0}\), the coefficient of \(Z_0^j\) in \(P_j\).
A common nonzero coefficient direction for the ordinary moments
at orders \(I\), over any characteristic-zero field, implies
\[
 jL+h\ \bigm|
 c_j^{D}\, R_h^{\,j}
            \prod_{q\in I'} B_{q,h}^{\,jD/q}.
 \tag{3}
\]
If \(j\) divides none of the other orders and \(R(a-h/j)\ne0\),
the integer on the right is nonzero. Its positive divisors give
a finite candidate list for the ordinary baseline \(L\).
Using several marks intersects their candidate lists.

To prove the theorem, evaluate the homogeneous resultant of all
the \(P_q(L,\cdot)\) in integers at the proposed baseline. It is
zero if the moment system vanishes in a nonzero direction over a
characteristic-zero field. Reduce its value modulo \(M=jL+h\).
For \(M=1\) the divisibility is automatic. Otherwise \(j\) is a
unit modulo \(M\). At order \(j\), every term other than
\(c_jZ_0^j\) contains a rising-factorial factor \(jL+h\),
because its positive exponent offset is at least \(a\).
Consequently \(P_j\equiv c_jZ_0^j\bmod M\).
The homogeneous resultant identity for replacing a form by a
variable power, and homogeneity in each form, give
\[
 0\equiv \pm c_j^D
       \operatorname{Res}_{Z'}(P_q(L,0,Z'):q\in I')^{\,j}
                  \pmod M.
\]
This is a universal integer polynomial identity and so also
holds for composite \(M\). See
[Busé–D'Andrea–Sombra–Weimann, Section 2 and Proposition 2.1](https://www.maia.ub.edu/~sombra/publications/Flexes/Flexes.pdf).

The degree in \(A\) of each coefficient of \(U_q\) is at most
\(qw\), so \(W_{q,h}\) is integral. By (2), evaluating
\(jL\equiv-h\bmod M\) gives the integral congruence
\[
 j^{q(a+w)}P_q(L,0,Z')\equiv B_{q,h}W_{q,h}(Z')\pmod M.
\]
Scaling the face forms uses the homogeneity exponents \(D/q\).
The resulting power of \(j\) is a unit modulo \(M\).
This proves (3). A factor \(j\nu-qh\) can vanish only if
\(j\mid q\), since \(\gcd(h,j)=1\). Finally,
\[
                   R_h=j^{(s-1)wD}\,R(a-h/j),
\]
again by homogeneity. This proves the claimed nonvanishing and
gives an exact finite test by evaluating the original resultant
at the remaining integer baselines.

In particular, if the face resultant \(R(A)\) is nonzero as a
polynomial and
\[
 \#\{1\le h\le a:\gcd(h,j)=1\}>\deg R ,
 \tag{4}
\]
some mark works. One always has \(\deg R\le(s-1)wD\) by the
coefficient degree bound and multihomogeneity. Thus, for a fixed
face profile with nonzero resultant polynomial at \(I'\), adding
one monomial at a sufficiently separated lower exponent in the
chosen row preserves a finite baseline decision at \(I\).
This step can be repeated if the marked row and separation
conditions hold at each deletion. A useful sufficient test for
the input \(R(A)\ne0\) is a single value \(A\) with projective
nonvanishing at \(I'\).

For a row with strictly increasing offsets, delete them in
increasing order, marking the largest remaining order each time.
At a stage with \(r\) remaining columns and order set \(I_r\),
let \(a_r\) be the gap from the least offset to the second least,
let \(w_r\) be the spread from the second least to the largest,
and let \(j_r=\max I_r\) and \(D_r=\prod_{q\in I_r\setminus\{j_r\}}q\).
Sufficient inequalities for a finite baseline decision are
\[
                  a_r\ge1+j_r(r-1)w_rD_r\qquad(2\le r\le s).
 \tag{5}
\]
Indeed the marks \(1+k j_r\) for \(0\le k\le(r-1)w_rD_r\)
all fit (4). The single-column moment resultant at the final
stage is nonzero, so induction gives nonzero resultant
polynomials for every support in this deletion order.
The other exponent rows are arbitrary throughout.

## 2. Flat faces and bivariate three-term windows

If all other columns have the same exponent offset \(a\), then
\(w=0\) and the face forms \(U_q\) are the ordinary factorial
moments for the deleted support in the remaining coordinates.
For a consecutive support-size window and \(j=\max I\),
\(I'\) is the shorter consecutive window starting at the same
order. A proved Strong Factorial theorem on that support
implies its face resultant is nonzero. This implication includes
coefficient directions with some zeros: their actual support
windows have length at most \(s-1\) and fit within \(I'\).
The resulting divisor certificate at \(h=1\) bounds the
ordinary baseline at every prescribed start, uniformly in the
phase coefficients and over every characteristic-zero field.

For example, consider the bivariate support
\[
             x^L y^{b_0},\quad x^{L+a}y^{b_1},
                 \quad x^{L+a}y^{b_2},\qquad b_1\ne b_2.
\]
For any three positive consecutive orders, the shorter pair of
orders has nonzero face resultant by
[ordinary univariate binomial Strong Factorial](ORDINARY_BINOMIAL_RESIDUAL_REGION.md).
Thus (3) gives a finite integer baseline list for every start
and arbitrary \(a>0\) and \(b_i\ge0\). The same use of the
[ordinary univariate quadratic window](UNIVARIATE_QUADRATIC_ALL_SUPPORTS_STRONG_FACTORIAL.md)
applies with up to three distinct remaining \(y\)-exponents
chosen from \(0,1,2\). For unequal positive offsets \(a,a+w\)
in a three-term support, the face is an ordinary binomial
family with one varying row. Its resultant polynomial is
nonzero, and its degree is at most \(w n(n+1)\), by the
[single-row specialization and weighted degree bound](MULTIVARIATE_SINGLE_ROW_ANGULAR_EXCEPTIONAL_FIBERS.md#the-smaller-bound-for-actual-common-roots).
Thus \(a\ge1+(n+2)w n(n+1)\) already suffices for the
window \(n,n+1,n+2\), with no restriction on the other
exponent coordinates. This sufficient inequality applies
also when the ordinary physical coprimality of the two-term
face at its original baseline is unknown.

For \(s=2\) and \(I=\{n,n+1\}\) with \(n\ge2\), either
order may be marked. With repeated copies of the same varied
row, (1) has the indicated rising factorial raised to the
multiplicity, and the argument still applies with the same
change in the powers of \(B_{q,h}\). In that case (3) provides
a divisor filter for the existing
[one unbounded row type finite reduction](MULTIVARIATE_ONE_UNBOUNDED_TYPE_FINITE_DECISION.md).
That result and the
[soft-type reduction](MULTIVARIATE_SOFT_TYPE_FINITE_REDUCTION.md)
already give the absolute finite box for ordinary bivariate
binomials. The new higher-support bound is for each prescribed
order set. At support size at least three no uniform terminal
cutoff is supplied by this argument.

### Initial windows with arbitrary deleted rows

**Theorem.** At the initial orders \(I=\{1,2,3\}\), every
admissible upper mark \(j=3\), \(1\le h\le a\),
\(\gcd(h,3)=1\), has \(R(a-h/3)>0\). The two deleted
columns may differ in any of the unmarked coordinates. There
is no separation condition on \(a\) or \(w\). Consequently
every fixed three-column profile with a unique row minimum
has an exact finite list of potentially failing translations
in that row, uniformly over characteristic-zero coefficient fields.

For the proof, put \(A=a-h/3\ge2a/3>0\), and write the
two face monomials on the positive orthant as
\[
 \psi_0(x)=x_1^A\prod_{i>1}x_i^{b_{i0}},\qquad
 \psi_1(x)=x_1^{A+w}\prod_{i>1}x_i^{b_{i1}}.
\]
Use the product-exponential measure \(d\mu=e^{-\sum x_i}dx\).
The real powers are unambiguous there, and their first and
second integrals converge. The two functions are linearly
independent because the exponent columns are distinct. Write
\[
 U_1=pZ_0+rZ_1,\qquad U_2=uZ_0^2+2vZ_0Z_1+zZ_1^2.
\]
The factorial normalization in (1) gives
\[
 p=\frac{\int\psi_0\,d\mu}{\Gamma(A+1)},\quad
 r=\frac{\int\psi_1\,d\mu}{\Gamma(A+1)},\quad
 (u,v,z)=\frac{\bigl(\int\psi_0^2,\int\psi_0\psi_1,
                         \int\psi_1^2\bigr)}{\Gamma(2A+1)},
\]
with the same measure in the three numerator integrals.
In particular \(p,r>0\), and
\[
 R(A)=p^2z-2prv+r^2u
     =\frac{1}{\Gamma(2A+1)}
          \int(r\psi_0-p\psi_1)^2\,d\mu>0.                 \tag{I1}
\]
Indeed \(U_1=0\) has only the projective direction \((-r,p)\),
and \(U_2\) is positive there. This proves nonvanishing over
\(\mathbb C\). At rational \(A\), the normalized forms have
rational coefficients, so their nonzero resultant excludes a
common direction over every characteristic-zero field as well.
This uses positivity of the first two Gamma moments directly;
it does not invoke an ordinary integer-exponent binomial theorem.

For an explicit divisor integer, let
\[
 R_h=3^{4w}(p^2z-2prv+r^2u)\in\mathbb Z_{>0}.
\]
Equation (3) specializes to
\[
 3L+h\mid c_3^2 R_h^3 B_{1,h}^{6}B_{2,h}^{3}.               \tag{I2}
\]
Every factor on the right is nonzero. Taking \(h=1\) alone
therefore gives a finite divisor list; using all admissible
\(h\) intersects those lists. For a homogeneous bivariate
trinomial one may use either exponent row. The other row,
including its common translation, is fixed when forming the
integer in (I2). Thus this is a finite decision for that fixed
row profile, not a bound independent of both translations or
the support. The candidate \(L=0\) always survives the
\(h=1\) divisibility test, so nonvanishing of the marked face
alone does not eliminate the physical residual.

The [exact face checker](verify_sfc_initial_beta_faces.py) checks
75 rational marked faces with arbitrary orientations in up to
two unmarked rows, a repeated-column failure of strictness,
and four specializations of the full ternary resultant to (I2).
These finite controls check the formulas and normalization;
the square integral proves the unrestricted statement. Replay
from the repository root with Python 3.12.3 and SymPy 1.14.0:

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_sfc_initial_beta_faces.py
~~~

At higher starts the deleted orders need not be \(1,2\).
The fractional/Gamma-binomial window and elimination of the
surviving integer translations remain separate obligations.

### A smaller initial-window divisor integer

**Proposition.** Fix a homogeneous bivariate profile
\((0,d),(a,d-a),(d,0)\), integers \(0<a<d\), and a
nonnegative integer height \(H\). Put \(b=d-a\) and use
the polynomial \(\Pi=S/(K_{a,d}(L+H+d)^3)\) from
[the total-degree factor](SFC_HOMOGENEOUS_TRINOMIAL_SHORT_SPANS.md#a-total-degree-factor-in-every-eliminated-moment).
Then
\[
 J_{a,d,H}=3^\nu\Pi(-1/3,H)\in\mathbb Z_{>0},\qquad
 \nu=4d-\lceil a/2\rceil-3.                              \tag{I3}
\]
If the ordinary moments at \(1,2,3\) have a common nonzero
coefficient direction at height \(L\ge0\), then
\(3L+1\mid J_{a,d,H}\). This applies to every support and
fixed height \(H\); the constant need not be expanded as a
polynomial to be computed.

Here is an exact product formula for the constant, including
the mixed cubic contribution. Normalize both exponent rows
by \((mL)!(mH)!\). Put \(\ell=-1/3\), \(A_*=a-1/3\), and
\[
 \begin{gathered}
 p=(H+1)_b,\quad r=(A_*+1)_b,\\
 u=(2H+1)_{2b},\quad
 v=(2A_*+1)_b(2H+1)_b,\quad z=(2A_*+1)_{2b},\\
 \mathcal B=p^2z-2prv+r^2u>0,\quad
 t_1=(\ell+1)_a,\quad t_2=(2\ell+1)_{2a},\\
 c_3=(3H+1)_{3d},\quad p_2=(\ell+1)_d,\quad
 C_*=t_1^2t_2\mathcal B,\quad D_*=c_3p_2^3.
 \end{gathered}
\]
The square-integral proof of (I1), with the common unmarked
factorials also removed, proves \(\mathcal B>0\).
At this mark the third form is exactly \(c_3Z_0^3\).
On the first zero plane its only nonzero coefficient is
\(D_*\), while the quadratic's \(Y^2\) coefficient is
\(C_*\). Therefore
\[
 J_{a,d,H}=
 \frac{3^\nu C_*^2D_*}
      {K_{a,d}(-1/3,H)(H+d-1/3)^3}.                       \tag{I4}
\]
All factors on the right are positive. Its integrality
follows from (I3)'s polynomial degree bound, not from an
unstated cancellation of rational Gamma values. For an
integer root \(L\) of \(\Pi\), evaluation modulo \(3L+1\)
allows substitution \(L\equiv-1/3\), because three is a
unit. This proves the divisibility.

One can also remove a forced fifth power from the old
resultant constant. With
\[
 R_*=3^{2b}\mathcal B\in\mathbb Z_{>0},\quad
 B_1=\prod_{j=1}^a(3j-1),\quad
 B_2=\prod_{j=1}^{2a}(3j-2),\quad
 M=c_3^2R_*^3B_1^6B_2^3,\quad \mu=3(H+d)-1,
\]
one has \(\mu^5\mid M\), and a failed physical translation
must satisfy
\[
            3L+1\mid\gcd\bigl(M/\mu^5,J_{a,d,H}\bigr).  \tag{I5}
\]
For integrality of \(R_*\), each term of \(\mathcal B\)
has degree at most \(2b\) in \(A_*\). Using this smaller
clearing exponent instead of the \(4b\) in (I2) only removes
powers of three, which are units in the divisor congruence.
The two-column case of (12)–(14) in the linked note gives
\(A_*+H+b\mid\mathcal B\) as an integer polynomial before
specialization. Hence \(\mu\mid R_*\); also \(\mu\mid c_3\).
This proves \(\mu^5\mid M\).

To justify removing that power in the divisor condition,
let \(F(L,H)\) be the full homogeneous resultant of the
three normalized forms. At the mark it equals
\(c_3^2C_*^3\), so
\(M=3^{12a+6b}F(-1/3,H)\). The exact total-degree factor
gives \(F=(L+H+d)^5F_0\) with \(F_0\) integral.
At an ordinary common zero, \(F_0=0\) since \(L+H+d>0\).
Thus reduction modulo \(3L+1\) gives the nonzero divisor
integer \(M/\mu^5\). Any remaining power of three in a
denominator is invertible in that congruence. Combining
this with (I3) proves (I5).

The [exact reduced-sieve checker](verify_sfc_trinomial_reduced_sieve.py)
compares (I4) with an independently expanded \(\Pi\) and
evaluates every candidate from (I5) in these six controls:

| \((a,d,H)\) | Mark \(M\) | After removing \(\mu^5\) | After (I5) |
|---|---:|---:|---:|
| \((1,2,0)\) | 108 | 48 | 15 |
| \((1,3,0)\) | 1,944 | 1,404 | 216 |
| \((2,3,0)\) | 1,836 | 1,296 | 144 |
| \((1,3,1)\) | 5,616 | 2,496 | 346 |
| \((2,4,0)\) | 30,780 | 13,680 | 1,404 |
| \((1,4,0)\) | 31,104 | 18,144 | 1,800 |

The entries count positive divisors congruent to one modulo
three. All 3,925 final candidates have \(\Pi(L,H)\ne0\).
These profiles already lie within the existing all-height
span theorem; they test the elimination mechanism rather
than enlarge its range. The universal conclusion is the
sharper necessary condition (I5). Uniform elimination of
its divisors remains open. In particular \(L=0\) still
passes a mark whose modulus is one. If \(\Pi=0\) at a
candidate, the second cubic remainder or the original
moment system must still be checked; \(\Pi=0\) alone
does not exhibit a counterexample.

Replay with Python 3.12.3, python-flint 0.9.0 and SymPy 1.14.0:

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_sfc_trinomial_reduced_sieve.py
~~~

## 3. Integer marks using an ordinary deleted face

Integer specializations cause all restricted higher forms to
vanish before normalization. Removing the exact common
resultant power gives another divisor test, this time using
an ordinary deleted-face resultant.

**Theorem.** Take \(n=\min I\), \(I^+=I\setminus\{n\}\),
\(D^+=\prod_{q\in I^+}q\), and let \(\ell\ge1\) satisfy
\(a\ge(n+1)\ell\). Use the face forms \(U_q\) of (2), now with
\(q\in I^+\), and suppose
\[
 R^+=\operatorname{Res}_{Z'}(U_q(a-\ell,Z'):q\in I^+)\ne0.
\]
If the original ordinary moment forms have a common nonzero
direction over a characteristic-zero field, then
\[
 L+\ell\ \bigm|\
 c_n^{D^+}\,(R^+)^n\,
     \prod_{q\in I^+}
      \left(q\prod_{\substack{1\le u\le qa\\u\ne q\ell}}
                       (u-q\ell)\right)^{nD^+/q}.
 \tag{6}
\]
The right hand side is nonzero. No bound on the width \(w\) of
the deleted exponent row is used.

To prove this, put \(L=-\ell+t^n\), \(Z_0=tV\), and
write \(b_q=q\prod_{1\le u\le qa,\ u\ne q\ell}(u-q\ell)\).
For \(q>n\), every term involving \(Z_0\) at exponents
\(1\le\alpha_0\le n\) has offset at least
\(a(q-\alpha_0)\ge(n+1)\ell(q-n)\ge q\ell\).
It thus contains \(qL+q\ell=qt^n\) and an extra \(t^{\alpha_0}\).
A term with \(\alpha_0>n\) has an extra power of \(t\) from
\(Z_0^{\alpha_0}\) alone. The terms on the face have limit
\(b_q U_q(a-\ell,Z')\) after dividing by \(t^n\), by (2).
In order \(n\) the pure anchor contributes \(c_n V^n\)
after division by \(t^n\). Every other term has offset at least
\(a>n\ell\), hence a factor \(nt^n\), and only face terms
survive beside that pure anchor. The limiting order-\(n\) form
is \(c_nV^n+b_n U_n(a-\ell,Z')\).

These limiting forms have nonzero resultant. In fact the
higher forms depend only on \(Z'\) and have no common nonzero
direction there. Adding a scalar multiple of the face term
\(b_n U_n(a-\ell,Z')\) to \(c_nV^n\) therefore leaves the
resultant nonzero for every complex scalar. As a polynomial
in that scalar it is constant. Set the scalar to zero and
apply dimension descent. This gives the integer on the right
in (6), up to sign, as the limiting resultant.

Let \(F(L)=\operatorname{Res}_Z(P_q(L,Z):q\in I)\) and
\(\Pi=n D^+\). Multihomogeneity multiplies the resultant by
\(t^{-n\Pi\sum_{q\in I}1/q}\) when every form is divided by
\(t^n\). The common change \(Z_0=tV\) multiplies it by
\(t^\Pi\), by the linear change-of-variables rule for homogeneous
resultants. For this diagonal change, the rule follows by
noting that the zero locus of the universal resultant is
unchanged for \(t\ne0\), and comparing the coefficient of the
resultant on the coordinate powers. The scaled resultant is
therefore
\[
      t^{-n k} F(-\ell+t^n),\qquad
      k=\sum_{q\in I^+}\Pi/q .
\]
Since the limit is nonzero, \(F(L)\) has a zero of exact order
\(k\) at \(L=-\ell\). The polynomial \(G(L)=F(L)/(L+\ell)^k\)
has integer coefficients and its value at \(-\ell\) is the
right hand side of (6) up to sign. At an ordinary failed
baseline \(G(L)=0\), so \(L+\ell\mid G(-\ell)\).

For a consecutive window, the input to (6) is the ordinary
shorter window starting one order later. For three bivariate
monomials the deleted family is binomial and has at most two
unequal rows. The existing absolute height reduction for
two-monomial profiles supplies a computable \(H\) such that \(a-\ell>H\) suffices when the
deleted columns differ in the marked row. If they agree in
that row, the remaining univariate binomial theorem suffices
at every baseline. With unequal deleted row exponents this
gives a ternary reduction for
\(a\ge\max\{n+1,H+2\}\), by taking \(\ell=1\), independently
of how large the upper offset is. By the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
for binomials, \(a\ge n+1\) suffices once the shorter window's
terminal is at or above that cutoff. Below it, for
\(a\ge(n+1)\ell\) the ordinary deleted resultant is itself
an exact finite test at the selected \(n,\ell\).

## 4. Fractional marks at primes of their order

Consider a three-term window \(I=\{n,n+1,n+2\}\). Write the offsets
in the marked row as \(0,a,a+w\) with \(a,w>0\). The other
exponent rows are arbitrary. To use (3), it suffices to prove that
the deleted pair of binomial moments at the two orders \(q\ne j\)
has no common root at baseline \(A=a-1/j\).
After dividing by its nonzero endpoint scalar, the polynomial
of degree \(q\) takes the form
\[
 V_q(Z)=\sum_{k=0}^q \binom qk (qA+1)_{wk}
         \prod_i \frac{(q b_i+(c_i-b_i)k)!}{(q b_i)!}\, Z^k .
 \tag{7}
\]
Here \((b_i,c_i)\) are the two face entries in each unmarked row.
Only unequal pairs need to be kept. The rational specialization
has positive Gamma parameters, so both homogeneous endpoints
are nonzero. Common roots can be tested in any algebraic closure
of a \(p\)-adic field because the coefficients are rational.

### Even starts with arbitrary unmarked rows

For a nonnegative integer \(t\) put \(B(t)=\lceil\log_2(t+1)\rceil\),
its binary length. Define
\[
 K_q=2 B(q)+\sum_i B(q\max(b_i,c_i)).
\]
If there are no unequal unmarked rows, use \(K_q=0\).
Set \(d=\sum_i(c_i-b_i)\).

**Theorem.** Suppose \(n\) is even and \(w>K_n+K_{n+1}\).
Then the upper mark \(j=n+2,\ h=1\) gives a nonzero integer
in (3). This yields a finite exact baseline decision even when
the lower gap \(a\) is one. If the two deleted columns agree
on the unmarked rows, every \(w\ge1\) works.

Write \(e=v_2(j)\). At the odd order \(q=n+1\) we have
\(v_2(qA)=-e\), so each factor in \((qA+1)_{wk}\) has that
valuation. At the even order \(q=n\), when \(e\ge2\) we have
\(v_2(q)=1\), and the corresponding factorial valuation is
\(-(e-1)wk\). When \(e=1\), \(qA\) is integral at 2, and the
valuation of this rising factorial is nondecreasing in \(k\).
We use only this monotonicity in the latter case.

For order \(q\), the other coefficient factors in (7) have
valuation \(dk\) plus a constant and an error of oscillation at
most \(K_q\). Indeed \(v_2(m!)=m-s_2(m)\), where \(s_2\) is the
binary digit sum. Each digit sum lies between zero and the
indicated binary length. Also \(v_2\binom qk\) lies between zero
and \(2B(q)\) by the same formula. When those rows are equal,
all the unmarked factorial ratios are one. In that case the
binomial valuations are nonnegative with zero endpoints, so
they do not affect the bounds for root valuations below.

If coefficient valuations are \(\sigma k+C+\epsilon_k\)
with error oscillation at most \(K\), all root valuations
belong to \([-\sigma-K,-\sigma+K]\). At a larger valuation the
constant term strictly dominates every other term, and at a
smaller one the leading term strictly dominates. Adding a
nondecreasing function of \(k\) to the coefficient valuations
preserves the upper bound on root valuations by the same
constant-term comparison. For the odd parent this gives the
lower bound \(ew-d-K_{n+1}\). For the even parent the upper
bound is \((e-1)w-d+K_n\), including \(e=1\). In the empty-core
case \(K_q=0\) follows directly by the endpoint observation.
The two bounds are disjoint under the stated inequality.
For an even \(n\ge4\), the same argument applies to the lower
mark \(j=n\) using \(q=n+1,n+2\) and the bound
\(w>K_{n+1}+K_{n+2}\). The two marks then give independent
divisor conditions on the same baseline.

### Odd starts when only the marked row differs on the face

The following elementary estimate uses the linked factorial
coefficients to retain half the pole depth. Let \(T\in\mathbb Q_p\)
with \(v_p(T)=-e<0\). For positive integers \(q,w\) put
\[
 J_q(U)=\sum_{k=0}^q \binom qk
           \prod_{\nu=1}^{wk}(1+\nu/T)\, U^k .
\]
Every root \(u\) is a unit, since the coefficients are integral
and both endpoints are units. Moreover
\[
                         v_p(1+u)\ \ge\ e/2.                 \tag{8}
\]
To see this, write \(E_\ell(x)\) for the polynomial giving
the \(\ell\)-th elementary symmetric function of \(1,\ldots,wx\)
at every nonnegative integer \(x\). Its degree is at most \(2\ell\):
for a set of length \(m\), the difference of the corresponding
\(\ell\)-th symmetric functions at \(m\) and \(m-1\) is \(m\)
times the \((\ell-1)\)-st one at \(m-1\). Induct and sum the
polynomial finite difference, then substitute \(m=wx\).
The expansion of \(E_\ell(x)\) in the binomial basis
\(\binom xr\) therefore has integral coefficients \(e_{\ell,r}\)
and \(r\le2\ell\), by the forward-difference formula at zero.
Expanding in \(1/T\) gives the exact identity
\[
 J_q(U)=(1+U)^q+
 \sum_{\ell=1}^{qw}T^{-\ell}\!
 \sum_{r=0}^{\min(q,2\ell)}
    e_{\ell,r}\binom qr U^r(1+U)^{q-r}.
 \tag{9}
\]
For a unit \(u\) with \(\theta=v_p(1+u)<e/2\), every summand
in the correction has valuation at least \(\ell e+(q-r)\theta\),
strictly greater than \(q\theta\). Thus it cannot cancel the
first term. This proves (8), with \(u=-1\) included.

**Theorem.** If the two deleted columns agree on every unmarked
row, then both even and odd starts admit a finite exact baseline
decision, with the only possible exceptions to this assertion
among odd \(n\ge3\) satisfying
\[
                         n(n+2)\mid(2^w-1)^2 .               \tag{10}
\]
For each fixed \(w\) the exceptional start list is finite and
computable independently of the lower gap \(a\) and all unmarked
exponents. The list is empty for \(w=1,2\) and for every
power of an odd prime.

The even case was proved above. At an odd endpoint \(j=n\) or
\(j=n+2\), both of the remaining orders are prime to any
\(p\mid j\), since \(p\) is odd and their differences from \(j\)
are one and two. For \(n\ge3\), neither remaining order is
a multiple of \(j\). Let \(e=v_p(j)\); both \(qA\) have
valuation \(-e\). A common root \(Z\) of the marked face would
give unit roots \(u_q=(q A)^w Z\) of the corresponding \(J_q\),
using \(T=qA\) in (8). Their ratio is the \(w\)-th power of the
ratio of the orders. By (8) its distance from one has valuation
at least \(e/2\). Modulo \(p^e\) the ratio of the order at distance
two from \(j\) to the one at distance one is \(2\).
It follows that \(2v_p(2^w-1)\ge e\), and hence the vanishing
of this marked face forces \(j\mid(2^w-1)^2\).
If neither endpoint works, then both \(n\) and \(n+2\)
divide this square. They are coprime for odd \(n\), giving (10).

For \(n=1\) use upper mark \(j=3\). The order-one root of the
face is real. The order-two face polynomial has negative
discriminant because \(\bigl((2A+1)_w\bigr)^2<(2A+1)_{2w}\)
for \(A>0,\ w>0\). For \(w=1,2\), the squares in (10) are \(1,9\) and neither has
the required pair of divisors. If \(w=\ell^u\) for an odd prime
\(\ell\) and \(u\ge1\), the order of \(2\) at every prime
dividing \(2^w-1\) is a nontrivial power of \(\ell\).
All such primes, and hence all divisors of the square in (10),
are \(1\bmod \ell\). Two of them cannot differ by two.

The integer mark of Section 3 also applies for every \(a\ge n+1\)
in this one-row face case. At \(\ell=1\) its ordinary
deleted face has a univariate binomial moment pair, coprime by
the [ordinary binomial theorem](ORDINARY_BINOMIAL_RESIDUAL_REGION.md).
Thus, for each \(w\), the starts and lower gaps missed by
these three marks lie in the finite set of odd \(n\ge3\)
satisfying (10) and \(1\le a\le n\).

## 5. A determinant bound for faces of arbitrary support size

Suppose the \(r=s-1\) columns on the deleted face agree in every
unmarked row. Their offsets \(f_i=e_i-a\) in the marked row are
then distinct nonnegative integers, including zero. Define the
nonzero integer
\[
 D_{\mathbf f}=\det(d^{f_i})_{1\le d,i\le r},
 \qquad
 C_{\mathbf f}=\operatorname{lcm}(1,\ldots,r)\,D_{\mathbf f}^{\,2}.
 \tag{11}
\]
The determinant is nonzero because a real linear combination
of the \(r\) distinct powers cannot have \(r\) distinct positive
zeros, by Descartes' rule of signs.

**Theorem.** At start \(n\), the upper fractional mark \(j=n+r\),
\(h=1\), gives a finite exact baseline decision whenever
\(j\nmid C_{\mathbf f}\). The lower mark \(j=n\) does the same
when \(n>r\) and \(n\nmid C_{\mathbf f}\). Thus the only starts
missed by both criteria satisfy \(n+r\mid C_{\mathbf f}\) and,
if \(n>r\), \(\operatorname{lcm}(n,n+r)\mid C_{\mathbf f}\).
This gives an effective finite list of exceptional starts for
each fixed face width profile, uniformly in \(a\) and the
other exponent entries, with no bound on the support size.

For the proof, first extend (8) to a homogeneous form with
any number of coefficient variables. Write
\[
 {\cal J}_q(Y)=\sum_{|\alpha|=q}\binom q\alpha
            \prod_{\nu=1}^{\sum_i f_i\alpha_i}(1+\nu/T)\,Y^\alpha ,
          \qquad v_p(T)=-b<0.
\]
If \({\cal J}_q(Y)=0\) and \(Y\ne0\), then
\[
             v_p\Bigl(\sum_i Y_i\Bigr)
                   \ \ge\ \min_i v_p(Y_i)+ b/2 .              \tag{12}
\]
To prove it, scale first so that every \(Y_i\) is integral and
some \(Y_i\) is a unit. At degree \(\ell\) in \(1/T\), the
elementary symmetric coefficient is an integer-valued polynomial
in \(\alpha\) of total degree at most \(2\ell\), by the same
finite difference as in (9) with \(m=\sum_i f_i\alpha_i\).
Its coefficients in the product basis
\(\prod_i\binom{\alpha_i}{\beta_i}\) are iterated forward
differences, hence integers, and \(|\beta|\le2\ell\).
For each such basis term the sum over \(|\alpha|=q\) is
\[
  \binom q\beta\,Y^\beta
     \Bigl(\sum_iY_i\Bigr)^{q-|\beta|},
 \qquad
 \binom q\beta=\frac{q!}{(q-|\beta|)!\prod_i\beta_i!}
\]
when \(|\beta|\le q\), and zero otherwise. If the sum of the
\(Y_i\) has valuation \(\theta<b/2\), the order-\(\ell\)
correction has valuation at least
\(\ell b+(q-|\beta|)\theta>q\theta\). The zeroth-order term
cannot cancel. Undoing the scale proves (12).

Take an endpoint mark \(j\) as in the theorem, and suppose
its deleted-face resultant vanishes. Its remaining orders are
\(q=j-d\) for the upper mark and \(q=j+d\) for the lower, with
\(1\le d\le r\). Fix \(p\mid j\), put \(e=v_p(j)\), and let
\(E=\max_{1\le d\le r}v_p(d)\). If \(e\le E\), no estimate is
needed at \(p\). Otherwise \(v_p(q)=v_p(d)\le E<e\) for every
\(q\), and \(A=a-1/j\) has valuation \(-e\). Divide out the
unmarked scalar moment factors. If \(Z'\) is a common nonzero
direction for the deleted-face equations, put \(c_i=A^{f_i}Z_i\)
and choose its homogeneous scale with \(\min_i v_p(c_i)=0\).
The order-\(q\) equation is \({\cal J}_q(Y)=0\) with
\(Y_i=q^{f_i}c_i\) and \(T=qA\). Each \(Y_i\) is integral and
\(b=e-v_p(q)\ge e-E\). Equation (12) gives
\[
 v_p\Bigl(\sum_i q^{f_i}c_i\Bigr)\ge(e-E)/2
\]
at every remaining order. If \(M=(q^{f_i})_{q,i}\), its adjugate
is integral. At a coordinate with unit \(c_i\), the equality
\(\operatorname{adj}(M)Mc=\det(M)c\) therefore gives
\(v_p(\det M)\ge(e-E)/2\).
Modulo \(p^e\), the determinant of \(M\) is \(D_{\mathbf f}\)
up to sign, since \(q\equiv\pm d\). The bound
\((e-E)/2\le e\) transfers to \(v_p(D_{\mathbf f})\).
Thus \(v_p(j)\le E+2v_p(D_{\mathbf f})\) at every prime of \(j\),
or \(j\mid C_{\mathbf f}\). The stated order conditions make
\(j\nmid q\) for every remaining \(q\), as required in (3).
For a nonvanishing mark, (3) bounds every possible ordinary
baseline and gives the exact finite test.

## 6. Simultaneous affine translations of the exponent rows

Let \(s=r+1\) distinct columns have entries
\[
                 \rho_i L+u_i+e_{ik}\qquad (0\le k\le r)
\]
in row \(i\), with positive integer speeds \(\rho_i\), integer
intercepts \(u_i\), and nonnegative offsets \(e_{ik}\) whose minimum
in each row is zero. The ordinary baselines are the nonnegative
integers \(L\) for which every \(\rho_i L+u_i\ge0\).
In one marked row take \(\rho_i=1,\ u_i=0,\ e_{i0}=0\), and
\(e_{ik}>0\) for \(k>0\).
On the deleted face put \(a_i^F=u_i+\min_{k>0}e_{ik}\),
\(f_{ik}=u_i+e_{ik}-a_i^F\), and \(g_k=\sum_i f_{ik}\).
Assume the \(g_k\) are distinct. This condition says that the
face columns have distinct total degrees, even if they differ
in every row. Put
\[
 C_{\mathbf g,\rho}=\operatorname{lcm}(1,\ldots,r)\,
      \operatorname{lcm}_i(\rho_i)\,
            \det(d^{g_k})_{1\le d,k\le r}^{\,2}.
\]

**Theorem.** For a consecutive \(s\)-window select the upper order
\(j=n+r\), or the lower order \(j=n\) when \(n>r\).
Suppose
\[
              \rho_i\notin (j u_i,\,j(u_i+e_{i0})]
                    \quad\text{for every row }i.              \tag{13}
\]
If \(j\nmid C_{\mathbf g,\rho}\), there is a computable finite
list of all ordinary baselines that can fail the window on this
affine family, over every characteristic-zero field.

For three terms with distinct total degrees on the deleted face,
an additional rule works at every even start. Choose the endpoint
\(j\) divisible by four, and assume (13) and that all the speeds
are odd. Then the same finite decision holds regardless of
\(C_{\mathbf g,\rho}\).

To prove these claims, strip the common factorial of each row at
each order. The resulting moment forms are polynomials in \(L\),
with the multiplier \((q(\rho_i L+u_i)+1)_{\sum_k e_{ik}\alpha_k}\)
for row \(i\). At \(L=-1/j\), the marked row kills every term of
the order-\(j\) form except the pure coefficient of column zero.
Condition (13) ensures this coefficient survives in all other
rows: it could acquire a zero factor there only at index
\(\rho_i-j u_i\), which must lie in \([1,j e_{i0}]\).
Thus projective nonvanishing reduces to the \(r\) other orders
on the deleted face.

For a prime \(p\mid j\), let \(e=v_p(j)\),
\(E=\max_{1\le d\le r}v_p(d)\), and \(R=\max_i v_p(\rho_i)\).
If \(e>E+R\), put \(A_i=a_i^F-\rho_i/j\). Then
\[
         v_p(q A_i)=v_p(d)+v_p(\rho_i)-e<0,\qquad q=j\pm d .
\]
In particular all the stripped scalars needed to rebase each row
at \(A_i\) are nonzero. Set \(c_k=Z_k\prod_i A_i^{f_{ik}}\),
and normalize \(\min_k v_p(c_k)=0\). Write \(Y_k=q^{g_k}c_k\).
After removing the face scalars, the form at order \(q\) is
\[
 \sum_{|\alpha|=q}\binom q\alpha Y^\alpha
   \prod_i\prod_{\nu=1}^{\sum_k f_{ik}\alpha_k}
                       (1+\nu/(q A_i)).
\]
The half-pole estimate (12) still holds with
\(b=\min_i(-v_p(q A_i))\). Indeed, at order \((\ell_i)_i\)
in the inverse parameters, the coefficient of \(Y^\alpha\)
has a product of integer-valued polynomials in \(\alpha\),
of total degree at most \(2\sum_i\ell_i\). The pole factor has
valuation at least \(b\sum_i\ell_i\), so the same binomial-basis
proof applies. A common face zero would therefore force
\[
 v_p\Bigl(\sum_k q^{g_k}c_k\Bigr)\ge(e-E-R)/2.
\]
The integral adjugate and reduction modulo \(p^e\) from
Section 5 give
\(e\le E+R+2v_p(\det(d^{g_k}))\). If \(e\le E+R\), this bound
already holds. Thus if the face fails, every prime of \(j\)
satisfies the corresponding valuation bound, giving the required
divisibility. To justify this argument when a face scalar at
the formal point vanishes, note that if \(j\nmid C_{\mathbf g,\rho}\)
there is a prime violating the bound and hence with \(e>E+R\).
At such a prime all the rebasing scalars are nonzero, as above.

For three terms and the stated even rule, at \(p=2\) one has
\(e\ge2\), \(R=0\), and the two remaining orders have valuations
zero and one. In each face polynomial the coefficients have
valuations on a line of slope \((v_2(q)-e)(g_2-g_1)\), plus
the nonnegative valuations of the binomial coefficients, with
zero binomial valuations at both endpoints. Every root has
valuation equal to the negative of that slope, by comparing the
constant or leading term. The two slopes differ because the
total face degrees differ. This rules out a common root.

Finally let \(F(L)\) be the homogeneous resultant of all stripped
moment forms. We have proved \(F(-1/j)\ne0\) under either criterion.
For any computable \(\nu\ge\deg F\), failure at an ordinary \(L\)
forces \(jL+1\mid j^\nu F(-1/j)\), by evaluation of the integer
polynomial modulo \(jL+1\). This gives an exact finite list by
testing its positive divisors. One can instead isolate the
nonnegative integer roots of \(F\) and impose the baseline
inequalities.

For example, in two variables with speeds \((1,1)\), take the
offset columns \((0,1),(1,0),(4,2)\) and intercepts \((0,1)\).
The deleted face varies in both coordinates and its total-degree
gap is five. Condition (13) holds at every eligible endpoint
and the marked decision works at every even start. At every
odd start \(n\ge3\), failure at both endpoints would force
\(n(n+2)\mid(2^5-1)^2\), which is impossible since \(2^5-1=31\).
At \(n=1\) the upper order \(3\) also does not divide
\(C_{\mathbf g,\rho}=2\cdot31^2\). This gives a finite baseline
list at every positive starting order. A one-coordinate pole test would leave an unmarked factorial
ratio in the deleted pair.
