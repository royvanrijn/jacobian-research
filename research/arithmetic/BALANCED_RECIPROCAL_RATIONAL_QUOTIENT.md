# A rational quotient for every balanced reciprocal sign class

The square and nonsquare reciprocal-parameter cases can be treated by one
explicit scheme over \(\mathbb Q\). No choice among quadratic fields is
needed in its defining equations.

For degree \(n=8k+2\), its raw scheme length is

\[
                         2k\binom{4k}{k}^{\,2}.
\]

In degree 26 this is a finite scheme of length 290400, defined by seven
rational equations in seven variables. Thirteen stronger equations cut
out the necessary square-divisibility locus. Every normalized rational
good balanced reciprocal model corresponds to exactly one rational point
of its good open subset, including models whose sign factors are defined
only over a quadratic field.

This gives a finite complete rational-candidate procedure. It does not
supply the degree-26 rational point list, and does not prove that this
list is empty.

## The labelled involution

Fix \(k\ge1\), put \(h=2k\), \(L=8k+1\), and use the
[raw equations](RECIPROCAL_PARAMETER_FINITENESS_AND_INTEGRALITY.md)
for the balanced split \(p=q=k\):

\[
 P=w^k+\sum_{i=0}^{k-1}p_iw^i,\quad
 Q=w^k+\sum_{i=0}^{k-1}q_iw^i,\qquad
 \operatorname{rem}_P X=\operatorname{rem}_Q Y=G=0,
 \quad G=F(r)/(2r).
\]

Let \(A\) be their coordinate algebra over \(\mathbb Q\). It is finite,
of length \(4k\binom{4k}{k}^2\). Since \(G|_{r=0}=L\), write
\(G=L+rH\). In \(A\), \(r(-H/L)=1\); thus \(r\) is a unit.

Changing the sign of the reciprocal square root acts by

\[
 \tau(r)=-r,\qquad
 \tau(p_i)=s_iq_i,\qquad \tau(q_i)=s_ip_i,\qquad
 s_i=(-1)^{k+i}.                                         \tag{1}
\]

Equivalently, \(\tau P(w)=(-1)^kQ(-w)\) and
\(\tau Q(w)=(-1)^kP(-w)\). Direct reconstruction gives

\[
 \tau F(v)=-F(-v),\qquad
 \tau X(w)=Y(-w),\qquad \tau Y(w)=X(-w),\qquad \tau G=G.
                                                               \tag{2}
\]

Thus \(\tau\) preserves both the raw ideal and the stronger ideal
\(P^2\mid X,\ Q^2\mid Y,\ G=0\).

## Theorem 1. Explicit quotient coordinates and equations

Define invariant coordinates

\[
 x_i=p_i+s_iq_i,\qquad
 z_i=\frac{p_i-s_iq_i}{r},\qquad
 \eta=r^2.                                                \tag{3}
\]

Their inverse formulas are

\[
 p_i=\frac{x_i+rz_i}{2},\qquad
 q_i=s_i\frac{x_i-rz_i}{2}.                               \tag{4}
\]

Substitute (4) into the coefficient of \(w^j\) of
\(\operatorname{rem}_P X\), then reduce powers of \(r\) using
\(r^2=\eta\). Write the result as

\[
                         E_j(x,z,\eta)+rO_j(x,z,\eta),
                         \qquad 0\le j<k.               \tag{5}
\]

The corresponding minus remainder coefficient is, by (2),

\[
                         (-1)^j(E_j-rO_j).              \tag{6}
\]

The same substitution in \(G\) gives a polynomial
\(\overline G(x,z,\eta)\) with no odd part. In particular,

\[
                         \overline G|_{\eta=0}=L.        \tag{7}
\]

The invariant algebra \(B=A^\tau\) therefore has the explicit presentation

\[
 B=
 \frac{\mathbb Q[x_0,\ldots,x_{k-1},z_0,\ldots,z_{k-1},\eta]}
 {(E_0,O_0,\ldots,E_{k-1},O_{k-1},\overline G)}.            \tag{8}
\]

There are \(2k+1\) equations and variables. No variable factor is divided
out. Equation (7) already makes \(\eta\) invertible in this algebra.

To verify the presentation, invert \(r\) and apply the linear change
(3)--(4) to the ambient polynomial algebra. Equations (5)--(6) are
equivalent to \(E_j=O_j=0\), since 2 and \(r\) are units. Together
with \(r^2=\eta\) they give the original raw equations exactly.
Equation (7) removes the need to retain an explicit localization
at \(\eta\).

There is a canonical isomorphism

\[
                    A=B[r]/(r^2-\eta).                  \tag{9}
\]

In particular the quotient is finite and its length is

\[
                       \dim_{\mathbb Q}B
                          =2k\binom{4k}{k}^{\,2}.        \tag{10}
\]

These assertions hold even with nilpotents and degenerate raw points.
Indeed, every \(a\in A\) has the unique decomposition

\[
 a=\frac{a+\tau a}{2}
      +r\,\frac{a-\tau a}{2r},
\]

with both coefficients in \(B\). Thus \(A=B\oplus rB\) as a
\(B\)-module. The element \(\eta=r^2\) is a unit, so the discriminant
\(4\eta\) of (9) is a unit as well: this is a quadratic finite
etale cover, with the involution free. The length division by two
does not presume a reduced scheme.

The stronger quotient is obtained by doing the same even/odd splitting
for the \(2k\) coefficients of \(\operatorname{rem}_{P^2}X\).
It has \(4k+1\) equations, including \(\overline G\), and is a
closed subscheme of (8). Its length can be smaller; (10) is only
the raw length.

## Theorem 2. Rational points include every quadratic descent

For a rational point of the quotient with \(\eta\ne0\), take
\(r\) in the quadratic algebra \(K=\mathbb Q[r]/(r^2-\eta)\) and
reconstruct \(P,Q,F\) by (4). When \(\eta\) is square, this algebra
splits; when it is nonsquare, it is the corresponding quadratic field.
In either case set

\[
                     N(u)=r^{-(4k+1)}F(ru).             \tag{11}
\]

Its coefficients are invariant under \(r\mapsto-r\), by (2),
and hence lie in \(\mathbb Q\). More explicitly, write
\(F(v)=\sum_j f_jv^j\). After the substitution (4),

\[
 f_j=
 \begin{cases}
 A_j(x,z,\eta), & j\text{ odd},\\
 rB_j(x,z,\eta), & j\text{ even}.
 \end{cases}
\]

Thus the coefficient of \(u^j\) in \(N\) is respectively

\[
 \eta^{(j-4k-1)/2}A_j
 \quad\text{or}\quad
 \eta^{(j-4k)/2}B_j.                                      \tag{12}
\]

Only powers of \(\eta\) occur in denominators. The reconstruction is
monic, with

\[
 N(0)=L\eta^{-2k},\qquad
 N(1)=2\eta^{-2k}\overline G,\qquad
 N'(1)=\eta^{-2k}\overline G.                             \tag{13}
\]

Restrict to the good open conditions in the
[mixed note](MIXED_RECIPROCAL_COUPLED_ODE.md):
\(P Q\) squarefree, no \(w=\pm2\) or fixed-zero collision, no vanishing
opposite sign factor, and distinct nonzero squared critical values.
Then its proved converse makes (11) a normalized good reciprocal model.

Conversely, every normalized rational good model in this balanced
class has its parameter \(\lambda\in\mathbb Q^*\). Take
\(\rho^2=\lambda,\ r=1/\rho\). If \(\lambda\) is square in
\(\mathbb Q\), the parameters \(r,P,Q\), and hence (3), are rational
directly. If \(\lambda\) is nonsquare, the plus and minus factors
lie in the quadratic field \(\mathbb Q(r)\) and are exchanged by
(1) under its nontrivial conjugation. The invariant coordinates (3)
are therefore rational in this case as well.
The choices \(r,-r\) give the same quotient point and the same \(N\).

The equal-value pairs and their aggregate signs are uniquely determined
on the good locus. Consequently the good rational quotient points
are in bijection with normalized rational balanced reciprocal numerators
\(N\). The relation between the arithmetic parameters is

\[
                            \lambda=\eta^{-1}.           \tag{14}
\]

In particular their global and local square classes agree.

The good open is explicitly algebraic. Its defining nonvanishing
conditions can be multiplied, denominators cleared using the unit \(r\),
and then multiplied by their \(\tau\)-conjugates. The resulting invariant
polynomial \(\Delta(x,z,\eta)\) is nonzero exactly on that open.
Thus no implicit choice of a quadratic field remains in the point search.

## A finite complete rational-candidate procedure

For each \(k\), the following exact procedure terminates:

1. Construct (8) by (1)--(7), or use the stronger quotient equations.
   Form its finite coordinate algebra and localize at \(\Delta\);
   equivalently, saturate the finite ideal by powers of \(\Delta\).
2. Remove nilpotents. The resulting reduced finite algebra over
   \(\mathbb Q\) is a product of number fields.
3. Enumerate integer linear combinations of its coordinate generators
   until a separating element \(\theta\) has power span of full dimension.
   Such an element exists: the combinations failing to separate two
   distinct geometric points lie in a finite union of proper hyperplanes.
4. Compute the squarefree minimal polynomial of \(\theta\).
   Its rational roots, found by exact rational factorization, are
   precisely the rational points. Recover every coordinate as a polynomial
   in \(\theta\), and reconstruct \(N\) by (12).

Each step is finite exact algebra. In particular the quadratic twists
are processed together by one rational scheme, rather than by separate
point searches over every possible quadratic field.

This is a complete procedure, not a claim that its largest exact
Groebner and factorization steps have been executed. The raw length is
an upper bound for the dimension of the reduced good algebra and
the degree of its separating polynomial. It is not a practical
runtime bound.

## Degree twenty-six

For \(k=3\), \(s_i=(-1,+1,-1)\). There are seven invariant variables,
seven raw equations and thirteen strong equations. Their full exact
coefficients and all numerator reconstruction formulas are retained in

[degree26_balanced_rational_quotient_system.json](degree26_balanced_rational_quotient_system.json).

The raw quotient length is 290400. The degree-26 rational-candidate
procedure above is therefore a single explicit finite problem.
It is not 1024 unrelated quadratic-field problems.

The uniform denominator theorem gives every \(x_i,z_i\) as an
\(S\)-integer and \(\eta\) as an \(S\)-unit, where \(S\) consists
of primes at most 25. The nine primes are
\(2,3,5,7,11,13,17,19,23\). Thus at most 1024 signed square classes
of \(\eta\) can occur. The local-square
[obstruction at 37](DEGREE_TWENTY_SIX_BALANCED_LOCAL_SQUARE_OBSTRUCTION.md)
requires the class to be nonsquare modulo 37. These are valid arithmetic
filters on (8), rather than substitutes for its rational-point computation.

A good rational point would give the desired rational paired-critical
degree-26 pencil. The existing passport and alternating-monodromy
arguments then apply. No such point is constructed here and no
nonexistence theorem for the nonsquare balanced case is proved.

## Replay and small exact control

~~~text
.venv/bin/python research/arithmetic/verify_balanced_reciprocal_rational_quotient.py --write-receipt research/arithmetic/balanced_reciprocal_rational_quotient_verification.json --write-degree26-system research/arithmetic/degree26_balanced_rational_quotient_system.json
~~~

The checker verifies the involution, every even/odd equation pairing,
and all numerator reconstruction identities in degrees 10, 18 and 26.
It writes both rational equations and primitive integer equations,
without discarding variable factors.

As a small complete control, it computes the strong \(k=1\) quotient
over \(\mathbb Q\). Its length is 18, and its \(\eta\)-eliminant is
exactly the degree-18 polynomial retained in the mixed degree-ten
classification. Thus the universal construction recovers the existing
projective modulo-11 arithmetic obstruction with \(\eta=r^2\).

The uniform quotient theorem is proved above; bounded symbolic checks
are supplementary. No novelty comparison, formal verification or
external human review is asserted.
