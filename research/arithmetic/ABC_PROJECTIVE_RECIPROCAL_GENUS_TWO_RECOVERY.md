# Genus-two recovery for projective reciprocal relations

An affine equation between two distinct values of
\(F(z)=z-1/z+1/(1-z)\) on \(0<z\le1/2\) has a
coefficient-uniform subpower count. A cyclic cubic quotient
of genus two has all its branch points over the fixed field
\(\mathbb Q(\sqrt{-3})\). Its split \(1-\zeta_3\) descent
and uniform Mordell–Lang bound the number of points in terms
of the bad prime support of the affine coefficients. In a
higher relation on primitive triples, reserve the pair
carrying the largest exclusive prime packet for this recovery.
The resulting saving holds at both odd and even orders.

## Affine pair recovery over a fixed quadratic field

Write \(H(t)\) for the maximum of the absolute numerator
and positive denominator of a reduced rational \(t\).
For every \(K,\varepsilon>0\), uniformly in rational
\(r\ne0,s\) with \(H(r),H(s)\le X^K\), the number of pairs
of distinct rational numbers \(x,y\in(0,1/2]\) satisfying
\[
                         F(x)=rF(y)+s
 \tag{GR1}
\]
is \(O_{K,\varepsilon}(X^\varepsilon)\). No height
condition on \(x,y\) is required.

Put \(E=\mathbb Q(\sqrt{-3})\) and
\[
 \eta=\frac{1+\sqrt{-3}}2,\qquad
 \alpha=\frac{1+3\sqrt{-3}}2,\qquad
 U(z)=\frac{z-\eta}{z-\bar\eta},\qquad
 T(t)=\frac{t-\alpha}{t-\bar\alpha}.
\]
Complex conjugation is denoted by a bar. Differentiation gives
\[
                F'(z)=\frac{(z^2-z+1)^2}{z^2(1-z)^2}.
\]
The points \(\eta,\bar\eta\) have ramification index three
under \(F\), and their images are \(\alpha,\bar\alpha\).
Since \(F\) has degree three, the zero and pole of \(T\circ F\)
are exactly \(3\eta\) and \(3\bar\eta\). Evaluation at
\(z=\infty\) gives
\[
                              U(z)^3=T(F(z)).
 \tag{GR2}
\]
For a solution of (GR1) set \(t=F(x)\) and \(v=U(x)U(y)\).
It gives an \(E\)-point on the smooth projective model \(C\) of
\[
                   v^3=T(t)T\big((t-s)/r\big).
 \tag{GR3}
\]
For a prescribed \(t\) there are at most three choices of
each original argument, and hence at most nine pairs over
a point of \(C\).

The branch sets of the two factors on the right of (GR3) are
\(\{\alpha,\bar\alpha\}\) and
\(\{r\alpha+s,r\bar\alpha+s\}\). An intersection implies
equality because \(r,s\) are rational and the two points
in either set are conjugate. Comparison of real and imaginary
parts then gives exactly \((r,s)=(1,0)\) or \((-1,1)\).
The first implies \(x=y\) because \(F\) is strictly
increasing on \((0,1)\). The second gives \(x=1-y\)
by \(F(1-y)=1-F(y)\) and the same injectivity. In the
prescribed domain this requires \(x=y=1/2\). Thus neither
exception contributes distinct pairs.

In every other case (GR3) is a geometrically integral
cubic cover with four distinct totally ramified points.
The branch orders of its defining rational function are
\(1,-1,1,-1\), and it is a unit at infinity. Riemann–Hurwitz
gives \(2g(C)-2=3(-2)+4(3-1)=2\), so \(g(C)=2\).

## Descent with four rational branch points

Write \(r=p/q\), \(s=a/b\) in reduced form with positive
denominators and set \(d=qb\). After the change \(t'=dt\),
the four branch positions \(b_i\) are algebraic integers
in \(E\). They are \(d\) times
\(\alpha,\bar\alpha,r\alpha+s,r\bar\alpha+s\), which are
distinct under the hypotheses just checked. Put
\[
           M=3d\prod_{i<j}
                       |N_{E/\mathbb Q}(b_i-b_j)|.
\]
This is a nonzero rational integer with \(M\ll_K X^{O(K)}\)
when \(H(r),H(s)\le X^K\), since the branch differences in
both embeddings are bounded by a fixed power of these
heights. Let \(S\) be the finite places of \(E\) over \(M\).
There are at most \(2\omega(M)\) of them.

At a finite place outside \(S\) the four branch sections
are disjoint and the residue characteristic is prime to
three. The curve has good reduction there. To see this
on the normalized cover over the valuation ring, away from
the branch sections the defining function is a unit and
the cubic cover is étale by the derivative in \(v\).
The same holds at infinity, where the leading quotient
is one. On a neighborhood of a simple zero use \(v\);
near a simple pole use \(v^{-1}\). The local equation
at the section becomes \(z^3=(t'-b_i)U_i(t')\), with
\(U_i\) a unit. Its derivative in \(t'\) is a unit at
\(t'=b_i,z=0\), so this chart is smooth there as well.
These finite local covers agree on overlaps and give a
smooth proper model by gluing over the projective line.
The deck action extends on all of the charts.

Let \(J\) be the Jacobian and let \(\sigma\) be its
endomorphism induced by \(v\mapsto\zeta_3 v\).
Put \(\phi=1-\sigma\). The norm \(1+\sigma+\sigma^2\)
factors through the Jacobian of the base line and vanishes.
Consequently \(\phi^2=-3\sigma\). In particular \(\phi\)
is an isogeny of degree \(9\).

The branch points \(P_i\) on \(C\) all lie in \(C(E)\):
above an \(E\)-rational branch location the ramification
index is three, exhausting the degree of the cover.
Their degree-zero divisor classes are killed by \(\phi\)
and by three. The latter follows by using the pullback of
\((t'-b_i)/(t'-b_j)\). The space of coefficient vectors
\((n_1,\ldots,n_4)\) modulo three with sum zero has
dimension three. Choose lifts of integer sum zero, which
give a well-defined map to \(J[\phi]\) by these divisors.
If the class of \(\sum_i n_iP_i\) is trivial, any function
with that divisor is a deck eigenfunction. In the cyclic
function field it has the form \(v^h R(t')\), with \(R\)
rational in \(t'\) over \(\overline E\). Its orders on
the branch points modulo three are thus proportional to
\((1,-1,1,-1)\); the function \(R(t')\) has orders
divisible by three there. Hence the kernel has dimension
at most one. The branch classes give at least nine
elements of \(J[\phi](E)\), so this is all of the
geometric kernel. The kernel is therefore constant
over \(E\). This is the explicit branch-divisor
description for these four split points.
[Poonen–Schaefer, Proposition 6.2 and Lemma 12.9](https://math.mit.edu/~poonen/papers/descent.pdf)

The Kummer injection
\[
 J(E)/\phi J(E)
       \hookrightarrow H^1(E,J[\phi])
       \cong (E^\times/E^{\times3})^2
\]
lands in classes unramified outside \(S\). Indeed the
Jacobian has good reduction there and \(\phi\) extends,
since the deck action does. Its square is multiplication
by three times an automorphism, so its derivative is an
isomorphism outside residue characteristic three. On
the proper abelian scheme the isogeny is finite étale.
A rational Jacobian point extends to the valuation ring
and its inverse image under \(\phi\) is finite étale,
giving the claimed unramified class. In the cube-class
description the valuations are divisible by three at
those places. The number of such classes is at most
\(3^{O(1+|S|)}\), by the \(S\)-unit rank and the fixed
class group of \(E\). For example \(E\) has class
number one, since its Eisenstein integer ring is
norm-Euclidean.
[Descent and support accounting](ABC_PRIME_POWER_FERMAT_DESCENT.md#uniform-rational-point-bound-from-the-splitting-field)

Since \(\phi^2\) is multiplication by three up to a unit,
\[
       3^{\operatorname{rank}J(E)}
              \ \le\ [J(E):3J(E)]\
              \le\ [J(E):\phi J(E)]^2.
\]
This bounds the rank by \(O(1+\omega(M))\). Embed \(C\)
in \(J\) using a branch point as base point. Its degree
for a fixed multiple of the theta polarization depends
only on the genus, and a genus-two Abel–Jacobi curve
contains no positive-dimensional abelian coset.
Uniform Mordell–Lang therefore gives
\[
               |C(E)|\le \exp(O(1+\omega(M))).
\]
The absolute constant in this bound comes from the fixed
genus, polarization degree and ground field. Since
\(M\) has polynomial height, this is \(O_{K,\varepsilon}
(X^\varepsilon)\) for every \(\varepsilon>0\).
Together with the bounded fiber of the map to \(C\) this
proves (GR1).
[Gao–Ge–Kühne, Theorem 1.1′](https://arxiv.org/html/2105.15085v4)

## A saving at every order of distinct active rows

Fix \(B,K,\gamma>0\) and an integer \(k\ge3\). Let the
\(k\) original rows be distinct primitive positive
\(a_i+b_i=c_i\) with \(a_i\le b_i\), \(X/2<c_i\le X\)
and \(\operatorname{rad}(a_ib_ic_i)\le n^B\), where
\(n=X^L\). For sufficiently small fixed \(L>0\), the
number of solutions of
\[
                      \sum_{i=1}^k w_i F(a_i/c_i)=t
\]
is
\[
     O_{B,K,k,\gamma,L,\varepsilon}
        \!\left(X^{\{kB/6-B/[3(k-1)]+\gamma\}L+\varepsilon}\right)
 \tag{GR4}
\]
for any fixed nonzero rational weights of reduced height
at most \(X^K\), uniformly in those weights and the
right-hand side \(t\). For \(k=2\) the count on distinct
rows is \(O_{K,\varepsilon}(X^\varepsilon)\) by (GR1).
For \(k=1\) a prescribed value determines the row.
Ordering the summands removes the reflection exception
while keeping one row from each unordered pair.

Here are the prime prices giving (GR4). If a solution
exists, the right-hand side has polynomial height in
\(X\) with exponent depending only on \(k,K\). Put its
denominator primes and the numerator and denominator
primes of the weights in a pool \(P_0\). For the
distinct active rows write \(b_i=\log_n\operatorname{rad}
(a_ib_ic_i)\) and \(D=kB-\sum_i b_i\). Let \(S\) be
the incidence surplus of (RS5), and \(q_{ij}\) the
normalized mass of primes outside \(P_0\) occurring on
just rows \(i,j\). No prime outside \(P_0\) occurs
on exactly one row, by the pole-order test in (RS1).
The pool incidence and the incidence from primes on at
least three rows is at most \(3S\), because
\(m_p\le3(m_p-2)\) for \(m_p\ge3\). Hence, writing
\(q=\max_{i<j}q_{ij}\),
\[
                 kB\le D+3S+2\binom{k}{2}q .
 \tag{GR5}
\]
Place a maximizing pair last in the active order. Every
prime outside \(P_0\) occurring for the first time at
the penultimate step is exclusive to that final pair.
The last step has no new prime. Therefore the new
radical masses \(v_i\) on the preceding \(k-2\) rows sum
to \((kB-D-S)/2-q\), by the exact prime incidence
identity of (RS5).

Count those preceding rows by (RS2) in the chosen order.
Bin each new mass on a fixed finite mesh and apply the
effective radical count after fixing the earlier rows.
The pool supported on the earlier rows and \(P_0\) has
bounded polynomial height. After these rows have been
exposed, the remaining two satisfy (GR1), with
coefficients of bounded polynomial height in \(X\).
They cost \(X^\varepsilon\) with an arbitrary exponent
allowance. Altogether the normalized exponent is at
most \((kB-D-S-2q)/6\) before the fixed mesh and counting
allowances. For \(k\ge3\) one has
\(\binom{k}{2}\ge3\), so (GR5) implies
\[
                    D+S+2q\ \ge\ \frac{2B}{k-1}.
\]
The finitely many choices of the maximizing pair and a
sufficiently fine fixed mesh give (GR4) with the stated
allowances. The budget is fixed before the height limit
at every stage.
[Prime occurrence and effective count](ABC_THREE_CAP_RATIO_SUMSETS.md#a-pointwise-bound-when-the-denominator-records-all-three-letters)
· [The incidence surplus identity](ABC_THREE_CAP_RATIO_SUMSETS.md#radical-surplus-and-separation-by-letter)

For the \(2m\)-term additive energy of \(F\) on a family
of sorted rows, cancel equal terms in opposite lists and
group repeated rows within each list. Nontrivial
solutions have \(k\ge2\) distinct active rows whose
merged weights are nonzero. A cancelled group costs at
most \(n^{B/3+\gamma}\) up to height allowances by (RS2)
with an empty pool, and each uses at least two original
slots. Formula (GR4), with the \(k=2\) pair recovery
when needed, bounds all nontrivial solutions by
\[
 O_{B,m,\gamma,L,\varepsilon}\!
     \left(X^{\{mB/3-B/[3(2m-1)]+\gamma\}L+\varepsilon}\right)
                          \qquad(m\ge2).
\]
Here the nontrivial solutions have unequal row
multisets on the two sides. The argument gives a
power gap for exact relations of either parity. Its
input is a relation with polynomial-height weights,
so using it for the unrestricted \(abc\) count still
requires a mechanism producing sufficiently many such
relations in a near-full family.

## Scalar relations on powerful representations of a fixed target

Fix \(C,K,\gamma>0\) and \(k\ge3\). Let \(p,q\ge2\) be
fixed with \(L=1/p+1/q\) sufficiently small, put \(n=N^L\),
and suppose \(\operatorname{rad}(N)\le n^C\). Consider
the ordered positive representations \(A+B=N\) with
\(A\le B\), where \(A\) is \(p\)-full and \(B\) is \(q\)-full.
There is no coprimality restriction. With fixed nonzero
rational weights of reduced height at most \(N^K\), the
number of distinct-row solutions to
\[
                        \sum_{i=1}^k w_i F(A_i/N)=t
\]
is
\[
 O_{C,K,k,\gamma,p,q,\varepsilon}
       \!\left(N^{\{k/8-1/[4(k-1)]+\gamma\}L+\varepsilon}\right)
 \tag{GR6}
\]
uniformly in the weights and target \(t\). At order two
the count is subpower even without fullness. Distinct
representations of fixed \(N\) have distinct fractions
\(A_i/N\), so (GR1) applies on their actual common target.

Fix \(g_i=\gcd(A_i,B_i)\mid N\) at divisor cost and set
\(W_i=N/g_i,\ U_i=A_i/g_i,\ V_i=B_i/g_i\). Outside the
prime support of \(N\), the factors on \(U_i,V_i\) are
still \(p\)-full and \(q\)-full, respectively. The reduced
denominator of \(F(A_i/N)=F(U_i/W_i)\) is \(U_iV_iW_i\).
Outside \(N\) this has exactly the prime support of the
two moving letters. Thus put the primes of \(N\) in the
initial pool, along with those of the weights and the
inherited rational target as in (GR4). The raw budget
outside the pool is at most one in base \(n\), since the
two outside radicals are at most \(W_i^{1/p}\) and
\(W_i^{1/q}\).

First suppose all \(W_i\ge N^{1/8}\). In the prime exposure
order, for assigned new outside masses \(d_U,d_V\) in
base \(n\), apply the fixed-letter estimate (FL1) at
height \(W_i\). The budget there is \(8L\), with caps
\(d_U\log N/(8\log W_i)\) and \(d_V\log N/(8\log W_i)\).
They are uniformly bounded and
\(\operatorname{rad}(W_i)\le W_i^{8CL}\), so the estimate
applies. The normalized cost in base \(n\) is
\(\min(d_U,d_V)/2\), up to an arbitrarily small fixed
allowance. In particular the cost is at most
\((d_U+d_V)/4\). Put a pair with largest exclusive
packet last. The prefix sum of new masses is
\((k-D-S)/2-q\) as in (GR5), with \(B=1\) applied to
the outside-pool incidences. The last pair costs subpower
by (GR1). Thus the normalized exponent is at most
\((k-D-S-2q)/8\) before allowances, which proves
(GR6) because \(D+S+2q\ge2/(k-1)\).

If some \(W_i<N^{1/8}\), expose those rows first.
For each, fix the two \(N\)-supported parts at
height-subpower cost. Enumerate the smaller outside
radical dictionary and then use the fixed \(W_i\) to
determine the other letter. The cost per such row is at
most \(N^{\min(1/p,1/q)/8+\varepsilon'}\le n^{1/16}
N^{\varepsilon'}\). The residual rational target has
polynomial height. If at least three tall-quotient rows
remain, apply their bound above. The saving relative
to their diagonal exponent is at least \(1/[4(k-1)]\),
and each exposed short-quotient row cost less than
\(1/8\) in normalized exponent. If at most two rows
remain, their count is subpower from (GR1) or injectivity.
Even \(n^{k/16}\) fits the exponent in (GR6) for \(k\ge3\).
All choices of the short positions are finite at fixed
\(k\). This proves the bound with the same short-quotient
split and prime prices as (JR8), using scalar recovery
instead of the coordinatewise reciprocal equations.
[Fixed-letter estimate with prime pools](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Joint reciprocal fixed-target proof](ABC_JOINT_RECIPROCAL_PAIR_RECOVERY.md#powerful-representations-of-one-target-with-arbitrary-gcd)

For the \(2m\)-row scalar energy, \(m\ge2\), the same
cancellation step bounds the contribution with unequal
row multisets by
\[
 O_{C,m,\gamma,p,q,\varepsilon}
     \!\left(N^{\{m/4-1/[4(2m-1)]+\gamma\}L+\varepsilon}\right).
\]
A single row on the original target costs at most
\(n^{1/4}\) before allowances, by the same tall and short
quotient counts. Each cancelled group uses at least two
slots, and for the nonzero groups (GR6) gives at least
the displayed gap. Applying (GR1) directly handles the
case of two distinct active rows.

## Relation to bounds for pencils of curves

Dimitrov, Gao and Habegger derive a subpolynomial bound
for \(y^2=(x-s)Q(x)\) with fixed squarefree \(Q\) and a
moving rational parameter \(s\). Their Corollary 6.1 and
its proof use a fixed splitting field for the
Weierstrass-point torsion, the support of the moving
discriminant, and a rank-dependent pencil bound. This is
the same type of support-to-rank mechanism. Here (GR3)
has a cyclic cubic deck transformation and a split
\(1-\zeta_3\) kernel over the fixed quadratic field. The
two affine parameters both vary, so we use uniform
Mordell–Lang for all curves of the relevant genus and
degree.
[Dimitrov–Gao–Habegger, Corollaries 1.3 and 6.1](https://arxiv.org/pdf/1904.07268v2)

## Fixed product recovery for the same projective coordinate

For every \(K,\varepsilon>0\), a nonzero rational \(h\) of
height at most \(X^K\) has \(O_{K,\varepsilon}(X^\varepsilon)\)
ordered representations as
\[
                         h=F(x)F(y),\qquad
                        x,y\in(0,1/2]\cap\mathbb Q .           \tag{GR7}
\]
Use \(t=F(x)\) in the same split cubic construction. The curve
is now \(v^3=T(t)T(h/t)\), with branch pairs
\(\{\alpha,\bar\alpha\}\) and
\(\{h/\alpha,h/\bar\alpha\}\). If they are disjoint, the
good reduction and descent proof of (GR1) applies. The
additional constant in the four-factor quotient is the
unit-at-infinity value \(T(0)=\alpha/\bar\alpha\); include
its numerator and denominator primes among the excluded
places. The branch positions still have polynomial height.
If the branch sets intersect they agree by conjugation.
Since \(\alpha^2\notin\mathbb Q\) and
\(\alpha\bar\alpha=7\), this requires \(h=7\).
In that case
\[
                   T(t)T(7/t)=\alpha/\bar\alpha .
\]
This constant cannot be a cube in \(E\). The ideal generated
by the algebraic integer \(\alpha\) has norm \(7\), and
its conjugate is a different prime ideal. Indeed a common
prime would divide \(\alpha+\bar\alpha=1\). The constant
therefore has valuations \(1,-1\) at the two primes over
\(7\). Equation (GR2) excludes every pair in this case.

The function \(F\) does not vanish at a rational argument:
its numerator equation at zero is \(z^3-z^2-2z+1=0\),
whose only possible rational roots are \(1,-1\), neither
a root. For any set \(\mathcal A\) of fractions
\(a/c\in(0,1/2]\) of height at most \(X\), every product from
\(F(\mathcal A)\) has polynomial height. Thus (GR7) gives
\[
                  |F(\mathcal A)\,F(\mathcal A)|
                             \gg_\varepsilon
                    |\mathcal A|^2 X^{-\varepsilon}.
\]
Together with (GR1), the same bound holds for two-fold sums.
There is no dependence on the radical budget in these two
pair conclusions. Sequential pricing in (GR4) depends on
the unique negative valuation in a linear equation. An
isolated denominator prime in a product may cancel against
a numerator prime, so the product bound by itself does
not give the corresponding higher-order radical pricing.


## Relation to the simplest cubic family

Clearing denominators in \(F(z)=t\) and setting \(u=-z\)
gives \(u^3+(t+1)u^2+(t-2)u-1=0\). This is Shanks'
simplest cubic polynomial at parameter \(-(t+1)\).
The cyclic symmetry used for (GR2) is thus the
simplest-cubic symmetry. The field-isomorphism and Thue
equation literature addresses coincidences of primitive
cyclic cubic fields for parameter values, including
irreducible specializations. The affine recovery here
counts pairs of completely split rational specializations
on a prescribed affine relation, with uniform dependence
on the coefficient height. The genus-two cover is handled
through its four split branch points and the
rank-dependent bound from uniform Mordell–Lang.
[Hoshi, simplest cubic fields and cubic Thue equations](https://arxiv.org/abs/0810.3374)

