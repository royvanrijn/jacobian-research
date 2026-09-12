# The Keller monodromy action spectrum: the \(\operatorname{PSL}_2(11)\) benchmark

## 1. Outcome and status

The inverse-Galois programme should be graded by permutation actions and by
the geometry of their cover presentations, not by an abstract list of groups.
This note makes that refinement concrete for

\[
 G=\operatorname{PSL}_2(\mathbb F_{11}).
\]

One group already supplies three materially different primitive actions:

1. the natural degree-twelve action on
   \(\mathbb P^1(\mathbb F_{11})\);
2. the degree-eleven action on the cosets of one conjugacy class of
   \(A_5<G\); and
3. the degree-eleven action on the cosets of the other conjugacy class of
   \(A_5<G\).

The two degree-eleven actions are nonisomorphic Gassmann twins.  Their
classical \((3,2,11)\) covers are genus-zero Shabat polynomials over
\(K=\mathbb Q(\sqrt{-11})\), and one derivative-unit coordinate gives two
determinant-one degree-eleven morphisms of smooth affine surfaces to the same
target.  This is a new relative/chart Keller realization of a Gassmann pair,
different from Fano point/line duality.

The natural degree-twelve action applied to the same regular
\((3,2,11)\)-cover has genus one.  Its rigid cover therefore has a stable
birational obstruction before any boundary-unit calculation begins.  The
degree-eleven charts are rational, but their source and target unit ranks are
six and two, and the derivative represents a nonzero class modulo the target
pullback lattice.  Thus the two actions of the same group teach different
affine-completion lessons.

The neighboring natural actions sharpen that point.  For every prime
\(p>3\), the standard modular \((2,3,p)\) triple has natural quotient genus
\[
 \bigl[p-6-3(-1\mid p)-4(-3\mid p)\bigr]/12.
\]
The values at \(p=7,11,13\) are \(0,1,0\): the elliptic obstruction at
\(11\) is a split/nonsplit arithmetic effect, not a monotone consequence of
larger degree.

Normalizing the two direct correspondences gives a further separation.  The
bidegree-\((5,5)\), total-degree-five component normalizes to the quotient by
an \(A_4\) intersection subgroup and has genus one.  The bidegree-\((6,6)\),
total-degree-six component normalizes to the quotient by a \(D_{10}\)
intersection subgroup and has genus two.  Their affine plane models have,
respectively, five and eight ordinary nodes, and the reduced node ideal is
the conductor.  The genus-one normalization has a \(K\)-point and is the
base change of the conductor-\(121\) curve
\(v^2+uv=u^3+u^2-2u-7\), with \(j=-121\).  Exact traces above \(23\) show
that it is not \(K\)-isogenous to the natural \(X_0(11)\) quotient.  In
particular neither correspondence hides a rational lower-color bridge
between the two degree-eleven charts.

The normalization step now goes all the way to the exact boundary-unit
lattices.  On \(X_0(11)\) the lattice has rank three and the evident
\(j,j-1728,f_T\) units generate an index-six sublattice.  On the genus-one
correspondence it is the rank-fourteen degree kernel.  On the genus-two
correspondence the two elliptic quotient traces give a rank-seventeen
lattice with quotient \(\mathbb Z^3\oplus\mathbb Z/5\); the cyclic cubic
Galois action on the \((2,2)\)-kernel proves that no rational two-torsion
class was lost.  These are presentation invariants with arithmetic content,
not merely genus calculations.

The two correspondence projections now have exact compact boundary
pullback matrices as well.  In both degrees their rank-six unit images meet
in exactly the two common triangle-base units, and their rank-ten sum is
primitive.  The remaining unit cokernels are free of ranks four and seven.
The minimal ledger containing the common base and the two pulled-back
derivative classes is also primitive, but two additional mask characters
cannot complete either full source-unit lattice.  This is an exact
obstruction to that declared two-mask architecture, not to arbitrary
nonlinear modifications of two polynomial outputs.

The residual symmetry is now explicit.  Projection exchange acts trivially
on the rank-four quotient and as five fixed characters plus one exchanged
pair on the rank-seven quotient.  Effective regular-mask bases exist.  On
\(C_5\), simple-pole masks have index two, so one double-pole mask is
unavoidable; on \(C_6\), five fixed mask classes and one exchanged pair give
an integral basis.  Every one of these masks is now constructed exactly in
the normalization algebra.  A new infinity-imbalance character on the
rank-seven quotient proves that the asymmetric pole profiles on three of
the \(C_6\) classes are intrinsic, not artifacts of the chosen lift.
Assigning one zero-section-preserving normal monomial to each class pins the
first factor-rich supports in normal degrees two and three.  These are
coefficient supports with explicit rational functions on the normalized
curves, not constructed polynomial Keller maps.

No absolute polynomial Keller map with either generic action is constructed
here.  The exact checker
[`verify_psl2_11_keller_action_spectrum.py`](../scripts/verify_psl2_11_keller_action_spectrum.py)
separates the symbolic polynomial calculation from the finite-group and
Riemann--Hurwitz calculation.  The identification of the displayed
polynomials with the two \(\operatorname{PSL}_2(11)\) dessins and their common
regular cover is the external Jones--Zvonkin input.

## 2. The spectrum to classify

Fix a characteristic-zero field \(k\).  Let
\(\mathfrak K_{\mathrm{abs}}(k)\) be the set of permutation-isomorphism
classes of faithful transitive actions \((G,\Omega)\) for which there are an
integer \(m\) and a noninvertible polynomial Keller map

\[
 F:\mathbb A^m_k\longrightarrow\mathbb A^m_k
\]

whose geometric generic inverse monodromy is \(G\curvearrowright\Omega\).
The arithmetic generic group is extra data and must be recorded separately.
Define \(\mathfrak K_{\mathrm{chart}}(k)\) similarly, allowing a
determinant-one finite etale morphism between smooth affine charts rather
than requiring both sides to be affine space.

This gives the open classification problem suggested by the current
examples.

> **Keller monodromy action-spectrum problem.**  Classify
> \(\mathfrak K_{\mathrm{abs}}(k)\), dimension by dimension and after stable
> polynomial equivalence.  Determine which action-theoretic, Hurwitz,
> birational, and boundary-divisor invariants obstruct promotion from
> \(\mathfrak K_{\mathrm{chart}}(k)\) to
> \(\mathfrak K_{\mathrm{abs}}(k)\).

The word *action* is indispensable.  If \(G\curvearrowright G/H\), then
regularity, the block interval \([H,G]\), the permutation character,
orientation, inertia cycle types, and the quotient geometry all depend on
\(H\), not only on \(G\).

There should also be a versal refinement
\(\mathfrak K_{\mathrm{versal}}(k)\).  It asks for a Keller family which is
versal for \(G\)-torsors, rather than for one cover whose function-field group
is \(G\).  Essential dimension constrains this versal spectrum.  It does not
give a lower bound on the parameter count of a single rigid cover.

## 3. What is already invariant

The first useful action card has four layers.

| layer | data | present force |
|---|---|---|
| permutation | \((G,H)\), degree, core, block interval, normalizer, permutation character | a nontrivial regular action is absolutely impossible; primitive monodromy makes a realizing map atomic |
| Hurwitz | inertia classes, Nielsen components, cycle indices, source genus, field of rationality | selects a cover input and can obstruct stable rationality of that input |
| affine boundary | source/target unit lattices, pullback matrix, derivative vector, class groups, log-canonical data | controls derivative suspension and affine-space completion |
| moduli | Hurwitz dimension versus \(\operatorname{ed}_k(G)\) and versality | separates one existential cover from a universal realization |

Two points are theorems rather than heuristics.

First, the Campbell--Razar--Wright Galois-case theorem gives

\[
 (G,G/1)\notin\mathfrak K_{\mathrm{abs}}(k)
 \qquad(G\ne1).                                      \tag{3.1}
\]

This is an obstruction to every nontrivial regular action, not only to cyclic
groups.  It says nothing against a core-free nonnormal stabilizer.

Second, if an absolute realization factors as two nonunits, its monodromy
action is imprimitive.  A primitive realization is therefore atomic.  The
converse fails without the polynomial-sandwich algebraization of an
intermediate field, so imprimitivity is not by itself a factorization
theorem.

The Hurwitz and boundary rows are presentation-dependent.  They become
action-level invariants only after minimizing over all covers and all Keller
presentations of the action.  That minimization is part of the proposed
classification, not a theorem already available.

## 4. Natural \(\operatorname{PSL}_2(q)\) actions

For \(q\ge4\), the natural action

\[
 \operatorname{PSL}_2(\mathbb F_q)
 \curvearrowright\mathbb P^1(\mathbb F_q)             \tag{4.1}
\]

has degree \(q+1\), is doubly transitive and primitive, and has a nontrivial
core-free Borel stabilizer.  It therefore passes the regular-action gate and
is an atomic absolute target.  Since the nonabelian simple group has no
nontrivial sign character, its natural permutation image is contained in
\(A_{q+1}\); any primitive polynomial or rational presentation must have
square discriminant over its geometric constant field.

The action also has a uniform cycle dictionary.

- A nontrivial unipotent element of order \(p\), for \(q=p^f\), fixes one
  projective point and has \(q/p\) cycles of length \(p\).
- A split semisimple element of order \(e\) fixes two points and has
  \((q-1)/e\) cycles of length \(e\).
- A nonsplit semisimple element of order \(e\) has no fixed point and has
  \((q+1)/e\) cycles of length \(e\).

Thus a single group-theoretic triangle triple can give quotient curves of
different genera in different actions.  This is exactly what happens at
\(q=11\).

For prime \(q=p>3\), the modular generators make this comparison uniform.
Let

\[
 S=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad
 R=\begin{pmatrix}0&-1\\1&1\end{pmatrix},\qquad
 T=S^{-1}R
                                                        \tag{4.2}
\]

in \(\operatorname{PSL}_2(\mathbb F_p)\).  Their orders are \(2,3,p\);
the branch triple \((R^{-1},S,T)\) has product one.  It generates because
\(T\) and its \(S\)-conjugate are the upper and lower elementary
unipotents.  Put

\[
 f_2=1+\left(\frac{-1}{p}\right),\qquad
 f_3=1+\left(\frac{-3}{p}\right).                    \tag{4.3}
\]

These are exactly the numbers of fixed projective points of \(S\) and
\(R\), since their fixed-point discriminants are \(-4\) and \(-3\).
The natural passport is therefore

\[
 3^{(p+1-f_3)/3}1^{f_3}\mid
 2^{(p+1-f_2)/2}1^{f_2}\mid p\,1,
\]

and Riemann--Hurwitz gives the closed formula

\[
 \boxed{\quad
 g_{\mathrm{nat}}(p)
 ={p-6-3\left(\frac{-1}{p}\right)
        -4\left(\frac{-3}{p}\right)\over12}.
 \quad}                                               \tag{4.4}
\]

In particular

\[
 g_{\mathrm{nat}}(7)=0,\qquad
 g_{\mathrm{nat}}(11)=1,\qquad
 g_{\mathrm{nat}}(13)=0.                             \tag{4.5}
\]

Thus the elliptic obstruction at \(11\) is arithmetic, not monotone in
\(p\): both adjacent prime natural actions return to genus zero.  This
closes the first genus comparison, but it does not construct absolute Keller
maps for the degree-eight or degree-fourteen actions.

## 5. The exact \(\operatorname{PSL}_2(11)\) action calculation

Construct \(G\) as

\[
 \operatorname{SL}_2(\mathbb F_{11})/\{\pm I\}.
\]

The checker enumerates all \(660\) elements.  Its conjugacy classes have the
following action fingerprints; the two order-five and two order-eleven rows
each denote two distinct conjugacy classes with the same displayed cycle
type.

| element order | class size | natural degree 12 | either exceptional degree 11 |
|---:|---:|---|---|
| 1 | 1 | \(1^{12}\) | \(1^{11}\) |
| 2 | 55 | \(2^6\) | \(2^4 1^3\) |
| 3 | 110 | \(3^4\) | \(3^3 1^2\) |
| 5 | 132, 132 | \(5^2 1^2\) | \(5^2 1\) |
| 6 | 110 | \(6^2\) | \(6\,3\,2\) |
| 11 | 60, 60 | \(11\,1\) | \(11\) |

### 5.1 The exceptional Gassmann pair

There are exactly \(22\) subgroups generated by a \((2,3,5)\) pair and of
order \(60\).  They split into two \(G\)-conjugacy classes of eleven
subgroups, all isomorphic to \(A_5\).  Choose representatives
\(H_+,H_-\).  Then

\[
 H_+\not\sim_G H_-,\qquad
 \mathbf 1_{H_+}^G=\mathbf 1_{H_-}^G.                \tag{5.1}
\]

The checker proves the stronger elementwise statement: every
\(g\in G\) has the same full cycle partition on \(G/H_+\) and
\(G/H_-\).  It also finds cross-action subdegrees \(5+6\), so neither
subgroup fixes a point in the other coset action.  This is not the Fano
point/hyperplane construction: the two stabilizers are exceptional
\(A_5\)'s in \(\operatorname{PSL}_2(11)\).

### 5.2 One triangle, two quotient genera

There are \(1320\) generating triples

\[
 (\sigma_2,\sigma_3,\sigma_{11}),\qquad
 \sigma_2\sigma_3\sigma_{11}=1,                     \tag{5.2}
\]

with respective orders \((2,3,11)\).  They form two simultaneous-conjugacy
orbits of size \(660\), distinguished by the two order-eleven classes.

In either exceptional degree-eleven action the passport is

\[
 3^3 1^2\mid 2^4 1^3\mid 11.                        \tag{5.3}
\]

The permutation indices are \(6,4,10\).  Hence

\[
 2g-2=-2\cdot11+(6+4+10)=-2,
 \qquad g=0.                                         \tag{5.4}
\]

In the natural degree-twelve action the same abstract triple has passport

\[
 3^4\mid2^6\mid11\,1.                               \tag{5.5}
\]

Now the indices are \(8,6,10\), so

\[
 2g-2=-2\cdot12+(8+6+10)=0,
 \qquad g=1.                                         \tag{5.6}
\]

This is the first lesson of the action spectrum: the same regular
\(G\)-closure gives a polynomial genus-zero quotient for \(H_\pm\) and an
elliptic quotient for the natural Borel stabilizer.

## 6. The explicit degree-eleven Shabat pair

Put \(K=\mathbb Q(r)\), \(r^2=-11\).  Following Jones--Zvonkin, define

\[
\begin{aligned}
 p_1(x)&=2x+11-3r,\\
 p_2(x)&=2x^2-(11-3r)x-(22+6r),\\
 p_3(x)&=x^2+11x+55+9r,
\end{aligned}
\]

and

\[
 P_r(x)=\frac{p_1(x)^3p_2(x)^3p_3(x)}{2^{12}3^{14}}. \tag{6.1}
\]

For the second branch value put

\[
\begin{aligned}
 q_1(x)&=2x+5+3r,\\
 q_2(x)&=2x^3+(15-3r)x^2-(12-12r)x+56+96r,\\
 q_3(x)&=2x^3-18x^2+(21+45r)x-(175+279r).
\end{aligned}
\]

Direct exact expansion gives

\[
 \boxed{
 P_r(x)-1=
 \frac{q_1(x)^2q_2(x)^2q_3(x)}{2^{11}3^{14}}.}       \tag{6.2}
\]

The published display has a minus sign and a misplaced parenthesis in this
line; (6.2) is the identity satisfied by the displayed \(p_i\)'s.  The
checker verifies (6.2) directly rather than silently inheriting that typo.

Differentiation gives the complete critical divisor

\[
 P_r'(x)=
 \frac{11}{2^{12}3^{14}}
 p_1(x)^2p_2(x)^2q_1(x)q_2(x).                       \tag{6.3}
\]

All six \(p_i,q_i\) are irreducible over \(K\).  The exact discriminant is

\[
 \operatorname{Disc}_x(P_r(x)-u)
 =-\frac{11^{11}}{2^{60}3^{140}}u^6(u-1)^4
 =\left(
 \frac{11^5r}{2^{30}3^{70}}u^3(u-1)^2
 \right)^2.                                         \tag{6.4}
\]

The conjugate polynomial \(P_{-r}\) gives the other dessin.  Jones and
Zvonkin identify their common regular closure with
\(\operatorname{PSL}_2(11)\) and the two quotients with the nonconjugate
\(A_5\) classes above.  Their paper also proves that the two dessins are
Galois conjugate over \(\mathbb Q\) and are defined over
\(\mathbb Q(\sqrt{-11})\).

The common closure is also visible directly in polynomial coordinates.  If
\(N_r=2^{12}3^{14}P_r\), exact factorization gives

\[
 N_r(X)-N_{-r}(Y)=64\,C_5(X,Y)C_6(X,Y),              \tag{6.5}
\]

where \(C_i\) has bidegree \((i,i)\).  The full coefficient lists are pinned
in the checker, which verifies (6.5) by direct expansion.  The bidegrees
match the two cross-action \(H_+\)-orbits of sizes five and six on
\(G/H_-\).  Together with the Jones--Zvonkin common-closure theorem, this
identifies \(C_5,C_6\) as the two irreducible double-coset correspondences.

Reference: G. A. Jones and A. K. Zvonkin,
[*Klein's ten planar dessins of degree 11, and beyond*](https://arxiv.org/abs/2104.12015),
especially Sections 4--6.

### 6.1 The normalizations are \(A_4\)- and \(D_{10}\)-quotients

Let \(\mathcal X\) be the common regular
\(\operatorname{PSL}_2(11)\)-cover.  For a representative \(g\) of a
double coset \(H_+gH_-\), the normalization of the corresponding fiber-product
component is

\[
 \mathcal X/(H_+\cap gH_-g^{-1}).                    \tag{6.6}
\]

The order-five and order-six cross-orbits give intersection subgroups

\[
 J_5\simeq A_4,\quad |J_5|=12,
 \qquad
 J_6\simeq D_{10},\quad |J_6|=10.                   \tag{6.7}
\]

The checker obtains these identifications internally from the element-order
censuses

\[
\begin{aligned}
 J_5 &: 1^1,2^3,3^8,\\
 J_6 &: 1^1,2^5,5^4.
\end{aligned}
\]

Consequently the two normalizations have degrees \(55\) and \(66\) over the
\(u\)-line.  Their passports and genera are

| curve | stabilizer | degree over \(u\) | passport \(u=0\mid1\mid\infty\) | genus |
|---|---:|---:|---|---:|
| \(\widetilde C_5\) | \(A_4\) | 55 | \(3^{17}1^4\mid2^{26}1^3\mid11^5\) | 1 |
| \(\widetilde C_6\) | \(D_{10}\) | 66 | \(3^{22}\mid2^{30}1^6\mid11^6\) | 2 |

For example, the total permutation indices are \(110\) and \(134\), so

\[
\begin{aligned}
 2g(\widetilde C_5)-2&=-2\cdot55+110=0,\\
 2g(\widetilde C_6)-2&=-2\cdot66+134=2.             \tag{6.8}
\end{aligned}
\]

Thus normalization closes the proposed lower-color shortcut:
\(\widetilde C_5\) has genus one and \(\widetilde C_6\) has genus two.
Neither becomes rational after adjoining independent variables.

### 6.2 Exact conductor support

Although \(C_d\) has bidegree \((d,d)\) in
\(\mathbb P^1\times\mathbb P^1\), its displayed affine polynomial has total
degree \(d\).  Its ordinary projective-plane closure has squarefree top form,
so all points at infinity are smooth.  Singular's exact normal.lib
calculation gives

| curve | arithmetic genus | affine nodes | conductor colength | geometric genus |
|---|---:|---:|---:|---:|
| \(C_5\) | 6 | 5 | 5 | 1 |
| \(C_6\) | 10 | 8 | 8 | 2 |

In both rows the conductor ideal equals the reduced singular ideal.  Put
\(p_i^-(y)=p_i(y)|_{r\mapsto-r}\), and similarly for \(q_i^-\).  The node
support of \(C_5\) is the disjoint union

\[
\begin{aligned}
 V(p_2(x),p_1^-(y))&\quad(\text{degree }2),\\
 V(p_1(x),p_2^-(y))&\quad(\text{degree }2),\\
 V(q_1(x),q_1^-(y))&\quad(\text{degree }1).          \tag{6.9}
\end{aligned}
\]

The node support of \(C_6\) is

\[
\begin{aligned}
 V(p_2(x),p_2^-(y))&\quad(\text{degree }4),\\
 V(p_1(x),p_1^-(y))&\quad(\text{degree }1),\\
 V\bigl(q_2^-(y),\,
 18x+(r+1)y^2+(8r-10)y+88-20r\bigr)
 &\quad(\text{degree }3).                            \tag{6.10}
\end{aligned}
\]

The Hessian determinant is a unit on each reduced singular scheme, so every
listed point is an ordinary node.  This also explains an apparent
compactification discrepancy.  The plane closure keeps the \(d\) smooth
directions at infinity separate.  The closure in
\(\mathbb P^1\times\mathbb P^1\) sends all of them to
\((\infty,\infty)\), creating an ordinary \(d\)-fold point of delta
\(\binom d2\).  That extra delta is a graph-compactification artifact, not
affine conductor.

### 6.3 Boundary pullback through the two projections

Each normalization has two projections

\[
 \widetilde C_d\longrightarrow X_+,\qquad
 \widetilde C_d\longrightarrow X_-                 \tag{6.11}
\]

of degree \(d\).  The ramification partitions are the same for both
projections.  In the table, a row applies to every geometric root of the
listed boundary factor.

| boundary point on \(X_\pm\) | inertia over \(u\) | degree-five projection | degree-six projection |
|---|---:|---|---|
| root of \(p_1\) or \(p_2\) | 3 | \(1^5\) | \(1^6\) |
| root of \(p_3\) | 1 | \(3,1,1\) | \(3,3\) |
| root of \(q_1\) or \(q_2\) | 2 | \(1^5\) | \(1^6\) |
| root of \(q_3\) | 1 | \(2,2,1\) | \(2,2,1,1\) |
| \(\infty\) | 11 | \(1^5\) | \(1^6\) |

Thus all projection ramification is concentrated above the unramified tails
\(p_3=0\) and \(q_3=0\).  Its total degrees are

\[
 R(\widetilde C_5/X_\pm)=10,\qquad
 R(\widetilde C_6/X_\pm)=14,                         \tag{6.12}
\]

which independently recover genera one and two from the degree-five and
degree-six maps to \(\mathbb P^1\).

### 6.4 The degree-five normalization is the conductor-\(121\) elliptic curve

The genus-one component has a \(K\)-point for a reason which is visible
before finding any Weierstrass coordinates.  The rational affine node

\[
 R=\left(-\frac{5+3r}{2},-\frac{5-3r}{2}\right)
\]

pulls back to a \(K\)-rational effective divisor \(D_2\) of degree two on
\(\widetilde C_5\).  The line at infinity pulls back to a \(K\)-rational
effective divisor \(D_5\) of degree five.  Hence \(3D_2-D_5\) is a
\(K\)-rational divisor of degree one.  Riemann--Roch on a genus-one curve
makes it linearly equivalent to an effective divisor of degree one, so

\[
 \widetilde C_5(K)\ne\varnothing.                  \tag{6.13}
\]

Thus there is no nontrivial torsor obstruction: after choosing a point,
\(\widetilde C_5\) is its own Jacobian.

There is also an exact plane-cubic reduction.  Put

\[
 s^2=2(11-r),\qquad t^2=2(11+r).
\]

A quadratic transformation centered at (R) and the conjugate pair
\(V(p_2(x),p_1^-(y))\) takes the quintic to a quartic.  The line through
that pair contracts to a \(K\)-rational smooth point.  A second quadratic
transformation centered at that point and the remaining conjugate pair
\(V(p_1(x),p_2^-(y))\) gives the sparse cubic

\[
 H(U,V,W)=\sum c_{ijk}U^iV^jW^k=0,                \tag{6.14}
\]

over (K(t)), with nonzero coefficients

\[
\begin{array}{c|l}
ijk&c_{ijk}\\ \hline
210&\frac{264627-29403t+8019r+5346rt}{8}\\
201&\frac{1234926-65043rt-421443t+753786r}{32}\\
120&\frac{264627-5346rt+29403t+8019r}{8}\\
111&78408+60588r\\
102&\frac{-32670-5049rt-35937t+196614r}{4}\\
021&\frac{1234926+421443t+65043rt+753786r}{32}\\
012&\frac{-32670+35937t+5049rt+196614r}{4}\\
003&-28314+16434r.
\end{array}                                       \tag{6.15}
\]

The checker replays both Cremona substitutions in the primitive degree-eight
field generated by (r,s,t).  Applying the exact ternary-cubic Jacobian
formula to (6.14) gives

\[
 j(\widetilde C_5)=-121.                           \tag{6.16}
\]

More precisely, the temporary cubic has twist parameter
\(2(11+r)(4455+891r)^2\) relative to the small curve

\[
 E_{121}:\qquad v^2+uv=u^3+u^2-2u-7.              \tag{6.17}
\]

The two temporary node fields leave four possible descent twists, represented
by (1,2(11-r),2(11+r),33).  They are separated by exact reduction at the
two primes above \(23\), corresponding to \(r=9,14\pmod {23}\).  At either
prime the plane quintic has (20) affine points and five points at infinity.
Its three rational nodes have nonsplit tangents, so normalization has

\[
 \#\widetilde C_5(\mathbb F_{23})=20+5-3=22,
 \qquad a_{\mathfrak p}=2.                         \tag{6.18}
\]

The quadratic characters of the four candidate twists at the ordered pair
of primes are

\[
 (1,1),\quad(1,-1),\quad(-1,1),\quad(-1,-1).
\]

Since (6.17) also has trace pair ((2,2)), (6.18) selects the first row:

\[
 \boxed{\ \widetilde C_5\simeq_K E_{121}\ }.       \tag{6.19}
\]

The model (6.17) has

\[
 c_4=121,\qquad c_6=5203,\qquad
 \Delta=-11^4,\qquad j=-121.                      \tag{6.20}
\]

Over \(\mathbb Q\) it is the conductor-(121\) curve with LMFDB label
[121.c2 (Cremona 121c1)](https://www.lmfdb.org/EllipticCurve/Q/121/c/2).
Its (-11)-twist is 121.a1 (Cremona 121a2), so those two rational
isogeny classes become isomorphic after base change to \(K\).

This resolves the comparison with the natural quotient decisively.  At
\(23\), the curve \(X_0(11)\) in (9.1) has (25) points and trace (-1),
whereas (6.18) has trace (2).  Therefore

\[
 \operatorname{Jac}(\widetilde C_5)
 \not\sim_K X_0(11).                               \tag{6.21}
\]

This is stronger than observing that the two curves came from different
subgroups: their elliptic isogeny types are genuinely different.

The boundary calculation on this curve is especially clean after
normalization.  Over the seven boundary blocks

\[
 (p_1,p_2,p_3,q_1,q_2,q_3,\infty)
\]

the degrees of the normalized \(K\)-prime punctures are

\[
 (1,4)\mid(4,4,2)\mid(4,2)\mid(2,3)\mid
 (3,6,6)\mid(6,3)\mid(5).                         \tag{6.21a}
\]

Thus there are fifteen primes, of total degree \(55\).  The node tangent
extensions are essential here: the three plane-node primes of degrees
\(2,2,1\) become normalized branch primes of degrees \(4,4,2\).

PARI proves that both \(E_{121}(\mathbb Q)\) and its \((-11)\)-twist have
rank zero.  Reduction at the split primes \(3\) and \(5\) gives group orders
\(2\) and \(5\), so \(E_{121}(K)\) also has trivial torsion.  Consequently

\[
 E_{121}(K)=0.                                      \tag{6.21b}
\]

If \(d=(d_1,\ldots,d_{15})\) is the degree vector in (6.21a), the full unit
lattice is therefore

\[
 \boxed{L_5=\ker\bigl(d:\mathbb Z^{15}\to\mathbb Z\bigr).} \tag{6.21c}
\]

It has rank \(14\).  Since \(d_1=1\), a primitive basis is
\(-d_i e_1+e_i\), \(2\le i\le15\), and
\(\mathbb Z^{15}/L_5\simeq\mathbb Z\).  This closes the full boundary-unit
problem for the degree-five normalization, not merely its Mordell--Weil
rank.

### 6.5 The degree-six normalization is bielliptic

The genus-two component also admits an exact canonical reduction.  For a
plane sextic with eight ordinary nodes, the canonical system is the pencil
of cubics through the node scheme.  In the present coordinates the checker
finds the following basis:

\[
\begin{aligned}
 A={}&4x^2y+(r-1)xy^2-2y^3+4x^2+2rxy+(r-7)y^2\\
    &-(8r+44)x+(20r+44)y-44r+220,\\
 B={}&x^3+3x^2-\frac12xy^2+\frac{11+3r}{4}xy
       +\frac{-55+33r}{2}x\\
    &+\frac{1+r}{4}y^3+\frac{1+3r}{4}y^2
       +(44-12r)y-198-44r.
                                                        \tag{6.22}
\end{aligned}
\]

Exact reduction modulo the three node ideals in (6.10) shows that both
cubics vanish on the full length-eight conductor scheme.  Conversely, the
restriction map from the ten-dimensional space of plane cubics to that
scheme has rank eight, so \(A,B\) is the whole canonical pencil.

Put \(t=B/A\).  Eliminating \(y\) from \(C_6=0\) and \(B-tA=0\), then removing the
node factors of multiplicities \(2,2,4\), leaves a quadratic in \(x\).  Its
discriminant is a cubic square times a squarefree sextic.  After setting
\(u=4t\), its square class gives

\[
 w^2=-(3r+1)h(u),                                  \tag{6.23}
\]

where

\[
\begin{aligned}
 h(u)={}&26u^6+(78r-138)u^5-(345r+659)u^4
 +(112r+2856)u^3\\
 &+(489r-2675)u^2+(-486r+1230)u+120r-142.
                                                        \tag{6.24}
\end{aligned}
\]

This form reveals an extra involution.  Besides the hyperelliptic
involution, (6.23) is preserved by

\[
 u\longmapsto \frac{-6u+1-7r}{(1-r)u+6}.           \tag{6.25}
\]

Its fixed points are \(\alpha=(3-r)/2\) and
\(\beta=(-5-r)/2\).  The coordinate

\[
 z=\frac{u-\alpha}{u-\beta},\qquad
 u=\frac{\beta z-\alpha}{z-1}                     \tag{6.26}
\]

sends (6.25) to \(z\mapsto-z\).  The exact polynomial identity

\[
 (z-1)^6[-(3r+1)h(u(z))]
 =512(1+3r)(1-11z^2-77z^4-121z^6)
\]

therefore gives the particularly small model

\[
 \boxed{\ \widetilde C_6:\quad
 Y^2=d(1-11z^2-77z^4-121z^6),\qquad
 d=2(1+3r).\ }                                    \tag{6.27}
\]

The two non-hyperelliptic involutions now give degree-two quotient maps

\[
\begin{array}{rcll}
 E_+:&Y^2=d(1-11X-77X^2-121X^3),
     &X=z^2,\quad Y=Y,&j(E_+)=-4096/11,\\
 E_-:&W^2=dX(1-11X-77X^2-121X^3),
     &X=z^2,\quad W=zY,&j(E_-)=-32768.
\end{array}                                       \tag{6.28}
\]

Both descents are already split over \(K\), with no further quadratic
extension.  For \(E_+\), the substitution that makes the cubic monic gives
invariants equal to those of

\[
 E_{11}:\quad y^2+y=x^3-x^2                       \tag{6.29}
\]

twisted by

\[
 44d=(22+6r)^2.
\]

Thus \(E_+\simeq_K E_{11}\).  Over \(\mathbb Q\), (6.29) is
[LMFDB 11.a3](https://www.lmfdb.org/EllipticCurve/Q/11/a/3), hence lies in
the conductor-\(11\) isogeny class of \(X_0(11)\).  For \(E_-\), the binary
quartic invariants are \(I=352,J=13552\).  Its Jacobian is the twist by

\[
 -36d=[6(3-r)]^2
\]

of

\[
 E_{\mathrm{CM}}:\quad y^2+y=x^3-x^2-7x+10.       \tag{6.30}
\]

Hence \(E_-\simeq_K E_{\mathrm{CM}}\).  The rational curve (6.30) is the
optimal curve [LMFDB 121.b2](https://www.lmfdb.org/EllipticCurve/Q/121/b/2)
and has complex multiplication by the order of discriminant \(-11\).

The two quotient pullbacks give the promised decomposition

\[
 \boxed{\ \operatorname{Jac}(\widetilde C_6)
 \sim_K^{(2,2)} E_{11}\times E_{\mathrm{CM}}.\ }   \tag{6.31}
\]

This splitting also closes the boundary-unit calculation.  In terms of the
canonical coordinates \(t,w\) in (6.23), put

\[
 z_n=8t-3+r,\qquad z_d=8t+5+r.
\]

The two quotient maps to the rational minimal models

\[
\begin{aligned}
 E_+&:y^2+y=x^3+2x^2+x,\\
 E_-&:y^2+y=x^3+2x^2-6x+3
\end{aligned}
\]

are

\[
\begin{aligned}
X_+&=\frac{-121d z_n^2-(110+330r)z_d^2}
 {z_d^2(-22-6r)^2},\\
Y_+&=\frac{3872dw-(7744-3168r)z_d^3}
 {z_d^3(-22-6r)^3},\\
X_-&=\frac{d z_d^2-(2+6r)z_n^2}
 {z_n^2(6-2r)^2},\\
Y_-&=\frac{-32dw+(288+64r)z_n^3}
 {z_n^3(6-2r)^3}.                                \tag{6.32}
\end{aligned}
\]

Singular separates the normalized punctures, and exact Riemann--Roch
interpolation in their residue fields computes the trace of every prime
through (6.32).  Order the twenty \(K\)-prime punctures as

\[
\begin{gathered}
p_{1n},p_{1a},p_{1b}\mid
p_{2a},p_{2n,a},p_{2n,b},p_{2b}\mid
p_{3a},p_{3b}\mid q_{1a},q_{1b}\mid\\
q_{2a},q_{2n},q_{2b},q_{2c}\mid
q_{3a},q_{3b},q_{3c}\mid\infty_1,\infty_5.
\end{gathered}                                      \tag{6.33}
\]

Their degree vector is

\[
\mathbf d=(2,2,2,2,4,4,2,2,2,3,3,3,6,3,6,3,3,6,1,5).
                                                               \tag{6.34}
\]

Let \(T=(0,0)\in E_+(K)\), of order five, and let

\[
 G_0=(3,5),\qquad
 G_2=\left(-8,\frac{-1+11r}{2}\right)\in E_-(K).   \tag{6.35}
\]

Here \(G_0\) is rational and \(G_2\) is anti-invariant.  PARI proves that the
rational curve and its \((-11)\)-twist both have rank one, so these points
are independent and \(E_-(K)\) has rank two.  Reduction at the split primes
\(3,5\) gives orders \(5,9\), proving that \(E_-(K)\) has trivial torsion.
For \(E_+\), both reduction orders are five; together with the PARI rank-zero
calculation this gives \(E_+(K)=\langle T\rangle\simeq\mathbb Z/5\).

In the order (6.33), the three class rows are

\[
\begin{aligned}
\mathbf c_T={}&(0,1,0,1,4,3,4,0,4,1,0,4,1,1,2,4,0,0,2,4),\\
\mathbf c_0={}&(4,-1,3,-1,0,12,1,3,1,2,4,0,2,2,14,0,4,10,1,5),\\
\mathbf c_2={}&(0,1,-1,-1,0,0,1,1,-1,-1,1,-1,0,1,0,1,-1,0,0,0).
                                                               \tag{6.36}
\end{aligned}
\]

Thus a boundary divisor \(\mathbf n\) maps to
\((\mathbf c_T\mathbf n)T\) on \(E_+\) and to
\((\mathbf c_0\mathbf n)G_0+(\mathbf c_2\mathbf n)G_2\) on \(E_-\).
The \(q_{2n}\) column can either be eliminated directly or recovered from
\(\operatorname{div}(q_2(x))\); both elliptic rows satisfy every weighted
principal-fiber relation, including the ramification-three \(p_3\) rows and
ramification-two \(q_3\) rows.

It remains to justify that the \((2,2)\)-isogeny has not lost a rational
two-torsion class.  Its geometric kernel is the even \(\mathbb F_2\)-subspace
on the three pairs of roots of

\[
 1-11X-77X^2-121X^3.                              \tag{6.37}
\]

This cubic is irreducible over \(K\), and its discriminant is
\((176r)^2\).  Its cyclic cubic Galois group therefore permutes the three
nonzero kernel points transitively.  The isogeny kernel has no nonzero
\(K\)-rational point.  Vanishing in both elliptic quotients is consequently
equivalent to vanishing in \(\operatorname{Jac}(\widetilde C_6)(K)\).

The full unit lattice is therefore

\[
\boxed{
L_6=\left\{\mathbf n\in\mathbb Z^{20}:
\begin{array}{l}
\mathbf d\mathbf n=0,\\
\mathbf c_0\mathbf n=\mathbf c_2\mathbf n=0,\\
\mathbf c_T\mathbf n\equiv0\pmod5
\end{array}\right\}.}                              \tag{6.38}
\]

It has rank \(17\).  The checker records a row-Hermite basis and computes

\[
 \mathbb Z^{20}/L_6\simeq\mathbb Z^3\oplus\mathbb Z/5. \tag{6.39}
\]

### 6.6 Both projection images and their cokernels

The full lattices make it possible to compare the two projections
integrally, rather than only through their ramification partitions.  Use the
compact target boundary

\[
 S_\pm=(p_1,p_2,p_3,q_1,q_2,q_3,\infty),\qquad
 \deg S_\pm=(1,2,2,1,3,3,1).                       \tag{6.40}
\]

Since \(X_\pm\simeq\mathbb P^1\), its full boundary-unit lattice is the
degree kernel \(U_\pm\simeq\mathbb Z^6\).  Let \(P_{d,x}\) and \(P_{d,y}\)
be the seven-row compact divisor-pullback matrices for the two projections.
The following table specifies them completely.  A superscript \(n\) denotes
the normalized prime above a node, and coefficients are ramification
indices.

| target prime | \(P_{5,x}\) | \(P_{5,y}\) | \(P_{6,x}\) | \(P_{6,y}\) |
|---|---|---|---|---|
| \(p_1\) | \(P_{11}+P_{12}^{n}\) | \(P_{11}+P_{21}^{n}\) | \(p_{1n}+p_{1a}+p_{1b}\) | \(p_{1n}+p_{2a}+p_{3a}\) |
| \(p_2\) | \(P_{21}^{n}+P_{22}+P_{23}\) | \(P_{12}^{n}+P_{22}+P_{32}\) | \(p_{2a}+p_{2n,a}+p_{2n,b}+p_{2b}\) | \(p_{1a}+p_{2n,a}+p_{2n,b}+p_{3b}\) |
| \(p_3\) | \(P_{33}+3P_{32}\) | \(3P_{23}+P_{33}\) | \(3p_{3a}+3p_{3b}\) | \(3p_{1b}+3p_{2b}\) |
| \(q_1\) | \(Q_{11}^{n}+Q_{12}\) | \(Q_{11}^{n}+Q_{21}\) | \(q_{1a}+q_{1b}\) | \(q_{2b}+q_{3b}\) |
| \(q_2\) | \(Q_{21}+Q_{22}+Q_{23}\) | \(Q_{12}+Q_{22}+Q_{32}\) | \(q_{2a}+q_{2n}+q_{2b}+q_{2c}\) | \(q_{1a}+q_{2n}+q_{2c}+q_{3a}\) |
| \(q_3\) | \(2Q_{32}+Q_{33}\) | \(2Q_{23}+Q_{33}\) | \(2q_{3a}+2q_{3b}+q_{3c}\) | \(2q_{1b}+2q_{2a}+q_{3c}\) |
| \(\infty\) | \(\infty_5\) | \(\infty_5\) | \(\infty_1+\infty_5\) | \(\infty_1+\infty_5\) |

For \(C_5\), the notation \(P_{ij}\) and \(Q_{ij}\) records the \(x\)- and
\(y\)-colors.  In the order used in (6.21a), the primes are

\[
\begin{gathered}
P_{11},P_{12}^{n}\mid P_{21}^{n},P_{22},P_{23}\mid
P_{33},P_{32}\mid Q_{11}^{n},Q_{12}\mid\\
Q_{21},Q_{22},Q_{23}\mid Q_{32},Q_{33}\mid\infty_5.
\end{gathered}                                      \tag{6.41}
\]

The degree identities

\[
 P_{5,*}\mathbf d_5=5\deg S_\pm,\qquad
 P_{6,*}\mathbf d_6=6\deg S_\pm                   \tag{6.42}
\]

check every row, including the ramification-three and ramification-two
tails.  Define

\[
 I_d=P_{d,x}^{*}U_+ + P_{d,y}^{*}U_-\subseteq L_d. \tag{6.43}
\]

Both summands have rank six.  Their intersection is exactly the common
rank-two pullback of the triangle-base units

\[
 (3,3,1,0,0,0,-11),\qquad
 (0,0,0,2,2,1,-11),                                \tag{6.44}
\]

so \(\operatorname{rank} I_5=\operatorname{rank} I_6=10\).  Exact Smith
reduction in bases of (6.21c) and (6.38) gives

\[
 \boxed{L_5/I_5\simeq\mathbb Z^4,
 \qquad L_6/I_6\simeq\mathbb Z^7.}                \tag{6.45}
\]

In particular both two-projection images are primitive.  The order-five
class of (6.39) has not disappeared: it lives in the prior quotient
\(\operatorname{Div}_{S_6}/L_6\), and is therefore an admission congruence
for a proposed boundary divisor.  Once a divisor is principal, its residual
class modulo the two projection images is torsion-free.  This distinguishes
the \(C_6\) arithmetic obstruction from the index-six saturation on
\(X_0(11)\).

That comparison can also be stated as a cokernel calculation.  If \(J_{12}\)
is generated by the pullbacks of \(j\) and \(j-1728\), then

\[
 L_{12}/J_{12}\simeq\mathbb Z,\qquad
 L_{12}/\langle J_{12},f_T\rangle\simeq\mathbb Z/6. \tag{6.46}
\]

Thus the order six is an internal saturation defect among evident units;
the order five on \(C_6\) is a divisor-class obstruction before unit
pullback is considered.

Finally, compactify the derivative row (8.3) as

\[
 \delta=(2,2,0,1,1,0,-10).                         \tag{6.47}
\]

The smallest two-output valuation ledger is generated by the two common
rows (6.44) and \(P_{d,x}^{*}\delta,P_{d,y}^{*}\delta\).  These four rows
are independent and primitive in both \(L_5\) and \(L_6\), giving

\[
 L_5/A_5\simeq\mathbb Z^{10},\qquad
 L_6/A_6\simeq\mathbb Z^{13}.                      \tag{6.48}
\]

So the weak two-output derivative ledger passes: there is no forced square,
cube, or fifth root at this level.  The stronger two-mask character
completion fails by rank.  Even starting from the whole rank-ten lattice
\(I_d\), two new divisor rows leave free rank at least two for \(C_5\) and
five for \(C_6\).  A monomial boundary architecture that requires all
characters therefore needs at least four and seven independent new mask
classes, respectively.  A nonlinear polynomial output may encode several
boundary factors, so this last statement does not exclude arbitrary
coupled two-output modifications.

### 6.7 Residual symmetry and the first factor-rich support

The free quotients in (6.45) can be made constructive.  Let \(\tau\) be the
semilinear involution which exchanges the two projections and sends
\(r\mapsto-r\).  Exact reduction modulo \(I_d\) gives

\[
 \boxed{\quad
 (L_5/I_5,\tau)\simeq\mathbb Z_{\mathrm{triv}}^4,\qquad
 (L_6/I_6,\tau)\simeq
 \mathbb Z_{\mathrm{triv}}^5\oplus\mathbb Z[\mathbb Z/2].\quad} \tag{6.49}
\]

The second equality is integral, not merely rational.  In an explicit
seven-element quotient basis, five rows are fixed and the final two are
exchanged; the change of basis is unimodular.  Thus the unique
permutation-basis signatures are

\[
\begin{array}{c|c|c|c}
 &\text{fixed classes}&\text{exchanged pairs}&\text{distinct classes}\\ \hline
 C_5&4&0&4\\
 C_6&5&1&7.
\end{array}                                         \tag{6.50}
\]

There are also representatives which are regular away from infinity.  In
the \(C_5\) notation (6.41), take

\[
\begin{aligned}
 \operatorname{div}(M_1)&=Q_{11}^{n}+Q_{33}-\infty_5,\\
 \operatorname{div}(M_2)&=P_{32}+Q_{33}-\infty_5,\\
 \operatorname{div}(M_3)&=P_{11}+2Q_{11}^{n}-\infty_5,\\
 \operatorname{div}(M_4)&=2Q_{11}^{n}+Q_{32}-2\infty_5.
                                                        \tag{6.51}
\end{aligned}
\]

All four rows lie in \(L_5\), and their images in \(L_5/I_5\) have
determinant minus one.  Moreover this pole profile is optimal.  There are exactly
twenty-six effective \(K\)-rational boundary divisors of degree five whose
finite part is supported on the fourteen finite punctures.  Their
simple-pole rows have quotient Smith diagonal \((1,1,1,2)\).  Hence
simple-pole masks span an index-two sublattice; any regular mask basis needs
at least one pole of order two, and (6.51) attains the minimum total pole
order \(1+1+1+2=5\).

For \(C_6\), one effective permutation-basis lift is recorded below.  The
pole pair \((a,b)\) means \(a\infty_1+b\infty_5\).

| mask | finite zero divisor | pole pair |
|---|---|---:|
| \(N_1\) | \(4p_{1a}+8p_{1b}+4p_{3a}\) | \((7,5)\) |
| \(N_2\) | \(8p_{1n}+p_{2a}+5p_{2b}+4p_{3b}+q_{2a}+q_{2b}\) | \((7,7)\) |
| \(N_3\) | \(5p_{1n}+p_{1a}+p_{2b}+2p_{3b}+q_{2n}\) | \((4,4)\) |
| \(N_4\) | \(5p_{1n}+5p_{2b}+5p_{3b}\) | \((5,5)\) |
| \(N_5\) | \(2p_{1n}+p_{1a}+p_{2a}+p_{2n,b}+3p_{2b}+3p_{3b}\) | \((4,4)\) |
| \(N_6\) | \(5p_{1a}+10p_{1b}+4p_{3a}+q_{1b}\) | \((6,7)\) |
| \(N_7\) | \(4p_{1b}+5p_{2a}+10p_{3a}+q_{3b}\) | \((6,7)\) |

The first five quotient classes are \(\tau\)-fixed and the last two are
exchanged.  Every displayed divisor satisfies all three Mordell--Weil rows
and the order-five congruence in (6.38).  Their quotient determinant is
minus one.

Now require a two-normal correction to vanish on the old zero section, and
give each independent residual class a separate nonconstant monomial slot.
Writing \(\tau(s)=t\), the first support is

\[
\begin{aligned}
 \Phi_5={}&M_1s+M_2t+M_3s^2+M_4st,\\
 \Psi_5={}&\tau(M_1)t+\tau(M_2)s+\tau(M_3)t^2+\tau(M_4)st, \tag{6.52}\\[2mm]
 \Phi_6={}&N_1s+N_2t+N_3s^2+N_4st+N_5t^2+N_6s^3,\\
 \Psi_6={}&\tau(N_1)t+\tau(N_2)s+\tau(N_3)t^2+\tau(N_4)st
             +\tau(N_5)s^2+N_7t^3.                 \tag{6.53}
\end{aligned}
\]

Projection-unit factors may be inserted without changing this residual
support.  The normal-degree bounds are minimal within the declared
one-class-per-monomial, zero-section-preserving architecture: two variables
have only two nonconstant monomials through degree one and five through
degree two.  Thus four slots first occur in degree two, while six slots first
occur in degree three.

### 6.8 Exact normalization-module interpolation

The divisor rows above can be lifted without first finding a smooth plane
model.  Suppose an affine plane curve has normalization module

\[
 \overline R=R+RT_1+\cdots+RT_s,\qquad T_i=u_i/c,
                                                        \tag{6.54}
\]

and the top homogeneous part of \(c\) is nonzero at every point of the
infinity divisor \(D_\infty\).  Give \(T_i\) pole weight
\(\deg u_i-\deg c\).  Multiples of \(T_i\) by plane polynomials of the
remaining degree then give a filtered subspace of
\(L(mD_\infty)\).  Once exact reduction modulo the normalization relations
shows that its dimension is
\(m\deg D_\infty-g+1\), Riemann--Roch proves equality.  Vanishing at a
finite prime to order \(e\) is ordinary linear evaluation in
\(\overline R/\mathfrak p^e\).  Thus a prescribed divisor becomes a kernel
calculation over \(K\), followed by a principal-ideal check which excludes
unintended affine zeros.  This is a reusable normalization-module
interpolation test, rather than a curve-specific search for a parametrization.

For \(C_5\), Singular returns

\[
\begin{aligned}
c_5={}&4x^2-2rxy-4y^2+(11-5r)x-(11+5r)y+10r,\\
u_{51}={}&2y^3-6xy-(9r+33)x-(18r+66)y-550,\\
u_{52}={}&2xy^2+16xy+(11-3r)y^2+(77-9r)x\\
 &\quad +(88-24r)y-165r+275,
\qquad T_i=u_{5i}/c_5 .                         \tag{6.55}
\end{aligned}
\]

The top denominator is coprime to the quintic at infinity, and
\(1,x,y,T_1,T_2\) is \(L(D_\infty)\).  Exact representatives of (6.51) are

\[
\begin{aligned}
M_1={}&39-5r-(2r+2)x+2y+(8-4r)T_1+14T_2,\\
M_2={}&259-101r-(16r+20)x+(2r+36)y
       +(84-40r)T_1+134T_2,\\
M_3={}&121-15r+(8-4r)x+(20-2r)y
       +(24-12r)T_1+42T_2,\\
M_4={}&101r-2105+(6r+84)x-(17r+445)y
       +(24r-768)T_1-(167r+475)T_2\\
 & +(18-2r)x^2+(45-5r)xy-(6r+38)y^2
       +(25r+51)xT_2+92yT_2.                    \tag{6.56}
\end{aligned}
\]

Reduction in the normalization ring proves that the four principal ideals
have exactly the finite parts in (6.51), of lengths \(5,5,5,10\).  The
exchange is also explicit:

\[
\begin{aligned}
 \tau(T_1)&={11-5r\over8}-{x\over2}-{r y\over4}
             -{rT_1\over2}+{7T_2\over4},\\
 \tau(T_2)&=-{11+3r\over4}-{y\over2}-T_1-{rT_2\over2}.
                                                        \tag{6.57}
\end{aligned}
\]

For \(C_6\), the normalization denominator is exactly the canonical adjoint
\(A\) in (6.23).  Its three numerator degrees are \(4,4,5\), against
\(\deg A=3\), so \(T_1,T_2,T_3\) have pole weights \(1,1,2\).  Consequently

\[
 V_m=K[x,y]_{\le m}+K[x,y]_{\le m-1}T_1
       +K[x,y]_{\le m-1}T_2+K[x,y]_{\le m-2}T_3. \tag{6.58}
\]

Exact normalization reduction gives

\[
 \dim V_4=23,\qquad \dim V_5=29,\qquad \dim V_7=41,
                                                        \tag{6.59}
\]

which equals \(6m-1=\ell(mD_\infty)\).  Evaluation modulo the boundary
prime powers in the table gives a one-dimensional kernel for each of
\(N_2,\ldots,N_7\), and the associated principal ideals have finite lengths
\(42,24,30,24,41,41\).  At the degree-eight \(p_2\)-node the rank test itself
separates the two degree-four normalized branches and selects
\(p_{2n,b}\).

The \(32\) finite conditions for \(N_1\) leave dimension nine in \(V_7\).
In the chart \(x=w/z,\ y=1/z\), write a candidate numerator and denominator
as \(H_0(w)+zH_1(w)+\cdots\) and \(f_0(w)+zf_1(w)+\cdots\).  Along the
degree-five infinity factor

\[
 g_5(w)=w^5+{r+1\over2}w^4-w^3+w^2+{r-1\over2}w-1,
\]

the two cancellation conditions are

\[
 H_0=0,\qquad f_0'H_1-f_1H_0'=0\pmod {g_5}.       \tag{6.60}
\]

Their ten exact rows raise the combined rank to \(40\) and select \(N_1\)
uniquely, with pole pair \((7,5)\).  Finally define

\[
 \delta_\infty(D)=\operatorname{coeff}_{\infty_1}(D)
                  -\operatorname{coeff}_{\infty_5}(D). \tag{6.61}
\]

Both projection images lie in its kernel, so it descends to \(L_6/I_6\).
On the ordered basis \(N_1,\ldots,N_7\),

\[
 \delta_\infty=(-2,0,0,0,0,1,1).                 \tag{6.62}
\]

Hence no change by projection units can give equal infinity orders to the
\(N_1,N_6,N_7\) classes.  This is an action-presentation invariant that the
divisor-class rank alone did not see.

Equations (6.52)--(6.53) are still a finite coefficient-support ansatz, not
a Keller map.  The Riemann--Roch representative gap is now closed in the
normalization algebra.  The next calculation must turn these rational
normalization functions into a complete ambient polynomial output block and
impose the full constant-Jacobian and inverse-adjugate divisibilities.
Polynomial descent, affine-space recognition, and preservation of both
degree-eleven fields remain open.

This comparison is externally meaningful.  One elliptic factor recovers
the isogeny type of the natural modular quotient, while the other is the
CM conductor-\(121\) class.  Together with (6.21), the three positive-genus
curves in the correspondence diagram therefore contribute exactly three
rational elliptic isogeny types over \(K\): the non-CM conductor-\(121\) class
with \(j=-121\), the conductor-\(11\) modular class, and the CM
conductor-\(121\) class.

## 7. A second global Gassmann Keller chart

Let

\[
 B=\operatorname{Spec}K[u,(u(u-1))^{-1}]
\]

and

\[
 X_\pm
 =\operatorname{Spec}
 K[x,(P_{\pm r}(x)(P_{\pm r}(x)-1))^{-1}].           \tag{7.1}
\]

Equations (6.2)--(6.3) show that \(P_{\pm r}'\) is a unit on \(X_\pm\).
Therefore

\[
\begin{aligned}
 \widehat\pi_\pm:X_\pm\times\mathbb A^1_z
 &\longrightarrow B\times\mathbb A^1_Z,\\
 (x,z)&\longmapsto
 \left(P_{\pm r}(x),\frac{z}{P_{\pm r}'(x)}\right)
\end{aligned}                                        \tag{7.2}
\]

is a finite etale degree-eleven morphism with determinant one.  The second
target coordinate reconstructs \(z=P_{\pm r}'(x)Z\), so no sheet is added or
lost.  The two maps have the exact same target, nonisomorphic Gassmann inverse
covers, and a common \(\operatorname{PSL}_2(11)\)-closure.

At every good finite-field fiber, the two zero-dimensional fibers have the
same zeta function.  Indeed their Frobenius elements lie in the common group,
and the two actions have equal cycle partitions element by element.  Thus
(7.2) is the degree-eleven analogue of the degree-seven global Sunada Keller
pair, but its subgroup geometry is exceptional \(A_5\), not projective
duality.

This proves

\[
 (G,G/H_+),(G,G/H_-)
 \in\mathfrak K_{\mathrm{chart}}(K).                 \tag{7.3}
\]

It does not put either action in \(\mathfrak K_{\mathrm{abs}}(K)\).

## 8. The exact affine-completion defect

Use the ordered source boundary basis

\[
 ([p_1],[p_2],[p_3],[q_1],[q_2],[q_3]).
\]

Both source and target rings in (7.1) are localizations of polynomial PIDs,
so their class groups vanish.  Their unit ranks and pullback matrix are

\[
 \operatorname{rank}\mathcal O(X_\pm)^\times/K^\times=6,
 \qquad
 \operatorname{rank}\mathcal O(B)^\times/K^\times=2, \tag{8.1}
\]

\[
 \pi^*[u]=(3,3,1,0,0,0),\qquad
 \pi^*[u-1]=(0,0,0,2,2,1).                           \tag{8.2}
\]

The derivative vector is

\[
 [P_r']=(2,2,0,1,1,0).                               \tag{8.3}
\]

It is not in the rational, hence not in the integral, span of the two target
rows: the \(p_3\) and \(q_3\) coordinates force both target coefficients to
vanish.  The quotient of source units by the target pullback has rank four,
and the derivative has a nonzero class in that quotient.

This gives three exact conclusions.

1. The rational derivative suspension (7.2) works precisely because all six
   boundary colors are inverted.
2. Neither the source nor the target chart becomes affine space after
   polynomial stabilization: nonconstant units survive adjoining polynomial
   variables.
3. A polynomial suspension which retains \(P_r(x)\) as an output and is
   triangular in one new coordinate has Jacobian divisible by the
   nonconstant \(P_r'(x)\).  It cannot be Keller.

These statements exclude stable straightening and the direct one-coordinate
completion.  They do not exclude a coupled affine modification which changes
at least two old outputs.

## 9. Why the natural action is a different obstruction

Let \(H_B<G\) be a natural point stabilizer.  It has order \(55\) and index
\(12\).  Quotienting the same regular \((3,2,11)\)-cover by \(H_B\) gives the
genus-one cover (5.5).  After removing the three branch values it is a finite
etale cover of the thrice-punctured line, and the universal volume-ratio
suspension again produces a determinant-one smooth-affine chart.

This quotient is the classical modular curve

\[
 E=X_0(11):\qquad Y^2+Y=X^3-X^2-10X-20,             \tag{9.1}
\]

whose discriminant is \(-11^5\).  The degree-twelve map can be made fully
explicit.  Put

\[
\begin{aligned}
 \mathbf b_{12}
 &=(1,X,Y,X^2,XY,X^3,X^2Y,X^4,X^3Y,X^5,X^4Y,X^6),\\
 \mathbf b_{13}&=(\mathbf b_{12},X^5Y),
\end{aligned}
\]

and define \(D=\mathbf d\cdot\mathbf b_{12}\),
\(N=\mathbf n\cdot\mathbf b_{13}\), where

\[
\begin{aligned}
\mathbf d={}&(8278096,-9028165,-3062640,799092,977047,288374,\\
             &\hspace{11mm}15778,-11099,-3511,-482,-34,-1),\\
\mathbf n={}&(122758525012,82027395739,5600100593,4120902300,\\
             &\hspace{11mm}19113285610,-12574009616,6795826151,
             -1879636243,\\
             &\hspace{11mm}262696343,-15351268,172103,-709,1).
                                                               \tag{9.2}
\end{aligned}
\]

For a function \(A(X)+B(X)Y\), write

\[
 \operatorname{Nm}(A+BY)=A(A-B)-B^2(X^3-X^2-10X-20).
\]

Direct exact calculation gives

\[
\begin{aligned}
 \operatorname{Nm}(D)&=(X-16)^{12},\\
 \operatorname{Nm}(N)&=-(X-16)A_4(X)^3,\\
 \operatorname{Nm}(N-1728D)&=-(X-16)B_6(X)^2,        \tag{9.3}
\end{aligned}
\]

with

\[
\begin{aligned}
A_4(X)={}&X^4-52820X^3+1333262X^2+4971236X+9789217,\\
B_6(X)={}&X^6-288318X^5+141521931X^4+169928888X^3\\
         &+8135691435X^2+30544230678X+28453700753.
                                                               \tag{9.4}
\end{aligned}
\]

Both displayed polynomials are squarefree.  Let \(O\) be the point at
infinity, \(P=(16,-61)\), and \(Q=(16,60)=-P\).  The denominator has divisor

\[
 \operatorname{div}(D)=Q+11P-12O.
\]

The numerator and \(N-1728D\) both vanish simply at \(Q\), so that common
point cancels.  Therefore

\[
 j=\frac ND:E\longrightarrow\mathbb P^1             \tag{9.5}
\]

has four triple zeros, six double points above \(1728\), an order-eleven
pole at \(P\), and a simple pole at \(O\).  This is the passport

\[
 3^4\mid2^6\mid11\,1
\]

in exact Weierstrass coordinates.  The identification of (9.1) with
\(X_0(11)\) and of (9.5) with its modular \(j\)-map is the classical modular
input; the norm factorizations and passport are replayed directly by the
checker.  See Noam Elkies,
[*Elliptic Curves in Nature*](https://people.math.harvard.edu/~elkies/nature.html),
the conductor-\(11\) entry.

The positive-genus boundary lattice is nevertheless completely explicit.
Let \(Z\) be the reduced degree-four zero fiber in (9.4), let \(W\) be the
reduced degree-six fiber above \(1728\), and retain \(P=(16,-61)\) and \(O\).
Both \(Z\) and \(W\) are \(K\)-prime.  PARI proves rank zero for (9.1) and
its \((-11)\)-twist; reduction at the split primes \(3\) and \(5\) has order
five in both cases.  Hence

\[
 E(K)=\langle P\rangle\simeq\mathbb Z/5.           \tag{9.6}
\]

The Miller function

\[
 f_T=\frac{(Y+6X-35)^2(Y+5X-19)}{(X-5)^2}
\]

has divisor \(5P-5O\).  The checker also records explicit primitive units
\(g_Z,g_W\), with

\[
\begin{aligned}
 \operatorname{div}(g_Z)&=Z-2P-2O,\\
 \operatorname{div}(g_W)&=W-3P-3O,\\
 g_Z^3&=j f_T,\qquad
 g_W^2&=-(j-1728)f_T.                               \tag{9.7}
\end{aligned}
\]

In the ordered basis \((Z,W,P,O)\), the full unit lattice is

\[
 \boxed{
 L_{12}=\left\langle
 (1,0,-2,-2),\ (0,1,-3,-3),\ (0,0,5,-5)
 \right\rangle.}                                   \tag{9.8}
\]

Equivalently it is the simultaneous kernel of the degree row
\((4,6,1,1)\) and the class row \((2,3,1,0)\bmod5\).  Its Smith form is

\[
 \mathbb Z^4/L_{12}\simeq\mathbb Z\oplus\mathbb Z/5. \tag{9.9}
\]

This exposes a genuine saturation phenomenon.  The three evident units
\(j,j-1728,f_T\) have rows

\[
 (3,0,-11,-1),\quad(0,2,-11,-1),\quad(0,0,5,-5),
\]

and generate an index-six sublattice of \(L_{12}\).  The missing cube and
square roots are exactly \(g_Z\) and \(g_W\); the boundary geometry detects
more than the two target units alone.

The comparison with (6.6) is now sharper.  Both the natural quotient and
\(\widetilde C_5\) have genus one, but they arise from different subgroups:
the Borel \(C_{11}\!:\!C_5\) of order \(55\) and the intersection subgroup
\(A_4\) of order \(12\).  Section 6.4 proves that the first is the
conductor-\(11\) curve (9.1), while the second is the conductor-\(121\)
curve (6.17), and their traces above \(23\) prove that they are not
\(K\)-isogenous.  The natural map has degree \(12\) over the \(j\)-line,
whereas \(\widetilde C_5\) has degree \(55\) over the same triangle base and
degree five over either exceptional rational quotient.  The other
correspondence normalization has genus two.

All three positive-genus function fields remain nonrational after adjoining
independent variables.  Hence neither the natural rigid input nor either
normalized correspondence can be turned into affine space by stable
straightening.  Equations (6.21c), (6.38), and (9.8) now add the complete
boundary-unit lattices to that stable-birational obstruction.

This is a cover-specific obstruction, not an absolute no-go theorem for the
action \(G\curvearrowright\mathbb P^1(\mathbb F_{11})\).  A different
higher-dimensional \(G\)-cover might have rational source and a different
boundary ledger.  The action-spectrum problem asks whether any such cover
admits absolute Keller completion.

## 10. Essential dimension and Hurwitz rigidity are different axes

The two \((3,2,11)\) Nielsen classes are rigid: after fixing the three branch
values, their Hurwitz dimension is zero.  Nevertheless the currently cited
birational literature records

\[
 3\le \operatorname{ed}_{\mathbb C}
 \bigl(\operatorname{PSL}_2(11)\bigr)\le4             \tag{10.1}
\]

and leaves the exact value open.  See Beauville,
[*Finite simple groups of small essential dimension*](https://arxiv.org/abs/1101.1372),
and Prokhorov,
[*Finite groups of birational transformations*](https://ems.press/content/book-chapter-files/28268),
Section 7.

There is no contradiction.  A rigid cover is one \(G\)-torsor over one
function field.  Essential dimension measures the least dimension of a
versal compression capable of representing all \(G\)-torsors after field
extension.  Thus (10.1) belongs on a versal Keller card, not as an
existential obstruction to (7.2).

For the wider family \(q=p^f\), the unconditional local essential-dimension
lower bound grows at least with \(f\).  Brosnan--Reichstein--Vistoli compute
the local essential dimension of \(\operatorname{PSL}_2(q)\) over
\(\mathbb C\) as \(f\) for even \(q\) and \(\max(2,f)\) for odd \(q\):
[*Essential Dimension in Mixed Characteristic*](https://doi.org/10.4171/DM/653),
Proposition 5.2.  This makes the natural \(\operatorname{PSL}_2(q)\) family a
meaningful test of whether a future Keller construction is merely
existential or genuinely versal.

## 11. A falsifiable next queue

The action spectrum changes the next targets.

1. **Pullback/cokernel ledger -- completed.**  Equations (6.40)--(6.48)
   compute both compact pullback matrices, prove that their rank-ten sums
   are primitive, and separate the index-six saturation on \(X_0(11)\) from
   the order-five divisor-class gate on \(C_6\).
2. **Normalization-module masks -- completed.**  Equations
   (6.49)--(6.62) determine the integral exchange modules, give effective
   divisor bases for every residual class, construct exact representatives
   in both normalization algebras, and isolate the descended \(C_6\)
   infinity-imbalance character.  They also force normal degrees two and
   three in the declared one-class-per-monomial architecture.  The next
   coefficient calculation is now an ambient one: clear the normalization
   denominators compatibly with the output block, then impose the
   constant-Jacobian and inverse-adjugate equations.  Any ansatz retaining
   \(P_r\) plus only one triangular normal coordinate is already closed.
3. **Compare \(q=7,11,13\) in the natural action -- genus audit
   completed.**  Formula (4.4) gives genera \(0,1,0\), controlled by the
   two Legendre symbols for the order-two and order-three fixed points.
   The next useful comparison is the boundary-unit and affine-completion
   ledger for explicit degree-eight and degree-fourteen rational Belyi
   presentations; genus zero alone supplies no Keller realization.
4. **Test Gassmann saturation.**  Ask whether
   \((G,G/H_1)\in\mathfrak K_{\mathrm{abs}}\) and
   \(\mathbf1_{H_1}^G=\mathbf1_{H_2}^G\) force
   \((G,G/H_2)\in\mathfrak K_{\mathrm{abs}}\).  The derivative suspension
   proves simultaneous chart realizability once a common closure is given;
   absolute affine-space descent need not respect Gassmann equivalence.

The first new invariant is therefore not a single number.  It is the filtered
action spectrum together with its minimal Hurwitz genus, stable birational
type, boundary-unit defect, and versal parameter cost.  The
\(\operatorname{PSL}_2(11)\) benchmark shows that all four coordinates are
independent enough to matter.

## 12. Reproduction

Run

```bash
.venv/bin/python scripts/verify_psl2_11_keller_action_spectrum.py
.venv/bin/python scripts/verify_psl2_11_keller_action_spectrum.py \
  --pari-mordell-weil
.venv/bin/python scripts/verify_psl2_11_keller_action_spectrum.py \
  --singular-normalization
.venv/bin/python scripts/verify_psl2_11_normalization_masks.py
Singular -q scripts/psl2_11_c5_boundary_pullbacks.sing
.venv/bin/python scripts/verify_psl2_11_c6_boundary_images.py
```

The first command is the dependency-light exact replay.  The second requires
PARI/GP and certifies the rational ranks and Heegner generators used by the
quadratic rank decompositions.  The third requires Singular \(4.4.1\) and
verifies normalization, conductor, all normalized boundary degrees, the
canonical elimination for \(C_6\), and the two quadratic Cremona
transformations for \(C_5\), followed by every normalization-module mask.
The fourth is the robust standalone mask replay.  It pins the Singular source
hash, rejects interpreter diagnostics, and constructs the filtered
Riemann--Roch spaces, the prime-power evaluation kernels, and the infinity
jets in (6.54)--(6.62).  The fifth labels every affine \(C_5\) boundary
component by both projection colors and separates its node branches.  The
sixth is the more expensive exact Singular/SymPy residue-field replay of
the class rows (6.36).  Together they verify:

1. the corrected Shabat factorization (6.2), derivative (6.3), and square
   discriminant (6.4);
2. the exact bidegree-five/six correspondence factorization (6.5);
3. irreducibility of all six boundary factors over
   \(\mathbb Q(\sqrt{-11})\);
4. the rank-six/rank-two unit pullback and non-descending derivative vector;
5. the order and conjugacy-class census of \(\operatorname{PSL}_2(11)\);
6. the two nonconjugate classes of eleven \(A_5\) subgroups;
7. equality of the two degree-eleven cycle partitions for every group
   element and the cross-subdegrees \(5+6\);
8. the two rigid \((2,3,11)\) Nielsen orbits;
9. the natural prime-triangle genus formula (4.4), including genera
   \(0,1,0\) at \(p=7,11,13\), and genus zero in the exceptional
   degree-eleven actions versus genus one in the natural degree-twelve action;
10. the \(A_4\) and \(D_{10}\) intersection subgroups, degree-\(55/66\)
    passports, genera \(1/2\), and both boundary-projection profiles;
11. the degree-two/degree-five divisor argument giving
    \(\widetilde C_5(K)\ne\varnothing\), its sparse cubic, \(j=-121\),
    exact descent to (6.17), and the trace obstruction to an isogeny with
    \(X_0(11)\);
12. the canonical adjoint pencil of \(C_6\), its exact even hyperelliptic
    model, its two elliptic quotients, and the conductor-\(11\) times
    CM-conductor-\(121\) decomposition (6.31);
13. the exact \(X_0(11)\) equation, norm factorizations, degree-twelve
    \(j\)-passport, primitive boundary units, and index-six saturation;
14. the full positive-genus unit lattices of ranks \(3,14,17\), including
    their Smith quotients and the rational-kernel audit for the
    \((2,2)\)-isogeny;
15. both compact correspondence pullback matrices, their primitive
    rank-ten images, the free cokernels of ranks four/seven, and the weak
    versus strong two-output ledger distinction;
16. the residual exchange modules
    \(\mathbb Z^4_{\mathrm{triv}}\) and
    \(\mathbb Z^5_{\mathrm{triv}}\oplus\mathbb Z[\mathbb Z/2]\), the
    effective mask bases, the index-two simple-pole lattice on \(C_5\), and
    the minimal normal-support degrees two and three;
17. the exact \(C_5\) formulas (6.55)--(6.57), all seven \(C_6\)
    normalization-module kernel representatives, the infinity-five jet
    selection of \(N_1\), and the descended infinity imbalance (6.62);
18. with the PARI flag, the rank-zero statements for \(X_0(11)\),
    \(E_{121}\), and \(E_+\), and the rank-two decomposition of \(E_-(K)\);
19. with the optional Singular flag, the five/eight ordinary nodes, the
    equality of conductor and reduced singular ideals, the canonical
    discriminant elimination for \(C_6\), and the two-step Cremona reduction
    of \(C_5\), together with the complete mask replay; and
20. with the boundary-image replay, the nineteen direct elliptic residue
    traces and the remaining \(q_2\)-node relation that produce (6.36).

The checker is an exact finite calculation and symbolic replay.  Jones and
Zvonkin's monodromy identification and common-regular-cover theorem remain
the cited external bridge from the displayed polynomial to the abstract
group calculation.
