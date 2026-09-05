# F2 geometric-degree-six Stein and cubic-germ reduction

> **Status.**  This note proves an unconditional localization theorem at
> geometric degree six and a complete classification under the additional
> hypothesis that the terminal Cartier slice in the local degree-three Stein
> factor is normal (equivalently here, conductor-free).
> Every affine nonproperness component punctures the target boundary at the
> unique terminal value `125/729`, and all of its local branching is carried
> by the single three-sheet packet there.  For a normal terminal slice the germ
> has the one-parameter cubic normal form
> `z=w^3+a(pi)w`.  It has exactly three possible reduced ramification
> patterns: one `k=1` three-cycle row, two `k=1` transposition rows, or one
> `k=2` transposition row.  In particular a mixed transposition/three-cycle
> pair cannot occur.  Every one of these patterns reproduces the already
> known terminal logarithmic cokernel `R/(w^3)`, so that cokernel does not
> distinguish them.  The fixed squarefree carrier further forces
> `ord_pi(a)>=5`.  Consequently the two-component normal row has only
> exponents `6,8,10`, while the one-component normal cusp has only exponents
> the numerical candidates `5,7,9`, with exact infinity semigroups and affine
> delta budgets below. The [sharp polynomial parametrization gap bound](F2_DEGREE_6_10_POLYNOMIAL_GAP.md)
> subsequently excludes `r=9`: a birational pair of degrees `(6,10)` has
> first odd gap at most 21, whereas that row requires 27. Thus only `r=5,7`
> remain among the normal odd rows.
> The two even target quintics satisfy an exact sparse-difference theorem:
> the `r=6,8,10` rows have respectively conic, line, and constant implicit
> differences.  The constant-difference locus is an explicit one-parameter
> family of disjoint four-node quintics.  Its affine branch divisor is SNC,
> so tame local monodromy contributes zero to Orevkov's local excess, whereas
> the degree-six identity requires excess one.  Thus the entire normal
> `r=10` row is excluded.  Exact Sage/SIROCCO braid monodromy and exhaustive
> transposition enumeration exclude the four algebraic `r=8` cusp pairs.
> The `r=6` cusp incidence is one irreducible rational surface: its ordinary
> cusp stratum reduces by tame-split regeneration to one certified rational
> witness, while its only locally admissible higher-cusp points are two
> certified `E_6` scale classes.  None has a transitive six-sheet
> transposition action.  Hence all three normal even rows `r=6,8,10` are
> excluded.
> Nonnormal terminal slices remain open, even when the
> ambient Stein surface is smooth; their positive normalization conductor is
> precisely the
> point information missed by the terminal determinant packet.  Thus this is
> a finite reduction, not an exclusion of degree six or of `(75,125)`.

The residue arithmetic, cubic discriminants, parity classification,
high-contact quintic saturation, endpoint blowup models, and three basic
conductor orders are replayed by
[`verify_f2_geometric_degree_six_stein_reduction.py`](../scripts/verify_f2_geometric_degree_six_stein_reduction.py).
The high-contact saturation step uses Singular; all displayed identities
around it are also checked directly in SymPy.

## 1. The terminal fiber exhausts degree six

Let

\[
 F:\mathbb A^2\longrightarrow\mathbb A^2
\]

be a hypothetical F2 `(75,125)` Keller map of geometric degree `d=6`, and
let

\[
 X\longrightarrow Z\mathop{\longrightarrow}^{\bar F}Y          \tag{1.1}
\]

be the Stein factorization on a proper resolved compactification.  The
finite map `bar F` has degree six.  The certified terminal divisor maps to
the target boundary component with residue cover

\[
 h(s)=\frac{125s(s+1)^5}{(9s^2+15s+5)^3}.                       \tag{1.2}
\]

Its finite nonzero branch value and its fiber partition are

\[
 \lambda_0=\frac{125}{729},\qquad
 h^{-1}(\lambda_0)=(3,1,1,1).                                  \tag{1.3}
\]

At every other finite nonzero value the partition is `(1^6)`.

Fix a target-boundary point `y_lambda` with residue `lambda`.  For every
point `z_i` of the finite Stein fiber above it, let `N_i` be the local degree
of `bar F`.  A terminal residue point of local index `q_i` contributes at
least `q_i` to `N_i`.  Therefore

\[
 6=\sum_iN_i\ge\sum_iq_i=6.                                   \tag{1.4}
\]

Every inequality in (1.4) is an equality.  There are no additional local
sheets over the point, and each terminal residue point has local surface
degree exactly its residue index.

There is one separation point hidden in this statement.  Two distinct
points of the terminal normalization cannot be identified by `X->Z`.  Such
an identification would put both points in one connected Stein fiber, hence
would require a connected chain of contracted curves joining them.  Every
curve over the target boundary lies in the source boundary.  But this
boundary is obtained from the line at infinity by point blowups, so its dual
graph is a tree; a chain joining two distinct points of the same horizontal
terminal vertex would create a cycle.  Therefore the terminal points remain
distinct in `Z`.

If `lambda` is finite, nonzero, and different from `lambda_0`, the six local
degrees are all one.  A finite flat local algebra of rank one is the base
ring, so all six germs are etale.  No affine branch divisor can meet the
target boundary there.  Consequently:

\[
 \boxed{\text{every affine nonproperness component punctures at }
        \lambda_0=125/729.}                                   \tag{1.5}
\]

At `lambda_0`, completion of the finite algebra splits it into three rank-one
etale factors and one local rank-three factor.  Hence every affine branch
germ lies in the branch locus of that same cubic factor.  The three simple
sheets cannot carry it.

This argument uses only finite-degree conservation and the terminal
passport.  It does not assume that the cubic Stein point is smooth.

## 2. Immediate global consequences

A meridian of an irreducible branch divisor in a three-sheet local cover is
either a transposition or a three-cycle.  Its moved-sheet cost is therefore
`2` or `3`.  Orevkov's residue-degree identity gives, at `d=6`,

\[
 \sum_j A_j\le d-1=5.                                          \tag{2.1}
\]

Thus there are at most two affine nonproperness components, and two
three-cycle components are impossible.  Before using the local algebra, the
only numerical multisets are

\[
 (2),\quad(3),\quad(2,2),\quad(2,3).                           \tag{2.2}

Section 3 removes the last multiset whenever the terminal Cartier slice of
the cubic Stein factor is normal.

## 3. Normal terminal Stein slice: complete local classification

Assume now that the terminal Cartier slice of the local degree-three Stein
factor is normal.  Its completed local ring is `k[[w]]`, with target residue
parameter `z=w^3`.  Since this regular Cartier section lies in a
two-dimensional Cohen--Macaulay local ring, the Stein point itself is
regular and `pi` is a source parameter.  The characteristic-zero
one-parameter `A_2` normal form then gives formal or analytic coordinates in
which

\[
 \boxed{(\pi,w)\longmapsto
        (\pi,z=w^3+a(\pi)w),\qquad a(0)=0.}                    \tag{3.1}
\]

A target translation has removed the constant term.  The ramification and
branch equations are

\[
 3w^2+a(\pi)=0,qquad
 \Delta(\pi,z)=-4a(\pi)^3-27z^2=0.                            \tag{3.2}

There are three and only three cases.

### 3.1 Cyclic specialization

If `a=0`, the reduced ramification curve is `w=0`, with generic index three.
Its image is the smooth branch `z=0`, transverse to `pi=0`.  It is one
`k=1` component with cubic inertia:

\[
 (A,k)=(3,1).                                                  \tag{3.3}

### 3.2 Even order

Suppose `r=ord_pi(a)<infinity` is even.  After adjoining a harmless unit
square root, (3.2) has two smooth ramification branches

\[
 w=\mathord\pm c\pi^{r/2}.                                    \tag{3.4}

Each is simply ramified and each has `pi` itself as normalization parameter.
Their images are two distinct `k=1` target branches.  Thus the pattern is

\[
 (A_1,k_1)+(A_2,k_2)=(2,1)+(2,1).                             \tag{3.5}

### 3.3 Odd order

If `r=ord_pi(a)` is odd, the reduced ramification curve is irreducible.  Its
normalization has

\[
 \pi=t^2,\qquad w=t^r\cdot\text{unit},\qquad
 z=t^{3r}\cdot\text{unit}.                                   \tag{3.6}

It is generically simply ramified, while its image has boundary contact
`k=2`.  The pattern is

\[
 (A,k)=(2,2).                                                  \tag{3.7}

For `r=1` the image is the ordinary cusp `(pi,z)=(t^2,t^3)`.
For larger odd `r` it is the higher unibranch germ `(2,3r)`.

Combining the three cases gives the promised exhaustive smooth-point list:

\[
 \boxed{
 (3;k=1),\qquad(2,2;k=1,1),\qquad(2;k=2).}                    \tag{3.8}

In particular the numerical possibility `(2,3)` in (2.2) is not the branch
divisor of one smooth cubic germ.

### 3.4 The carrier Newton polygon forces `r>=5`

At degree six the double-carrier row is already unavailable: its two
terminal degree-six packets force geometric degree at least twelve.  Hence
the carrier is the certified squarefree row.  The divisorial valuation over
the target ray `(5,36)` has

\[
 v(\pi)=5,\qquad v(z)=36,                                  \tag{3.9}
\]

and the cubic factor has `(e,f)=(1,3)` there.  Apply `v` to the three terms
of

\[
 w^3+\pi^r w-z.
\]

The segment joining `(0,36)` and `(3,0)` in its Newton polygon has height
`24` at horizontal coordinate one.  The middle term has height `5r`.
An unramified residue-degree-three extension requires the one-segment
polygon, with root value `12`; otherwise the cubic valuation splits into
smaller packets.  Therefore

\[
 5r>24,
 \qquad\boxed{r\ge5}.                                       \tag{3.10}
\]

For `r>=5`, after scaling `w` by value `12`, the residual cubic is the
irreducible carrier cubic; the middle term has strictly larger value.

### 3.5 Exact genus and semigroup ledgers

Suppose first that `r` is odd.  On the normalization of the target branch,

\[
 (\pi,z)=(t^2,t^{3r}\cdot\text{unit}).
\]

Returning through the `(5,2)` toric chart gives Puiseux characteristic
exponents

\[
 (4;10,10+3r)
\]

at the unique projective point at infinity.  Equivalently its local value
semigroup is

\[
 \langle4,10,20+3r\rangle.                                  \tag{3.11}
\]

Its conductor and delta invariant are therefore

\[
 c_\infty=27+3r,
 \qquad\delta_\infty=\frac{27+3r}{2}.                       \tag{3.12}
\]

The curve has degree ten and rational normalization, so its total affine
delta is `36-delta_infinity`.  Polynomiality supplies one more restriction.
The Abhyankar delta sequence at infinity starts with `(10,6)` and, because
(3.11) has two characteristic pairs, has one last value `s`.  Its conductor
is `15+s`.  This equals twice the genus of a smooth member of the one-place
pencil, namely `2(36-delta_infinity)`.  Hence `s=30-3r`.

A delta-sequence value is positive.  Combining this with (3.10), parity,
and Chau's forced affine singularity leaves exactly

\[
\begin{array}{c|c|c|c|c}
r&\Gamma_\infty^{\rm local}&\text{Abhyankar delta sequence}
 &\delta_\infty&\sum_{p\in C\cap\mathbb A^2}\delta_p\\ \hline
5&\langle4,10,35\rangle&(10,6,15)&21&15\\
7&\langle4,10,41\rangle&(10,6,9)&24&12\\
9&\langle4,10,47\rangle&(10,6,3)&27&9.
\end{array}                                                   \tag{3.13}
\]

Thus the apparent odd exponents `11,13` allowed by the genus inequality are
not polynomial one-place curves of degrees `(6,10)`.

The last row of (3.13) is also impossible, by the stronger
[polynomial parametrization theorem PF2D6O1](F2_DEGREE_6_10_POLYNOMIAL_GAP.md).
The numerical delta-sequence conditions do not guarantee parametrizability
in the prescribed degrees. An exact ideal-membership certificate forces
every degree-`(6,10)` pair with odd Puiseux gap greater than 21 to factor
through a quadratic polynomial. The `r=9` gap is 27, contradicting a
normalization parametrization. The `r=7` gap 21 is attained by an explicit
birational polynomial pair, so this bound leaves precisely `r=5,7` here.
The original numerical table is retained as an input ledger to that proof.

Now let `r` be even and put `m=3r/2`.  The two target components are
rational quintics.  Each has the fixed `(2,5)` infinity cusp, so each has
`delta_infinity=2` and total affine delta four.  In the toric chart the two
smooth strict transforms have mutual contact `m`.  The inverse monomials

\[
 a=\pi^5h^{-3},\qquad b=\pi^2h^{-1}
\]

show, by the two conjugate Puiseux substitutions, that their original
intersection multiplicity at infinity is

\[
 I_\infty(C_+,C_-)=10+m.                                    \tag{3.14}
\]

Bezout gives `I_infinity<=25`.  Together with (3.10), this leaves

\[
\begin{array}{c|c|c|c}
r&m&I_\infty(C_+,C_-)&I_{\mathbb A^2}(C_+,C_-)\\ \hline
6&9&19&6\\
8&12&22&3\\
10&15&25&0.
\end{array}                                                   \tag{3.15}
\]

Equations (3.13) and (3.15) are complete numerical classifications of the
finite normal-form rows.  They are not existence statements: the listed
curves must still support the required connected six-sheet cover and the
Keller filling.

### 3.6 The odd rows are Orevkov's unique clean cubic germs

The odd normal form is not an unclassified analytic cover.  Let
`a=pi^r` after absorbing a unit.  Pulling the target discriminant back to
the source gives the exact factorization

\[
 -4\pi^{3r}-27z^2
 =-(3w^2+\pi^r)^2(3w^2+4\pi^r).                \tag{3.16}
\]

Thus the ramification curve

\[
 R:(3w^2+\pi^r=0)
\]

and the residual unramified curve

\[
 S:(3w^2+4\pi^r=0)
\]

are both irreducible `(2,r)` cusps.  Each maps bijectively to the target
`(2,3r)` cusp, and

\[
 I_0(R,S)
 =\ell k[[\pi,w]]/(w^2,\pi^r)=2r.             \tag{3.17}
\]

This is Case (b) of Orevkov's 2026 germ theorem with parameters

\[
 (k_1,k_2;l_1,l_2)=(2,r;1,0),
 \qquad(N,n)=(3,2),                             \tag{3.18}
\]

after swapping the two target exponents.  Since quadratic ramification is
unique in that theorem, (3.16) exhausts the smooth, bijective, irreducibly
ramified local cover; there is no second clean analytic escape.

Orevkov's extra-property calculation also assigns this germ the resolved
ramification--residual intersection charge

\[
 Q_r=2\delta(R)+\operatorname{mult}(C)-1=r.    \tag{3.19}
\]

For the three numerical rows, `Q_r` is therefore `5,7,9`; PF2D6O1 separately
excludes the last. Proving that
this positive charge injects into the root-subtracted localized
`ch_2`/`Fitt_1` point budget would exclude all three rows at once.  Such a
comparison is not presently proved: (3.19) is a resolved fiber-product
intersection, whereas the global logarithmic remainder is a virtual class.
Equation (3.19) identifies the exact missing bridge rather than assuming it.

### 3.7 Exact logarithmic module before SNC resolution

There is a useful warning for that bridge.  Put

\[
 h=3w^2+\pi^r,qquad D_{\rm src}=(\pi h=0).
\]

The reduced plane divisor `D_src` is free.  A Saito basis of logarithmic
derivations is

\[
 E=2\pi\partial_\pi+rw\partial_w,
 \qquad
 H=6\pi w\partial_\pi-
   (3w^2+(r+1)\pi^r)\partial_w.                \tag{3.20}
\]

Indeed their coefficient determinant is `-2(r+1)pi*h`.  Evaluate the
target logarithmic forms `(dlog(pi),dz)` on `(E,H)` for
`z=w^3+pi^r w`.  The resulting matrix has first row `(2,6w)` and

\[
 \det\Theta=-2(r+1)h^2.                        \tag{3.21}
\]

Because the entry `2` is a unit, elementary operations give the exact
module

\[
 \boxed{
 \operatorname{coker}\Theta\simeq
 k[[\pi,w]]/(3w^2+\pi^r)^2,
 \qquad\operatorname{Fitt}_1(\operatorname{coker}\Theta)=k[[\pi,w]].}
                                                               \tag{3.22}
\]

Thus the unresolved clean cusp has no isolated `Fitt_1` correction at all.
For odd `r`, the conductor length `(r-1)/2`, the intersection length `2r`
in (3.17), and the resolved charge `r` in (3.19) are normalization and
resolution data attached to the singular support of this cyclic module.
They are lost if (3.22) is replaced by its scalar determinant alone.

This also explains why the endpoint packet remains `R/(w^3)` for every
`r`: SNC resolution transports one cyclic divisor packet, but the scalar
endpoint does not remember how its singular support was normalized.  A
global proof must compare the resolved filtration with (3.19); searching
for an isolated raw `Fitt_1` point cannot close these rows.

### 3.8 Forced SNC chain and cyclic Smith atlas

For odd `r=2s+1`, resolve the reduced source boundary
`pi(3w^2+pi^r)=0`.  The minimal embedded SNC resolution has exceptional
valuation rays

\[
 (1,1),(1,2),\ldots,(1,s),(1,s+1),(2,r).       \tag{3.23}
\]

The last ray is produced by blowing up the final triple point.  Its divisor
is trivalent and has self-intersection `-1`; the preceding `(1,s+1)`
divisor is a `(-2)` leaf, the `(1,s)` divisor has weight `-3`, and the
earlier chain divisors have weight `-2`.  Its final vertex has the same
trivalent `(-1)` combinatorial type as the incidence slot left open by the
compiled carrier graph.  Equality of those vertices is a global attachment
identification and is not assumed here.

Track the logarithmic determinant divisor through these blowups.  At a
center of reduced boundary multiplicity `m`,

\[
 K_{X'}+D_{X'}=g^*(K_X+D_X)+(2-m)E.
\]

Starting from the double ramification support `2R` in (3.22), the
exceptional determinant coefficients are

\[
 3,6,\ldots,3s,\quad3s+2,\quad3r.              \tag{3.24}
\]

Every exceptional divisor is contracted to the target point **on the
logarithmic boundary**.  This differs from contraction to an ordinary
affine target point: if `E` has `ord_E(pi)=a>0`, then the pullback of
`dlog(pi)` has normal coefficient `a`, a unit in characteristic zero.
Consequently the logarithmic matrix has a unimodular entry on every
exceptional divisor, and the generic Smith atlas is

\[
\begin{array}{c|c|c}
\text{ray}&\text{determinant order}&
 \text{generic Smith exponents}\\ \hline
(1,j),\ 1\le j\le s&3j&(0,3j)\\
(1,s+1)&3s+2&(0,3s+2)\\
(2,r)&3r&(0,3r).
\end{array}                                                    \tag{3.25}
\]

Thus SNC resolution preserves generic cyclicity and the generic
`Fitt_1` ideal remains the unit ideal.  For the three live rows the cyclic
determinant-order lists are

\[
 r=5:(3,6,8,15),\quad
 r=7:(3,6,9,11,21),\quad
 r=9:(3,6,9,12,14,27).                         \tag{3.26}
\]

This removes the source-graph ambiguity in the normal cusp case.  What
remains is global: identify this forced chain with the compiled carrier
vertices and retain the normalization/conductor filtration of its singular
determinant support.  A rank-two contracted-divisor correction would be the
wrong local model here.

### 3.9 Why the scalar global `ch_2` number cancels the exponent

The same resolution gives a decisive limitation of the numerical Chern
route.  Among its `s+2` blowups, exactly `s+1` centers have reduced boundary
multiplicity three; the remaining center is the tangency of two reduced
branches and has multiplicity two.  Therefore

\[
 \Delta L_X^2=-(s+1)=-\frac{r+1}{2}.           \tag{3.27}
\]

Since the target model is unchanged and
`D_log=L_X-f^*L_Y`, the determinant square changes by the same amount:

\[
 \Delta D_{\log}^2=-\frac{r+1}{2}.             \tag{3.28}
\]

Consequently the two changes cancel in the complete-chain combination:

\[
 \Delta\left(B_f-\frac12D_{\log}^2\right)=0.  \tag{3.29}
\]

The losses in `L_X^2` are respectively `3,4,5` for `r=5,7,9`, but none
changes the scalar residual `d-1-A=3`.  Thus neither more Laurent
coefficients nor the unfiltered virtual `ch_2` degree can distinguish the
three normal cusp rows.  A closure must use the normalization filtration
behind (3.19), a global cover/filling obstruction, or the remaining
coefficient equations.

### 3.10 The even rows have the parallel tangency chain

For even `r=2s`, the ramification equation factors into two smooth branches

\[
 R_\pm:\quad w=\mathord\pm c\pi^s,
 \qquad I_0(R_+,R_-)=s.                         \tag{3.30}
\]

Resolving `T+R_++R_-` requires exactly the rays

\[
 (1,1),(1,2),\ldots,(1,s).                     \tag{3.31}
\]

Every center has reduced boundary multiplicity three.  The exceptional
chain has weights `(-2,...,-2,-1)`; the final `(-1)` divisor is trivalent,
meeting the preceding chain and the two separated ramification branches at
three distinct points.  The determinant orders and generic log Smith types
are

\[
 \operatorname{ord}_{E_j}\det\Theta=3j,
 \qquad\operatorname{Smith}_{\eta_{E_j}}(\Theta)=(0,3j).
                                                               \tag{3.32}
\]

Thus the three rows have determinant lists

\[
 r=6:(3,6,9),\quad
 r=8:(3,6,9,12),\quad
 r=10:(3,6,9,12,15).                          \tag{3.33}
\]

Here `Delta L_X^2=Delta D_log^2=-s`, so the scalar cancellation (3.29)
holds again.  The normal-slice classification is therefore complete on
both parity branches: every unresolved datum is now global gluing,
normalization filtration, or affine complement topology.

### 3.11 High contact forces conic, line, and constant differences

The even rows admit a second exact classification which does not use any
uncompiled Laurent coefficient.  Normalize the first rational quintic as

\[
 p_1=t^3+at,
 \qquad q_1=t^5+bt^4+ct^2+dt,                 \tag{3.34}
\]

and write its implicit equation as `F_1(P,Q)=0`.  After an independent
affine change of the parameter on the second component, equality of the
leading infinity coefficient puts it in the form

\[
\begin{aligned}
 p_2&=t^3+At+p_0,\\
 q_2&=t^5+Bt^4+Et^3+Ct^2+Dt+q_0.
\end{aligned}                                  \tag{3.35}
\]

The equality of the leading coefficient follows already from contact
greater than ten at the common `(2,5)` infinity branch.  The coprimality of
three and five then permits the simultaneous monic parameter normalization
used in (3.35).

Pull `F_1` back along `(p_2,q_2)`.  Its coefficients of
`t^14,t^13,t^12,t^11,t^10` vanish successively and uniquely.  They give

\[
\begin{aligned}
 B&=b, & E&={5(A-a)\over3},\\
 C&={4Ab-4ab+3c+5p_0\over3},
 &D&={5A^2-15Aa+10a^2+12bp_0+9d\over9},\\
 q_0&={2A^2b-8Aab+6Ac+10Ap_0+6a^2b-6ac-15ap_0\over9}.
\end{aligned}                                  \tag{3.36}
\]

Let `F_2` be the exact implicit quintic of (3.35), after removing its
`Et^3` term by the target shear `Q-E(P-p_0)`.  Substitution of (3.36) gives
the sparse identity

\[
 \boxed{
 F_2-F_1=h_{30}P^3+h_{20}P^2+h_{11}PQ
          +h_{10}P+h_{01}Q+h_{00}.}            \tag{3.37}
\]

There are no other monomials.  In the projective chart `Q=1` at infinity,
the normalization has `ord(P,W)=(2,5)`.  Homogenizing (3.37) to degree five
therefore gives the exact order table

\[
\begin{array}{c|c|c}
\text{term}&I_\infty&\deg_t\text{ of its affine pullback}\\ \hline
P^3W^2&16&9\\
PQW^3&17&8\\
P^2W^3&19&6\\
QW^4&20&5\\
PW^4&22&3\\
W^5&25&0.
\end{array}                                                   \tag{3.38}
\]

The orders are all distinct, so cancellation between rows of this table is
impossible.  Comparing (3.38) with (3.15) proves

\[
\begin{array}{c|c|c}
r&I_\infty&\text{forced implicit difference}\\ \hline
6&19&F_2-F_1=\alpha P^2+\beta P+\gamma Q+\delta,
             \quad\alpha\ne0,\\
8&22&F_2-F_1=\beta P+\delta,\quad\beta\ne0,\\
10&25&F_2-F_1=\delta,\quad\delta\ne0.
\end{array}                                                   \tag{3.39}
\]

Thus all three affine intersections in the middle row lie on one target
line.  In the last row the two curves are distinct fibers of one polynomial
pencil.

The last row can be solved completely.  Put `x=A-a`.  Saturating
`h_30,h_11,h_20,h_01,h_10` by `h_00` gives the prime ideal

\[
\begin{aligned}
 125x^2&=81b^4,&p_0&=-{bx\over5},\\
 a&=-{3b^2+5x\over10},
 &c&=-{b(b^2+5x)\over10},\\
 d&=-{x(81b^2+100x)\over405}.                 \tag{3.40}
\end{aligned}
\]

Moreover

\[
 h_{00}={81b^{13}x\over5^{11}},                \tag{3.41}
\]

so `h_00!=0` forces `b x!=0`; no omitted component lives on the saturation
boundary.  Equivalently, write `b=beta`, `a=kappa beta^2`.  Then

\[
 125\kappa^2+75\kappa-9=0,                     \tag{3.42}
\]

and the first curve has

\[
 c=(\kappa+1/5)\beta^3,
 \qquad d={(10\kappa-1)\beta^4\over25}.        \tag{3.43}
\]

The second coefficients are

\[
\begin{aligned}
 A&=-{(5\kappa+3)\beta^2\over5},&B&=\beta,
 &E&=-{(10\kappa+3)\beta^2\over3},\\
 C&=-{(5\kappa+2)\beta^3\over5},
 &D&={(15\kappa+14)\beta^4\over25},\\
 p_0&={(10\kappa+3)\beta^3\over25},
 &q_0&=-{(25\kappa+48)\beta^5\over375}.
\end{aligned}                                  \tag{3.44}
\]

After the translation and shear used above, its normalized parameter is
`kappa'=-3/5-kappa`, the other root of (3.42), and

\[
 \boxed{F_2-F_1=-{81\beta^{15}(10\kappa+3)\over5^{12}}\ne0.} \tag{3.45}
\]

This family is everywhere nodal on the open (3.45), not merely generically
nodal.  Normalize `beta=1` and let

\[
 R(u)=u^4+u^3+\kappa u^2+(\kappa-1/5)u
      -\kappa^2-(10\kappa-1)/25                \tag{3.46}
\]

be the collision quartic.  Reduction modulo (3.42) gives

\[
\begin{aligned}
 \operatorname{Disc}R&=-{81(200\kappa-21)\over390625},\\
 \operatorname{Res}(R,3u^2+4\kappa)
   &={81(25\kappa+21)\over3125},\\
 \operatorname{Res}(R,J_{\rm tan})
   &=-{177147(25\kappa-3)\over1220703125}.
\end{aligned}                                  \tag{3.47}
\]

Every right-hand side is coprime to (3.42).  Hence the four collisions are
distinct, off the diagonal, and tangent-separated.  Their common target
values are

\[
 P(u)=-u(\kappa+u^2),\qquad
 Q(u)={(\kappa+u^2)(10\kappa u+5u^3-1)\over5}. \tag{3.48}
\]

The exact saturation of `R(u)=R(v)=0`, equality of both values in (3.48),
and `u!=v` is the unit ideal.  Thus each curve has exactly four distinct
ordinary affine nodes.  Equation (3.45) makes the two curves disjoint.

### 3.12 The entire normal `r=10` row is impossible

Use the residue-normalized form of Orevkov's identity.  For every affine
dicritical row `l`, with transverse and residue degrees `(e_l,f_l)`, and
for every finite point `x` on its affine normalization, put

\[
 \epsilon_{l,x}=\mu_x-e_lq_x\ge0,
\]

where `q_x` is the local degree of the residue map.  At geometric degree
six,

\[
 \boxed{5=\sum_l e_lf_l+\sum_{l,x}\epsilon_{l,x}.}           \tag{3.49}
\]

The following tame-local observation is the decisive point.

**SNC excess lemma.**  At a smooth or ordinary-SNC point of the reduced
branch divisor of a finite normal characteristic-zero surface cover, a
simple-inertia, residue-degree-one row has `epsilon=0`.

Indeed, pass to a strict henselization of the regular target surface.  The
tame fundamental group of the complement of one or two coordinate branches
is respectively procyclic or a product of two procyclic groups.  Its local
meridian permutations therefore commute.  Two commuting transpositions are
equal or disjoint, so every local orbit touched by a given transposition has
exactly two sheets.  Consequently `mu_x=2=e_l q_x`.  Normality causes no
flatness escape here: a finite normal surface algebra over a regular
two-dimensional local ring is Cohen--Macaulay and hence locally free.  The
same argument applies separately to the two normalization points of a
self-node.

In the even normal row there are exactly two dicritical rows and both have
`(e,f)=(2,1)`.  Thus (3.49) requires

\[
 \sum_{l,x}\epsilon_{l,x}=5-2-2=1.            \tag{3.50}
\]

For `r=10`, Section 3.11 proves that the complete affine branch divisor is
the disjoint union of two four-node quintics.  Every one of its points is
smooth or ordinary SNC, so the SNC excess lemma makes the left side of
(3.50) zero.  This contradiction proves

\[
 \boxed{\text{the normal terminal-slice row }r=10
        \text{ cannot occur in a degree-six Keller map}.}    \tag{3.51}
\]

The proof covers the whole saturated locus (3.40), including both conjugate
values of `kappa`; it is not a generic-fiber argument.  For `r=8`, the same
budget shows that a surviving line-difference pair must contain a non-SNC
point carrying exactly one unit of excess.  For `r=6`, the analogous
survivor must lie on the non-SNC part of the conic-difference locus.  Those
are the next finite classifications.

### 3.13 The `r=8` row reduces to one ordinary-cusp quartic

There is a useful parity refinement of the SNC excess lemma.  Suppose a
target singular point has at least two normalization branches and some
local orbit touched by their transpositions has size `N>=3`.  Every one of
the `s` normalization points in that orbit contributes `N-2`, so the point
contributes

\[
 s(N-2)\ge2.                                      \tag{3.52}
\]

Orbits of size two contribute zero.  Thus the total budget one in (3.50)
can only be carried by one **unibranch** target singularity whose local
permutation orbit has size three.  Every multibranch singularity must be
tame-split.

This necessary locus is also exactly computable.  Keep the notation of
Section 3.11 and suppose first that the unibranch point lies on `C_1`.  The
normalization is nonimmersive exactly when

\[
 \operatorname{Crit}_1=
 25a^4+48a^3b^2-144a^2bc+90a^2d+108ac^2+81d^2=0. \tag{3.53}
\]

Saturate

\[
 h_{30}=h_{11}=h_{20}=h_{01}=\operatorname{Crit}_1=0
\]

by `h_10`.  The part with `b=0` is empty.  After the weighted normalization
`b=1`, the saturated ideal is a reduced quartic field.  Its primitive
equation is

\[
\boxed{
196000000a^4+260940000a^3+82362825a^2
-2390688a+20736=0.}                              \tag{3.54}
\]

All other coefficients are linear over this field:

\[
\begin{aligned}
c={}&-{51940000a^3+65374350a^2+15840099a-898128\over4361202},\\
d={}&-{53410000a^3+72989725a^2+21185232a-301824\over7268670},\\
x={}&-{7840000a^3+22668000a^2+25088409a+4472496\over5814936},\\
p_0={}&{42140000a^3+55884050a^2+19349013a+295344\over2422890}.
\end{aligned}                                                   \tag{3.55}
\]

The quartic is squarefree and has two real roots

\[
 -0.7452332634\ldots,qquad -0.6134611257\ldots,                \tag{3.56}
\]

and one conjugate nonreal pair.  Exact subresultants over (3.54) show that
`p_1'` and `q_1'` have exactly one common root.  The determinant of their
second and third derivative vectors is nonzero there, so this point is an
ordinary `A_2` cusp.  The collision quartic remains squarefree.  Since every
`k=1` quintic has affine delta four, the other three collision pairs are
three distinct ordinary nodes.  The critical resultant of `C_2` is coprime
to (3.54), and its collision quartic is squarefree, so `C_2` has four
ordinary nodes and no cusp.

Finally, restriction of the line `h_10P+h_00=0` to either normalization is
a squarefree cubic.  Its exact resultants with both collision quartics are
coprime to (3.54).  Therefore the three mutual intersections are distinct,
transverse, and avoid every self-singularity.  We have proved the uniform
configuration

\[
 \boxed{
 C_1:A_2+3A_1,qquad C_2:4A_1,qquad
 C_1\cap C_2=3\text{ transverse points}.}        \tag{3.57}
\]

Interchanging the two components gives the only other locus.  In
particular, the `A_4,E_6,E_8` critical degenerations do not occur on the
line-difference open.  As a separate finite check, the torus-knot groups of
`A_2,A_4,E_6,E_8` have respectively `6,0,6,0` labelled transitive actions
on three letters with geometric meridian a transposition; (3.54) selects
the first row.

The exploratory script
[`explore_f2_r8_cusp_braid.py`](../scripts/explore_f2_r8_cusp_braid.py)
uses the larger real root in (3.56), tracks the six points of a generic
vertical fiber around all twelve discriminant values, and ports Sage's
piecewise-linear crossing algorithm.  Exhaustive enumeration of the
`15^6` transposition tuples, reduced to `15^5=759375` by simultaneous
conjugacy, leaves exactly seven tuples.  They have the form

\[
 (\tau,\sigma,\tau,\sigma,\sigma,\tau),
 \qquad [\tau,\sigma]=1,                        \tag{3.58}
\]

in the base-fiber component order `(2,1,2,1,1,2)`.  Hence none is
transitive on six sheets and none gives the cusp a three-sheet orbit.
Repeating the computation with the other real embedding gives the identical
seven tuples.

The final replay is exact.  The Sage/SIROCCO verifier
[`verify_f2_r8_cusp_braid.sage`](../scripts/verify_f2_r8_cusp_braid.sage)
constructs the two implicit quintics over the quartic field (3.54), fixes
each of its four complex embeddings in turn, and computes twelve certified
Zariski--van Kampen braids.  SIROCCO returns disjoint tubes containing the
true continued roots, so the crossing words do not rely on floating root
matching.  For each embedding, exhaustive enumeration of the
`15^5=759375` transposition tuples leaves exactly the seven tuples (3.58).
All are intransitive.  A connected degree-six normal cover has transitive
monodromy, so no one of the four algebraic pairs can be its branch divisor.
Consequently

\[
 \boxed{\text{the entire normal terminal-slice row }r=8
        \text{ is impossible}.}                              \tag{3.59}
\]

### 3.14 The `r=6` cusp surface and its two `E_6` points

The conic-difference row initially looks larger, but its complete
excess-one incidence is rational.  Let `s` be the normalization parameter of
the unibranch point on `C_1`.  Criticality gives

\[
 a=-3s^2,qquad d=-5s^4-4bs^3-2cs.             \tag{3.60}
\]

Put `y=p_0-sx`.  After (3.60), `h_30=h_11=0` has an irreducible elimination
resultant.  On the dense chart `yU!=0` it is the graph

\[
 b={5V\over12U},\qquad
 c=-{N\over54(sx-p_0)},                         \tag{3.61}
\]

where

\[
\begin{aligned}
U={}&36s^2x^2y+3sx^4+36sxy^2+x^3y+9y^3,\\
V={}&-216s^3x^2y-18s^2x^4-216s^2xy^2+6sx^3y
     -54sy^3+x^5+9x^2y^2,\\
N={}&-108bp_0s^2-36bp_0x+108bs^3x
     -45p_0^2-45s^2x^2+5x^3.
\end{aligned}                                                   \tag{3.62}
\]

The charts omitted by (3.61) are not extra components.  Direct elimination
gives the following useful audit.  Write

\[
 J=20s^3+6bs^2-c,                              \tag{3.63}
\]

so `J!=0` is the ordinary-`A_2` determinant.

\[
\begin{array}{c|c|c}
\text{chart}&h_{20}&J\\ \hline
x=0&-5p_0^3/54&5p_0/6\\
y=0&5x^4(18s^2+x)/(8748s)&5x(18s^2+x)/(81s)\\
U=V=0&x^4(27bp_0-10x^2)/(2187p_0)&-c
\end{array}                                                     \tag{3.64}
\]

In the last row necessarily `s=0`, `9p_0^2+x^3=0`, and

\[
 c={x(-18bp_0+5x^2)\over27p_0}.                \tag{3.65}
\]

Thus every point of the `h_20!=0` incidence has an ordinary `A_2` cusp
unless its local type is `A_4`, `E_6`, or `E_8`.  The finite local
transposition check recorded after (3.57) removes `A_4` and `E_8`.  Imposing
the remaining `E_6` condition `s=c=0`, and normalizing `b=1`, gives exactly

\[
 (x,p_0)=\left(-{36\over25},{72\over125}\right),\qquad
 (x,p_0)=\left(-{32\over25},{256\over375}\right).              \tag{3.66}
\]

For both rows `C_1` has coefficients `(a,b,c,d)=(0,1,0,0)`.  The respective
second rows `(A,E,C,D,p_0,q_0)` are

\[
\begin{aligned}
&\left(-{36\over25},-{12\over5},-{24\over25},{48\over25},
        {72\over125},-{288\over625}\right),\\
&\left(-{32\over25},-{32\over15},-{128\over225},{2048\over1125},
        {256\over375},-{2048\over3375}\right).                 \tag{3.67}
\end{aligned}
\]

One rational ordinary-cusp witness on the dense surface is

\[
 (a,b,c,d,x,p_0)=
 \left(-3,{35\over12},-{325\over54},-{125\over27},1,-1\right).
                                                               \tag{3.68}
\]

At (3.68), the second coefficients are

\[
 (A,E,C,D,q_0)=
 \left(-2,{5\over3},-{205\over54},-{170\over27},-{365\over162}\right).
                                                               \tag{3.69}
\]

The certified verifier
[`verify_f2_r6_cusp_braid.sage`](../scripts/verify_f2_r6_cusp_braid.sage)
computes respectively `17,14,14` Zariski--van Kampen braids for (3.68) and
the two rows (3.67).  In every case the finite transposition enumeration
leaves seven assignments and no transitive assignment.

It remains to justify why one ordinary witness covers the special
multibranch strata of the surface.  This is the **tame-split regeneration
lemma**.  By the total excess-one argument (3.52), every multibranch point
has zero excess.  Its branch transpositions are therefore equal or disjoint,
hence pairwise commuting.  A pure braid fixes a tuple of commuting
transpositions: this follows first for the standard pure generator
`sigma_i^2` from the Hurwitz rule, and then for all of the pure braid group
by conjugation.  Regenerating a labeled multibranch singularity inside the
irreducible incidence surface adds exactly such pure local braid relations.
Consequently an admissible transposition representation on any special
ordinary-cusp member would extend to the dense member (3.68).  The certified
nonexistence there excludes it.  The only non-ordinary unibranch members
which pass the local test are the two `E_6` rows, already checked directly.
Component exchange gives no new case.  Therefore

\[
 \boxed{\text{the entire normal terminal-slice row }r=6
        \text{ is impossible}.}                              \tag{3.70}
\]

## 4. Why the terminal `R/(w^3)` packet cannot choose a case

The exact terminal calculation gives a cyclic logarithmic cokernel
`R/(w^3)` at the endpoint over `lambda_0`.  All three cases in (3.8) have
that same resolved endpoint.

For `a=0`, the terminal and ramification curves are already transverse.  In
logarithmic bases the map is diagonal with last entry `3w^3`.

For `ord(a)=1`, the ramification curve is tangent to the terminal curve.
After the two point blowups that separate them, a chart at the
terminal--exceptional node is

\[
 \pi=vw^2,qquad z=(1+v)w^3.                                  \tag{4.1}

The logarithmic determinant is

\[
 w^3(3+v)\cdot\text{unit}.                                   \tag{4.2}

For `r=ord(a)>=2`, one blowup separates the terminal tangent direction from
the ramification direction.  At the terminal--exceptional node,

\[
 \pi=vw,qquad
 z=w^3+v^rw^{r+1},                                            \tag{4.3}

and the logarithmic determinant is

\[
 w^3\left(3+v^rw^{r-2}\right).                               \tag{4.4}

Each parenthesis is a unit.  Since the first logarithmic row is unimodular,
all cases reduce to

\[
 \boxed{\operatorname{coker}(d\bar F^{\log})\simeq R/(w^3).} \tag{4.5}

This explains structurally why the terminal determinant and its generic
Smith data could not reveal the missing affine branch.

On the current common fan, the immediate component on the carrier side of
the terminal divisor is a `(-1)` component of valency two.  A smooth cubic
cusp attachment would meet it at a second point and make it trivalent.  This
is exactly the sole local shape not removed by the existing contracted
divisor gate.  This incidence match is a necessary compatibility, not a
construction of the global cover.

## 5. Nonnormal terminal slices are a conductor problem

The finite local cubic algebra is flat: a normal two-dimensional local
domain finite over the regular target is Cohen--Macaulay, hence free of rank
three.  Its closed length-three local fiber has one of the two algebra
types

\[
 k[\epsilon]/(\epsilon^3),
 \qquad
 k[u,v]/(u,v)^2.                                               \tag{5.1}

That dichotomy does **not** decide normality of the terminal slice.  The first
two rows below both have curvilinear closed fiber, although only the first is
normal.

The restriction to the target boundary makes the normalization parameter
visible.  Put `z=w^3`.  Three basic rank-three orders illustrate the missing
information:

\[
\begin{array}{c|c|c|c}
\text{boundary order}&\text{basis over }k[[z]]&
 \delta=\ell(\widetilde A/A)&\operatorname{ord}_z\operatorname{disc}\\ \hline
k[[w]]&(1,w,w^2)&0&2\\
k[[w^2,w^3]]&(1,w^2,w^4)&1&4\\
k[[w^3,w^4,w^5]]&(1,w^4,w^5)&2&6.
\end{array}                                                     \tag{5.2}

The last closed fiber is square-zero.  The discriminant order rises by
`2delta`, while the normalized residue map remains `z=w^3` in every row.
There are higher-conductor orders as well; (5.2) is diagnostic, not an
exhaustive classification of singular surface sections.

### 5.1 Exact conductor--contact identity

The same discriminant calculation gives an exact classification without
classifying the order itself.  Let `A_0` be the terminal Cartier slice of
the local cubic algebra and let

\[
 A_0\subset \widetilde A_0=k[[w]],\qquad z=w^3,
 \qquad \delta_T=\ell_k(\widetilde A_0/A_0).                  \tag{5.3}
\]

Both rings are free rank-three lattices over the DVR `k[[z]]`.  The
normalized lattice has discriminant `-27z^2`, while passage to a sublattice
of colength `delta_T` multiplies its discriminant by the square of the
index.  Therefore

\[
 \operatorname{ord}_z\operatorname{disc}(A_0/k[[z]])
 =2+2\delta_T.                                                \tag{5.4}
\]

The boundary itself is not a branch component: its generic normal index is
one.  Hence intersecting the discriminant divisor of the finite cubic map
with the target boundary counts precisely the affine branch germs through
the terminal value.  If `C_j` has boundary contact `k_j` and generic tame
inertia index `e_j` inside the cubic packet, its discriminant coefficient is
`e_j-1`.  Consequently

\[
 \boxed{\sum_j(e_j-1)k_j=2+2\delta_T.}                        \tag{5.5}
\]

Combining (5.5) with the four inertia multisets in (2.2) gives the complete
degree-six conductor--contact atlas:

\[
\begin{array}{c|c}
\text{generic inertia}&\text{exact contact equation}\\ \hline
(3)&k=1+\delta_T\\
(2)&k=2+2\delta_T\\
(2,2)&k_1+k_2=2+2\delta_T\\
(3,2)&2k_3+k_2=2+2\delta_T.
\end{array}                                                    \tag{5.6}
\]

Here every `k_j` lies in `1,...,24` by affine purity.  In particular the
first family has `0<=delta_T<=23`, the second has
`0<=delta_T<=11`, and the transposition contact `k_2` in every mixed row is
even.  The mixed row also spends the entire Orevkov budget
`3+2=d-1`, so every residue-defect term in that identity vanishes.

At `delta_T=0`, (5.6) gives exactly the three normal-slice patterns (3.8):
`(3;k=1)`, `(2;k=2)`, and `(2,2;k_1=k_2=1)`; the mixed family is empty.
Thus positive conductor is not merely one possible explanation of an extra
contact.  Equation (5.5) measures *all* of the missing boundary contact
exactly.

Ambient smoothness does not remove this regime.  The exact cyclic cubic map

\[
 (x,y)\longmapsto(u=x^3-y^2,\ v=y)                            \tag{5.7}
\]

has smooth source and target.  The target boundary `u=0` pulls back to the
cusp `x^3=y^2`, whose normalization is `(x,y)=(t^2,t^3)`.  Hence its
boundary order is `k[[t^2,t^3]]`, its terminal residue map is still
`v=t^3`, and its normalization conductor has length one.  The ramification
curve `x=0` has generic index three and maps to

\[
 u=-v^2,                                                       \tag{5.8}
\]

a smooth branch of contact `k=2` with the target boundary.  Thus even a
smooth Stein point can carry the extra pattern `(A,k)=(3,2)` once its
terminal slice is nonnormal.  More generally `x^3=u+v^m` produces higher
contact orders when the boundary slice is irreducible.  These models are
local compatibility examples, not Keller maps.

This is exactly the normalization mismatch anticipated by the logarithmic
perfect-complex programme.  The terminal `R/(w^3)` packet sees the normalized
degree-three branch but not the conductor of the contracted Stein slice.
Any exclusion of the conductor regime must either:

1. recover and bound this conductor through the localized `ch_2`/`Fitt_1`
   filtration; or
2. classify the normal rank-three surface algebras compatible with the
   compiled carrier chain and exclude their resolution graphs.

## 6. Exact remaining degree-six ledger

The degree-six problem is now divided as follows.

| local cubic regime | target pattern | present disposition |
| --- | --- | --- |
| normal slice, `a=0` | one `k=1` cubic-inertia component | covered by the existing one-component degree-six cubic-row exclusion |
| normal slice, `ord(a)=6` | irreducible rational ordinary-cusp surface plus the two `E_6` rows (3.66) | excluded by certified braid monodromy and tame-split regeneration, (3.70) |
| normal slice, `ord(a)=8` | the two component-swapped copies of the exact `A_2+3A_1` / `4A_1` quartic (3.54)--(3.57) | excluded for all four embeddings by certified braid monodromy, (3.59) |
| normal slice, `ord(a)=10` | two disjoint four-node quintics with constant difference | excluded by (3.51) |
| normal slice, `ord(a)=9` | one `k=2` simple-inertia component; third semigroup in (3.13) | excluded by the sharp polynomial parametrization gap theorem [PF2D6O1](F2_DEGREE_6_10_POLYNOMIAL_GAP.md) |
| normal slice, `ord(a)=5,7` | one `k=2` simple-inertia component; first two semigroups in (3.13) | unresolved two-row cusp-at-infinity case |
| nonnormal terminal slice | one of the four finite conductor--contact families (5.6) | unresolved normalization/contracted-fiber case |

Thus a degree-six counterexample, if one exists, must lie in one of the last
two rows.  Conversely, excluding those two rows closes geometric degree six
completely.  The last row includes smooth ambient Stein points such as
(5.7); it must not be shortened to “singular Stein point.”

## Sources

- Rick Miranda,
  [*Triple Covers in Algebraic Geometry*](https://www.math.colostate.edu/~miranda/preprints/TripleCoversInAG.pdf),
  Sections 2--5, for finite flat cubic algebras, the discriminant, and the
  smooth total-ramification cusp model.
- S. Yu. Orevkov,
  [*On three-sheeted polynomial mappings of C2*](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
  Lemma 4.2, for the global local-degree budget.
- Abdallah Assi and Pedro A. Garcia-Sanchez,
  [*On curves with one place at infinity*](https://arxiv.org/abs/1407.0490),
  Proposition 2 and Section 4, for the conductor formula and the
  delta-sequence polynomiality criterion used in (3.13).
- S. Yu. Orevkov,
  [*On germs of mappings C2 to C2*](https://www.math.univ-toulouse.fr/~orevkov/k-en.pdf),
  Theorem 2 and Section 6, for the uniqueness and resolved intersection
  charge in (3.18)--(3.19).
- Miguel Marco-Buzunariz,
  [Sage's Zariski--van Kampen implementation](https://doc.sagemath.org/html/en/reference/curves/sage/schemes/curves/zariski_vankampen.html),
  especially `braid_monodromy` and `braid_from_piecewise`, for the certified
  braid computation and the crossing convention ported by the exploratory
  script.
- Miguel Marco-Buzunariz and Marcos Rodriguez,
  [*SIROCCO: a library for certified polynomial root continuation*](https://zaguan.unizar.es/record/131386),
  for the disjoint certified root tubes used by the two Sage verifiers.

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_geometric_degree_six_stein_reduction.py
.venv/bin/python scripts/verify_f2_degree_6_10_gap.py
.venv/bin/python scripts/explore_f2_r8_cusp_braid.py
.venv/bin/python scripts/explore_f2_r8_cusp_braid.py --real-root-index 0
# Requires SageMath with its optional SIROCCO package:
sage -python scripts/verify_f2_r6_cusp_braid.sage
sage -python scripts/verify_f2_r8_cusp_braid.sage
```
