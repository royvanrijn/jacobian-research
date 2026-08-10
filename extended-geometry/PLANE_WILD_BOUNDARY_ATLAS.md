# A wild-boundary atlas for finite plane covers

> **Status.** Sections 2--5 give, for every characteristic `p` and every
> `N` divisible by `p`, finite flat generically separable plane covers and an
> exact normalization theorem for the one-retained-sheet subfamily.  Writing
> `N=p^n*m` separates its boundary residue degree as `(f_sep,f_insep)=(m,p^n)`;
> the characteristic-two hidden cubic is the first row.  Section 7 proves a
> stronger reconstruction gate: for every `N>2`, deleting only the proposed
> fierce boundary leaves a tame ramification
> branch of index `N-1`.  Exact Fitting-ideal certificates are supplied for
> `N=3,5,7` and the controls `N=2,4,6`, including mixed residue rows in
> characteristics two and three.  Thus these rows do **not** construct
> new polynomial Keller maps, even if their remaining open surface happened to
> be an affine plane.  Section 7.3 then changes the boundary presentation: the
> balanced gluing `P^(N-1)QT` removes the companion different uniformly, but
> its natural affine-plane chart has Jacobian `-u^(2N-4)`.  The divisor
> localization sequence then computes the normalized complement's class group
> as `Z/(N-1)`, excluding every alternative affine-plane reconstruction for
> `N>2`.  Together these endpoints exclude the complete monomial-gluing band
> `P^aQT` for every `a>=0`: `a=0` has a free unit and `Z/(N+1)` class
> torsion, `1<=a<N-1` retains a tame different, `a=N-1` retains
> boundary-class torsion, and `a>=N` has a wild index-`N` branch over `P=0`.
> More generally, a coefficient `C(P)` with a factor away from `P=0` creates
> an additional different divisor, so the monomials exhaust all one-variable
> gluings having no extra target support.  For arbitrary `C(P,Q)`, the same
> support theorem reduces the search to `cP^aQ^b`.  The thickened balanced
> rows are finite base changes of the `b=0` surface.  Their prime-to-`p` part
> is excluded by compactly supported Euler characteristic, their pure-`p`
> part for `N>2` preserves the exact boundary-class order by push--pull, and
> the remaining `p=N=2` tower has full core class group `(Z/2)^2` at every
> Frobenius height.  Together with the unchanged `P=0` Newton rows, this
> closes the entire two-parameter monomial quadrant: the characteristic-two
> cubic `N=2,a=1,b=0` is the unique affine-plane Keller row under this
> one-omitted-boundary and target-support hypothesis.  The
> retained-polynomial extension is closed as well: its support identity
> first forces `A=a0+T*B(T^p)`, and the normalized complement then has
> `chi_c=deg(A)` after the roots split, excluding every nonlinear `A`.
> The linear case is already excluded by `Cl=Z/(N-1)` for `N>2`.  Thus
> characteristic two is exceptional throughout this balanced,
> squarefree-retained architecture, not only in the monomial subfamily.
> Section 8 records the determinant ledger, and Section 9 gives a characteristic-zero
> local-cohomology receptacle, not a new proof of `JC_2`.

The starting row is the hidden cubic in the
[characteristic-two plane normalization](../verified/HUQ_KURUVILLA_PLANE_BOUNDARY_NORMALIZATION.md):

\[
 T^3+T^2+(PQ+P^3)T+P^3
 =(T-1)(T^2-P^3)+PQT
 \quad(\operatorname{char}k=2).                 \tag{1.1}
\]

The right side, rather than the displayed source formula, is the part which
generalizes.  It separates four ingredients:

1. a retained factor `T-1`;
2. a fierce factor `T^2-P^3`;
3. the gluing term `PQT`; and
4. one transverse different supplied by ordinary differentiation in `T`.

This viewpoint is compatible with the direct modulo-four result.  The
[de Rham--Cartier obstruction](../verified/HUQ_KURUVILLA_PLANE_W2_OBSTRUCTION.md)
now rules out not only the literal characteristic-two source formula but
every polynomially left--right equivalent plane representative.  A
mixed-characteristic continuation must therefore change the finite boundary
presentation, not merely change coordinates on that formula.

## 1. Ledger conventions

Let `D` be a target prime and `E` a prime of a normalized finite cover above
it.  The row recorded here is

\[
 (e_E,f_{E,\mathrm{sep}},f_{E,\mathrm{insep}},d_E,
   \mathfrak c_E,s_E).                            \tag{1.2}
\]

Here `d_E` is the Kähler-different exponent and

\[
 s_E=e_Ef_{E,\mathrm{sep}}f_{E,\mathrm{insep}}   \tag{1.3}
\]

is the number of generic sheets lost if `E` is deleted.  The conductor entry
must name its object.  Three different notions occur and are not
interchangeable:

- the conductor of a primitive surface order in its normalization;
- the conductor of a reduced boundary curve in its own normalization; and
- the Artin or Swan conductor of a local Galois character.

For a node `k[[a,b]]/(ab)`, the boundary-normalization conductor is `(a,b)`.
For a monomial curve `k[[z^u,z^v,...]]`, it is `z^c k[[z]]`, where `c` is the
numerical-semigroup conductor.

## 2. The characteristic-divisible hidden-order family

Let `k` have characteristic `p>0`, let

\[
 N=p^n m\qquad(n\geq1,\ (m,p)=1),                 \tag{2.1}
\]

and define

\[
 \boxed{H_N=(T-1)(T^N-P^{N+1})+PQT.}              \tag{2.2}
\]

Put

\[
 B_N=k[P,Q,T]/(H_N).                               \tag{2.3}
\]

### Proposition 2.1 -- finite separable degree prime to `p`

The map `Spec(B_N) -> A^2_(P,Q)` is finite flat of degree `N+1`.  Its
function field has exact separable degree `N+1`, which is prime to `p`.

Indeed, `H_N` is monic of degree `N+1` in `T`.  Viewed as a polynomial in
`Q`, it is primitive and linear, with coefficient `PT`.  The other
coefficient `(T-1)(T^N-P^(N+1))` has no common factor `P` or `T`, so `H_N` is
irreducible.  Finally,

\[
 (H_N)_T=T^N-P^{N+1}+PQ\ne0.                      \tag{2.4}
\]

Thus this is already a uniform all-characteristic family of finite plane
covers of prime-to-characteristic degree.  It is not yet a uniform family of
Keller endomorphisms of the affine plane.

### Proposition 2.2 -- the fierce row

At `Q=0`,

\[
 H_N=(T-1)(T^N-P^{N+1}).                           \tag{2.5}
\]

The retained component is

\[
 A:\quad Q=0,\quad T=1.                            \tag{2.6}
\]

The other component is

\[
 E_0:\quad Q=0,\quad T^N=P^{N+1}.                 \tag{2.7}
\]

At its generic point, `P` and `T` are units and `Q` is a uniformizer.  The
residue extension has normalization parameter `z` with

\[
 P=z^N,\qquad T=z^{N+1}.                           \tag{2.8}
\]

The extension `k(P) subset k(z)`, `P=z^N`, has separable degree `m` and
inseparable degree `p^n`.  Hence the generic row is

\[
 \boxed{(e,f_{\rm sep},f_{\rm insep},d,s)
       =(1,m,p^n,1,N).}                            \tag{2.9}
\]

For the different, write `a=T-1` and `b=T^N-P^(N+1)`.  The two identities

\[
 H_N=ab+PQT,
 \qquad
 a(H_N)_T-H_N=-PQ                                  \tag{2.10}
\]

give

\[
 (H_N)_T=-\frac{PQ}{T-1}\quad\text{at the generic point of }E_0. \tag{2.11}
\]

Thus the monogenic different has valuation one.  The retained component is
generically ordinary: `(H_N)_T=1-P^(N+1)` on `A`.

### Proposition 2.3 -- `N+1` conductor nodes

The two reduced boundary components meet where

\[
 T=1,\qquad P^{N+1}=1.                             \tag{2.12}
\]

There are `N+1` reduced points because `p` does not divide `N+1`; its
derivative is `P^N` in `k[P]`.  At any such point,
`a=T-1` and `b=T^N-P^(N+1)` are regular parameters and

\[
 Q=-\frac{ab}{PT}.                                 \tag{2.13}
\]

The reduced boundary is therefore the node `k[[a,b]]/(ab)`, with conductor
`(a,b)`.  Since `dP` is a unit multiple of `db`,

\[
 dP\wedge dQ=(\text{unit})\,b\,db\wedge da.        \tag{2.14}
\]

The different remains `(b)` through every node.  Here `E_0=(b)` and
`A=(a)`.

## 3. Exact surface normalization

The primitive order is not normal along `(P,T)`.  In its fraction field set

\[
 \boxed{W=\frac{P^N}{T}.}                          \tag{3.1}
\]

It is integral: besides `WT=P^N`, it satisfies

\[
 W^2-WP^N+WQ+P^{N-1}T^{N-2}(T-1)=0.               \tag{3.2}
\]

Let `C_N=B_N[W]`.

### Theorem 3.1 -- uniform normalization and primitive conductor

For every `N` divisible by `p`,

\[
 \boxed{C_N=\operatorname{Norm}(B_N),\qquad
        \mathfrak c_{B_N\subset C_N}=(P,T).}       \tag{3.3}
\]

A presentation of `C_N` is given by `H_N` and

\[
\begin{aligned}
 WT-P^N&=0,\\
 (T-1)(T^{N-1}-PW)+PQ&=0,\\
 W^2-WP^N+WQ+P^{N-1}T^{N-2}(T-1)&=0.              \tag{3.4}
\end{aligned}
\]

For completeness, the normality argument is short.  As a `B_N`-module,

\[
 C_N=B_N+B_NW=\frac{(T,P^N)}{T}.                  \tag{3.5}
\]

The quotient by `(T,P^N)` is one-dimensional Cohen--Macaulay, so (3.5) is a
maximal Cohen--Macaulay module over the hypersurface `B_N`.  Hence `C_N`
satisfies `S_2`.  The singular locus of `B_N` is the line `(P,T)`.  Above its
generic point, (3.2) reduces to `W(W+Q)=0`, giving two height-one primes.
On `W=0`, `Q` is a unit and (3.2) plus the second relation in (3.4) gives a
regular parameter; on `W=-Q`, `W` is a unit and `WT=P^N` gives a regular
parameter.  Thus `C_N` is regular in codimension one and is normal by
Serre's criterion.

Moreover,

\[
 TW=P^N,
 \qquad
 PW=P^{N+1}+T^{N-1}-T^N-PQ                       \tag{3.6}
\]

belong to `B_N`, so `(P,T)` is contained in the conductor.  The class of `W`
along the singular line shows that no element outside `(P,T)` annihilates
`C_N/B_N`; hence equality holds.

The two upstairs conductor branches are

\[
 (P,T,W),\qquad(P,T,W+Q).                          \tag{3.7}
\]

This exactly recovers the two characteristic-two branches.

## 4. The normalized missing boundary and the exceptional row

The closure `E` of the fierce component in `Spec(C_N)` has parameterization

\[
 P=z^N,\qquad T=z^{N+1},\qquad
 W=z^{N^2-N-1},\qquad Q=0.                         \tag{4.1}
\]

Therefore

\[
 \mathcal O(E)=k[z^N,z^{N+1},z^{N^2-N-1}].        \tag{4.2}
\]

The third exponent is the Frobenius number of the semigroup generated by
`N,N+1`.  Adjoining it fills exactly the largest gap.  Consequently the
boundary-normalization invariants are

\[
\begin{array}{c|c|c}
 &\text{conductor exponent}&\delta\\ \hline
 N=2&0&0\\
 N>2&N(N-2)&\dfrac{N(N-1)}2-1.
\end{array}                                        \tag{4.3}
\]

### Corollary 4.1 -- characteristic two is exceptional under a smooth-edge
hypothesis

In the characteristic-divisible hidden-order family, the normalized fierce
boundary `E`
is a smooth affine line if and only if `N=2`.  Thus, under the natural
hypothesis

> the cover has one retained sheet, one fierce prime, different exponent
> one, and the omitted fierce prime is a smooth boundary line after surface
> normalization,

the characteristic-two row is the unique nontrivial row of this atlas.

This is an actual exception theorem for the boundary presentation, not a
theorem that odd-characteristic Keller counterexamples do not exist.  For
`N>2`, the remaining boundary defect has positive length (4.3); it is the
first extra module which any odd-characteristic reconstruction must absorb.

## 5. Prescribing the prime-to-characteristic degree

The one-retained-sheet family can be enlarged without changing the fierce
row.  Let `A(T)` be monic, squarefree, of degree `r`, with `A(0) != 0`, and
put

\[
 \boxed{H_{N,A}=A(T)(T^N-P^{N+1})+PQT.}            \tag{5.1}
\]

The same linear-in-`Q` argument proves irreducibility, and the cover is
finite flat of exact separable degree

\[
 d=N+r.                                            \tag{5.2}
\]

The `Q=0` fiber has the fierce component of sheet contribution `N` and `r`
ordinary retained components.  At the fierce generic point,

\[
 (H_{N,A})_T
 =PQ\left(1-\frac{TA'(T)}{A(T)}\right),            \tag{5.3}
\]

modulo the equation.  The parenthesis is not identically zero because its
constant coefficient is `A(0)`.  Hence the different exponent is still one.
Each root of `A` contributes `N+1` reduced nodes, for a total of
`r(N+1)`.

Given any prescribed `d>N` with `p` not dividing `d`, take `r=d-N`; then
`p` does not divide `r` and (5.1) realizes a finite separable plane cover of
degree `d` with one fierce `N`-sheet.  What is not supplied is a distinguished
open isomorphic to `A^2` after deleting that sheet.

## 6. Artin--Schreier and Artin--Schreier--Witt comparison rows

The simplest fierce additive block is

\[
 Z^N+Q^mZ-P=0,\qquad N=p^n.                        \tag{6.1}
\]

It is irreducible and generically separable, with derivative `Q^m`.  At
`Q=0` its row is

\[
 \boxed{(1,1,N,m,\text{no algebra conductor},N).}  \tag{6.2}
\]

After adjoining an `(N-1)`st root of `Q^m`, (6.1) becomes a generalized
Artin--Schreier equation.  This is the cleanest way to compare fierce rows
having the same `(e,f)` but different differents.

For the usual Artin--Schreier extension

\[
 Y^p-Y=cQ^{-m},\qquad p\nmid m,                    \tag{6.3}
\]

the totally ramified row is

\[
 \boxed{(p,1,1,(p-1)(m+1),m+1,p),}                \tag{6.4}
\]

where `m+1` is the Artin conductor, not an algebra-normalization conductor.

For a reduced Artin--Schreier--Witt vector of length `n`, let
`u_1<...<u_n` be the upper jumps of the cyclic degree-`p^n` extension.  Its
row is

\[
\boxed{
 \left(p^n,1,1,
 \sum_{j=1}^n(p^j-p^{j-1})(u_j+1),
 (u_1+1,\ldots,u_n+1),p^n\right).}                \tag{6.5}
\]

The formula follows by summing the conductor of the `p^j-p^(j-1)` characters
of exact order `p^j`.  The ramification-break input is standard
Artin--Schreier--Witt theory; useful primary references are
[Thomas (2005)](https://doi.org/10.5802/jtnb.514) and
[Elder--Keating (2025)](https://arxiv.org/abs/2503.16830).
For two-dimensional Artin--Schreier ramification along normal-crossing
boundaries, see [Zhukov (2002)](https://arxiv.org/abs/math/0209183).

These comparison rows show why `(e,f)` is insufficient.  A tame Kummer row
of index `r` prime to `p` has different `r-1`; (6.3) has wild index `p` and
different `(p-1)(m+1)`; (6.2) has index one, inseparable residue degree `N`,
and freely varying different `m`.

## 7. The one-boundary reconstruction gate in characteristics `3,5,7`

The fierce ledger alone does not detect all ramification.  Let

\[
 \mathfrak D_N=\operatorname{Fitt}_0
 \Omega_{C_N/k[P,Q]}                              \tag{7.1}
\]

be the relative Kähler different.  In the presentation (3.4), it is generated
by the two-by-two minors of the four-by-two Jacobian matrix with columns
`dW,dT`.  Write

\[
 I_E=(Q,T^N-P^{N+1}),\qquad I_L=(P,T,W)            \tag{7.2}
\]

inside `C_N`.  Saturating `\mathfrak D_N` by `I_E` removes the proposed fierce
boundary from its support.

### Proposition 7.1 -- the companion tame branch

For every `N` divisible by `p`, the reduced inverse image of `P=0` has exactly
three components:

\[
 A_0=(P,T-1,W),\qquad
 L_0=(P,T,W),\qquad
 L_1=(P,T,W+Q).                                   \tag{7.3}
\]

Their generic sheet contributions are respectively

\[
 1,\qquad N-1,\qquad1.                            \tag{7.4}
\]

In particular, they exhaust the degree `N+1`.  The middle component `L_0` is
the first upstairs conductor branch from (3.7).  At its generic point `Q` is
a unit.  Its localized maximal ideal is generated by `P,T,W`; the first
relation in (3.4) gives

\[
 P=-\frac{(T-1)T^{N-1}}{Q-(T-1)W}
   =({\rm unit})T^{N-1},
 \qquad
 W=\frac{P^N}{T}=({\rm unit})T^{N^2-N-1}.         \tag{7.5}
\]

Thus `T` has the least positive valuation among the displayed maximal-ideal
generators and is a uniformizer.
Since `p` does not divide `N-1`, this is tame ramification.  Its full generic
row, with the conductor object named, is

\[
 \boxed{(e,f_{\rm sep},f_{\rm insep},d,
          \mathfrak c,\text{sheet loss})
       =(N-1,1,1,N-2,(P,T,W),N-1).}               \tag{7.6}
\]

Here `(P,T,W)` means the corresponding branch of the primitive-order
conductor, not an Artin conductor.  Formula (7.5) proves the ramification and
different assertions for every `N`; it does not depend on a bounded search.
The decomposition (7.3) follows directly by setting `P=0` in (3.4): either
`T=1,W=0`, or `T=0` and `W(W+Q)=0`.

### Exact Fitting certificates

For the requested odd-characteristic rows, exact Singular elimination gives

\[
\begin{array}{c|c|c|c}
p=N&\mathfrak D_N:I_E^\infty&\sqrt{\mathfrak D_N}
 &\operatorname{length}_{k(Q)[[T]]}
   k(Q)[[T]]/(\mathfrak D_N:I_E^\infty)\\ \hline
3&(P,W,T)&I_E\cap I_L&1\\
5&(P,W,T^3)&I_E\cap I_L&3\\
7&(P,W,T^5)&I_E\cap I_L&5.
\end{array}                                       \tag{7.7}
\]

The characteristic-two control `N=2` has
`\mathfrak D_2:I_E^\infty=C_2` and
`\sqrt{\mathfrak D_2}=I_E`; the `N=4` control gives `(P,W,T^2)` and support
`I_E\cap I_L`.  The mixed-residue controls `(p,N)=(2,6),(3,6)` both give
`(P,W,T^4)` with the same support `I_E\cap I_L`.  Thus (7.7) and these
controls are exact global certificates, while Proposition 7.1 supplies the
uniform local theorem.

The residual determinant module is

\[
 \mathcal R_N=C_N/(\mathfrak D_N:I_E^\infty).     \tag{7.8}
\]

At the generic point of `L_0` it is
`k(Q)[[T]]/(T^{N-2})`, of lengths `1,3,5` in characteristics `3,5,7`.
This support is a whole divisorial component over `P=0`, not the finite set of
conductor nodes.  It therefore cannot be removed by a local-cohomology class
supported only at those nodes.

### Corollary 7.2 -- characteristic two is exceptional for this one-boundary
presentation

Deleting only `E` leaves `L_0`.  For every `N>2`, the restricted map is still
ramified there with different exponent `N-2`; no change of source coordinates
can turn its Jacobian into a unit.  Hence the canonical hidden-order family
cannot yield a Keller map by the stipulated one-boundary reconstruction in
characteristics `3,5,7`, or indeed for any `N>2`.  The row `N=2` is the unique
one for which the companion different vanishes.

This is the third proposed win under a precise natural hypothesis: the finite
cover is (2.2), and the reconstruction omits exactly its fierce boundary.  It
is not a theorem that characteristic two is exceptional among all possible
plane boundary presentations.  Removing both `E` and `L_0` would be a
different, two-boundary problem and would lose sheets over two distinct target
divisors.

### 7.3 Balanced gluing removes the companion different

The obstruction in Proposition 7.1 is caused by the valuation-one gluing
coefficient `PQT`.  Keep the same `Q=0` boundary but replace the primitive
order by

\[
 \boxed{\widehat H_N=(T-1)(T^N-P^{N+1})+P^{N-1}QT,}
 \qquad
 \widehat B_N=k[P,Q,T]/(\widehat H_N).             \tag{7.9}
\]

The linear-in-`Q` argument from Proposition 2.1 again gives an irreducible,
finite flat, generically separable cover of exact degree `N+1`.  Its `Q=0`
fiber, fierce row, `N+1` nodes, and sheet loss `N` are unchanged.

Let `\widehat C_N` be its normalization and let `\widehat E` be the closure
of the fierce component.  In characteristic `p`, direct differentiation gives

\[
 (\widehat H_N)_T=T^N-P^{N+1}+P^{N-1}Q,
 \qquad
 (T-1)(\widehat H_N)_T-\widehat H_N=-P^{N-1}Q.    \tag{7.10}
\]

Consequently the primitive order is already étale wherever `PQ` is a unit.
At `P=0`, the Newton polygon of the `T=0` cluster has vertices

\[
 (0,N+1),\qquad(1,N-1),\qquad(N,0),               \tag{7.11}
\]

and slopes `-2,-1`.  After normalization its two rows are

\[
 (1,1,1,0,\text{none},1),qquad
 (1,N-1,1,0,\text{none},N-1).                    \tag{7.12}
\]

The second residue extension is separable because `p` does not divide `N-1`.
Together with the ordinary `T=1` sheet, these contributions sum to `N+1`.
Thus no height-one different remains over `P=0`.  Over `Q=0`, the retained
sheet is ordinary and `\widehat E` has the same fierce row `(1,1,N,1,N)` as
before when `m=1`; in general its row is `(1,m,p^n,1,N)` from (2.9).
Purity of the branch locus therefore gives

\[
 \boxed{\widehat C_N\setminus\widehat E
       \longrightarrow\mathbb A^2_{P,Q}
       \text{ is étale}.}                         \tag{7.13}
\]

Here purity is applied to the finite normalization over the regular target;
see [Stacks Project, Tag 0BMB](https://stacks.math.columbia.edu/tag/0BMB).
Equation (7.13) is a uniform finite-cover determinant theorem, not an
affine-plane reconstruction.

The exact normalization diagnostics for the requested rows give

\[
\begin{array}{c|c|c|c|c}
p=N&\mathfrak c_{\widehat B_N\subset\widehat C_N}
 &\#\widehat C_N(\mathbb F_p)&\#\widehat E(\mathbb F_p)
 &\#(\widehat C_N\setminus\widehat E)(\mathbb F_p)\\ \hline
3&(P,T)^2&12&3&9\\
5&(P,T)^4&30&5&25\\
7&(P,T)^6&56&7&49.
\end{array}                                       \tag{7.14}
\]

For `p=3`, the complete relative Fitting computation also gives

\[
 \mathfrak D_{\widehat C_3}:I_{\widehat E}^{\infty}=(1),
 \qquad
 \sqrt{\mathfrak D_{\widehat C_3}}=I_{\widehat E}. \tag{7.15}
\]

The point counts in (7.14) match `\mathbb A^2`; they are diagnostics, not an
isomorphism proof.  For `p=5,7`, the conductor and point counts are exact
normalization computations, while (7.13) is the all-`N` valuation-and-purity
proof and does not depend on computing a large global Fitting ideal.
The mixed controls `(p,N)=(2,6),(3,6)` likewise have conductor `(P,T)^5`
and counts `(6,2,4)` and `(12,3,9)` respectively for
`(#C,#E,#(C-E))` over the prime field.

There is a natural birational affine-plane chart.  Put

\[
\begin{aligned}
 r&=1+x^{N-1}y,\\
 u&=1+x^{N+1}r,\\
 P&=xru^{N-1},\qquad T=ru^N,\\
 Q&=-\frac{ru^N-1}{x^{N-1}}u^{N-2}.
\end{aligned}                                      \tag{7.16}
\]

The quotient in (7.16) is polynomial: `r-1` is divisible by `x^(N-1)`
and `u^N-1` is divisible by `u-1`, hence by `x^(N+1)`.  When `m=1`, the
Frobenius identity compresses it to the earlier two-term formula

\[
 \frac{ru^N-1}{x^{N-1}}
 =y+x^{N^2+1}r^{N+1}.                             \tag{7.16a}
\]

Substitution gives `\widehat H_N=0` for every `p|N`.  Moreover, if
`F=T^N-P^(N+1)`, then

\[
 x=\frac{PT^{N-1}}F,qquad
 u=\frac{T^N}F,                                   \tag{7.17}
\]

and `r,y` follow rationally, so this chart has the correct function field and
generic degree `N+1`.  Its Jacobian is nevertheless

\[
 \boxed{\det\frac{\partial(P,Q)}{\partial(x,y)}
       =-u^{2N-4}.}                                \tag{7.18}
\]

This identity is uniform, not an extrapolation from the bounded rows.  With
`x,r` as coordinates and using `N=0` in `k`, differentiation gives

\[
 dP=u^{N-2}(r\,dx+x\,dr),\qquad
 dT=u^N\,dr,
\]

while `dP wedge du=0`.  For
`S=x^{-(N-1)}(T-1)` and `Q=-S*u^(N-2)`, one obtains

\[
 dP\wedge dQ=-u^{2N-4}x^{-(N-1)}dx\wedge dr.
\]

Finally `dx wedge dr=x^(N-1) dx wedge dy`, proving (7.18) for every
`p|N`.

Thus (7.16) is Keller exactly for `N=2`, when it recovers the known
characteristic-two formulas.  For every `N>2` divisible by `p`, including
the mixed rows `(p,N)=(2,6),(3,6)`, it exhibits the remaining reconstruction
defect explicitly rather than producing a counterexample.  The following
calculation excludes every other affine-plane chart of this balanced
complement.

### Proposition 7.4 -- exact class group of the balanced complement

Put

\[
 U_N=\widehat C_N\setminus\widehat E,
 \qquad F=T^N-P^{N+1}.                            \tag{7.19}
\]

Let `A_0` be the ordinary `T=1` component over `P=0`, and let `L_0,L_1`
be the slope `-1,-2` components from (7.11).  Their valuation rows give the
complete principal-divisor ledger

\[
\begin{aligned}
 \operatorname{div}(P)&=A_0+L_0+L_1,\\
 \operatorname{div}(T)&=L_0+2L_1,\\
 \operatorname{div}(F)&=\widehat E+NL_0+(N+1)L_1.
\end{aligned}                                      \tag{7.20}
\]

Using (7.17), this becomes

\[
\begin{aligned}
 \operatorname{div}_{\widehat C_N}(x)
   &=A_0+(N-2)L_1-\widehat E,\\
 \operatorname{div}_{\widehat C_N}(u)
   &=(N-1)L_1-\widehat E.
\end{aligned}                                      \tag{7.21}
\]

Hence `x,u` are regular on `U_N`, and

\[
 \operatorname{div}_{U_N}(x)=A_0+(N-2)L_1,
 \qquad
 \operatorname{div}_{U_N}(u)=(N-1)L_1.            \tag{7.22}
\]

Assume `N>2`.  Then `D_{U_N}(x)` removes exactly `A_0` and `L_1`, and `u`
is a unit there.  Conversely, on this chart

\[
 r=\frac{u-1}{x^{N+1}},qquad
 y=\frac{u-1-x^{N+1}}{x^{2N}},                    \tag{7.23}
\]

so the formulas (7.16) and the universal property of normalization give

\[
 \boxed{D_{U_N}(x)\simeq
        \operatorname{Spec}k[x^{\pm1},u^{\pm1}]
        =\mathbb G_m^2.}                          \tag{7.24}
\]

The divisor-class localization sequence is therefore

\[
 k[x^{\pm1},u^{\pm1}]^*/\Gamma(U_N,\mathcal O)^*
 \longrightarrow \mathbb ZA_0\oplus\mathbb ZL_1
 \longrightarrow\operatorname{Cl}(U_N)\longrightarrow0.       \tag{7.25}
\]

This is an instance of the general
[class-trivial-core Smith criterion](../plane-jc/BOUNDARY_LATTICE_PREFILTER.md#dual-torus-core-localization):
`U_N` is normal, (7.24) identifies the torus chart, and (7.22) lists the
complete codimension-one complement.  These hypotheses are what make the
kernel and cokernel calculation conclusive; the same matrix without them
would be only a numerical candidate.

Every Laurent-chart unit is `c*x^a*u^b`.  By (7.22), such a unit extends to a
unit of `U_N` only if

\[
 a=0,\qquad a(N-2)+b(N-1)=0,                     \tag{7.26}
\]

hence only if `a=b=0`.  Thus `\Gamma(U_N,\mathcal O)^*=k^*`, and the two
principal-divisor columns in (7.25) are

\[
 \begin{pmatrix}1\\N-2\end{pmatrix},
 \qquad
 \begin{pmatrix}0\\N-1\end{pmatrix}.             \tag{7.27}
\]

Their Smith normal form has invariant factors `1,N-1`.  Consequently

\[
 \boxed{\operatorname{Cl}(U_N)\simeq\mathbb Z/(N-1).}          \tag{7.28}
\]

The stronger individual-class formula from the class-trivial-core theorem gives

\[
 \boxed{\operatorname{ord}[L_1]=N-1,\qquad
        \operatorname{ord}[A_0+(N-2)L_1]=1.}       \tag{7.28a}
\]

Indeed, augmenting (7.27) by `(0,1)` changes its top determinantal divisor
from `N-1` to `1`, whereas `(1,N-2)` is already the first
principal-divisor column.  Thus the obstruction is carried by the named
boundary class, not just by an unspecified nonzero element of `Cl(U_N)`.

In particular,

\[
 \operatorname{Cl}(U_3)=\mathbb Z/2,\qquad
 \operatorname{Cl}(U_5)=\mathbb Z/4,\qquad
 \operatorname{Cl}(U_7)=\mathbb Z/6.             \tag{7.29}
\]

The additional composite control is `Cl(U_6)=Z/5`.
Since the affine plane has trivial class group, `U_N` is not an affine plane
for any `N>2`.  For `N=2`, both divisor columns become the standard basis and
the known reconstruction has trivial class group.  Thus balanced gluing
passes the entire different ledger but still singles out the
characteristic-two cubic through a global boundary-class obstruction.

### Corollary 7.5 -- the monomial-gluing dichotomy

The original and balanced orders are the endpoints of the natural band

\[
 H_{N,a}=(T-1)(T^N-P^{N+1})+P^aQT,
 \qquad 1\leq a\leq N-1.                         \tag{7.30}
\]

Put `g=gcd(a,N-1)`.  The `T=0` Newton polygon over `P=0` has vertices

\[
 (0,N+1),\qquad(1,a),\qquad(N,0).                \tag{7.31}
\]

Its long segment has reduced slope
`-(a/g)/((N-1)/g)`.  The corresponding generic row is

\[
 \boxed{\left(\frac{N-1}{g},g,1,
               \frac{N-1}{g}-1,
               \text{none},N-1\right).}          \tag{7.32}
\]

The ramification is tame because `p` does not divide `N-1`; the residual
degree-`g` binomial is separable.  Thus the different in (7.32) vanishes if
and only if `a=N-1`.  At that unique endpoint, Proposition 7.4 gives instead
the nonzero global obstruction `Cl(U_N)=Z/(N-1)`.

Consequently, for every `N>2`, no monomial gluing in the complete band
(7.30) can produce a Keller affine-plane reconstruction after omitting only
the fierce boundary: `a<N-1` fails the different gate, while `a=N-1` fails
the class-group gate.  When `N=2`, the band consists only of `a=1`, both
obstructions vanish, and the known cubic reconstruction is recovered.  This
is a characteristic-two exception theorem under a broader and explicit
one-boundary hypothesis, not merely a comparison of two isolated formulas.

### Theorem 7.6 -- all nonnegative powers and one-variable coefficients

The restriction `1<=a<=N-1` can be removed.  First let

\[
 H_{N,C}=(T-1)(T^N-P^{N+1})+C(P)QT,\qquad
 0\ne C(P)\in k[P].                              \tag{7.33}
\]

In characteristic `p|N`, put `F=T^N-P^(N+1)`.  The universal identities are

\[
 (H_{N,C})_T=F+C(P)Q,\qquad
 (T-1)(H_{N,C})_T-H_{N,C}=-C(P)Q.                \tag{7.34}
\]

Suppose an irreducible factor `R(P) != P` occurs in `C` with multiplicity
`h`.  The locus `R=F=0` is a height-one locus on the cover and is contained
in `H=(H_{N,C})_T=0`, so every normalized prime above it belongs to the
different support.  At any such prime `D` where `Q,T,T-1` are units, the
cover equation gives

\[
 F=-\frac{C(P)QT}{T-1},\qquad
 (H_{N,C})_T=-\frac{C(P)Q}{T-1}.                 \tag{7.35}
\]

and therefore the exact normalized coefficient is

\[
 \operatorname{ord}_D((H_{N,C})_T)
   =h\operatorname{ord}_D(R).                    \tag{7.35a}
\]

It need not equal `h`: normalization can increase `ord_D(R)`.  For example,
write `N=p^n m` and set `q=p^n`.  Over a vertical factor `R=P-a` with
`a!=0`, away from the collision `a^(N+1)=1`, the polynomial `F` has `m`
distinct residue roots, each of multiplicity `q`; on each corresponding
normalized branch `ord_D(R)=q`, so (7.35a) equals `h q`.  At a collision
with `T=1`, division by `T-1` is unavailable and the coefficient requires
its own Newton polygon, but the inclusion in the different support remains.

Consequently, if the permitted target different support is contained in
`P Q=0`, unique factorization in `k[P]` forces

\[
 C(P)=cP^a,\qquad c\in k^*,\quad a\geq0.          \tag{7.36}
\]

The scalar `c` is absorbed by rescaling `Q`.  It remains to classify every
nonnegative integer `a`.

For `a=0`, write `B_{N,0}=k[P,Q,T]/(H_{N,0})`.  The equation is irreducible
because it is primitive and linear in `Q`.  Its Jacobian singular locus is
the single point `(P,Q,T)=(0,0,0)`: `H_Q=T`, then the equation forces `P=0`,
and `H_T=F+Q` forces `Q=0`.  The hypersurface is `S_2` and regular in
codimension one, hence normal.  Let `E=(Q,F)` and

\[
 U_{N,0}=\operatorname{Spec}(B_{N,0})\setminus E,
 \qquad L=(P,T).
\]

The affine UFD core

\[
 W=D_{U_{N,0}}(TF)
   =\operatorname{Spec}k[P,T^{\pm1},F^{-1}]       \tag{7.37}
\]

has unit lattice generated by `T,F`, and `L` is its complete codimension-one
complement in `U_{N,0}`.  At the generic point of `L`, `P` is a uniformizer
and

\[
 \operatorname{ord}_L(T)=N+1,qquad
 \operatorname{ord}_L(F)=N+1.                    \tag{7.38}
\]

The class-trivial-core valuation matrix is therefore the single row

\[
 V_{N,0}=\begin{pmatrix}N+1&N+1\end{pmatrix}.     \tag{7.39}
\]

Its kernel is generated by `(1,-1)` and its cokernel is `Z/(N+1)`.  Hence

\[
 \boxed{\Gamma(U_{N,0},\mathcal O)^*/k^*\simeq\mathbb Z,
 \qquad \operatorname{Cl}(U_{N,0})\simeq\mathbb Z/(N+1).}    \tag{7.40}
\]

Thus the previously omitted endpoint `a=0` fails both affine-plane gates.

For `a>=N`, the Newton point `(1,a)` lies strictly above the segment joining
`(0,N+1)` to `(N,0)`.  The lower polygon is the single reduced-slope segment

\[
 -\frac{N+1}{N}.                                  \tag{7.41}
\]

It gives a branch of ramification index `N` over `P=0`.  Since `p|N`, this is
a wild affine ramification component and cannot survive in a Keller source.
Combining this with Corollary 7.5 gives the complete classification

\[
\begin{array}{c|c}
a=0&\text{free unit and }\mathbb Z/(N+1),\\
1\leq a<N-1&\text{companion tame different},\\
a=N-1&\mathbb Z/(N-1)\text{ generated by }[L_1],\\
a\geq N&\text{wild index-}N\text{ branch over }P=0.
\end{array}                                       \tag{7.42}
\]

Therefore among all covers (7.33) with no target different divisor away from
`P Q=0`, the unique affine-plane Keller reconstruction in this architecture
is `N=2,C(P)=cP`, up to rescaling `Q`.  This is stronger than the finite
monomial-band statement: it exhausts every polynomial coefficient depending
only on `P` under the same one-boundary support hypothesis.

### Corollary 7.7 -- arbitrary target-polynomial coefficients

The support reduction is not one-variable.  Let

\[
 H=(T-1)F+C(P,Q)QT,\qquad 0\ne C(P,Q)\in k[P,Q].  \tag{7.43}
\]

The same identity `(T-1)H_T-H=-C(P,Q)Q` shows that an irreducible factor
`R(P,Q)` of `C`, not associated to `P` or `Q`, creates a height-one different
component on `R=F=0`.  At every normalized prime with `T-1` a unit its exact
coefficient is again `h ord_D(R)`, where `h` is the multiplicity of `R` in
`C`; collision primes remain support-positive and require their local Newton
polygon for the coefficient.  Hence

\[
 \operatorname{Supp}(\operatorname{Diff})_{\rm target}
 \subseteq V(PQ)
 \quad\Longrightarrow\quad
 C(P,Q)=cP^aQ^b                              \tag{7.44}
\]

for some `c in k^*` and `a,b>=0`.  Thus the arbitrary target-polynomial
coefficient search collapses to a two-parameter monomial quadrant.  Theorem
7.6 completely closes its `b=0` edge.  Proposition 7.8 compiles the local
`b>0` ledger, and Theorem 7.9 plus Corollary 7.10 close the remaining
monomial quadrant globally.

### Proposition 7.8 -- the first `b>0` ledger

Put `c=b+1` in the remaining monomial family

\[
 H_{N,a,b}=(T-1)F+P^aQ^cT.                      \tag{7.45}
\]

At the generic fierce component over `Q=0`, the quantities `P,T,T-1` are
units and `Q` remains a uniformizer.  The residue equation is still
`T^N=P^(N+1)`, so its mixed residue degree and sheet loss are unchanged.
On the cover,

\[
 (H_{N,a,b})_T=-\frac{P^aQ^c}{T-1},
\]

and the exact generic row is therefore

\[
 \boxed{(e,f_{\rm sep},f_{\rm insep},d,
          \text{boundary conductor},\text{sheet loss})
  =(1,m,p^n,b+1,N(N-2),N).}                     \tag{7.46}
\]

The conductor entry here is the conductor exponent of the reduced fierce
boundary curve; that curve is still `T^N=P^(N+1)`, so its normalization
semigroup and conductor do not depend on `b`.

The collision model changes.  At any of the `N+1` points with `T=1`,
`Q=0`, and `P^(N+1)=1`, put `u=T-1`.  Since `F_P=-P^N` is a unit, `F` is a
regular transverse coordinate.  The exact rearrangement

\[
 uF+P^aQ^c(1+u)
 =u(F+P^aQ^c)+P^aQ^c                         \tag{7.47}
\]

gives, after unit coordinate changes, the completed local equation

\[
 uv+Q^c=0.                                      \tag{7.48}
\]

Thus `b=0` has a smooth total-space crossing, while `b>0` has an isolated
`A_b`-type hypersurface singularity (with the usual characteristic-dependent
qualification on the rational-double-point label).  Its reduced special
fiber remains the same node.  These singular points lie on the omitted
fierce component, so their existence alone does not exclude the reconstructed
open.  Equations (7.46)--(7.48) reduce the `b>0` programme to exact global
normalization and boundary-class questions.  The following theorem resolves
them without requiring a closed formula for every normalization ring.

### Theorem 7.9 -- base-change and multiple-fibre obstruction

Consider the balanced thickening

\[
 \widehat H_{N,c}=(T-1)F+P^{N-1}Q^cT,
 \qquad c=b+1.                                    \tag{7.49}
\]

Let `U_(N,c)` be its normalization with the fierce prime removed.  Every
`c>=2` gives a non-affine-plane open.  Thus no `b>0` balanced row can be a
Keller affine-plane reconstruction.

For the proof, write `c=p^s d`, with `(d,p)=1`.  Whether `U_(N,c)` is an
affine plane can be tested after extending the ground field, so assume that
`k` is algebraically closed.  Substitution `S=Q^c` gives a finite map

\[
 \pi_c:U_{N,c}\longrightarrow U_{N,1}             \tag{7.50}
\]

of degree `c`.  The degree is exact because the target function `S` has
valuation one on the retained divisor over `S=0`, so `Z^c-S` is Eisenstein
there.  Its `p^s` part is radicial.  Its degree-`d` part is Kummer and is
finite etale away from that retained divisor.

The base surface has

\[
 \chi_{\mathrm c}(U_{N,1})=1.                    \tag{7.51}
\]

For `N>2`, Proposition 7.4 decomposes it into `G_m^2`, the ordinary fill
`A_0=A^1`, and `L_1=G_m`, whose compactly supported Euler characteristics
are `0,1,0`.  For `N=2`, this is the known affine-plane reconstruction, so
the same value follows directly.  The retained branch over `S=0` is

\[
 B\simeq\mathbb A^1\setminus\mu_{N+1},
 \qquad \chi_{\mathrm c}(B)=-N,                  \tag{7.52}
\]

because its `N+1` meetings with the deleted fierce divisor are absent.
Universal homeomorphisms preserve the ell-adic Euler characteristic by
[topological invariance of the etale site](https://stacks.math.columbia.edu/tag/04DY),
while finite etale degree-`d` covers multiply it.  Splitting (7.50) into the
branch and its complement therefore gives

\[
 \boxed{\chi_{\mathrm c}(U_{N,c})
   =d(\chi_{\mathrm c}(U_{N,1})-\chi_{\mathrm c}(B))
      +\chi_{\mathrm c}(B)
   =(N+1)d-N.}                                    \tag{7.53}
\]

If `d>1`, this is not one, so `U_(N,c)` is not `A^2`.

It remains to consider `d=1`, so `c` is a pure `p`-power.  If `N>2`, finite
push--pull on Weil classes gives

\[
 \pi_{c*}\pi_c^*[L_1]=c[L_1].                    \tag{7.54}
\]

The class `[L_1]` has exact order `N-1`, and `gcd(c,N-1)=1`.  Hence its
pullback still has exact order `N-1`; in particular `Cl(U_(N,c))` is nonzero.

The only rows not covered by (7.53) or (7.54) have `p=N=2` and
`c=2^s`, with `s>=1`.  They admit a uniform multiple-fibre calculation.
Pull back the functions `x,u` from `U_(2,1)` and put

\[
 W_c=D_{U_{2,c}}(xu),\qquad t=x/u,
 \qquad A=1-t^3,\qquad r=c/2.                    \tag{7.55}
\]

The base-change equation on this open is

\[
 \boxed{W_c=\operatorname{Spec}
 k[t^{\pm1},u^{\pm1},Q]/
 (Q^ct^4u^2+Au+1).}                               \tag{7.56}
\]

Indeed, `P=tT`, `F=T^2/u`, and
`T=(u-1)/(t^3u)` on `D(xu)`, which gives (7.56).  Its `t`-derivative is the
unit `t^2u`, so this base-changed open is already smooth and normal.  Define

\[
 z=Q^rt^2u+1.
\]

Then

\[
 z^2=Au.                                           \tag{7.57}
\]

For the three cube roots `zeta in mu_3`, let
`D_zeta=(t-zeta,z)`.  Equation (7.57) gives the exact vertical relations

\[
 \operatorname{div}(t-\zeta)=2D_\zeta,
 \qquad
 \operatorname{div}(z)=\sum_{\zeta\in\mu_3}D_\zeta.         \tag{7.58}
\]

These are all relations among the three vertical classes.  To see this, put
`K=k(t)`.  The generic fibre has coordinate ring

\[
 B_r=K[z^{\pm1},Q]/
 \left(Q^r-\frac{A(z+1)}{t^2z^2}\right).          \tag{7.59}
\]

The finite map `Spec(B_r)->G_(m,z)` is radicial.  After geometric scalar
extension it has the same two points at infinity as `G_m`; taking divisors
on the normal projective completion therefore shows that `B_r^*/K^*` is
free of rank one.  Let `g` be a primitive generator and write `z=a g^h`,
with `a in K^*`.  Since the pullback of every base valuation at `t=zeta` is
even by (7.58), while `ord_(D_zeta)(z)=1`, both `h` and every
`ord_(D_zeta)(g)` are odd.  Hence the valuation vector modulo two of any
principal divisor supported on the `D_zeta` is either zero or `(1,1,1)`.
Together with (7.58), this proves that the kernel of

\[
 \bigoplus_{\zeta\in\mu_3}\mathbb ZD_\zeta
 \longrightarrow\operatorname{Cl}(W_c)
\]

is generated by `2D_zeta` and `sum D_zeta`.  Consequently the vertical
subgroup is

\[
\boxed{\langle[D_\zeta]:\zeta\in\mu_3\rangle
        \simeq(\mathbb Z/2)^2
        \subseteq\operatorname{Cl}(W_c).}         \tag{7.60}
\]

In fact this is the full class group.  Over `K=k(t)`, put `y=z^(-1)` and
write `r=2^h`.  After adjoining an `r`-th root of the nonzero constant
`A/t^2`, the generic-fibre equation (7.59) becomes

\[
 R^r=y^2+y.                                       \tag{7.60b}
\]

It has the explicit geometric parameterization

\[
 R=s^2+s,\qquad y=s^r,
 \qquad
 s=y+\sum_{i=0}^{h-1}R^{2^i}.                    \tag{7.60c}
\]

Thus the smooth projective completion of the generic fibre is geometrically
rational.  The affine point `(z,Q)=(1,0)` is `K`-rational, so the
[genus-zero characterization](https://stacks.math.columbia.edu/tag/0C6L)
identifies the completion with `P^1_K`.  The affine curve omits two points,
one of which is the `K`-rational point `(y,Q)=(0,0)`; the divisor sequence on
`P^1_K` therefore gives `Cl(B_r)=0`.  Localization to the generic
fibre now shows that `Cl(W_c)` is generated by vertical prime divisors.
For `t` away from `mu_3`, the same equation is integral because `z+1` has
valuation one and is not a square; the fibre is a single reduced principal
divisor.  The only nontrivial vertical classes are therefore the three
`D_zeta`.  Hence

\[
 \boxed{\operatorname{Cl}(W_c)\simeq(\mathbb Z/2)^2
        \quad\text{for every }c=2^s,\ s\geq1.}    \tag{7.60d}
\]

The vertical-relation argument has a reusable form.  Let a normal variety
map to a smooth curve, let `D_1,...,D_r` be prime fibres with
`div(s-s_i)=m_iD_i`, and suppose a rational function regular along these
fibres has `div(z)=sum D_i`.  If the generic-fibre unit group modulo base
constants is free of rank one, then the same valuation argument gives

\[
 \boxed{\langle[D_1],\ldots,[D_r]\rangle
   \simeq
   \left(\bigoplus_{i=1}^r\mathbb Z/m_i\right)
      /\langle(1,\ldots,1)\rangle.}              \tag{7.60a}
\]

Indeed, if `g` is a primitive generic unit and `z=a g^h`, then
`1=m_i ord_(s_i)(a)+h ord_(D_i)(g)`.  Thus `h` is invertible modulo every
`m_i`, and all the displayed valuations of `g` are congruent to the same
inverse of `h` modulo their respective `m_i`.  Every principal vertical
relation is therefore diagonal.  The known relations `m_iD_i` and
`sum D_i` are the complete relation lattice.  For a prime `ell`, sort the
numbers `v_ell(m_i)`; the `ell`-primary invariant factors in (7.60a) have
exactly the same exponents with the largest one deleted: after choosing an
index of largest exponent, subtracting the diagonal kills that coordinate
and leaves the other cyclic summands uniquely.  In particular the
group has order `prod_i(m_i)/lcm_i(m_i)` and is nonzero exactly when some
prime divides at least two of the multiplicities.  When all `m_i=m`, this
recovers `(Z/m)^(r-1)`.  This multiple-fibre gate needs only the generic unit
rank, not factoriality or the full class group of a large core.
If, in addition, the generic fibre has trivial class group and every other
vertical prime has trivial divisor class, localization promotes (7.60a) to
the full class group of the total space.  The latter condition holds, for
example, when the base is an open of `A^1` with global parameter `s` and
every other fibre is irreducible and reduced, since such a fibre is
`div(s-s_0)`.  Irreducibility alone over an arbitrary base curve is not
enough.

The class group also lifts back exactly to the reconstructed open.  The
complement of `W_c=D(xu)` in `U_(2,c)` consists of the prime
`A_0=(x=0)` and the prime `L_c=(u=0)`.  They are reduced: on `A_0` the
base-change equation is `Q^c=y`, while on `L_c` it is `Q^c=x^(-1)`.
Here the `N=2` chart has
`S=y+x^5(1+xy)^3`; the two displayed restrictions follow directly from
`u=1+x^3(1+xy)`.  Thus no hidden multiplicity enters the localization row.
The two core units `t=x/u` and `u` have boundary-valuation matrix

\[
 \begin{pmatrix}1&0\\-1&1\end{pmatrix}.           \tag{7.60e}
\]

This matrix is unimodular.  The
[exact sequence for an open](https://stacks.math.columbia.edu/tag/0B5Z)
therefore makes restriction an isomorphism, and

\[
 \boxed{\operatorname{Cl}(U_{(2,c)})
        \simeq\operatorname{Cl}(W_c)
        \simeq(\mathbb Z/2)^2
        \quad(c=2^s,\ s\geq1).}                  \tag{7.60f}
\]

Thus every pure Frobenius row is excluded on the full source open.  For
`c=2`, putting
`w=z/u=Qt^2+u^(-1)` recovers the simpler UFD localization

\[
 (W_2)_w\simeq
 \operatorname{Spec}k[t^{\pm1},w^{\pm1},A^{-1}], \tag{7.61}
\]

but factoriality of this localization is not needed for the uniform
obstruction.

### Corollary 7.10 -- complete arbitrary-coefficient exception theorem

Every monomial in (7.45) is now classified.  At `a=0`, the nonconstant unit
from (7.40) pulls back under the finite dominant substitution `S=Q^c`.  For
`1<=a<N-1`, the generic `P=0` Newton row (7.32) is unchanged because `Q` is
a unit there.  Theorem 7.9 handles `a=N-1`.  For `a>=N`, the wild index-`N`
row (7.41) is likewise unchanged.  Thus

\[
 \boxed{N=2,\qquad a=1,\qquad b=0}               \tag{7.62}
\]

is the unique affine-plane Keller row among all `H_(N,a,b)` with `p|N`.
Combining this with Corollary 7.7 proves the same statement for every
polynomial coefficient `C(P,Q)` whose target different support is contained
in `V(PQ)`: up to rescaling `Q`, the characteristic-two cubic is the unique
row under the one-omitted-fierce-boundary hypothesis.  This is an exception
theorem for this natural boundary architecture, not for arbitrary finite
plane covers.

### Theorem 7.11 -- balanced retained-sheet support theorem

The prescribed-degree construction of Section 5 also has an exact balanced
support sieve.  Let `A(T)` be monic and squarefree with `A(0)!=0`, and put

\[
 H_{N,A}^{\rm bal}
 =A(T)(T^N-P^{N+1})+P^{N-1}QT.                  \tag{7.63}
\]

This remains an irreducible finite flat generically separable cover of degree
`N+r`, where `r=deg(A)`.  Its fierce row is unchanged.  In characteristic
`p|N`, if `F=T^N-P^(N+1)`, direct differentiation gives the exact identity

\[
 \boxed{A(H_{N,A}^{\rm bal})_T-A'H_{N,A}^{\rm bal}
   =P^{N-1}Q(A-TA').}                            \tag{7.64}
\]

On the cover above `D(PQ)`, the function `A` is a unit.  Indeed, `A=0`
would force `T=0`, contrary to `A(0)!=0`.  Therefore the relative derivative
on this open is a unit multiple of `A-TA'`.  The function `T` is also a unit
there, and `H_Q=P^(N-1)T` is a unit.  Hence the total space is regular on
this locus and no normalization can remove the detected relative
ramification.  If `alpha` is a root of
`A-TA'`, then `alpha!=0`; squarefreeness also gives `A(alpha)!=0`.
Setting `T=alpha` in (7.63) determines `Q` rationally from the still-free
parameter `P`, and gives a height-one different component meeting `D(PQ)`.
Consequently

\[
 \operatorname{Supp}(\operatorname{Diff})_{\rm target}\subseteq V(PQ)
 \quad\Longleftrightarrow\quad A-TA'\in k^*.     \tag{7.65}
\]

Write `A=sum_j a_jT^j`.  Since the coefficient of `T^j` in `A-TA'` is
`(1-j)a_j`, condition (7.65) is equivalent to

\[
 \boxed{A(T)=a_0+T B(T^p),\qquad a_0\ne0.}       \tag{7.66}
\]

In particular a monic retained polynomial can pass this support gate only if
`r congruent 1 mod p`, and then the cover degree `N+r` is also congruent to
one modulo `p`.  Conversely, for every such `r`, the squarefree control
`A=1+T^r` realizes (7.66).  Thus balanced gluing does not preserve arbitrary
prime-to-characteristic cover degree: it preserves exactly the degree-one
congruence class within this retained-polynomial architecture.

For `r=1` this is the balanced family already excluded by
`Cl(U_N)=Z/(N-1)` when `N>2`.  The first genuinely new support-admissible
odd-characteristic rows with `N=p` therefore have

\[
 (p,r,\deg H)=(3,4,7),\quad(5,6,11),\quad(7,8,15). \tag{7.67}
\]

These are the rows surviving the support sieve alone.  The next theorem
excludes their proposed affine-plane complements uniformly.

### Theorem 7.12 -- geometric retained-root obstruction

Let `A` have degree `r>=1`, let `C_(N,A)` be the normalization of the surface
(7.63), and let `E_(N,A)` be
the reduced fierce divisor above

\[
 Q=0,\qquad T^N=P^{N+1},
\]

and put `U_(N,A)=C_(N,A)-E_(N,A)`.  The following calculation needs only that
`A` is squarefree and `A(0)!=0`; it does not require the support condition
(7.65).

On `D(A)`, put

\[
 P=xu,\qquad T=x^2u,\qquad
 Q=A(x^2u)\bigl(u-x^{N-1}\bigr).                 \tag{7.68}
\]

These formulas annihilate (7.63).  Conversely, in the function field put
`x=T/P`.  After inverting `A`, it satisfies the monic equation

\[
 x^N+\frac QA x-P=0,                             \tag{7.69}
\]

and `u=x^(N-1)+Q/A`; hence `P=xu`, `T=x^2u`.  Adjoining `x` therefore gives
the finite birational normal algebra

\[
 k[x,u,A(x^2u)^{-1}].
\]

It is exactly the normalization over `D(A)`.

Let the ground field be `F_q`, and write `n_q(A)` for the number of roots of
`A` in `F_q`.  Since every root `alpha` is simple and nonzero, the chart
(7.68) removes `q-1` points for each equation `x^2u=alpha`.  On the actual
surface, the fibre `T=alpha` has reduced support

\[
 \{P=0\}\ \cup\ \{Q=0\},                         \tag{7.70}
\]

two affine lines meeting once, hence `2q-1` points.  The original surface is
already regular along this fibre: at `P=0`, the `T`-derivative contains the
unit `A'(alpha)alpha^N`, while at `Q=0,P!=0`, the `Q`-derivative
`P^(N-1)alpha` is a unit.  Normalization consequently changes nothing there.
Finally `E_(N,A)` is the affine line parametrized by

\[
 P=s^N,\qquad T=s^{N+1},\qquad Q=0.
\]

Thus the exact counts are

\[
\boxed{
\begin{aligned}
 \#C_{N,A}(\mathbb F_q)&=q^2+n_q(A)q,\\
 \#E_{N,A}(\mathbb F_q)&=q,\\
 \#U_{N,A}(\mathbb F_q)&=q^2+\bigl(n_q(A)-1\bigr)q.
\end{aligned}}                                   \tag{7.71}
\]

After a finite extension splitting `A`, the same stratification gives, over
an algebraic closure and with `L=[A^1]`,

\[
 [C_{N,A}]=L^2+rL,\qquad
 [U_{N,A}]=L^2+(r-1)L,\qquad
 \chi_c(U_{N,A})=r,                               \tag{7.72}
\]

where `r=deg(A)`.  Since `chi_c(A^2)=1`, every `r>1` complement is
geometrically non-affine-plane.  This is stronger than a base-field point
count and avoids a separate unit or class-group calculation.

Combining (7.72) with Theorem 7.11 leaves only `r=1`.  Proposition 7.4
excludes that row for every `N>2` by
`Cl(U_N)=Z/(N-1)`.  Because `p|N`, the sole affine-plane possibility in the
entire balanced retained-polynomial architecture is therefore

\[
 \boxed{p=N=2,\qquad r=1,}                        \tag{7.73}
\]

the known characteristic-two cubic.  This proves that characteristic two is
exceptional under the natural hypotheses of squarefree retained boundary,
one omitted fierce divisor, balanced gluing, and target different supported
in `V(PQ)`.  It is not a theorem about arbitrary plane-cover presentations.

## 8. Atlas table and determinant ledger

For the canonical `n=1,r=1` rows, the two decisive components are:

| `p` | cover degree | fierce six-tuple `(e,f_sep,f_insep,d,c,s)` | companion six-tuple `(e,f_sep,f_insep,d,c,s)` | result |
|---:|---:|---|---|---|
| `2` | `3` | `(1,1,2,1,c_E=0,2)` | `(1,1,1,0,(P,T,W),1)` | actual Keller reconstruction |
| `2`, `N=6` | `7` | `(1,3,2,1,c_E=24,6)` | `(5,1,1,4,(P,T,W),5)` | mixed residue row; both one-boundary gates fail |
| `3` | `4` | `(1,1,3,1,c_E=3,3)` | `(2,1,1,1,(P,T,W),2)` | one-boundary Keller gate fails |
| `3`, `N=6` | `7` | `(1,2,3,1,c_E=24,6)` | `(5,1,1,4,(P,T,W),5)` | mixed residue row; both one-boundary gates fail |
| `5` | `6` | `(1,1,5,1,c_E=15,5)` | `(4,1,1,3,(P,T,W),4)` | one-boundary Keller gate fails |
| `7` | `8` | `(1,1,7,1,c_E=35,7)` | `(6,1,1,5,(P,T,W),6)` | one-boundary Keller gate fails |
| general `p` | `p+1` | `(1,1,p,1,c_E=p(p-2),p)` | `(p-1,1,1,p-2,(P,T,W),p-1)` | uniform local gate |

In the fierce tuple, `c_E` is the exponent of the boundary-normalization
conductor; its corresponding `delta` values for `p=2,3,5,7` are `0,2,9,20`.
The fierce component has `N+1` conductor nodes.  The two `N=6` rows make
the split between separable and inseparable residue degree explicit while
holding total sheet loss, conductor, and companion row fixed.

### 8.1 Machine-compiled survivor report

The exact
[survivor artifact](../artifacts/generated-results/plane_wild_boundary_survivor_atlas.json)
keeps proved reconstructions, finite covers without reconstructed sources,
and local comparison rows in separate status classes.  Through cover degree
`15` it records:

| packet | rows | obstructed | still needs reconstruction | known Keller |
|---|---:|---:|---:|---:|
| proved monomial hidden-order controls | `46` | `45` | `0` | `1` |
| original prescribed-degree covers | `23` | `3` | `20` | `0` |
| balanced prescribed-degree covers | `23` | `23` | `0` | `0` |
| additive/AS/AS--Witt/Kummer comparisons | `25` | `0` | `6` | `0` |

The remaining nineteen comparison rows are marked local-only, not as global
survivors.  The five rows that pass the support sieve alone are exactly

\[
 (p,d)=(3,7),(3,10),(3,13),(5,11),(7,15).       \tag{8.3}
\]

Theorem 7.12 now rejects all five by geometric point count, so the balanced
retained-polynomial reconstruction queue is empty, not merely bounded through
degree `15`.

The same compiler runs the unequal multiple-fibre gate on every unordered
packet of length two or three, multiplicities at most `12`, containing a
`p`-divisible entry.  It rejects `186/240`, `103/142`, and `43/77` packets
in characteristics `3,5,7`; the remaining `54,39,34` are exactly the
pairwise-coprime packets.  The minimal abstract templates are `(2,3)`,
`(2,5)`, and `(2,7)`.  They remain abstract design targets for other source
fibrations, not fibre packets derived from the now-excluded covers in (8.3).

### 8.2 First odd row: complete prime-field coefficient scan

For the first row `(p,N,r,d)=(3,3,4,7)`, condition (7.66) leaves exactly

\[
 A=T^4+bT+a_0,qquad b\in\mathbb F_3,quad
 a_0\in\mathbb F_3^*.                            \tag{8.4}
\]

All six polynomials are squarefree.  The exact
[degree-seven scan](../artifacts/generated-results/plane_wild_boundary_p3_degree7_scan.json)
uses Singular normalization and relative Fitting ideals.  In every row the
normalization `C_A` is smooth, its primitive conductor is `(P,T)^2`, its
fierce boundary `E_A` is an affine line, and

\[
 \operatorname{Diff}(C_A/\mathbb A^2)_{\rm red}=E_A. \tag{8.5}
\]

Thus all six pass the complete different-away-from-the-boundary gate.  Their
exact prime-field counts are:

| `(a_0,b)` | `#C_A(F_3)` | `#E_A(F_3)` | `#(C_A-E_A)(F_3)` | result |
|---|---:|---:|---:|---|
| `(1,0)` | `9` | `3` | `6` | point-count obstruction |
| `(1,1)` | `12` | `3` | `9` | survives over `F_3` |
| `(1,2)` | `12` | `3` | `9` | survives over `F_3` |
| `(2,0)` | `15` | `3` | `12` | point-count obstruction |
| `(2,1)` | `9` | `3` | `6` | point-count obstruction |
| `(2,2)` | `9` | `3` | `6` | point-count obstruction |

Because `#A^2(F_3)=9`, four rows cannot be the stipulated affine-plane open
over `F_3`.  The remaining quartics factor respectively as a linear factor
times an irreducible cubic.  Singular normalization over the extension
fields gives

| `A` | field | `#C_A` | `#E_A` | `#(C_A-E_A)` | `#A^{-1}(0)` |
|---|---:|---:|---:|---:|---:|
| `T^4+T+1` | `F_9` | `90` | `9` | `81` | `1` |
| `T^4+2T+1` | `F_9` | `90` | `9` | `81` | `1` |
| `T^4+T+1` | `F_27` | `837` | `27` | `810` | `4` |
| `T^4+2T+1` | `F_27` | `837` | `27` | `810` | `4` |

Thus both apparent survivors fail over `F_27`, since `810!=27^2`.  Formula
(7.71) explains all six prime-field counts and both extension-field jumps,
and (7.72) excludes every retained quartic geometrically.  The degree-seven
row is therefore the first exact odd row to pass the different ledger and
then fail the source gate; it is not an odd-characteristic Jacobian
counterexample.

The local determinant cancellation required of a Keller realization is now
unambiguous.  If `(z,q)` are regular normalization parameters at the generic
point of `E`, then

\[
 v_E\det\frac{\partial(P,Q)}{\partial(z,q)}=1.     \tag{8.1}
\]

If `x,y` are polynomial coordinates on a reconstructed affine-plane open and
the Jacobian in those coordinates is one, the chain rule forces

\[
 v_E\det\frac{\partial(z,q)}{\partial(x,y)}=-1.    \tag{8.2}
\]

Equation (8.2) is the missing reciprocal entry.  The finite cover and its
different do not manufacture it.  For `N>2`, there is already a prior
failure: the positive valuation `N-2` on `L_0` remains after `E` is removed.
In characteristic two with `N=2`, that residual valuation is zero and the
reciprocal fierce entry is supplied by the three compatible reconstruction
charts.  Treating any other finite-cover row itself as a counterexample would
confuse a boundary cover with a polynomial Keller map.

## 9. Characteristic-zero `JC_2` boundary module

The same atlas gives a precise receptacle for a Case-1 residue in
characteristic zero.  This is a formulation, not a completed computation.

Let `S=Spec(A)` be a coefficient chart carrying the complete Case-1 Newton
bands, let `B -> S` be the compiled reduced boundary family, and let

\[
 \nu:\widetilde B\longrightarrow B               \tag{9.1}
\]

be its normalization.  Write

\[
 \mathcal Q_B=\nu_*\mathcal O_{\widetilde B}/\mathcal O_B. \tag{9.2}
\]

This is coherent and finite over the conductor locus.  Locally it is exactly
the module exposed by the atlas:

\[
\begin{array}{c|c}
\text{node}&(k[[a]]\oplus k[[b]])/k[[a,b]]/(ab)\simeq k,\\
\text{cusp}&k[[z]]/k[[z^{m_1},\ldots,z^{m_s}]],
\end{array}                                        \tag{9.3}
\]

with the cusp basis indexed by the semigroup gaps.

Let

\[
 \mathcal L=\omega_{\overline X/S}\otimes
             \pi^*\omega_{\overline Y/S}^{-1}(-\Delta)          \tag{9.4}
\]

be the determinant line after subtracting the compiled different `Delta`.
The source reconstruction jets, target-coordinate jets, and allowed gauge
changes give a finite coherent map

\[
 \Phi:\mathcal J_{\rm src}\oplus\mathcal J_{\rm tgt}
          \oplus\mathcal J_{\rm gauge}
 \longrightarrow
 \mathcal Q_B\otimes\mathcal L.                  \tag{9.5}
\]

The jet orders are not arbitrary: they are the conductor exponents and pole
orders supplied by the boundary compiler.  Define

\[
 \boxed{M=\operatorname{coker}(\Phi).}             \tag{9.6}
\]

The omitted Case-1 bands are now restored as exact universal expressions by
the eight-layer continuation in
[`CASE1_FULL_BAND_CONTINUATION.md`](../plane-jc/CASE1_FULL_BAND_CONTINUATION.md).
After transporting those expressions to the alternate chart, the determinant
residue is a well-defined section `rho in M`.  If it vanishes off the finite
surviving boundary stratum `Z=V(I_Z)`, then canonically

\[
 \boxed{\rho\in H_Z^0(M)=0:_M I_Z^\infty.}         \tag{9.7}
\]

For a presentation `F_1 -> F_0 -> M -> 0` with image `R`, this becomes the
checkable saturation quotient

\[
 H_Z^0(M)=\frac{R:I_Z^\infty}{R}.                 \tag{9.8}
\]

This is the promised coherent local-cohomology class.  The divisorial module
`\mathcal R_N` in (7.8) also shows what this finite-support receptacle cannot
do: a residual determinant defect must first be proved absent away from `Z`.
The balanced calculation supplies the complementary global gate.  There the
residual different vanishes, but

\[
 [L_1]\ne0\in\operatorname{Cl}(U_N),
 \qquad (N-1)[L_1]=0.                            \tag{9.9}
\]

Thus the compiled determinant line must first be trivial in the Picard/class
group of the reconstructed open.  Only after both the divisorial-different
gate and this boundary-class gate vanish can its remaining mismatch be a
finite-support class in (9.7).  In characteristic zero this means that the
Case-1 compiler must output the boundary divisor lattice and the class of
`\mathcal L|_U`, not just the conductor quotient `\mathcal Q_B`.
When a `G_m^2` core is available, the output can be made canonical: list all
codimension-one fill divisors, evaluate two character generators on them,
and compute the full Smith form.  A nonzero kernel is a unit obstruction and
a nonzero cokernel is a Weil-class obstruction.  Only a unimodular matrix
passes this second gate, and even then affineness and coordinate
reconstruction remain separate obligations.
Theorem 7.12 adds a third global check for retained-polynomial charts.  Once
the retained roots split, the proposed open contributes
`[U]=L^2+(r-1)L` and `chi_c(U)=r`.  A node-supported class in (9.7) cannot
turn this into the affine-plane value `1`.  A characteristic-zero Case-1
presentation using an `r>1` retained module must therefore compile additional
global boundary pieces whose scissor/Euler contribution cancels
`(r-1)L` before the finite-support residue is formed.
The survivor compiler now emits the exact relation matrices `R` for the
minimal coprime packets `(2,3),(2,5),(2,7)`.  Their Smith diagonals are all
unimodular, so their vertical packet class group is zero.  This supplies the
lower-right block of `[[V,A],[0,R]]`; it does not supply the source-fill
matrix `V`, the relation-lift corrections `A`, or the jet map `Phi`.  Those
are precisely the data that must come from an actual characteristic-zero
Case-1 boundary presentation before (9.8) can be evaluated.  The remaining
proof obligations now use the common `P0`--`P4`
[support-saturation gate](../plane-jc/CONDUCTOR_JET_TRUNCATION.md#42-integrated-support-saturation-gate):

1. `P0`: compile the branchwise conductor/contact-loss ledger and build
   (9.5) functorially from the conductor square; if the bound fails, restore
   only the displayed lower-jet deficit;
2. `P1`: verify that the archived residue represents the resulting section
   and is zero away from `Z`, after the divisorial-different, boundary-class,
   and retained-Euler gates above pass;
3. `P2`: prove that `Z` contains no minimal Fitting component of `M`;
4. `P3`: prove that `M` is `S_1`, or directly exclude associated primes over
   `I_Z`; and
5. `P4`: apply
   [`SST1`](../verified/SUPPORT_SATURATION_PRINCIPLE.md#plane-jc-conductor-residue)
   to prove (9.8) vanishes, or retain an explicit nonzero torsion witness.

<!-- status-consumer: SST1 12c5cb15e8b6de26 -->

The specific `(72,108)` Case-1 coefficient systems are already empty by the
[determinantal closure](../plane-jc/JC2_72_108_DETERMINANTAL_CLOSURE.md).
Thus (9.6)--(9.8) is a geometric repackaging target and a module for future
Case-1 strata, not a missing step in that numerical exclusion.  The older
open-problem wording which calls those three rows surviving is superseded by
the status ledger.

The abstract truncation theorem is now proved: on branch (i), an expression
with contact loss (lambda_i) is determined in the conductor quotient by
inputs modulo (t_i^{c_i+\lambda_i}), and this bound is sharp in general.  The
stronger dependency ledger tests one conductor-plus-path-loss inequality only
on each actual input-to-output path, while a valuation adapter computes the
available order from a certified frontier of omitted Newton support.  It does
not yet
certify the archived Case-1 truncation because the map `Phi`, its expression
paths, and the band-to-normal valuation data have not been compiled.  Its gain
is that the missing input is now a finite exact ledger with per-input deficits
rather than an unspecified full Newton tail.

## 10. Mixed-characteristic decision gate

The boundary module changes the order of operations:

1. lift the finite order, its normalization module, conductor square, and
   determinant line;
2. compare their derived base changes and Bockstein classes;
3. search for a lifted reciprocal determinant section in (9.5);
4. only after that section exists, ask whether the reconstructed open is an
   affine plane and whether its coordinates algebraize polynomially.

For the literal characteristic-two map, step 3 fails already over `W_2`, and
the failure is invariant under polynomial left--right equivalence.  The
characteristic-divisible boundary presentation (2.2), or a different conductor
presentation with the same coarse row, is therefore the correct object to
lift first.

## 11. Reproduction and claim boundary

Run the dependency-light identity and semigroup audit with

```bash
.venv/bin/python scripts/verify_wild_boundary_atlas.py
```

Compile the status-separated survivor atlas, balanced retained-sheet support
and geometric root-count sieves, bounded odd-characteristic packet scan, and
characteristic-zero module template with

```bash
.venv/bin/python scripts/compile_plane_wild_boundary_survivors.py
```

Regenerate its pinned JSON only after an intentional theorem or bound change:

```bash
.venv/bin/python scripts/compile_plane_wild_boundary_survivors.py --write
```

Run the complete `F_3` coefficient scan and the decisive `F_9/F_27`
extension counts for the first nonlinear balanced retained row with

```bash
.venv/bin/python scripts/verify_plane_wild_boundary_p3_degree7.py
```

Run the bounded exact normalization regressions with

```bash
.venv/bin/python scripts/verify_wild_boundary_atlas.py --singular
```

Run the balanced-gluing normalization and point-count diagnostics with

```bash
.venv/bin/python scripts/verify_wild_boundary_atlas.py --balanced-singular
```

The Singular run independently normalizes the rows

\[
 (p,N)=(2,2),(2,4),(2,6),(3,3),(3,6),(5,5),(7,7), \tag{11.1}
\]

including two non-prime-power controls.  It checks that the normalization
module is generated by `1,P^N/T`, verifies the
presentation (3.4), computes the primitive-order conductor `(P,T)`, decomposes
the reduced pullback of `P=0`, and verifies the relative-different saturation
and radical in (7.7), including both characteristic-two controls.  It is a
bounded regression for Theorem 3.1 and an exact certificate for (7.7); the
all-`N` normalization proof is the `S_2+R_1` argument above, and the all-`N`
companion-branch calculation is (7.5).

The balanced command normalizes
`(p,N)=(2,2),(2,6),(3,3),(3,6),(5,5),(7,7)`, computes the primitive-order
conductors and base-field point counts in (7.14), and for `(p,N)=(3,3)` computes the
complete relative Fitting saturation (7.15) and absolute smoothness.  It
independently replays the normalization for `N=2,3`; the larger rows are
direct algorithm outputs with exact conductor and
quotient-ring checks, not independent `norTest` replays.  The default
dependency-light command verifies the balanced hidden equation, the
polynomial quotient in (7.16), rational reconstruction identities, Jacobian
formula (7.18), and Smith determinant in (7.27)--(7.28) for all seven test
rows.

No command in this note proves any of the following:

- an odd-characteristic polynomial Keller counterexample;
- `Spec(C_N) minus E` is an affine plane for `N>2`;
- the reciprocal determinant entry (8.2) is polynomially realizable;
- the actual Case-1 module (9.6): the bands are now supplied, but the
  alternate-chart matching map and residue class are not yet constructed; or
- a mixed-characteristic lift of a boundary presentation.

<!-- status-consumer: CJT1 afb70f90ff10f3d7 -->

<!-- status-consumer: C1FBC1 0f14ef01fff25097 -->
