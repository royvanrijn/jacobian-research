# Shifted-resultant positivity and the exact limits of contiguity

Internal research audit, 21 November 2026. Base repository commit
`19eb9912fb60142f8b620dac24896931d4ff3e46`. This note proposes research
targets; it does not change the status of unrestricted SFC or any canonical
proof. No novelty or external-review claim is made.

## Recommendation

The most concrete surviving target in this direction is a uniform positive
expansion of the **shifted monomial-linear resultant in its exponent
parameters**. This is different from Gamma-dispersion coefficient
positivity already disproved in the repository. It would close the whole
monomial-times-linear family, including primitive unequal-exponent cases,
but would still need a separate bridge to arbitrary supports. The proposed
Euler-width certificate with a constants-only polynomial kernel has a
valid conditional implication, but the parallel arithmetic review has
now refuted its universal premise; it is not a recommended target.

For a full-support route, investigate the interaction of a distinguished
factorial functional with a finite ordinary multiplication complex. Generic
GKZ regularity, toric Cohen--Macaulayness, invertible contiguity, and absence
of Bernstein--Sato roots at positive orders are already available or automatic.
They do not supply that interaction. Section 2 below makes this distinction
explicit for the actual factorial generating series, rather than an analogy.

The arithmetic alternative is an explicit bound on valuations of maximal
Macaulay minors along the retained digit-marker approximation. Its exact
noncollapse condition is useful as a proof interface, but without a new
bound it is equivalent to the desired rank statement. The currently proved
amplified approximation also preserves the balanced-profile restriction.

## 1. Shifted resultants in the common-monomial exponents

For support size `s`, starting order `n>=1`, and indeterminates
`r=(r_1,...,r_s)`, set

\[
 g_m(c;r)=[z^m]\prod_{i=1}^s(1-c_i z)^{-mr_i-1},
 \qquad
 R_{s,n}(r)=\operatorname{Res}_c(g_n,\ldots,g_{n+s-1}).                 \tag{1}
\]

Each coefficient is a polynomial over `Q` in the `r_i`. For nonnegative
integer exponents these are the actual normalized factorial moments of
`F=(prod_i X_i^(r_i)) sum_i c_i X_i`.

The canonical [initial-resultant proof](SFC_INITIAL_RESULTANT_AND_SHIFTED_WINDOWS.md)
establishes

\[
 R_{s,1}(r)=\prod_{\varnothing\ne I\subseteq[s]}
 \left(1+\sum_{i\in I}r_i\right)^{(|I|-1)!(s-|I|)!}.                 \tag{2}
\]

The value `R_(s,n)(0)` is `+1` or `-1`. One way to see its absolute value
is to use the complete homogeneous functions `g_m(c;0)=h_m(c)` over an
arbitrary field. If the actual active support has size `t`, write
`P(z)=prod_active(1-c_i z)`, so `sum h_m z^m=1/P(z)`. A block of `t`
successive zeros propagates backwards under the recurrence for `1/P`,
because its last coefficient is a nonzero product of active `c_i`. This
would force `h_0=0`. Thus the resultant is nonzero modulo every prime and
is a unit in `Z`. Normalize its orientation by this value, writing
`R^+_(s,n)(0)=1`.

### A sufficient, genuinely stronger positivity conjecture

Expand in the product Newton basis

\[
 R^+_{s,n}(r)=\sum_{\alpha\in\mathbb N^s}
       A_{s,n,\alpha}\prod_i\binom{r_i}{\alpha_i},
 \qquad A_{s,n,\alpha}=\Delta^\alpha R^+_{s,n}(0).                  \tag{3}
\]

**Candidate N.** Every coefficient in (3) is nonnegative, uniformly in
`s,n,alpha`.

Since the constant coefficient is one, Candidate N implies
`R^+_(s,n)(r)>=1` at every nonnegative integer exponent tuple. Hence it
implies arbitrary-start SFC for every common-monomial linear family.
It is stronger than nonvanishing and can be falsified by one exact forward
difference. Ordinary monomial-coefficient positivity is stronger still;
Schur positivity of this symmetric polynomial is another stronger possible
structure. Monomial positivity also gives nonvanishing for all nonnegative
real exponent parameters, which is beyond the ordinary polynomial claim.

This does **not** establish full SFC: arbitrary exponent matrices do not
become `r+I` by a proved support-preserving moment transformation.

### Exact small falsification screen

For `s=2`, use `c=(t,1)` and compute the rational Sylvester resultant
directly. For `n=1,2,3,4,5`, the exact results are:

| n | Total degree of R | Nonzero ordinary coefficients | Negative coefficients | Nonzero Schur coefficients after the factor below | Negative Schur coefficients |
|---|---:|---:|---:|---:|---:|
| 1 | 3 | 8 | 0 | 1 | 0 |
| 2 | 9 | 43 | 0 | 9 | 0 |
| 3 | 18 | 148 | 0 | 41 | 0 |
| 4 | 30 | 386 | 0 | 124 | 0 |
| 5 | 45 | 841 | 0 | 295 | 0 |

In all five cases the removed factor is
`(r+1)(u+1)(r+u+1)^n`; the quotient has respectively 1, 15, 75, 235, and 570
nonzero ordinary coefficients, all positive. All five full resultants
have constant coefficient one. The `n=5` quotient has total degree 38;
that exact calculation took 53.5 seconds with a 120-second cutoff.
For each fixed `n` in the table, these complete symbolic identities prove
strict positivity for every nonnegative real pair `(r,u)`, not just for
sampled exponent values. Consequently each listed two-monomial window
is regular for every nonnegative integer exponent pair over every
characteristic-zero field. The order remains restricted to `1<=n<=5`;
no all-order positivity conclusion follows.

For two variables, the Schur calculation requires no symmetric-function
library: for each homogeneous degree `d`, if `a_k` is the coefficient of
`r^(d-k)u^k`, the coefficient of `s_(d-k,k)` is `a_k-a_(k-1)`, with
`a_(-1)=0`, for `0<=k<=floor(d/2)`.

The repository's existing Gamma counterexample is a different deformation:
its parameter changes the normalized law of a fixed exponent row, rather
than the common-monomial exponents in (1). That counterexample excludes
raw Gamma-resultant positivity, not Candidate N. Conversely these five
binary controls provide no evidence that the excluded Gamma statement can
be repaired merely by calling its coefficients connected diagrams.

### A precise limitation on arrangement-determinant identifications

Already for `s=2,n=2`, the residual factor in the canonical proof is

\[
\begin{split}
Q(r,u)={}&81r^3u^2+74r^3u+17r^3+81r^2u^3+257r^2u^2\\
&+181r^2u+37r^2+74ru^3+181ru^2+115ru+22r\\
&+17u^3+37u^2+22u+4.
\end{split}                                                        \tag{4}
\]

It cannot be a product of affine linear polynomials even over `C`.
Its highest homogeneous part is `81 r^2 u^2(r+u)`. Thus such a product
would have a vertical and a horizontal linear factor. But the gcd in
`Q[r]` of its coefficients as a polynomial in `u` is one, and the
corresponding gcd in `Q[u]` is one. The Euclidean identities for those
gcds remain valid over `C`, excluding both kinds of factors.

Varchenko's [SIGMA 7 (2011), 032](https://sigma-journal.com/2011/032/sigma11-032.pdf),
Theorem 2.1, printed p.11, states nondegeneracy of the arrangement
contravariant form when all dense-edge weight sums are nonzero. For a
fixed arrangement with weights affine in `(r,u)` and no dense-edge weight
identically zero, its nonzero determinant therefore has zero locus
contained in finitely many affine hyperplanes.
Equation (4) prevents identifying the shifted resultant with that full
determinant up to nonzero constants and affine weight factors. A useful
arrangement connection would have to involve a more selective minor,
subquotient, or a different parameter map. It cannot follow from the
initial subset factors alone.

There is also an early obstruction to a positive-semidefinite **affine**
determinant explanation. The one-variable restriction

\[
                    Q(r,0)=17r^3+37r^2+22r+4
\]

has discriminant `-428`, hence a nonreal conjugate root pair. A pencil
`det(B_0+r B_1)` with `B_0` Hermitian positive definite and `B_1`
Hermitian has only real finite roots: conjugate by `B_0^(-1/2)` and
diagonalize the resulting Hermitian matrix. Thus the shifted resultant
cannot be such a pencil, even though all the coefficients examined are
positive. This does not exclude polynomial Gram representations or a
combinatorial positive expansion. As a calibration, the diagonal slice
does factor: `Q(r,r)=(2r+1)(3r+1)^2(3r+2)^2`.

**Work package.** Try to construct a uniform Newton/Schur expansion or a
sign-preserving recurrence for its coefficients. Falsify first at general
ternary exponents and small shifted order; the already proved repeated-
exponent ternary slices are calibration cases. A negative Newton
coefficient stops Candidate N as stated. A positive finite table does not
advance its quantifiers.

### A finite-operator implication whose universal premise is false

The following certificate was proposed in the parallel differential
review. Its implication is elementary, but the premise is now refuted
as a universal route. It remains useful to locate exactly where a
selected-period argument must do more than find an annihilator.

For fixed nonnegative integer `r` and a nonzero coefficient vector `c`,
write `G(t)=sum_(m>=0) g_m(c;r)t^m`, so `G(0)=1`. Seek a nonzero operator

\[
 L=\sum_{j=0}^{w}t^j P_j(\theta),\qquad
 \theta=t\partial_t,\qquad w\le s,\qquad LG=0,                       \tag{5a}
\]

with `P_j` polynomials and

\[
                       \ker L\cap\mathbb C[t]\subseteq\mathbb C.
                                                                    \tag{5b}
\]

There is no upper bound on the degrees of the `P_j` in this target.
Condition (5b) allows the polynomial kernel to be zero or the constants;
it forbids only nonconstant polynomial solutions.

**Proved implication.** Conditions (5a)--(5b), together with the canonical
initial resultant, imply every `s`-term window for this `r,c` is nonzero.
If orders `n,...,n+s-1` vanished, let
`Q=sum_(m<n)g_m t^m`. Then `ord(G-Q)>=n+s`. Each summand of `L` preserves
or increases order, so `ord L(G-Q)>=n+s`, whereas
`deg LQ<=n-1+w<=n+s-1`. Since `LQ=-L(G-Q)`, both sides are zero.
Condition (5b) makes `Q=1`. Combining its coefficients with the proposed
gap gives `g_1=...=g_s=0`, contradicting (2). QED.

This is a genuinely additional certificate requirement, not a theorem
following from generic holonomicity. It is also not the rejected
``homogeneous recurrence order at most s'' statement: the differential
order can be large, and a constant solution is explicitly permitted.
Polynomial-kernel control is finite once `L` is known. A polynomial
solution of degree `d` must satisfy `P_w(d)=0` when `P_w!=0`; the finitely
many nonnegative integer roots give candidate degree bounds, after which
one solves the exact coefficient system. This finite test can certify individual parameter strata; the following
counterexample rules out its universal premise.

The [parallel differential audit](SFC_DIFFERENTIAL_BRANCH_OBSTRUCTIONS.md),
Section 4, gives the actual profile `r=(1,2)`, `c=(1,u)`, with
`12u^2-5u+5=0`. Its nonselected residue-period control has `s=2`,
`ell_2=ell_3=0`, and `ell_1!=0`. It is a constant linear combination of
conjugates of the actual `rho`, so every rational homogeneous annihilator
of `rho` also annihilates it. For any such operator of width at most two,
the same truncation argument forces the **nonconstant** polynomial
`ell_0+ell_1 t` into its kernel. Thus (5b) fails at every differential
order in this control. The period is not the selected factorial series;
this is not an SFC counterexample. It excludes this universal operator
certificate even after constant polynomial solutions are allowed.

## 2. The actual GKZ lift is already Cohen--Macaulay

Let `A` be a nonnegative integer `d by s` matrix with nonzero columns
`a_j`. For shapes `beta_i>0`, define the formal series

\[
 G_\beta(c)=\sum_{k\in\mathbb N^s}
       \frac{\prod_i(\beta_i)_{(Ak)_i}}{\prod_j k_j!}c^k.             \tag{5}
\]

No analytic convergence near `c=0` is asserted. Its homogeneous pieces
are Gamma moment forms; at `beta=(1,...,1)` they are ordinary factorial
moments divided by `m!`.

Introduce coefficient variables `u_1,...,u_d` for the linear damping
monomials. The corresponding formal exponential integral is

\[
 \mathcal G_\beta(u,c)=\frac1{\prod_i\Gamma(\beta_i)}
 \int_{\mathbb R_+^d}x^{\beta-1}
       \exp\left(\sum_i u_i x_i+\sum_jc_jx^{a_j}\right)dx,           \tag{6}
\]

interpreted coefficientwise in `c` near `u=(-1,...,-1)`. Every coefficient
in that interpretation is a convergent Gamma integral, so the following
identities can alternatively be checked term by term.

The confluent GKZ matrix is

\[
                      B=[I_d\mid A].                               \tag{7}
\]

Its toric ideal is exactly

\[
 I_B=(\partial_{c_j}-\prod_i\partial_{u_i}^{a_{ij}}:1\le j\le s),
 \qquad S_B\simeq\mathbb C[t_1,\ldots,t_d].                         \tag{8}
\]

Indeed the quotient substitutes each `partial_(c_j)` by its indicated
monomial in the independent `partial_(u_i)`, proving both inclusions and
the isomorphism. Its semigroup is `N^d`, and its lattice is `Z^d`.
The Euler equations of (6) have GKZ parameter `alpha=-beta`; the sign
comes from integration by parts in `x_i`, or directly from the Gamma
coefficient formula. A constant support monomial factors out as
`exp(c_0)` and is omitted from (7) to keep the positive grading hypotheses.

Thus the natural toric ring is normal and Cohen--Macaulay for **every**
ordinary support. Theorem 9.1 and Corollary 9.2 of
[Matusevich--Miller--Walther, *Homological methods for hypergeometric families*](https://arxiv.org/pdf/math/0406383),
printed pp.19--20, identify GKZ rank jumps with exceptional local-
cohomology quasidegrees, and exclude them for a Cohen--Macaulay toric
ring. Their hypotheses include a pointed integer configuration and the
Euler--Koszul construction. Both apply to (7). This controls the
hypergeometric module's holonomic rank; it says nothing by itself about
zeros of (5) at a specified complex coefficient vector.

### Positive shape contiguity is also already invertible

Use the grading `deg partial_j=-b_j` of
[Schulze--Walther, *Hypergeometric D-modules and twisted Gauss--Manin systems*](https://agag-mschulze.math.rptu.de/download/ekdi.pdf).
Remark 3.7, printed p.7, gives the exact criterion for the **right**
contiguity map
`K(S_B,alpha) -> K(S_B,alpha+v)` induced by a column `v`:
it is a quasi-isomorphism precisely when

\[
                  -\alpha\notin v+\operatorname{qdeg}(S_B/t^v).
                                                                    \tag{9}
\]

This is distinct from their Theorem 3.6 about left multiplication on one
fixed complex. For the polynomial semigroup ring in (8), the quasidegrees
in (9) are the union, over `i` with `v_i>0` and `0<=h<v_i`, of the
hyperplanes `z_i=-h`. Set `alpha=-beta-v`. The forbidden coordinates of
`-alpha=beta+v` would be `v_i-h`; none occurs when `beta_i>0`.
Consequently positive shape increases `beta -> beta+v` are already
invertible contiguity operations, with the solution-space direction
opposite to the displayed module map.

This does **not** say that positive Gamma parameters are strongly
nonresonant, or that localization by every coordinate is an isomorphism.
For an identity column `v=e_i`, the source `alpha=-beta` is forbidden
for the right shift `alpha -> alpha+e_i` when `beta_i=1`. The positive
shape increment instead uses the source `alpha=-beta-e_i`, for which
the same excluded coordinate would require `beta_i+1=1`. The direction
of the parameter shift is essential. The numbering quoted here is that
of the linked author PDF; the corresponding remark in arXiv v3 is 3.6.

At the normalized formal-series level, the same fact is elementary:

\[
 G_{\beta+e_i}=\frac{\beta_i+\sum_j a_{ij}\theta_j}{\beta_i}G_\beta,
 \qquad \theta_j=c_j\partial_{c_j}.                                 \tag{10}
\]

The operator has a positive nonzero eigenvalue on every monomial when
`beta_i>0`. The actual Horn equations are

\[
 \partial_{c_j}G_\beta=
 \prod_i\prod_{\ell=0}^{a_{ij}-1}
   (\beta_i+\sum_h a_{ih}\theta_h+\ell)G_\beta.                      \tag{11}
\]

Neither (10) nor (11) asserts that evaluating a polynomial at a point
commutes with its derivatives. Nor does an invertible linear operator on
homogeneous pieces preserve ordinary multiplication ideals.

### An exact countercontrol using genuine positive shape contiguity

Take two independent linear Gamma features, initially of shapes `(1,1)`.
The relevant forms, with the `m!` normalization in (5), are

\[
 p_2=x^2+xy+y^2,\qquad p_3=x^3+x^2y+xy^2+y^3.
\]

Raise the first shape to two. The contiguity operator is `T=1+theta_x`,
giving

\[
 q_2=3x^2+2xy+y^2,\qquad q_3=4x^3+3x^2y+2xy^2+y^3.
\]

Both pairs are regular: their dehomogenized resultants are respectively
`1` and `16`. Nevertheless

\[
 T^{-1}(xq_2)=\tfrac34x^3+\tfrac23x^2y+\tfrac12xy^2
                         \notin(p_2,p_3).                           \tag{12}
\]

The degree-three functional `[x^2y]-[xy^2]` annihilates `xp_2`, `yp_2`
and `p_3`, but evaluates to `1/6` on (12). Thus even an exact positive-
shape GKZ isomorphism does not give the proposed ordinary ideal map.
This does not disprove regularity preservation; both endpoints here are
regular, and their saturated ideals are both unit ideals.

Shape contiguity also keeps the integer matrix fixed. Adding a private
exponent row, merging rows, or changing a common-monomial exponent changes
that matrix. These operations are not merely positive parameter shifts
to which (9) automatically applies.

### What Bernstein--Sato theory does and does not add

Kashiwara, [*B-functions and holonomic systems*, Invent. Math.38 (1976),33--53](https://www.kurims.kyoto-u.ac.jp/~kenkyubu/kashiwara/Bfunct.pdf),
Corollary 5.2, printed p.48, proves that the roots of the Bernstein--Sato
polynomial of a nonzero holomorphic function are negative rational
numbers. The proof uses the resolution/divisibility theorem 5.1 and the
explicit normal-crossing product. Thus its shift scalar is nonzero at
nonnegative integer powers. This is a statement about differential
operators applied to `F^(m+1)`, not about their factorial integrals.

For instance `F=1+x^2` has the exact identity

\[
 [(m+1)-\tfrac{x}{2}\partial_x]F^{m+1}=(m+1)F^m.                    \tag{13}
\]

Writing `M_m=Lambda(F^m)`, factorial integration gives

\[
 (m+1)M_m=(m+\tfrac32)M_{m+1}
                   -\tfrac12\Lambda(xF^{m+1}).                     \tag{14}
\]

The additional inserted moment cannot be dropped. More generally
`Lambda(partial_i h)=Lambda(h)-Lambda_face_i(h|_(x_i=0))`; higher
operators bring both insertions and boundary terms. Every nonconstant
squarefree univariate polynomial already has `b_F(s)=s+1`: choose `UF+VF'=1`
and use `((s+1)U+V partial)F^(s+1)=(s+1)F^s`. Thus nonresonance of
this scalar alone cannot capture the unresolved support-length bound.

Oaku's [*Algorithms for integrals of holonomic functions over domains
defined by polynomial inequalities*](https://www.lab.twcu.ac.jp/oaku/Oaku2013.pdf),
Theorem 1, printed p.14 in the author preprint, gives an integration ideal
annihilating a well-defined distributional integral when the input ideal
is holonomic. Its proof retains total derivatives and the correct domain
distribution. It guarantees a holonomic system, not a support-size
recurrence or a nonvanishing scalar period. Its boundary-aware method is
an appropriate tool for constructing certificates, with those limits.

**Work package.** Start with the actual finite Macaulay multiplication map,
not a generic GKZ rank statement. Attempt a boundary-aware contiguity or
integration certificate for its maximal minors. An ordinary unit-row
adjoining theorem plus an ordinary integer-row-merging theorem would
give full SFC, as shown in the retained unrestricted-mechanisms note.
A universal private-row *removal* rule is already equivalent to full
SFC. These are proof targets, not consequences of (8)--(11). Stop a
proposed construction if it identifies an Euler--Koszul complex with the
ordinary moment Koszul complex without an explicit map, replaces the
selected factorial functional by all periods, or suppresses terms such as
the last term of (14).

## 3. Arithmetic noncollapse: exact condition and scope

Let `d_i=n+i`, `0<=i<s`, and

\[
 N=1+\sum_i(d_i-1)=1+s(n-1)+s(s-1)/2.
\]

Let `M_0` be the integer coefficient matrix of the ordinary multiplication
map `direct_sum R_(N-d_i) -> R_N` for the unnormalized factorial moment
forms. Full row rank is equivalent to absence of a common projective
zero. This equivalence is proved in Part B.4 of the canonical
[de-amplification audit](SFC_DEAMPLIFICATION_BRIDGES.md).

Suppose integer matrices `M_r` of the same size satisfy
`M_r == M_0 mod p^r`. Let `v_r` be the least `p`-adic valuation among
their maximal row-size minors, with infinity if all are zero.

**Elementary noncollapse lemma.** If `M_0` has full row rank and `v_0`
is its least minor valuation, then `v_r=v_0` for every `r>v_0`.
If `M_0` has rank deficit `delta>=1`, then `v_r>=delta*r`.

**Proof.** Determinants preserve the congruence, giving the first claim.
For the second, expand each determinant by multilinearity in its columns
using `M_r=M_0+p^r E_r`. A term with fewer than `delta` columns from `E_r`
contains more columns from `M_0` than its rank and vanishes. Every surviving
term contains at least `delta` factors `p^r`. QED.

Thus `liminf v_r/r<1` is sufficient for regularity, and eventual boundedness
of `v_r` is actually equivalent to it. The useful new target is an
**explicit arithmetic upper bound**, derived before knowing the target
rank, not the assertion that boundedness must hold.

For a balanced base support, the retained digit construction supplies
such matrices. Put `M=n+s-1`, choose a prime `p>M`, and for precision `r`
let

\[
 L=(p-1)p^{r-1},\qquad D_r=\sum_{i=1}^{L}p^{ri}.
\]

The exact retained theorem says
`K_(D_r)(k)=(D_r k)!/(D_r!)^k == 1 mod p^r` for `0<=k<=M`.
Decorating each column with its private `X_j^(D_r)` and replacing `c_j`
by `c_j/D_r!` therefore yields integer forms congruent to the originals.
Their Macaulay matrices are the required `M_r`. For all sufficiently
large `r`, balanced private amplification proves their full rank.

This alone does not bound `v_r`: full-rank matrices can converge to a
rank-deficient one, as `diag(1,p^r)` shows. Crude archimedean heights also
do not give a useful precision-independent bound, because these marker
degrees grow enormously with `r`.

The balanced hypothesis cannot be discarded here. Profile completion adds
fixed factors which this marker limit does not remove. In fact every
balanced coefficient family, after arbitrary coefficient rescaling and
separate nonzero normalization at each degree, satisfies the closed pure-
coefficient identities

\[
 a_{m,j}^{\ell}a_{\ell,k}^{m}
                 =a_{m,k}^{\ell}a_{\ell,j}^{m}.                     \tag{15}
\]

An arbitrary unbalanced factorial family need not satisfy (15). Thus it
cannot be recovered as a coefficientwise limit of those balanced frames
with the same normalizations. The canonical audit already records this
obstruction. The digit approximation is an exact potential bridge for
balanced supports; extending it to arbitrary supports requires another
construction or a genuine integer removal theorem.

**Work package.** Seek a valuation bound for a specified maximal minor or
the ideal of all such minors, exploiting digit structure and exact
cancellation. A generic claim that a nonzero amplified determinant stays
nonzero in the limit is invalid. A bound depending on precision as `r` or
faster is insufficient. Do not market the unrestricted bounded-valuation
assertion as an easier theorem: its nonquantitative content is the missing
regularity itself.

## 4. Priority and finite checks

1. **Shifted exponent-resultant structure:** the most concrete surviving
   route here to an
   all-order, arbitrary-support-size theorem inside the monomial-linear
   class. The initial resultant is an exact anchor and finite differences
   provide direct falsification tests. Seek a uniform identity, not more
   isolated exponent specializations. The constants-only Euler-width
   proposal is a failed shortcut, as recorded above.
2. **Distinguished-functional descent:** highest possible full-SFC payoff,
   but requires genuinely new compatibility with ordinary multiplication.
   The natural GKZ and Bernstein--Sato facts above eliminate several
   tempting prerequisite projects that would not address it.
3. **Arithmetic minor noncollapse:** a precise alternative for balanced
   frames. Continue only with a candidate valuation estimate; the bare
   limit argument and the unbalanced extension are already blocked.

The exact checks use rational SymPy arithmetic: the five
binary Sylvester resultants and their factorizations/Schur expansions;
the two coefficient gcds, leading form, and cubic discriminant of (4); the resultants and
nonmembership witness in (12); and the specialization of (13) at
`m=0,...,5`. The all-parameter identities and the noncollapse lemma have
the written proofs above. The checker provides finite exact replay and source-pin validation; it does
not verify the universal conjecture, formalize the written proofs, or imply
external review. Registry and publication status are unchanged.

For reproducing the main finite screen, the entire construction is:

```python
import sympy as S
r, u, t = S.symbols('r u t')
def g(m):
    return S.expand(sum(
        S.rf(m*r+1,k)*S.rf(m*u+1,m-k)*t**k
        / (S.factorial(k)*S.factorial(m-k)) for k in range(m+1)))
for n in (1,2,3,4,5):
    R = S.Poly(S.resultant(g(n),g(n+1),t),r,u)
    Q = R.exquo(S.Poly((r+1)*(u+1)*(r+u+1)**n,r,u))
    ordinary_negative = sum(bool(a < 0) for a in R.coeffs())
    schur = []
    for d in range(Q.total_degree()+1):
        previous = 0
        for k in range(d//2+1):
            current = Q.coeff_monomial(r**(d-k)*u**k)
            schur.append(current-previous)
            previous = current
    print(n, R.total_degree(), len(R.terms()), ordinary_negative,
          sum(bool(a != 0) for a in schur), sum(bool(a < 0) for a in schur))
```

## Exact replay

The companion checker requires Python 3 and the pinned SymPy version 1.14.0.
It checks all five symbolic resultants, exact common-factor division, full
Schur reconstruction, independently computed numerical Sylvester determinants,
the cubic discriminant, and the positive-shape ideal nonmembership witness.
Its note and existing source files are pinned by SHA-256. It rejects optimized
Python execution and missing or altered pinned inputs before writing a receipt.
The cross-linked differential obstruction is a separately reviewed result;
its checker is not silently included in this replay.

```sh
python3 research/factorial/image-rigidity/verify_sfc_shifted_resultant_positivity.py
```
