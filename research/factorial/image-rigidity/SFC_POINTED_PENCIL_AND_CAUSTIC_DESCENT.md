# Pointed polynomial pencils, Burgers inversion, and the sharp contact gap

**Status.** The pencil identities and the caustic-field proposition below are
proved. The closing lemmas in Section 6 are proposed research targets, not
theorems. They do not establish unrestricted SFC. The proved caustic certificate
is distinct from the proposed finite-contact-to-composition step.

## 1. The actual positive-exponent frame

Work over an arbitrary field `K` of characteristic zero; pass to an algebraic
closure only to discuss individual critical points. For an ordinary
monomial-linear phase with `S` active terms and positive integral exponents, put

\[
 P(z)=\prod_{i=1}^{S}(1-c_i z),\quad
 A(z)=\prod_{i=1}^{S}(1-c_i z)^{r_i},\quad
 a(z)=zA(z),\quad H(0)=0,\quad H'(z)=A(z)/P(z).
\]

The nonzero coefficients `c_i` may repeat. More generally this polynomial setup
applies to grouped data whenever every group total exponent is at least its
number of variables. Set

\[
 W=H',\qquad D=a'/H',\qquad d=\deg a,\quad h=\deg H.
\]

Then `a(0)=H(0)=0`, `a'(0)=H'(0)=D(0)=1`, `d-h=S`,
`deg D=S`, and `W` divides both `a` and `a'`. The distinguished inverse branch
gives

\[
 B(t)=H(a^{-1}(t)),\qquad \rho(t)=B'(t)=1/D(a^{-1}(t)).
\]

The desired bound is

\[
 \operatorname{ord}_0(H-b(a))\le n+S
 \quad(\deg b\le n,\ n\ge1).
 \tag{1}
\]

The repository already proves uniform but larger contact and window bounds in
[the algebraic-window note](SFC_ALGEBRAIC_UNIFORM_WINDOWS.md).
The point here is to retain structure lost by a general homogeneous linear
annihilator. We do not treat the cases where `H` has logarithmic terms.

## 2. The exact polynomial pencil and its moving caustic

Consider

\[
 a_\varepsilon(z)=a(z)+\varepsilon H(z),\qquad
 a_\varepsilon'(z)=W(z)(D(z)+\varepsilon).
 \tag{2}
\]

For generic `epsilon` the moving critical points are exactly the `S` roots of
`D+epsilon`, counted with multiplicity. A root `alpha` of `W` is a fixed critical
point. Because `a(alpha)=0`, its critical-value line is

\[
 t=\varepsilon H(\alpha).
 \tag{3}
\]

All these lines pass through `(epsilon,t)=(0,0)`. This concurrence is an exact
consequence of the factorial root alignment. Merely knowing that `a'/H'` is a
polynomial of degree `S` does not give it.

The derivative resultant factors exactly as

\[
 \operatorname{Res}_z(a+\varepsilon H-t,a_\varepsilon')
 =\operatorname{Res}_z(a+\varepsilon H-t,W)\,
  \mathcal C(t,\varepsilon),
 \qquad
 \mathcal C=\operatorname{Res}_z(a+\varepsilon H-t,D+\varepsilon).
 \tag{4}
\]

The first factor is a nonzero constant times the product of the lines in (3),
with multiplicities `ord_alpha W`. The moving factor `C` has degree exactly `S`
in `t`: use the norm over the `S` roots of `D+epsilon`; its leading coefficient
in `t` is a nonzero constant. The moving critical-value curve is parametrized by

\[
 \varepsilon=-D(z),\qquad t=T(z):=a(z)-D(z)H(z).
 \tag{5}
\]

Its distinguished point is `(-1,0)`, corresponding to `z=0`. It is different
from the concurrence point `(0,0)` of the fixed lines. Branches or singularities
at `(-1,0)` must be handled on the normalization, retaining the marked preimage
`z=0`.

## 3. A proved nonlinear certificate for genuine pullback descent

**Proposition.** Let `K` be any field of characteristic zero and let
`a,H` be nonconstant polynomials in `K[z]`, with `a'=H'D` for a polynomial
`D` of degree `S>=1`. Thus `d=deg a=deg H+S`. Let `delta` be the nonzero
leading coefficient of `D`, let `T=a-DH`, and put

\[
 e=[K(z):K(a,H)].
\]

Then

\[
 \boxed{K(D,T)=K(a,H)},\qquad
 \boxed{\mathcal C(t,\varepsilon)=c\,F(t,\varepsilon)^e},
 \qquad \deg_t F=S/e,
 \tag{6}
\]

where `F` belongs to `K[t,epsilon]`, is absolutely irreducible and monic in `t`,
and, with the standard ordered Sylvester resultant convention,

\[
             c=(-1)^{S(d+1)}\delta^d\in K^\times.
\]

In particular `e` divides `S`, `deg a`, and `deg H`. There are polynomials
`v,a_0,H_0` over `K` such that `K(a,H)=K(v)`, `deg v=e`,
`a=a_0(v)` and `H=H_0(v)`. Thus `e` is exactly the degree of a maximal common
polynomial right factor, over the original field as well as its algebraic
closure. For an actual positive-exponent frame, `e>1` gives the genuine
ordinary-frame descent already established in the
[ramified-pullback note](SFC_RAMIFIED_PULLBACK_CLOSURE.md).

**Proof.** The identity

\[
 T'=-D'H
 \tag{7}
\]

implies `H=-dT/dD`. Here `D'` and `H'` are nonzero by the degree hypotheses and
characteristic zero. A derivation on a characteristic-zero rational function field
extends uniquely to a finite algebraic extension and preserves intermediate
fields: differentiating a separable minimal polynomial gives its value inside
that extension. Hence `H` belongs to `K(D,T)`, and so does `a=T+DH`. Conversely the
derivation `d/dH` on `K(H)` gives `D=da/dH` in `K(a,H)`, and then `T` belongs
there. This proves the field identity. Equivalently, differentiate a separable
minimal polynomial for `T` over `K(D)`, or for `a` over `K(H)`.

The field degree `[K(z):K(D)]` is `S`. Since `z` is integral over `K[D]`, `T`
is integral too. Its monic minimal polynomial over `K(D)` has coefficients in
the integrally closed ring `K[D]` and degree `S/e`. Its norm from `K(z)` is its
`e`-th power. Substituting
`D=-epsilon` gives the moving resultant up to a nonzero constant. The factor is
constant because `deg a=deg H+S`, so the leading coefficient and degree in `z`
of `a+epsilon H-t` do not change with `epsilon`. Explicitly, using roots `alpha`
of `D+epsilon` over `K(epsilon)`, the ordered resultant is
`(-1)^(dS) delta^d prod_alpha(T(alpha)-t)`, proving the displayed scalar.

The coordinate ring of the image curve is `K[D,T]`, a subring of `K[z]`.
After extension of constants to an algebraic closure this inclusion stays
injective by flatness, and its target remains a domain. Hence its defining
polynomial `F` stays irreducible. In particular the multiplicity `e` and the
reduced degree `S/e` do not conceal a constant-field splitting.

Finally apply Lüroth's theorem over `K` to write the intermediate field as
`K(q)` with `q` rational over `K`; see [L]. Write `a=U(q)`, `H=V(q)` with
rational outer functions over `K`. Since `a` has its only pole at infinity,
every pole of `U` has all its preimages under `q` at infinity. There is exactly
one such value, namely `beta=q(infinity)` in `P^1(K)`. Choose a Möbius map
over `K` sending `beta` to infinity and set `v` to its composition with `q`.
Then `v` has no finite poles, so is a polynomial in `K[z]`. Both transformed
outer functions are polynomials as well: a finite pole would have a finite
preimage under the nonconstant polynomial `v`. Its degree is the field degree
`e`. This proves the right-factor statement over `K`, and degree multiplication
gives divisibility of `deg a` and `deg H`. No scalar period or residue-module
rank assertion enters. QED.

On the normalization of the reduced caustic, one can reconstruct

\[
 H=\frac{dt}{d\varepsilon},\qquad
 a=t-\varepsilon\frac{dt}{d\varepsilon}.
 \tag{8}
\]

Thus the reduced moving caustic keeps the pair up to its genuine common right
factor. If a defining equation is `F(t,epsilon)=0`, the derivative is
`-F_epsilon/F_t` in its function field; this does not require either partial
derivative to be nonzero at the marked point itself.

This is a concrete plane-curve duality certificate, not a novelty claim for
biduality itself. The proposition received separate internal mathematical
reviews by the root and factorial-review agents. The latter also reconstructed
the arbitrary-field scope, absolute irreducibility, exact norm scalar,
polynomial generator over `K`, and the split grouped-descent corollary below.
These are not external or formal verification claims.

### Corollary: all common right factors give ordinary-frame descent

This corollary is the retained pullback note's Section 13, now coupled with the
caustic detector; it is not a newly claimed descent theorem. Assume the
normalized actual frame of Section 1. Normalize its maximal factor over `K`
so that `v(0)=0` and `v'(0)=1`. This is possible because
`1=a'(0)=a_0'(v(0))v'(0)`. Then also `a_0(0)=H_0(0)=0` and
`a_0'(0)=H_0'(0)=1`.

Let `beta` be a root of `H_0'`, and let `alpha` above it have local degree `l`
under `v`. Writing `h_beta=ord_beta H_0'` gives

\[
 \operatorname{ord}_\alpha H'=l h_\beta+l-1.
\]

It is positive, so actual root alignment forces `a(alpha)=0`; hence
`a_0(beta)=0`. If `R_beta=ord_beta a_0`, the original inequality
`ord_alpha H'<=ord_alpha a-1` yields

\[
             l h_\beta+l-1\le l R_\beta-1,
             \qquad h_\beta\le R_\beta-1.
\]

Thus every nonzero root of `a_0` has
`p_beta=R_beta-ord_beta H_0'>=1`, and `R_beta>=p_beta`. The zero at the marked
origin is simple. At each other root, assign `p_beta` variables the shared
coefficient `1/beta`, allocating `R_beta` as a sum of `p_beta` positive
integers. This constructs a grouped ordinary factorial frame, with

\[
             S_0=\sum_{\beta\ne0}p_\beta
                 =\deg a_0-\deg H_0=S/e.
\]

For an actual frame whose coefficients lie in `K`, the seed coefficients also
lie in `K`: every root of `a` lies in `K`, every root `beta` of `a_0` is the
image under `v` of such a root, and `v` has coefficients in `K`. For a general
unsplit polynomial pair the grouped interpretation is over a splitting field.

Finally the marked inverse gives exactly

\[
 H(a^{-1}(t))=H_0(a_0^{-1}(t)).
\]

All normalized moments are preserved. If `e>1`, actual support strictly drops;
with a maximal factor the seed pair has no further common right factor. No
cycle-only restriction is present. This is why the exponent in (6) is useful
for the full pullback reduction rather than only its cyclic examples.

## 4. Burgers inversion and the precise Jacobian/tree bridge

Let `z_epsilon(t)` be the inverse of the pencil at zero, in
`C[[epsilon,t]]`, and set `Y(t,epsilon)=H(z_epsilon(t))`. Direct differentiation
gives

\[
 Y_t=\frac1{D(z_\varepsilon)+\varepsilon},\qquad
 Y_\varepsilon=-\frac{Y}{D(z_\varepsilon)+\varepsilon},\qquad
 \boxed{Y_\varepsilon+Y Y_t=0}.
 \tag{9}
\]

Equivalently `Y=B(t-epsilon Y)`. Formal Lagrange inversion gives

\[
 Y(t,\varepsilon)=\sum_{k\ge0}
 \frac{(-\varepsilon)^k}{(k+1)!}\,
 \partial_t^k\bigl(B(t)^{k+1}\bigr).
 \tag{10}
\]

The deformation is algebraic and retains the distinguished branch. Formula
(9) alone is universal for initial formal series and supplies no gap bound.

There is an exact normalization into Zhao's formal inverse equation. Put
`C(t)=B(t)-t`, which has order at least two, and

\[
 x=\frac{t}{1+\varepsilon},\qquad
 \eta=-\frac{\varepsilon}{1+\varepsilon}.
\]

If the inverse of `u-eta C(u)` is `G_eta(x)=x+eta N_eta(x)`, then

\[
 Y(t,\varepsilon)
 =\frac{t+N_\eta(t/(1+\varepsilon))}{1+\varepsilon}.
 \tag{11}
\]

Indeed `u+epsilon B(u)=(1+epsilon)(u-eta C(u))` and
`N_eta(x)=C(G_eta(x))`. In one variable the gradient condition is automatic.

Zhao's Theorem 2.4 identifies formal inverse deformation with
`partial_eta N=N partial_x N`; Theorem 5.2 expands its potential over binary
rooted trees. The source allows formal series over any Q-algebra. Its use here
is legitimate for `ord C>=2`. However, `C` is generally algebraic and not
polynomial. Neither polynomial termination nor Jacobian/Hessian nilpotence
follows. The useful import is the exact deformation and tree recursion, not a
Jacobian-conjecture theorem. See [Z], PDF pp.7–8 and 19–20.

## 5. Exact positive and negative controls

The companion checker verifies four pencils over the rationals, resultant
factorization, the dual identity, and a bounded formal Lagrange/Burgers replay.
These checks verify the stated examples, not the proposed closing lemmas.

### 5.1 Primitive positive frame

For `r=(1,1,2)` and `c=(1,2,-2)`,

\[
 a=8z^5-4z^4-6z^3+z^2+z,\quad H=z+z^2,\quad
 D=20z^3-18z^2+1.
\]

The fixed factor is `8(epsilon+4t)`, and the moving factor has `t`-degree three.
The pair is primitive because a common factor degree divides both `5` and `2`.
Therefore (6) proves that its moving factor is geometrically irreducible.
At the marked point `D'(0)=0`: the parameter expansions are
`epsilon+1=18z^2-20z^3` and `t=12z^3-6z^4-12z^5`. Thus the moving curve has a
cusp there despite being primitive. A singular caustic is not by itself a
composition certificate. The exact polynomial is retained in the receipt.

### 5.2 Genuine descent, not an arbitrary rank drop

For the retained `t=2` descent example,

\[
 a=z(1-z)(1-2z)^2,\quad H=z-z^2,\quad
 D=1-8H,\quad T=4H^2,
\]

the moving resultant is exactly

\[
 16\bigl(16t-(\varepsilon+1)^2\bigr)^2.
\]

The exponent two agrees with the common factor `H`; the reduced caustic is the
parabola `16t=(epsilon+1)^2`. This gives a nonlinear explanation of the actual
descent without treating all residue-determinant zeros as compositions.

### 5.3 Moving degree and Burgers are insufficient

The retained misaligned example is

\[
 H=z-z^2/2,\quad
 D=1+3z/5+6z^2/25,\quad
 a=z-z^2/5-3z^3/25-3z^4/50.
\]

It still satisfies `a'=H'D` and Burgers. The moving factor has `t`-degree two.
But with `b(u)=u-3u^2/10`,

\[
 H-b(a)=-27z^5/1250+O(z^6),
\]

so the putative `n+S=4` bound fails. Its fixed critical-value line is
`t=31/50+epsilon/2`, because `a(1)=31/50`. It misses the concurrence point.
This is a required negative control for any degree-only caustic argument.

The additional boundary control `a=z-z^2,H=z` has `W=1` and `S=1`.
There are no fixed critical lines; its moving curve is linear in `t`.

## 6. Two precise closing targets and what is actually missing

### Target A: a primitive pointed-pencil finite-jet lemma

For every normalized polynomial pair with `H'` dividing `gcd(a,a')`,
`deg a>deg H>=1`, and `K(a,H)=K(z)`, prove (1) for every `n>=1` and every
polynomial `b` of degree at most `n`.

Equivalently, in the general aligned class prove that a failure of (1) forces
the moving resultant in (4) to be a proper power. Proposition (6) would then
turn this into actual lower-support descent. The smallest-support hypothetical
counterexample is primitive, so that implication would close the entire
positive-exponent monomial-linear class.

This is still the sharp unresolved statement in geometric form, not an
additional proved reduction of its difficulty. The new leverage is that its
data comprise a degree-`S` moving covering, concurrent fixed critical-value
lines, and a marked normalization point. A viable proof must use all three.
One concrete attack is a local intersection estimate on the normalized
caustic that charges contact exceeding `n+S` to nontrivial covering
multiplicity. Ordinary Bezout or the degree in `t` alone cannot do that.

The tree version should work with the exact deformation (10) or (11), retaining
the aligned pencil relation. Expanding an arbitrary initial series into trees
does not produce the needed estimate. Complex coefficients also prevent an
unjustified termwise positivity argument. A successful tree identity must
prove the contact-to-composition implication uniformly in `n,S,r_i,c_i`.

### Target B: finite contact to a persistent zero-cycle center

Over the complex numbers, an alternative sufficient bridge is the following
unproved statement. From a
failure of (1) in an actual aligned pair, construct a nonzero integral
zero-cycle `C` of `a` whose full displacement under `a+epsilon H` vanishes
identically for every regular `t` and all sufficiently small `epsilon`.

The existing zero-cycle center theorem [ABC-M], Theorem A, would then imply
`a=a_0 o v`, `H=H_0 o v` with a common right factor, and the pushforward cycle
trivial as a chain. Since the cycle is nonzero, that factor has degree greater
than one. Retained ordinary-frame descent again reduces support.

The missing step is substantial: a finite Taylor gap does not presently give
such a cycle or an identity in two variables. First-order/tangential vanishing
is not enough; the same source explicitly distinguishes that failure. This
target should be pursued only if the pencil/tree expansion produces an
actual integral cycle and proves all deformation coefficients vanish. An
arbitrary complex linear functional on residue coordinates is not that cycle.

## 7. Logarithmic and twisted-period guardrails

For distinct `c_i`, set `alpha_i=1/c_i` and

\[
 \beta_i=\operatorname{Res}_{z=\alpha_i}\frac{dz}{zP(z)}.
\]

If `ell_i` is the Lagrange cardinal polynomial at these `S` nodes, then

\[
 -\beta_i=\ell_i(0),\qquad
 -\sum_i\beta_i f(\alpha_i)=f(0)\quad(\deg f<S).
 \tag{12}
\]

To check it, write `P(z)=prod(z-alpha_j)/prod(-alpha_j)` and compare
`1/(alpha_i P'(alpha_i))` with the cardinal product. Thus the retained
logarithmic residues `beta_i-r_i` combine integer multiplicities with a
specific evaluation functional at the marked point. They are not arbitrary
positive electrostatic charges. Formula (12) is a useful constraint to retain
in any Padé or Heine–Stieltjes formulation.

At integral moment order the twist `a^{-m}` is rational and its local
monodromy is trivial. Generic nonresonant twisted-period results therefore
cannot simply identify the retained residue observable. Critical-algebra
nondegeneracy is also different from observability of a chosen scalar period.
In particular, do not conflate an arrangement's unbalanced weights with
nonresonance of every local-system exponent.

The [retained bad scalar periods](SFC_MONOMIAL_LINEAR_CONTACT_BOUND.md#6-exact-obstruction-to-consecutive-scalar-observability)
show why branch information matters. They are
linear combinations of conjugates of `rho` and can have forbidden zero windows
that the actual branch does not. The separate
[differential-branch obstruction note](SFC_DIFFERENTIAL_BRANCH_OBSTRUCTIONS.md)
additionally obstructs a universal narrow homogeneous annihilator whose only
polynomial solutions are constants. A linear route would need a different
additional branch constraint;
normalizing its constant coefficient to one is insufficient.

## 8. Primary-source checks and limits

**[Z]** W. Zhao, *Inversion Problem, Legendre Transform and Inviscid Burgers'
Equations*, arXiv:math/0403020v2 (2004), Theorems 2.4 and 5.2.
<https://arxiv.org/pdf/math/0403020>.
The exact formal inverse/PDE and binary-tree statements are used in Section 4.
Their polynomial/nilpotence consequences are not imported.

**[ABC-M]** A. Álvarez, J. L. Bravo, C. Christopher, P. Mardešić,
*Infinitesimal Center Problem on Zero Cycles and the Composition Conjecture*,
arXiv:2006.07600, Theorem A and Definition 1.4, PDF pp.3–4.
<https://arxiv.org/pdf/2006.07600>.
The theorem concerns an identically zero full displacement and an integral
zero-cycle, not a finite local jet. Its composition condition includes trivial
cycle pushforward. It is used only for the conditional implication in Target B.

**[BMNP]** J. L. Bravo, P. Mardešić, D. Novikov, J. Pontigo-Herrera,
*Infinitesimal and tangential 16-th Hilbert problem on zero-cycles*,
Bulletin des Sciences Mathématiques 202 (2025), 103634.
<https://doi.org/10.1016/j.bulsci.2025.103634>;
inspected primary preprint <https://arxiv.org/pdf/2312.03081>.
Theorems A/B give degree-dependent maxima for isolated regular zeros of
zero-cycle integrals/displacements. For degrees `d>h`, `d>2`, the general
maximum is `h(d-1)!`. Theorem D disproves a broad asymptotic Chebyshev
expectation, even for simple cycles. These results neither count the desired
marked Taylor gap at a possibly globally critical value nor give an `S`-bound.

**[GP]** L. Gavrilov, F. Pakovich, *Moments on Riemann surfaces and
hyperelliptic Abelian integrals*, Comment. Math. Helv. 89 (2014), 125–155,
Theorems 2.1–2.2. <https://www.math.bgu.ac.il/~pakovich/Publications/abelian.pdf>.
The all-moment/zero-cycle reduction requires every moment to vanish; its
stronger composition conclusion has monodromy hypotheses. It does not bridge
a finite shifted block to an identity.

**[DLMF]** Section 31.15, Heine–Stieltjes polynomials.
<https://dlmf.nist.gov/31.15>.
The real-root/electrostatic conclusions use positive residues and ordered real
singular points. The complex coefficients and barycentric residues in (12)
do not satisfy those hypotheses uniformly.

**[L]** The standard arbitrary-field Lüroth theorem is also stated precisely
in the official Mathlib documentation, `RatFunc.Luroth.eq_adjoin_generator`:
<https://leanprover-community.github.io/mathlib4_docs/Mathlib/FieldTheory/RatFunc/Luroth.html>.
It assumes only that `K` is a field. The passage from its rational generator
to a polynomial generator over `K` is proved above. This reference is not a
claim that this note or its caustic proposition has been formalized.

## Replay

From the repository root:

```text
python research/factorial/image-rigidity/verify_sfc_pointed_pencil.py
```

The checker uses exact rational SymPy arithmetic. It rejects disabled assertions,
an altered or missing note, and altered or missing pinned retained sources before
writing a receipt. It is a bounded example replay, not a proof engine for either
proposed closing target.
