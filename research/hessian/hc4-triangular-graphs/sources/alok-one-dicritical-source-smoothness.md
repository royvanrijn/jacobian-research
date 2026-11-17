# Smooth finite normalization in the one-dicritical degree-six passports

## Claim boundary

Assume that a hypothetical topological-degree-six plane Keller map has one
dicritical component and one of the two surviving passports

\[
(G,e,d)=(A_6,3,1)
\qquad\text{or}\qquad
(G,e,d)=(S_6,2,1).
\]

Orevkov's minimal finite-boundary model is a linear chain

\[
L_C-E,
\]

where `E` is the unique dicritical and the connected, possibly empty, chain
`L_C` is constant.  In Section 5 Orevkov makes this model relatively minimal:
every point produced by contracting `L_C` may be assumed singular, since a
smooth one would make the corresponding blowups unnecessary.  This note
proves the strict strengthening, in that minimal model,

\[
\boxed{L_C=\varnothing.} \tag{0.1}
\]

Consequently the affine finite normalization is smooth in both passports.
For `A6` this eliminates contracted source points of local degree three as
well as the previously excluded degree-five point.  For `S6` it removes the
only caveat in the local knot classification: every positive jump is exactly
`T(2,3)`, `T(3,4)`, or `T(4,5)`.

This does **not** eliminate either one-dicritical passport, generic degree six,
or `JC(2)`.  The proof is conditional on the established one-dicritical
finite-normalization setup.  In particular, the endpoint-curvette and
sole-local-ramification hypotheses below are essential.

## 1. The cyclic endpoint obstruction

Let

\[
\rho:(X,x)\longrightarrow(\mathbb C^2,0)
\]

be a finite germ of local degree `mu`.  Assume:

1. `(X,x)` is the nonsmooth normal germ obtained by minimally contracting a
   nonempty negative-definite rational linear chain;
2. the ramification prime `D` is a curvette at an endpoint of that chain;
3. `D` is the only ramification prime of `rho` through `x`; and
4. its generic normal ramification index is `e`.

The local theorem is

\[
\boxed{
e=2\ \Longrightarrow\ \text{no such germ exists},
\qquad
e=3\ \Longrightarrow\ \mu\text{ is even}.
}
\tag{1.1}
\]

### 1.1 The Hirzebruch--Jung cover

The Hirzebruch--Jung classification makes `(X,x)` a nontrivial cyclic
quotient singularity.  Let

\[
q:(\mathbb C^2,0)\longrightarrow(X,x)
\tag{1.2}
\]

be its universal quasi-etale cyclic cover, of order `n>1`.  The endpoint
condition matters here.  In the toric model, an endpoint curvette lifts to a
coordinate axis.  Equivalently, the endpoint knot generates the fundamental
group of the lens-space link, so its inverse image in the universal `S3` link
is connected.  Thus the reduced inverse image of `D` is one smooth branch.
Straighten it equivariantly to

\[
q^{-1}(D)_{\mathrm{red}}=\{x=0\}. \tag{1.3}
\]

Put

\[
h=\rho\circ q=(u,v).
\]

The local degree multiplies:

\[
\deg_0h=\mu n. \tag{1.4}
\]

The cover `q` is etale in codimension one.  The different of `rho` has
coefficient `e-1` along `D`, and by hypothesis it has no other local prime.
Pulling the different back through `q`, rather than writing a coordinate
Jacobian on the singular germ `X`, gives

\[
J_h=\varepsilon x^{e-1} \tag{1.5}
\]

for a unit `epsilon`.

The cyclic deck representation is small.  If it had an invariant linear form,
its generator would have eigenvalue one and would be a pseudoreflection,
contrary to quasi-etaleness.  The components `u,v` are deck invariant, hence

\[
du_0=dv_0=0. \tag{1.6}
\]

### 1.2 Index two is impossible

If `e=2`, equation (1.5) has order one.  Equation (1.6) gives

\[
\operatorname{ord}u,\operatorname{ord}v\ge2,
\qquad
\operatorname{ord}J_h\ge
\operatorname{ord}u+\operatorname{ord}v-2\ge2,
\]

a contradiction.  This proves the first half of (1.1), independently of
`mu` and `n`.

### 1.3 Index three forces the `A1` quotient

Now let `e=3`.  The Jacobian in (1.5) has exact order two.  The order bound
forces both `u` and `v` to have nonzero quadratic initial forms `U,V`, with

\[
J(U,V)=c x^2\ne0. \tag{1.7}
\]

The quadratics cannot be proportional.  If they were coprime, they would
define a degree-two map

\[
[U:V]:\mathbb P^1\longrightarrow\mathbb P^1
\]

whose entire degree-two ramification divisor was concentrated at one point.
That is impossible: every ramification index of a degree-two map is at most
two, so every local contribution is at most one, while Riemann--Hurwitz gives
total contribution two.  Thus `U,V` have a common linear factor.  After
linear source and target changes,

\[
U=x^2,\qquad V=xy. \tag{1.8}
\]

Deck invariance now determines the quotient.  If `g^*x=chi_g x`, write
`U=xA` and `V=xB`.  The independent residual forms `A,B` both have character
`chi_g^{-1}`, so `g` acts on the whole cotangent space by that scalar.  But
`x` lies in the cotangent space and also has character `chi_g`; hence

\[
\chi_g^2=1.
\]

Faithfulness and `n>1` give

\[
\boxed{n=2,\qquad g=-I.} \tag{1.9}
\]

In particular, `u` and `v` are even power series and

\[
\deg_0h=2\mu. \tag{1.10}
\]

### 1.4 The intersection-parity endpoint

The tangent cone `v_2=xy` is reduced, so `v=0` has exactly two smooth
transverse branches.  Parametrize the branch tangent to `y=0` by

\[
\gamma_y(t)=(t,\alpha(t)).
\]

Along it, `v_y=t+O(t^2)`.  Differentiating
`v\circ\gamma_y=0` and using (1.5) gives

\[
\frac d{dt}(u\circ\gamma_y)
=\frac{J_h\circ\gamma_y}{v_y\circ\gamma_y}.
\]

The numerator has order two and the denominator order one, so this branch
contributes exactly two to `I_0(u,v)`.

Write the branch tangent to `x=0` as

\[
\gamma_x(t)=(\beta(t),t),
\qquad r=\operatorname{ord}\beta.
\]

Here `beta` is not identically zero.  Otherwise `x=0` would be a component of
`v=0`, the same chain-rule identity would make `u` constant on it, and the
finite map `h` would contract a curve.  Since
`v_x\circ\gamma_x=t+O(t^2)`,

\[
\frac d{dt}(u\circ\gamma_x)
=-\frac{J_h\circ\gamma_x}{v_x\circ\gamma_x}.
\]

This branch contributes exactly `2r`.  Equations (1.10) and the additivity of
intersection multiplicity give

\[
2\mu=I_0(u,v)=2+2r,
\qquad r=\mu-1. \tag{1.11}
\]

The involution `-I` preserves the unique branch tangent to `x=0`, so its graph
obeys

\[
\beta(-t)=-\beta(t).
\]

Therefore `beta` is odd and its finite order `r` is odd.  Equation (1.11)
forces `mu` to be even, proving the second half of (1.1).

## 2. Application to the two surviving passports

In the one-dicritical `A6` passport, `e=3`.  The exact jump theorem gives one
local-degree-five point and local degree three everywhere else on `D`.  If
`L_C` were nonempty, its contraction point would therefore have

\[
\mu\in\{3,5\},
\]

contradicting the evenness conclusion in (1.1).

In the one-dicritical `S6` passport, `e=2`.  Any nonempty contracted chain is
already impossible by the order-one Jacobian contradiction, without needing
its local degree.

Thus in both cases

\[
\boxed{L_C=\varnothing.} \tag{2.1}
\]

The application uses the relative-minimality statement just recalled.  If
`L_C` were nonempty, its contraction point would be nonsmooth and the local
theorem would apply.  Thus the conclusion is not a claim that an arbitrary
nonminimal resolution contains no artificially inserted constant curves.

Orevkov's minimal contraction is consequently an isomorphism near every finite point
of `E`.  The finite normalization is smooth along the affine dicritical
`D=E^\circ`, as well as on the original affine source.  In particular:

- the `A6` jump point is smooth, rank one, and has target knot `T(2,5)`;
- every `S6` jump-one point has knot `T(2,3)`;
- every `S6` jump-two point has knot `T(3,4)`; and
- every `S6` jump-three point has knot `T(4,5)`.

## 3. The resulting finite-flat surface

Let `W` be the normalization of the affine target in `C(x,y)`, and let

\[
\rho:W\longrightarrow\mathbb A^2
\]

be the finite normalization map.  The established one-dicritical boundary
identification and (2.1) now give

\[
W\setminus D\cong\mathbb A^2,
\qquad
D=E^\circ\cong\mathbb A^1,
\qquad
W\text{ smooth affine}. \tag{3.1}
\]

This has several exact global consequences.

1. **Finite flat rank six.**  The smooth surface `W` is Cohen--Macaulay, the
   target is regular, and a finite morphism has zero-dimensional fibers.
   Miracle flatness makes `rho` finite flat of rank six.

2. **A free algebra module.**  The module `rho_* O_W` is finite locally free
   of rank six over `C[P,Q]`.  Quillen--Suslin makes its underlying module
   free.  This does not trivialize its algebra structure.

3. **The Picard group.**  Divisor localization and `Pic(A2)=0` show that
   `Pic(W)=Cl(W)` is generated by `[D]`.  If `nD=div(f)`, then `f` restricts
   to a unit on `W-D=A2`, hence to a nonzero constant; equality in the common
   function field makes `f` constant and `n=0`.  Therefore

   \[
   \boxed{\operatorname{Pic}(W)\cong\mathbb Z[D].} \tag{3.2}
   \]

4. **Canonical and discriminant divisors.**  The original affine map is
   etale and `D` is the only boundary prime, so `D` is the whole ramification
   divisor support.  Fixing the target form `dP\wedge dQ`, finite-map
   Riemann--Hurwitz gives

   \[
   \operatorname{div}\rho^*(dP\wedge dQ)=(e-1)D,
   \qquad
   \omega_W\cong\mathcal O_W((e-1)D). \tag{3.3}
   \]

   Since `d=1`, the branch discriminant has generic multiplicity `e-1`.

These conclusions concern the affine finite normalization only.  They do not
say that the compactification is smooth over target infinity, and they do not
make `W` isomorphic to `A2`; indeed (3.2) rules out that isomorphism.  A free
rank-six pushforward is also not a split cover.  The remaining obstruction
must use the multiplication/trace form, braid data, or the compactification
graph, not module freeness alone.

## 4. Exact executable endpoint

The finite certificate
  [`OneDicriticalContractedSourceEndpoint`](../scripts/six_sheet_monodromy.py)
records only the arithmetic endpoint of the geometric proof:

- for `e=2`, Jacobian order one is below the invariant lower bound two;
- for `e=3`, the forced double cover gives contact `r=mu-1`, which conflicts
  with deck-oddness exactly when `mu` is odd.

The tests enumerate `mu=3,5` for `A6` and `mu=2,3,4,5` for `S6`.  They do not
pretend to certify the Hirzebruch--Jung classification, connected endpoint
lift, or different formula; those are the mathematical premises above.

## Sources

- S. Yu. Orevkov,
  [“On three-sheeted polynomial mappings of `C^2`”](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
  *Mathematics of the USSR-Izvestiya* 29 (1987), 587--596.  Lemma 2.1 gives
  the linear endpoint chain, Lemma 3.1 the generic local ramification form,
  Lemma 4.2 the defect identity, Section 5 the relative-minimality statement
  for points in `pi(L_C)`, and Lemma 5.2 the zero-jump embedding.
- Egbert Brieskorn,
  [“Rationale Singularitaten komplexer Flachen”](https://doi.org/10.1007/BF01425318),
  *Inventiones Mathematicae* 4 (1968), 336--358.  The
  Hirzebruch--Jung classification identifies a rational linear-chain
  singularity with a cyclic quotient; the endpoint lift is its standard
  toric curvette.
- [Stacks Project, Lemma 10.128.1](https://stacks.math.columbia.edu/tag/00R4)
  is the miracle-flatness criterion used in Section 3.

The cyclic-character reduction, index-two obstruction, intersection-parity
argument, and finite-flat consequences are derived in this repository.  No
claim of historical priority is made.
