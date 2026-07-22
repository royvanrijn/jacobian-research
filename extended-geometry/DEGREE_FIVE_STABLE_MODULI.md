# Infinitely many stable classes in generic degree five

This note gives an explicit positive-dimensional stable-moduli family for the
weighted marked-root construction.  Generic affine-mark faithfulness proves
exact parameter recovery on the working open.  The earlier
Hessian/Fitting-divisor argument is retained
because it is an independent, weaker fallback using only the ramification
divisor of the discriminant normalization and affine three-point geometry.

The resulting invariant is reusable: for split weighted seeds, stable
polynomial left--right equivalence forces the root divisor of `H''` on the
normalization line to be affinely equivalent.  In degree five this is an
unordered triple of points, hence carries one cross-ratio modulus.  In degree
four it is only an unordered pair, which explains why this invariant cannot
separate the two split quartic islands.

## 1. The explicit family

Put

\[
 H_\lambda(W)=
 \frac{W^2(W-1)\bigl(3W^2-(5\lambda+1)W+3\lambda\bigr)}{60}
\]

and

\[
 c_\lambda=-H_\lambda'(1)=\frac{\lambda-1}{30}.
\]

The decisive identity is

\[
 H_\lambda''(W)=
 \frac{(W-\lambda)(10W^2-8W+1)}{10}.                 \tag{1.1}
\]

Let

\[
 p=\frac{4-\sqrt6}{10},\qquad
 q=\frac{4+\sqrt6}{10}.                              \tag{1.2}
\]

Thus the three roots of `H_lambda''` are exactly

\[
 \{p,q,\lambda\}.                                    \tag{1.3}
\]

Define the explicit exceptional polynomial

\[
 \mathcal E(\lambda)=
 \lambda(\lambda-1)(25\lambda-1)
 (10\lambda^2-8\lambda+1)
 (100\lambda^2-29\lambda+10)                         \tag{1.4}
\]

and work on the nonempty Zariski-open parameter set

\[
 \Lambda=\{\lambda\in\mathbb C:\mathcal E(\lambda)\ne0\}.
\]

For `lambda in Lambda`, put

\[
 v=xy,\qquad S=x^2z,\qquad u=1+v,
\]

\[
 \gamma=1-\frac87v+S,\qquad W=u\gamma,               \tag{1.5}
\]

and write `p_lambda=H_lambda'` and

\[
 q_\lambda(W)=\frac{Wp_\lambda(W)-H_\lambda(W)}{c_\lambda}.
\]

Define

\[
 F_\lambda=(A_\lambda,B_\lambda,C_\lambda)
 :\mathbb A^3\longrightarrow\mathbb A^3
\]

by

\[
 C_\lambda=x\gamma,
\]

\[
 B_\lambda=
 \frac{c_\lambda+p_\lambda(W)/\gamma}{x},
\]

\[
 A_\lambda=
 \frac{u+q_\lambda(W)/\gamma^2}{x^2}.                 \tag{1.6}
\]

The displayed quotients are polynomials in `x,y,z`; their coefficients are
regular functions of `lambda` on `Lambda`.

The family varies only through `A_lambda` and `B_lambda`.  Its weighted
coordinate and third output are fixed:

\[
 \gamma=1-\frac87xy+x^2z,qquad C_\lambda=x\gamma.
\]

Equivalently, the weighted parameter is constantly `a_0=-8/7`.  This fixed
component is the input that permits the common rank-two adapted-coordinate
system constructed in the
[degree-five symplectic descent](DEGREE_FIVE_RANK_TWO_DESCENT.md).

## 2. Main theorem

### Theorem 2.1 — exact degree-five parameter recovery

For every `lambda in Lambda`, the map `F_lambda` is a polynomial Keller map
with

\[
 \det DF_\lambda=c_\lambda=\frac{\lambda-1}{30}\ne0, \tag{2.1}
\]

generic fiber degree five, and a fiber containing five distinct points.

If `F_lambda` and `F_mu` are polynomially left--right equivalent after
adjoining any number of identity variables, then

\[
 \boxed{\lambda=\mu.}                               \tag{2.2}
\]

Consequently the maps `F_lambda`, `lambda in Lambda`, are an uncountable
family of pairwise stably polynomially left--right inequivalent Keller maps
of generic degree five.

For the independent coarse invariant, define

\[
 \chi(\lambda)=\frac{\lambda-p}{q-p}                 \tag{2.3}
\]

and the six-element anharmonic orbit

\[
 \mathfrak A(x)=
 \left\{
 x,\ 1-x,\ \frac1x,\ \frac1{1-x},\
 \frac{x}{x-1},\ \frac{x-1}{x}
 \right\}.                                           \tag{2.4}
\]

The Hessian/Fitting-divisor argument alone gives

\[
 \chi(\mu)\in\mathfrak A(\chi(\lambda)).             \tag{2.5}
\]

Equivalently, `mu` belongs to the explicitly described set

\[
 \mathfrak O(\lambda)=
 \left\{p+(q-p)a:a\in\mathfrak A(\chi(\lambda))\right\}. \tag{2.6}
\]

Thus it bounds a stable class by six parameters even without invoking the
stronger full-cover theorem.  It is not the sharp separation statement.

The orbit obstruction can also be written without `sqrt(6)`.  Put

\[
 J(\lambda)=
 \frac{64(50\lambda^2-40\lambda+17)^3}
 {75(10\lambda^2-8\lambda+1)^2}.                     \tag{2.7}
\]

Then stable equivalence implies

\[
 J(\lambda)=J(\mu).                                  \tag{2.8}
\]

The rational map `J:P^1->P^1` has generic degree six, and its fibers are the
orbits (2.5).

## 3. Weighted admissibility and a five-point fiber

Direct differentiation gives

\[
 H_\lambda(0)=H_\lambda'(0)=H_\lambda(1)=0,
\]

\[
 H_\lambda'(1)=-\frac{\lambda-1}{30},\qquad
 \frac{H_\lambda''(1)}{c_\lambda}=-9.                \tag{3.1}
\]

Thus the weighted construction parameter is

\[
 a_0=-\frac{1+(-9)}{2+(-9)}=-\frac87,
\]

independent of `lambda`.  The weighted polynomiality and triangular-Jacobian
calculation therefore apply verbatim and give (2.1).  The marked inverse
equation is

\[
 E_{\lambda;A,B,C}(W)
 =H_\lambda(W)-BCW+c_\lambda AC^2=0.                 \tag{3.2}
\]

Its degree in `W` is five.

For completeness, fix `lambda` and choose a complex number `alpha`
transcendental over `overline(Q)(lambda)`.  At

\[
 (A,B,C)=(\alpha,0,1)
\]

the inverse equation is

\[
 H_\lambda(W)+c_\lambda\alpha=0.                    \tag{3.3}
\]

If a root `w` of (3.3) also satisfied `H_lambda'(w)=0`, then `w` and
`H_lambda(w)` would be algebraic over `overline(Q)(lambda)`, contradicting
`H_lambda(w)=-c_lambda alpha`.  Hence (3.3) has five distinct roots and all
reconstruction denominators are nonzero.  Every root reconstructs a distinct
source point.  This proves both exact generic degree five and the existence of
a five-point fiber.

## 4. The two canonical boundary vertices

Write

\[
 R_\lambda(W)=3W^2-(5\lambda+1)W+3\lambda.
\]

The exclusions in (1.4) give

\[
 \operatorname{Disc}_W(R_\lambda)
 =(25\lambda-1)(\lambda-1)\ne0,
\]

\[
 R_\lambda(0)=3\lambda\ne0,
 \qquad R_\lambda(1)=2(1-\lambda)\ne0.               \tag{4.1}
\]

Thus `H_lambda` has a double root at zero, the distinguished simple root one,
and two further distinct simple roots.

Over the target divisor `C=0`, the zero cluster has scaled equation

\[
 -\frac{\lambda}{20}u^2-Bu+c_\lambda A+O(C)=0.       \tag{4.2}
\]

It supplies the affine residue-degree-two branch.  The root at one supplies the
other affine branch.  For an additional simple root `rho`, reconstruction has

\[
 \gamma\longrightarrow-\frac{H_\lambda'(\rho)}{c_\lambda},
 \qquad
 y=\frac{W-\gamma}{C}.
\]

The branch is visibly boundary when

\[
 c_\lambda\rho+H_\lambda'(\rho)\ne0.                \tag{4.3}
\]

Up to a nonzero constant, the resultant of (4.3) with `R_lambda(rho)=0` is

\[
 \lambda(\lambda-1)^2(100\lambda^2-29\lambda+10),   \tag{4.4}
\]

which is nonzero on `Lambda`.  Hence both additional-root branches are
unramified boundary primes over `C=0`.  The local degree sum is

\[
 2+1+1+1=5.
\]

Over the reduced discriminant there is one ramified boundary prime with
`(e,f)=(2,1)` and three affine simple branches, again exhausting degree five.
The reconstruction denominators show that no third divisorial boundary image
can occur.  Consequently the two target boundary vertices are exactly

\[
 Z_\Delta=\text{the unique vertex receiving ramification},
 \qquad Z_0=V(C).                                    \tag{4.5}
\]

They are intrinsic and cannot be exchanged by a left--right equivalence.

As a regression, exact discriminant saturation gives

\[
 \left.C^{-2}\operatorname{Disc}_W(E_{\lambda;A,B,C})
 \right|_{C=0}
 =u_\lambda\bigl(\lambda(\lambda-1)A+150B^2\bigr),   \tag{4.6}
\]

with `u_lambda != 0` on `Lambda`.  Thus the target boundary intersection is
reduced.  This numerical fact is not the modulus used below.

## 5. The intrinsic Hessian-root divisor

Let

\[
 D_H\subset\mathbb A^2_{s,t}
\]

be the reduced discriminant curve of

\[
 H(W)-sW+t.
\]

Its normalization is

\[
 \nu_H:\mathbb A^1_r\longrightarrow D_H,
 \qquad
 r\longmapsto\bigl(H'(r),rH'(r)-H(r)\bigr).          \tag{5.1}
\]

Indeed, `r` is integral over the discriminant coordinate ring because it
satisfies `H'(r)-s=0`, and birationality follows from

\[
 \frac{d(rH'(r)-H(r))}{dH'(r)}=r
\]

in the function field.

Let

\[
 R_H=\mathbb C[H'(r),rH'(r)-H(r)]\subset\mathbb C[r].
\]

The relative differential module of the normalization is

\[
 \Omega_{\mathbb C[r]/R_H}
 \simeq
 \frac{\mathbb C[r]}{(H''(r))}\,dr.                 \tag{5.2}
\]

To see this, the image of the two generators of `R_H` in
`Omega_{C[r]/C}` is generated by

\[
 dH'(r)=H''(r)dr,
\]

\[
 d(rH'(r)-H(r))=rH''(r)dr.
\]

Their generated ideal is exactly `(H'')`.  Therefore the ramification divisor
of the normalization morphism is intrinsically the zero divisor of `H''` on
the normalization line.

### Proposition 5.1 — scheme-theoretic stable Hessian-root invariant

Let `F_H` and `F_G` be split weighted marked-root maps for which the complete
canonical target boundary consists of the intrinsically ordered pair
`(Z_Delta,Z_0)`.  If `F_H` and `F_G` are stably polynomially left--right
equivalent, then their full zeroth-Fitting divisors are isomorphic after
stabilization.  In particular every Hessian-root multiplicity is preserved.
If each divisor has at least two distinct support points, the effective
divisors `div(H'')` and `div(G'')` are affinely equivalent on `A^1`.

#### Proof

The normalization, relative-differential, and Fitting-ideal base-change steps
are instances of the general
[stable normalization functoriality theorem](../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md).
This proof only identifies the intrinsic degree-five boundary vertex and its
Hessian divisor.

Remove the intrinsic second boundary vertex.  On `C!=0`, the coordinate change

\[
 (A,B,C)\longmapsto(s,t,C)=(BC,cAC^2,C)              \tag{5.3}
\]

is an isomorphism, and it identifies

\[
 Z_\Delta\setminus Z_0\simeq D_H\times\mathbb G_m. \tag{5.4}
\]

After adjoining `s` identity variables, the general theorem identifies its normalization as

\[
 \mathbb A^1_r\times\mathbb G_m\times\mathbb A^s,
\]

with coordinate ring

\[
 B_H=\mathbb C[r,\xi,\xi^{-1},z_1,\ldots,z_s].       \tag{5.5}
\]

Stable boundary functoriality identifies the two intrinsic vertices, restricts
to an isomorphism of the opens (5.4), and lifts uniquely to an isomorphism of
the normalizations.  The lifted isomorphism preserves the relative
differential module and therefore its zeroth Fitting ideal.  Fitting ideals
commute with the flat base change used in stabilization.  By (5.2), the full
effective divisor, with coefficients, is

\[
 \sum_{H''(\alpha)=0}\operatorname{ord}_\alpha(H'')[r=\alpha].
                                                               \tag{5.6}
\]

Its irreducible support components are

\[
 V(r-\alpha)\times\mathbb G_m\times\mathbb A^s,
 \qquad H''(\alpha)=0.                               \tag{5.7}
\]

The units of (5.5) are exactly

\[
 B_H^*=\mathbb C^*\xi^{\mathbb Z}.                   \tag{5.7}
\]

Let `r'` be the normalization coordinate for `G`.  For two distinct roots
`beta_1,beta_2` of `G''`, the induced ring isomorphism has

\[
 r'-\beta_i\longmapsto
 u_i(r-\alpha_i),                                    \tag{5.8}
\]

where `u_i` is a unit and `alpha_i` is a root of `H''`.  Subtracting the two
relations gives

\[
 \beta_2-\beta_1
 =u_1(r-\alpha_1)-u_2(r-\alpha_2).                  \tag{5.9}
\]

Comparing coefficients of the polynomial variable `r` over
`C[xi,xi^{-1},z_1,...,z_s]` gives `u_1=u_2=u`.  Equation (5.9) then makes `u`
a nonzero complex constant.  Therefore

\[
 r'\longmapsto ur+v,
 \qquad u\in\mathbb C^*,\ v\in\mathbb C,            \tag{5.10}
\]

and the same affine transformation carries every component in (5.7).  Since
the Fitting ideals, not only their radicals, are identified, the valuation
along each component also agrees.  This is the asserted affine equivalence of
effective divisors.  QED.

## 6. Exact recovery and the coarse orbit statement

Normalize the seed by

\[
 \widetilde H_\lambda=H_\lambda/c_\lambda,
 \qquad \widetilde H_\lambda'(1)=-1.
\]

The [faithfulness theorem for the full marked cover](DECORATED_NORMALIZATION_INVARIANT.md)
applies because the seed has a double zero at `0`, a distinguished simple
root at `1`, and the complete reconstruction open is part of the stable
boundary object.  Stable equivalence therefore recovers the normalized seed
itself:

\[
 \widetilde H_\lambda=\widetilde H_\mu.
\]

Their second derivatives have the fixed roots `p,q` and the remaining roots
`lambda,mu`, respectively.  Equality of the normalized polynomials therefore
forces `lambda=mu`, proving the sharp statement (2.2).

For completeness, the independent Fitting-divisor argument gives the
following coarser conclusion.

For the present family, (1.1) and the exclusion
`10 lambda^2-8 lambda+1 != 0` make the intrinsic Hessian-root divisor reduced
with support

\[
 S_\lambda=\{p,q,\lambda\}.                          \tag{6.1}
\]

Proposition 5.1 shows that stable equivalence implies an affine transformation
of `A^1` carrying `S_lambda` to `S_mu`.

Normalize the ordered pair `(p,q)` to `(0,1)`.  The third point becomes
`chi(lambda)`.  Permuting the three finite points produces exactly the six
values in (2.3).  Hence (2.4) and (2.5) follow.

The standard symmetric function of this six-element orbit is

\[
 256\frac{(1-x+x^2)^3}{x^2(1-x)^2}.
\]

Substitution of `x=chi(lambda)` simplifies to (2.7), proving (2.8).  This
six-element orbit and `J(lambda)` remain useful explicit coarse invariants,
but exact parameter recovery is stronger.

## 7. Why the next boundary jet is not the modulus

Along `C=0`, equation (4.6) makes the reduced intersection a smooth conic.  At
a smooth point of that intersection, `Z_Delta` is a graph

\[
 A=f_\lambda(B,C)
\]

in the completed target ring.  A formal shear

\[
 A\longmapsto A-f_\lambda(B,C)+f_\mu(B,C)
\]

fixes `C=0` and identifies the two completed divisor pairs.  Consequently raw
coefficients in

\[
 C^{-2}\operatorname{Disc}_W(E)
 =D_0(A,B)+CD_1(A,B)+\cdots
\]

are not invariants before quotienting by this large shear group.  The local
completed pair at the reduced boundary contact is too flexible to retain the
seed parameter.

The normalization ramification divisor is global on
`Z_Delta\setminus Z_0`, intrinsic, and unaffected by those shears.  It is the
first place where degree five has enough marked points to carry a modulus.

## 8. Exact reproduction

Run

```bash
.venv/bin/python scripts/verify_degree_five_stable_moduli.py
.venv/bin/python scripts/verify_affine_branch_mark_audit.py
.venv/bin/python scripts/verify_stable_generator_rigidity.py
```

The first checker verifies the seed identities, admissibility, polynomial
weighted map, constant Jacobian, inverse equation, split-root and boundary
exceptional polynomials, exact `C^2` discriminant saturation, the reduced
boundary conic, and the rational orbit invariant (2.7).  The latter two audit
the affine root-one mark and reconstruction-pole calculation used by generic
affine-mark faithfulness.

No external specialist review of the Hessian/Fitting-divisor theorem or this
stable-moduli consequence is currently recorded.  Long's external consequence
papers concern different witnesses and do not supply such a review.
