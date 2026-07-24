# Same-degree stable inequivalence before moduli

This note answers the first classification question after the existence of
noninvertible Keller maps:

\[
 \boxed{\text{two explicit Keller maps of the same geometric degree need
 not be stably polynomially left--right equivalent.}}
\]

In fact, degree four already contains three pairwise distinct stable classes
coming from the weighted, cancellation, and root-engineered quadratic-gauge
constructions.  The separation uses only intrinsic strata of the canonical
Zariski--Main normalization package.  It therefore survives arbitrary
polynomial source and target automorphisms and arbitrary identity
stabilization.

The explicit theorem comes first.  The positive-dimensional weighted-moduli
consequence is recorded only in Section 6.

Work over an algebraically closed characteristic-zero field.  The weighted
and quadratic witnesses are defined over `Q`.  The cancellation witness is
defined over

\[
 K=\mathbb Q(\theta),\qquad \theta^2-4\theta+6=0.
\]

All geometric comparisons may be made after extending the common ground
field to its algebraic closure.

## 1. Three explicit quartic Keller maps

### 1.1 A weighted map

Put

\[
 u=1+3xy,\qquad \gamma=1-4xy-x^2z
\]

and define

\[
 F^{\rm wt}=(A_{\rm wt},B_{\rm wt},C_{\rm wt})
\]

by

\[
 \begin{aligned}
 A_{\rm wt}={2u+u^2-3u^4\gamma^2\over x^2},\qquad
 B_{\rm wt}={1+u-2u^3\gamma^2\over x},\qquad
 C_{\rm wt}=x\gamma.                                  \tag{1.1}
 \end{aligned}
\]

Both displayed quotients are polynomials.  This is the split weighted seed

\[
 H(W)=\frac{W^2-W^4}{2},
\]

equivalently the normalized quartic parameter `alpha=-1/2`.  Its marked
inverse equation, after multiplication by two, is

\[
 E_{\rm wt}(W)=W^2-W^4-2B_{\rm wt}C_{\rm wt}W
                  +A_{\rm wt}C_{\rm wt}^2.            \tag{1.2}
\]

Direct calculation gives

\[
 \det DF^{\rm wt}=-6,
\]

and (1.2), together with regular reconstruction away from its derivative,
gives geometric degree four.

### 1.2 A cancellation map

Put

\[
 A=1+xy^2,\qquad
 h(A)=\theta+(4\theta-6)A,
\]

\[
 B=A^2z+y^3h(A),\qquad P=AB,\qquad Q=y+xB.             \tag{1.3}
\]

Define the following polynomial:

\[
\begin{aligned}
R=-{x\over12}\Big[&
\theta\bigl(
 8x^4y^5z+16x^3y^6+10x^3y^3z+20x^2y^4
\bigr)\\
&+x^5y^4z^2-12x^4y^5z+2x^4y^2z^2-60x^3y^6
 -12x^3y^3z+x^3z^2\\
&\hspace{35mm}-48x^2y^4+4x^2yz+18xy^2-12
\Big].
                                                               \tag{1.4}
\end{aligned}
\]

Then

\[
 F^{\rm can}=(P,Q,R):\mathbb A^3_K\longrightarrow\mathbb A^3_K
                                                               \tag{1.5}
\]

is the cancellation map of type `(m,r)=(2,1)` associated with the parameter
root `theta`.  Equation (1.4) is the denominator-free expansion of

\[
 {x\over A}
 -\frac{Q^2}{2}\left({x\over A}\right)^2
 +\frac{2PQ}{3}\left({x\over A}\right)^3
 -\frac{P^2}{4}\left({x\over A}\right)^4.
                                                               \tag{1.6}
\]

The relation `theta^2-4theta+6=0` is exactly the finite cancellation
condition which removes the apparent denominator.  Direct calculation gives

\[
 \det DF^{\rm can}=-1.
\]

The inverse polynomial is

\[
 \Psi_{\rm can}(T)
 =T-\frac{Q^2T^2}{2}+\frac{2PQT^3}{3}
       -\frac{P^2T^4}{4}-R.                            \tag{1.7}
\]

Its derivative is

\[
 \partial_T\Psi_{\rm can}=1-T(Q-PT)^2.
\]

The cancellation reconstruction theorem therefore gives exact geometric
degree four.  For example, at `(P,Q,R)=(1,0,0)`, equation (1.7) is
`T(1-T^3/4)`, with four distinct geometric roots and nonzero reconstruction
denominator.

### 1.3 A quadratic-gauge map

Put

\[
 t=1+xy,\qquad q=t^2z-y^2(1+3t)
\]

and define

\[
\begin{aligned}
 F^{\rm quad}=\bigl(&tq,\\
 &-y+3xq+tq-2t^2x^2q^4,\\
 &-2x+3x^2y-x^3z+x^4q^4\bigr).                        \tag{1.8}
\end{aligned}
\]

This is the root-engineered quadratic gauge attached to

\[
 G(S)=S(S-1)(S+1)(S-2)
     =S^4-2S^3-S^2+2S.
\]

It satisfies

\[
 \det DF^{\rm quad}=-2.
\]

After changing the signs of the last two displayed target coordinates and
calling the resulting coordinates `(B,C)`, its inverse equation is

\[
 E_{\rm quad}(S)
 =2S-PS^2-2PS^3+P^4S^4-BS^2-C.                       \tag{1.9}
\]

At `(P,B,C)=(1,0,0)` this is exactly `G(S)`, so the four simple roots
reconstruct the complete fiber.  Thus its geometric degree is also four.

## 2. The stable invariants

For a dominant quasi-finite polynomial map `F:U->Y`, let

\[
 \mathcal B(F)=
 \left(
 \overline X_F=\operatorname{Norm}_Y k(U)\longrightarrow Y,\;
 U\hookrightarrow\overline X_F,\;
 \partial_F
 \right)                                               \tag{2.1}
\]

be its canonical finite-normalization package.  For all three displayed
maps, the complete boundary has two intrinsically ordered target images:

- `Z_Delta`, the unique vertex receiving a ramified boundary prime;
- `Z_0`, the other boundary vertex.

The order, the finite stratum maps, their relative Fitting ideals, and the
scheme-theoretic intersection `Z_Delta cap Z_0` are reconstructed from
`mathcal B(F)`.  Stable normalization functoriality shows that stabilization
only takes their product with affine space.  In particular, unit rank,
Laurent-support rank of a principal Fitting divisor, and nilradical
nilpotency index are stable polynomial left--right invariants.

## 3. Weighted versus the two reciprocal maps: `A1` against `Gm`

Delete the intrinsic second vertex `Z_0` and normalize the remaining
ramified target stratum.  The three answers are

\[
\begin{array}{c|c|c}
\text{map}&
\operatorname{Norm}(Z_\Delta\setminus Z_0)&
\operatorname{rank}(\mathcal O^\times/k^\times)\\ \hline
F^{\rm wt}&\mathbb A^1\times\mathbb G_m&1\\
F^{\rm can}&\mathbb G_m^2&2\\
F^{\rm quad}&\mathbb G_m^2&2.
\end{array}                                             \tag{3.1}
\]

The first factors are the promoted versions of the one-puncture
`\mathbb A^1` weighted critical normalization and the two-puncture
`\mathbb G_m` reciprocal critical normalizations.  The extra common
`\mathbb G_m` is forced by the three-dimensional intrinsic stratum; a
fixed-coordinate plane slice alone would not be invariant under arbitrary
target automorphisms.

Polynomial affine-space factors introduce no new units.  Hence unit rank
separates `F^wt` stably from both `F^can` and `F^quad`.

## 4. Cancellation versus quadratic gauge: Fitting-support geometry

Both remaining ramified strata are two-dimensional tori, so unit rank is
neutral.  Their relative zeroth-Fitting divisors are not isomorphic.

For cancellation type `(2,1)`, use torus coordinates `(Y,u)` with
`Q=(u+1)Y`.  Up to a Laurent monomial unit, the Fitting generator is

\[
 J_{\rm can}=2u-1.                                    \tag{4.1}
\]

Its Laurent support has two points and affine rank one.

For the quartic quadratic gauge, shear away the quadratic seed coefficient
and divide by the linear coefficient.  The normalized seed has
`a_3=-1`, `a_4=1/2`, so the Fitting generator is

\[
 J_{\rm quad}=-1-3PS^2+4P^4S^3.                      \tag{4.2}
\]

Its support

\[
 (0,0),\qquad(1,2),\qquad(4,3)
\]

has affine rank two.  Every automorphism of a two-dimensional torus acts on
Laurent exponents by an affine `GL_2(Z)` transformation.  Stabilization
introduces no additional units, so it cannot change either affine support
rank.  Thus (4.1) cannot be carried to (4.2), proving that `F^can` and
`F^quad` are stably inequivalent.

## 5. Independent boundary-contact check

The scheme-theoretic intersection of the two intrinsic target boundary
images gives a second complete separation:

\[
\begin{array}{c|c|c}
\text{map}&\text{reduced intersection pieces}&
\text{nilradical index }\mu\\ \hline
F^{\rm wt}&\mathbb A^1&1\\
F^{\rm quad}&\mathbb A^1\sqcup\mathbb G_m&2\\
F^{\rm can}&\mathbb A^1\sqcup\mathbb G_m&6.
\end{array}                                             \tag{5.1}
\]

For the weighted seed `alpha=-1/2`, the saturated boundary contact is a
reduced conic.  For every quartic quadratic gauge it has the form

\[
 k[B,C]/\bigl(B^2(BC-1)\bigr),
\]

so its nilradical has exact index two.  For cancellation type `(2,1)` it
has the form

\[
 k[Q,R]\Big/\left(
 Q^6(2RQ^2-c_0)
 \right),\qquad c_0\ne0,
\]

so its nilradical has exact index six.  Polynomial extension preserves the
least nilpotency exponent.  Thus (5.1) independently proves all three
pairwise stable inequivalences.

Combining Sections 3--5 gives the promised small theorem:

> **Quartic stable-inequivalence theorem.**  The explicit Keller maps
> `F^wt`, `F^can`, and `F^quad` have the same geometric degree four and are
> pairwise inequivalent under polynomial left--right equivalence, even after
> adjoining arbitrarily many identity variables.

## 6. Only then: positive-dimensional stable moduli

Fix `N>=4`.  Let `A_N` be the normalized degree-`N` weighted seed space

\[
 H(0)=H'(0)=H(1)=0,\qquad H'(1)=-1.
\]

It is an affine space of dimension `N-3`.  On its nonempty Hessian-clean,
ordinary boundary-clean open, the canonical ramified stratum has
normalization

\[
 \mathbb A^1_r\times\mathbb G_m
\]

and relative Fitting divisor

\[
 \operatorname{Fitt}_0
 \Omega_{\widetilde Z_\Delta/Z_\Delta}
 =\bigl(H''(r)\bigr).                                  \tag{6.1}
\]

The resulting marked Hessian-divisor map

\[
 H\longmapsto
 \bigl(\mathbb P^1;\operatorname{div}(H''),0,\infty\bigr)
                                                               \tag{6.2}
\]

is generically etale of degree `N-2` onto an image of dimension `N-3`.
The finite ambiguity consists exactly of rerooting at the `N-2` nonzero
simple roots of `H`.  Thus the Hessian divisor alone already proves that
the stable-class image has dimension `N-3`.

The normalized cover contains one further intrinsic mark: over `Z_0`, the
root-one sheet is the unique unramified residue-degree-one sheet contained
in the affine reconstruction open and outside the double-zero cluster.
Adding that affine sheet to (6.2) kills the rerooting ambiguity.  After
shrinking the open, stable polynomial left--right equivalence then forces
equality of the normalized seeds.

Consequently:

> **Degreewise stable-moduli theorem.**  For every `N>=4`, the weighted
> construction contains a nonempty `(N-3)`-dimensional algebraic family of
> pairwise stably polynomially left--right inequivalent Keller maps of
> geometric degree `N`.

In degree five this is the explicit one-parameter family

\[
 H_\lambda(W)=
 {W^2(W-1)\bigl(3W^2-(5\lambda+1)W+3\lambda\bigr)\over60},
\]

for which

\[
 H_\lambda''(W)
 ={(W-\lambda)(10W^2-8W+1)\over10}.
\]

The three Hessian roots carry a cross-ratio modulus; the affine root sheet
then refines the generic six-element anharmonic ambiguity to exact recovery
of `lambda`.

This answers the uniqueness question in two distinct senses:

1. already in geometric degree four there are at least three pairwise
   stable classes arising from three different construction skeletons;
2. in every degree `N>=4` there is an `(N-3)`-dimensional continuum of
   stable classes inside the weighted skeleton alone.

The second statement is a theorem about continuously many inequivalent
maps, not a claim that there are continuously many construction skeletons
or irreducible components of a globally constructed moduli space.

## 7. Dependencies and exact regression

Stable functoriality is proved in
[Stable normalization functoriality](STABLE_NORMALIZATION_FUNCTORIALITY.md).
The family-specific boundary calculations are in
[Stable separation of the quadratic and weighted families](QUADRATIC_WEIGHTED_STABLE_SEPARATION.md),
[Stable intersection of the quadratic-gauge and cancellation families](QUADRATIC_CANCELLATION_STABLE_INTERSECTION.md),
and
[Boundary geometry](../cancellation/BOUNDARY_GEOMETRY.md).
The weighted Hessian-divisor and affine-sheet results are in
[The decorated-normalization invariant](../extended-geometry/DECORATED_NORMALIZATION_INVARIANT.md).

Run

```bash
.venv/bin/python scripts/verify_same_degree_stable_inequivalence.py
```

The regression checks the three displayed maps, their constant Jacobians,
quartic inverse equations or a four-root control fiber, the two Fitting
support ranks, and the three boundary-contact indices.  Stable functoriality
itself is a theorem, not a finite symbolic calculation.
