# The foundational incidence construction

This note isolates the coordinate-free core of the cubic counterexample.  It
separates three questions that are easy to conflate:

1. why normalized multiplication is étale for **every** target hyperplane;
2. why the normalization is the marked-factor projective open and has generic
   degree three;
3. why the tangent nonosculating hyperplane is exceptional: its normalized
   source is affine three-space.

Work over a characteristic-zero field `k`.  Put

\[
 V_i=\operatorname{Sym}^i(k^2).
\]

Points of `V_1` and `V_2` are binary linear and quadratic forms.  Write
`R(L,Q)=Res(L,Q)`.

## Lemma 1: universal coefficient--resultant étaleness

Define

\[
 \Phi:V_1\times V_2\longrightarrow V_3\times\mathbb A^1,
 \qquad (L,Q)\longmapsto(LQ,R(L,Q)).                 \tag{1}
\]

Then `Phi` is étale on `R!=0`.

### Proof

Suppose the differential of multiplication kills a tangent vector:

\[
 \dot LQ+L\dot Q=0.                                  \tag{2}
\]

On `R!=0`, the factors are coprime.  Equation (2) implies that `L` divides
`dot L` and `Q` divides `dot Q`.  The degree bounds force

\[
 (\dot L,\dot Q)=\lambda(L,-Q).                      \tag{3}
\]

Thus relative scaling is the entire kernel of infinitesimal multiplication.
The resultant is homogeneous of degree two in the coefficients of `L` and
degree one in the coefficients of `Q`.  Therefore

\[
 dR_{(L,Q)}(\lambda L,-\lambda Q)
 =\lambda(2-1)R(L,Q)=\lambda R(L,Q).                 \tag{4}
\]

This is nonzero when `lambda!=0` and `R!=0`, so the differential of (1) is
injective.  Its source and target both have dimension five.  Consequently

\[
 \boxed{d\Phi_{(L,Q)}\text{ is an isomorphism},}
\]

and `Phi` is étale at `(L,Q)`.  This conclusion concerns the differential;
it does not assert that `Phi` is globally an isomorphism.  The coordinate
identity `det D Phi=-R^2` is a consequence, not the mechanism of the proof.

## Lemma 2: every normalized hyperplane slice is étale

Let `ell in V_3^*` be nonzero and define

\[
 X_\ell=\{(L,Q):R(L,Q)=1,\ \ell(LQ)=1\}.             \tag{5}
\]

Multiplication induces

\[
 \mu_\ell:X_\ell\longrightarrow
 H_\ell:=\{C\in V_3:\ell(C)=1\}\simeq\mathbb A^3.
                                                                    \tag{6}
\]

Then `mu_ell` is étale.  In particular, invertibility of its differential does
not depend on the contact type of `ell`.

### Proof

The diagram

\[
\begin{array}{ccc}
X_\ell&\longrightarrow&V_1\times V_2\\
\downarrow\mu_\ell&&\downarrow\Phi\\
H_\ell&\longrightarrow&V_3\times\mathbb A^1,
\qquad C\longmapsto(C,1)
\end{array}                                           \tag{7}
\]

is Cartesian.  Its upper-left corner automatically lies in `R!=0`.
Étaleness is preserved by base change, so Lemma 1 proves the assertion.

There is also a canonical projective description.  Put

\[
 U_\ell=
 \bigl(\mathbb P(V_1)\times\mathbb P(V_2)\bigr)
 \setminus\bigl(\{R=0\}\cup\{\ell(LQ)=0\}\bigr).   \tag{8}
\]

Projectivization gives an isomorphism `X_ell -> U_ell`.  Indeed, for any
representatives of a point of (8), set

\[
 m=\ell(LQ),\qquad r=R(L,Q).
\]

Both are nonzero.  Under `(L,Q)->(alpha L,beta Q)`, they transform as

\[
 m\longmapsto\alpha\beta m,qquad
 r\longmapsto\alpha^2\beta r.
\]

The unique rescaling for which both values equal one is

\[
 \alpha={m\over r},\qquad \beta={r\over m^2}.         \tag{9}
\]

Formula (9) is independent of the initial representatives and is regular on
`U_ell`.  It is the inverse to projectivization.

For a squarefree cubic `C`, the fiber of (6) consists of its three choices of
a linear factor `L`; the residual quadratic is `Q=C/L`, and (9) fixes the
scalings uniquely.  Hence `mu_ell` is generically three-to-one.

## Proposition 3: hyperplane orbits and the exceptional affine slice

Assume now that `k` is algebraically closed.  Up to projective change of
variables and rescaling `ell`, cubic hyperplanes have exactly three contact
types with the twisted cubic:

\[
 (1,1,1),\qquad(2,1),\qquad(3).                      \tag{10}
\]

Their normalized sources have classes

\[
\begin{array}{c|c}
\text{contact type}&[X_\ell]\ \\ \hline
(1,1,1)&\mathbb L^3-\mathbb L\\
(2,1)&\mathbb L^3\\
(3)&\mathbb L^3-\mathbb L^2.
\end{array}                                          \tag{11}
\]

Over `C`, it follows that

\[
 \boxed{X_\ell\simeq\mathbb A^3
 \quad\Longleftrightarrow\quad
 \ell\text{ has contact type }(2,1).}                \tag{12}
\]

The positive case is the tangent but nonosculating hyperplane.  Thus the
three mechanisms in the construction are logically separate:

\[
\begin{array}{rcl}
\text{étaleness}&\Longleftarrow&\text{Lemma 1 and base change, for every }\ell,\\
\text{generic noninjectivity}&\Longleftarrow&\text{three choices of a linear factor},\\
X_\ell\simeq\mathbb A^3&\Longleftarrow&\text{the exceptional }(2,1)\text{ hyperplane geometry}.
\end{array}                                           \tag{13}
\]

This is an internal repository proposition with the complete proof below.  It
is not attributed to the available public material and has no recorded
external review.

### Proof of the orbit list and the two negative cases

The restriction of `ell` to the twisted cubic is a binary cubic.  Its zero
divisor is invariant under the projective `SL_2` action.  Conversely, over an
algebraically closed field, `PGL_2` is transitive on ordered triples of
distinct points and on ordered pairs of distinct points.  Projectively the
`SL_2` and `PGL_2` orbits agree, so the three multiplicity partitions in
(10) are exactly the three hyperplane orbits.

Compute the class of `U_ell`, hence of `X_ell` by Lemma 2, by projecting

\[
 U_\ell\subset\mathbb P(V_1)\times\mathbb P(V_2)
\]

to the marked linear factor `[L] in P(V_1)=P^1`.  In the residual `P^2`, the
resultant condition removes the line of quadratics divisible by `L`.  The
hyperplane condition removes the kernel of

\[
 Q\longmapsto\ell(LQ).                                \tag{14}
\]

There are three possibilities.

- At every point which is not a multiple contact point, the two lines are
  distinct, so the fiber has
  class
  \[
  [\mathbb P^2]-2[\mathbb P^1]+[\mathrm{pt}]
  =\mathbb L^2-\mathbb L.                             \tag{15}
  \]
- At a double contact point, (14) vanishes precisely on the resultant line.
  The two removed lines coincide and the fiber is `A^2`, of class
  `L^2`.
- At a triple contact point, (14) is identically zero.  The pulled-back
  hyperplane contains the whole residual fiber, so the open fiber is empty.

For type `(1,1,1)`, every marked point has the first fiber type.  Hence

\[
 [X_\ell]=(1+\mathbb L)(\mathbb L^2-\mathbb L)
 =\mathbb L^3-\mathbb L.                              \tag{16}
\]

For type `(2,1)`, the double point contributes `L^2`, while its affine-line
complement in the base contributes fibers of class `L^2-L`:

\[
 [X_\ell]=\mathbb L^2+\mathbb L(\mathbb L^2-\mathbb L)
 =\mathbb L^3.                                        \tag{17}
\]

For type `(3)`, the fiber over the triple point is empty and the remaining
base is `A^1`, so

\[
 [X_\ell]=\mathbb L(\mathbb L^2-\mathbb L)
 =\mathbb L^3-\mathbb L^2.                            \tag{18}
\]

This proves (11).  Applying the Hodge--Deligne realization over `C` separates
the first and third classes from that of `A^3`.

### Proof of the positive case

Suppose the contact divisor is

\[
 2p+q,\qquad p\ne q.                                  \tag{19}
\]

Then

\[
 \boxed{X_\ell\simeq\mathbb A^3.}                    \tag{20}
\]

Every hyperplane of contact type `(2,1)` gives the same normalized source up
to an isomorphism induced by `PGL_2` and factor rescaling.

`PGL_2` acts transitively on ordered pairs `(p,q)` with `p!=q`, and binary-form
multiplication commutes with that action.  If the action changes the
resultant and hyperplane functional by nonzero scalars `delta` and `eta`, the
additional factor rescaling

\[
 L\longmapsto {\eta\over\delta}L,
 \qquad Q\longmapsto {\delta\over\eta^2}Q
\]

restores both equations in (5).  It is therefore enough to treat one
tangent-nonosculating representative.

For the representative

\[
 L=aT+bS,\qquad Q=cT^2+dTS+eS^2,qquad
 \ell(LQ)=ad+bc,
\]

the normalized incidence is the complete intersection

\[
 a^2e-abd+b^2c=1,\qquad ad+bc=1.                     \tag{21}
\]

It has the global polynomial coordinates

\[
\begin{aligned}
b&=1+ay,\\
c&=1-\frac32ay+a^2z,\\
d&=\frac12y-az+\frac32ay^2-a^2yz,\\
e&=-2z+4y^2-4ayz+3ay^3-2a^2y^2z,
\end{aligned}                                        \tag{22}
\]

with inverse

\[
 y=2bd-ae,\qquad
 z=2d^2+ce+6bd^2+3bce-\frac92e.                      \tag{23}
\]

Both compositions are polynomial identities modulo (21), including on
`a=0`.  Thus (22)--(23) prove (20), completing the proof of (12).  The detailed ideal reductions and the
residual-torus equivariance are recorded in the
[explicit normalized model](NORMALIZED_FACTORIZATION_MODEL.md).

## Foundational consequence

On the representative (21), discard the fixed coefficient `ad+bc=1` and put

\[
 G=(ac,ae+bd,be):\mathbb A^3_{a,y,z}\longrightarrow\mathbb A^3.
\]

Lemma 2 makes `G` étale, Lemma 3 makes its source affine three-space, and the
generic fiber has three points.  Moreover `det DG=-1`, and the originally
announced map satisfies

\[
 F_{\rm original}=B\circ G\circ A,
\]

where

\[
 A(z_1,z_2,z_3)=(z_1,z_2,-z_3/2),\qquad
 B(u_1,u_2,u_3)=(u_3,2u_2,2u_1).
\]

The construction is therefore the short chain

\[
\boxed{
\text{coprime factorization}
\Longrightarrow\text{étale normalized multiplication}
\Longrightarrow\text{exceptional affine slice}
\Longrightarrow\text{announced polynomial}.}
\]

The exact checker is
[`verify_normalized_factorization_slice.py`](../scripts/verify_normalized_factorization_slice.py).
