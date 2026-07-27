# Intrinsic algebraic-torus exclusion for the quartic Keller map

## 1. Result and scope

Let \(k\) be an algebraically closed field of characteristic zero and let
\(F:\mathbb A^3\to\mathbb A^3\) be the determinant-one quartic map

\[
\begin{aligned}
t&=1+xy,\qquad q=t^2z-y^2(1+3t),\\
F_1&=-\frac12tq,\\
F_2&=y-3xq-tq+2t^2x^2q^4,\\
F_3&=x(5-3t)+x^3z-(xq)^4.
\end{aligned}
\tag{1}
\]

The canonical normalization-boundary package of \(F\) intrinsically selects
the ramified target divisor \(Z_\Delta\), the second boundary divisor
\(Z_0\), and, after deleting \(Z_0\), the finite normalization

\[
\operatorname{Spec}k[P^{\pm1},r^{\pm1}]
 \longrightarrow Z_\Delta\setminus Z_0.
\tag{2}
\]

Retain its relative-differential Fitting divisor, the intrinsically ordered
toric punctures, and the base character determined up to scalar by \(Z_0\).
The algebraic automorphism group of this decorated ramified-normalization
stratum is

\[
\boxed{\mu _5,}
\tag{3}
\]

and its infinitesimal automorphism algebra is zero. More importantly,

\[
\boxed{\text{\(F\) admits no nontrivial algebraic
\(\mathbb G_m\) source--target equivariance.}}
\tag{4}
\]

Consequently \(F\) is not polynomially left--right equivalent to an
algebraic-torus-equivariant map. The conclusion is geometric over every
characteristic-zero ground field, by extension to its algebraic closure.
In the terminology of Shaska's
[*Graded Keller maps and the Jacobian
Conjecture*](https://arxiv.org/abs/2607.20210), this says that no polynomial
left--right representative of \(F\) is graded for a nontrivial weight
signature.

Statement (4) is stronger than the affine-linear calculation in
[the original counterexample note](NO_LINEAR_TORUS_COUNTEREXAMPLE.md).
The \(18\times18\) determinant \(-5\) and the \(785\times24\) rank
calculation themselves still exclude only affine-linear actions. The
upgrade here uses the canonical boundary package, a separate pointwise-fixed
hypersurface lemma, and generic deck rigidity.

This note does **not** compute every discrete or unipotent polynomial
left--right self-equivalence of \(F\). It computes the exact automorphism
group of the selected decorated ramified stratum and proves that the kernel
of restriction from polynomial left--right self-equivalences contains no
positive-dimensional torus. Identity stabilization is outside the claim:
an added identity factor carries tautological torus actions.

This last qualification is logically necessary.  For example,
\(F\times\operatorname{id}_{\mathbb A^1}\) is equivariant for the action
which fixes the original three coordinates and scales the added coordinate
on both source and target.  Thus finiteness of the unstabilized decorated
automorphism group cannot imply literal absence of positive-dimensional
symmetry after identity stabilization.  What does survive stabilization is
the following narrower intrinsic statement: every connected torus acting on
the stabilized decorated stratum acts trivially on its two-dimensional unit
lattice.  Such an action may still move the polynomial stabilization
coordinates.

### A sparser rational representative

There is a sparser representative of the same geometric quartic
left--right class.  Take

\[
 G_{\rm sp}(S)
 =S(S-1)(S-2)(3S+2)
 =3S^4-7S^3+4S.
\tag{4a}
\]

Thus \(g_2=0\), while \(g_1g_3g_4\ne0\).  Put

\[
 t=1+xy,\qquad
 q=t^2z-\frac47y^2(1+3t)
\]

and define

\[
\boxed{
\begin{aligned}
F_{{\rm sp},1}&=-\frac12tq,\\
F_{{\rm sp},2}&=y-\frac{21}{4}xq+3t^2x^2q^4,\\
F_{{\rm sp},3}&=x(5-3t)+\frac74x^3z-\frac32(xq)^4.
\end{aligned}}
\tag{4b}
\]

Then \(\det DF_{\rm sp}=1\), and the four distinct rational points

\[
\begin{aligned}
(0,0,1),\qquad&
\left(-\frac45,\frac94,-\frac{265}{32}\right),\\
\left(\frac12,-\frac32,100\right),\qquad&
\left(\frac3{10},-\frac{29}{6},-\frac{24820}{729}\right)
\end{aligned}
\tag{4c}
\]

all map to \((-1/2,0,0)\).

In expanded affine coordinates the three components have respectively

\[
\boxed{(7,51,38)}
\tag{4d}
\]

nonzero monomials and ordinary degrees \((7,26,24)\).  This support count is
minimal in the displayed determinant-one quartic quadratic-gauge normal
form.  Indeed, after writing \(a_i=g_i/g_1\), every expanded coefficient is
a nonzero rational Laurent monomial in \(a_3,a_4\), except for exactly seven
coefficients of the second component which are nonzero Laurent monomials
times \(a_2\).  Hence \(a_2=0\) deletes exactly those seven terms, and no
specialization with \(a_3a_4\ne0\) deletes any other term.

This is a family-relative sparsity statement, not an absolute minimum over
all polynomial coordinate systems.  The
[degree-four stable-moduli theorem](../verified/QUADRATIC_GAUGE_STABLE_MODULI.md)
places (4b) and (1) in the same polynomial left--right class over the
algebraic closure.  Therefore (4) applies equally to \(F_{\rm sp}\): no
polynomial left--right representative of (4b) is equivariant for a
nontrivial algebraic torus.  Its visible diagonal symmetry is only the
finite remnant

\[
(x,y,z)\longmapsto(\zeta x,\zeta^{-1}y,\zeta^{-2}z),
\qquad \zeta^5=1,
\tag{4e}
\]

with target weights \((-2,-1,1)\) modulo \(5\).

The geometric-degree qualifier should not be overstated.  Degree two is
excluded for every noninvertible Keller map, while the foundational
geometric-degree-three example is graded.  The repository does not yet
classify every geometric-degree-three counterexample: the cubic
package-extraction theorem in
[the minimal-boundary program](MINIMAL_BOUNDARY_CLASSIFICATION.md#5-the-degree-three-proof-target)
is open.  Consequently (4b) is a geometric-degree-four genuinely ungraded
example and the smallest certified upper bound, but unconditional global
minimality in geometric degree remains open.

## 2. The decorated stabilizer

Divide the seed

\[
G(S)=2S-S^2-2S^3+S^4
\]

by \(g_1=2\), and remove its quadratic coefficient by the intrinsic target
shear used in the
[quadratic-gauge stable-moduli theorem](../verified/QUADRATIC_GAUGE_STABLE_MODULI.md).
Thus

\[
a_3=-1,\qquad a_4=\frac12.
\]

On the normalization (2), the relative-differential Fitting generator is

\[
\boxed{
J(P,r)=-1+3a_3Pr^2+8a_4P^4r^3
      =-1-3Pr^2+4P^4r^3.
}
\tag{5}
\]

Here is a finite lattice certificate before imposing the ordered boundary.
Write exponent vectors as rows and let a Laurent automorphism have exponent
matrix \(M\in\operatorname{GL}_2(\mathbb Z)\).  An automorphism of the
divisor of \(J\) must carry its support

\[
 \mathcal S=\{(0,0),(1,2),(4,3)\}
\tag{6}
\]

to a translate of itself.  Since the support is affinely independent, every
such affine map is determined by one of the six permutations of
\(\mathcal S\).  Solving those six \(2\times2\) systems gives

\[
\begin{array}{c|c|c|c}
\text{permutation}&M&\text{support translation}&\text{verdict}\\ \hline
(0,1,2)&\begin{pmatrix}1&0\\0&1\end{pmatrix}
 &(0,0)&\text{integral unimodular}\\[3pt]
(0,2,1)&\begin{pmatrix}-2&-1\\3&2\end{pmatrix}
 &(0,0)&\text{integral unimodular}\\[3pt]
(1,0,2)&\frac15\begin{pmatrix}9&8\\-7&-9\end{pmatrix}
 &(1,2)&\text{nonintegral}\\[3pt]
(1,2,0)&\frac15\begin{pmatrix}-11&-7\\13&6\end{pmatrix}
 &(1,2)&\text{nonintegral}\\[3pt]
(2,0,1)&\frac15\begin{pmatrix}6&7\\-13&-11\end{pmatrix}
 &(4,3)&\text{nonintegral}\\[3pt]
(2,1,0)&\frac15\begin{pmatrix}1&-3\\-8&-1\end{pmatrix}
 &(4,3)&\text{nonintegral}.
\end{array}
\tag{7}
\]

The second matrix is a genuine involutive symmetry of the bare Fitting
divisor.  One representative is

\[
 P\longmapsto-\frac34P^{-2}r^{-1},\qquad
 r\longmapsto-\frac43P^3r^2,
\tag{8}
\]

which interchanges the \(Pr^2\) and \(P^4r^3\) terms.  Together with the
\(\mu _5\) translations computed below, it gives the order-ten dihedral
stabilizer of the bare divisor.  But (8) does not preserve the intrinsic
base character \(P\) up to scalar, so the ordered second-boundary image
rejects this sole nonidentity lattice matrix.

The ordered boundary package first forces an automorphism of the normalized
stratum to have

\[
P\longmapsto\beta P,\qquad
r\longmapsto\alpha P^m r^\varepsilon,
\quad
\alpha,\beta\in k^\times,\quad
m\in\mathbb Z,\quad
\varepsilon\in\{1,-1\}.
\tag{9}
\]

This is a unit-lattice statement on
\(k[P^{\pm1},r^{\pm1}]\), not a choice of affine coordinates. The support
of (5) intrinsically orders the two \(r\)-punctures, so
\(\varepsilon=1\), and its \(Pr^2\) term forces \(m=0\). These are exactly
the support-rigidity steps proved in the stable-moduli theorem.

Preservation of the divisor of \(J\) gives

\[
J(\beta P,\alpha r)=uJ(P,r)
\tag{10}
\]

for a Laurent unit \(u\). The unique constant term and the extremal support
force \(u=1\). Comparing the other two coefficients gives

\[
\beta\alpha^2=1,\qquad
\beta^4\alpha^3=1.
\tag{11}
\]

Therefore

\[
\beta=\alpha^{-2},\qquad \alpha^5=1,
\tag{12}
\]

scheme-theoretically. Conversely these transformations are induced by
actual left--right symmetries of \(F\), so (3) is exact rather than an upper
bound.

For \(\zeta^5=1\), the explicit symmetries in the coordinates of (1) are

\[
\sigma_\zeta(x,y,z)
 =(\zeta x,\zeta^{-1}y,\zeta^{-2}z)
\tag{13}
\]

on the source and

\[
\tau_\zeta(U,V,W)
=
\left(
\zeta^{-2}U,\,
\zeta^{-1}V-2(\zeta^{-1}-\zeta^{-2})U,\,
\zeta W
\right)
\tag{14}
\]

on the target. Direct substitution gives

\[
F\circ\sigma_\zeta=\tau_\zeta\circ F.
\tag{15}
\]

The shear in (14) is the conjugate of the diagonal action after restoring
the quadratic seed coefficient.

To compute the infinitesimal algebra, put
\(\alpha=1+\epsilon u\), \(\beta=1+\epsilon v\), and
\(\epsilon^2=0\). Linearizing (11) gives

\[
\begin{pmatrix}2&1\\3&4\end{pmatrix}
\binom uv=0.
\tag{16}
\]

The determinant is \(5\), so \(u=v=0\) in characteristic zero. Thus the
decorated automorphism algebra vanishes even though the finite residual
group \(\mu _5\) is nontrivial.

## 3. A pointwise-fixed hypersurface lemma

The remaining issue is faithfulness: a torus might act trivially on the
decoration while acting nontrivially away from it. The following elementary
lemma excludes that possibility here.

**Lemma.** Let \(A=k[x_1,\ldots,x_n]\) carry an algebraic
\(\mathbb G_m\)-action, let \(f\in A\) be irreducible, and suppose the action
on \(A/(f)\) is trivial. If \(A/(f)\) is not normal, then the action on
\(A\) is trivial.

**Proof.** The stable prime ideal \((f)\) has a homogeneous generator; write
its weight as \(d\). Decompose

\[
A=\bigoplus_{j\in\mathbb Z}A_j.
\]

Because the quotient action is trivial, every \(A_j\) with \(j\ne0\) lies
in \((f)\).

If \(d=0\), a nonzero element of \(A_j\), \(j\ne0\), is divisible by every
power of \(f\), which is impossible. Hence all weights are zero and the
action is trivial.

Suppose \(d\ne0\), and reverse the action if necessary so \(d>0\). A
negative-weight element is successively divisible by arbitrarily high
powers of \(f\), so there are no negative weights. There are likewise no
weights strictly between \(0\) and \(d\). Induction gives

\[
A_{md}=f^mA_0\quad(m\geq0),
\qquad
A_j=0\quad(d\nmid j).
\]

Distinct weights show that \(f\) is algebraically independent over \(A_0\);
therefore

\[
A=A_0[f],\qquad A/(f)=A_0.
\]

But normality of the polynomial ring \(A=A_0[f]\) implies normality of
\(A_0\): an element of \(\operatorname{Frac}(A_0)\) integral over \(A_0\)
is integral over \(A_0[f]\), and its membership in \(A_0[f]\) then places it
in \(A_0\). This contradicts the hypothesis. Thus \(d\ne0\) is impossible,
and the action is trivial. \(\square\)

## 4. Faithfulness for the quartic cover

Let a \(\mathbb G_m\)-action on the source and target make (1) equivariant.
[Stable-normalization functoriality](../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md)
transports the intrinsically ordered boundary package and induces a
homomorphism

\[
\mathbb G_m\longrightarrow\operatorname{Aut}(\mathcal D_F)=\mu _5,
\tag{17}
\]

where \(\mathcal D_F\) is the decorated stratum of Section 2. A connected
torus has no nontrivial morphism to the finite étale group \(\mu _5\) in
characteristic zero. Hence it acts trivially on the normalization (2), and
therefore fixes the dense open \(Z_\Delta\setminus Z_0\), and then all of
\(Z_\Delta\), pointwise.

The divisor \(Z_\Delta\) is prime and nonnormal. Indeed its finite
birational normalization on \(P\ne0\) has

\[
\operatorname{Fitt}_0
\Omega_{k[P^{\pm1},r^{\pm1}]/R_\Delta}=(J).
\tag{18}
\]

The Laurent polynomial (5) is not a unit, so the differential module is
nonzero and the normalization is not an isomorphism. The lemma now forces
the target \(\mathbb G_m\)-action on \(\mathbb A^3\) to be trivial.

The source action is consequently a target-fixed deck action. The inverse
quartic has geometric monodromy \(S_4\) by the
[quadratic-gauge monodromy theorem](ROOT_ENGINEERED_QUADRATIC_GAUGE.md#7-discriminant-normalization-and-monodromy).
Its sheet stabilizer is \(S_3\), which is self-normalizing in \(S_4\), so the
generic deck group is

\[
N_{S_4}(S_3)/S_3=1.
\tag{19}
\]

Thus the source action is also trivial. This proves (4).
A positive-dimensional algebraic torus over \(k\) contains a nontrivial
one-parameter subtorus, so the same argument excludes every
positive-dimensional torus action.

Finally, if a polynomial left--right equivalent map carried a nontrivial
algebraic torus equivariance, conjugating both actions by the polynomial
source and target changes would give such an equivariance of \(F\).
Therefore no polynomial left--right representative of (1) is
algebraic-torus-equivariant.

## 5. The exact statement after stabilization

Although literal symmetry-freeness after stabilization is impossible, the
unit-lattice conclusion has a stable form.  Let

\[
 A_s=k[P^{\pm1},r^{\pm1},t_1,\ldots,t_s].
\]

Its units are exactly \(k^\times P^{\mathbb Z}r^{\mathbb Z}\).  A connected
algebraic torus action preserving the pulled-back ordered decoration cannot
act nontrivially on this discrete exponent lattice.  For a one-parameter
subtorus its coaction on the two primitive units therefore has the form

\[
 P\longmapsto\lambda^pP,\qquad
 r\longmapsto\lambda^qr
\tag{20}
\]

for integers \(p,q\).  Preservation of the pulled-back Fitting divisor gives

\[
 p+2q=0,\qquad4p+3q=0.
\tag{21}
\]

The determinant is \(-5\), so \(p=q=0\).  Consequently every connected
torus action on the stabilized decorated normalization is vertical over
\(\operatorname{Spec}k[P^{\pm1},r^{\pm1}]\).  This proves that stabilization
creates no positive-dimensional symmetry visible to the canonical unit
lattice or Fitting data.  It does **not** prove that the vertical action is
conjugate to an action on the added affine coordinates alone; such a
splitting theorem would require additional input not supplied by the
boundary certificate.

## 6. Exact reproduction

Run

```bash
make verify-algebraic-torus-free
```

The checker reconstructs (5), enumerates all six permutations of its support,
certifies the two lattice matrices in (7), rejects the involution by the
intrinsic base character, verifies the stabilizer equations (11), computes
the tangent determinant \(5\), and checks (15) modulo
\(\zeta^5-1\). The pointwise-fixed hypersurface lemma and the monodromy
faithfulness argument are exact proofs in Sections 3--4; they are not
inferred from a bounded search or a computer rank calculation.
