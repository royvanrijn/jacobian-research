# New counterexample searches: first exact pass

This note records the first active pass on the five proposed searches.  It
separates exact conclusions from computational evidence and from new search
targets.

## A. Positive-characteristic infinitesimal synchronization

Let \(k\) be a field of characteristic \(p>0\), let \(R\) be monic and
original of degree \(r>1\), and let

\[
 R_1=R+\epsilon U,\qquad R_2=R,\qquad
 U\in xk[x],\quad \deg U<r.
\]

For an outer polynomial \(H\),

\[
 H(R_1)-H(R_2)=\epsilon H'(R)U.
\]

More generally, let the projection forget every coefficient of degree at
most \(q\), and define

\[
 V_q(H,R)=
 \{U\in xk[x]_{<r}:\deg(H'(R)U)\le q\}.
\]

### Theorem A.1 (complete one-sided tangent classification)

If \(H'=0\), then
\[
 V_q(H,R)=xk[x]_{<r},\qquad \dim V_q=r-1.
\]
If \(H'\ne0\) and \(d=\deg H'\), then
\[
 V_q(H,R)=
 \operatorname{span}\{x,\ldots,x^s\},\qquad
 s=\max(0,\min(r-1,q-rd)),
\]
and in particular
\[
 \dim V_q=\max(0,\min(r-1,q-rd)).
\]

Indeed, for nonzero \(U\), degree additivity in \(k[x]\) gives
\[
 \deg(H'(R)U)=r\deg H'+\deg U.
\]

For the Hessian cutoff \(q=1\), this becomes the trichotomy

\[
\begin{array}{c|c|c}
H'&\text{form of }H&V_1(H,R)\\ \hline
0&H=G(x^p)&xk[x]_{<r}\\
a\in k^*&H=ax+G(x^p)&kx\\
\deg H'\ge1&\text{all remaining cases}&0.
\end{array}
\]

Thus the characteristic-two example is the primitive separable row:
\[
 H=z^2+z,\qquad U=x,\qquad
 H(R+\epsilon x)-H(R)=\epsilon x.
\]
The purely inseparable specialization \(H=z^2\) is larger: every normalized
right tangent is killed.

This sharpens the proposed criterion.  Divisibility \(p\mid m\) is necessary
and sufficient for a **uniform** failure somewhere in the full monic
outer-degree-\(m\) family, because \(H=x^m\) then has zero derivative.  It is
not a pointwise classification: many degree-\(m\) polynomials with \(p\mid m\)
have \(\deg H'\ge1\) and no one-sided invisible tangent.  Also, separability
of \(H\) alone is not sufficient after Hessian projection: \(H=ax+G(x^p)\)
is separable but retains the one-dimensional \(kx\) defect.

The exact checker
[`scripts/verify_positive_characteristic_ritt_infinitesimals.py`](scripts/verify_positive_characteristic_ritt_infinitesimals.py)
exhausts small finite-field cases and verifies both characteristic-two
models.

## B. Degree-forty-two support saturation

The residual problem is the eleven-variable ring with five normal variables
and six base variables from the transported \(\{2,7\}\) power chart.  With
residual ideal \(I\), the remaining support is

\[
 \mathfrak k=(w_0,w_1w_2,ABw_2).
\]

The exact target remains
\[
 (I:\mathfrak k^\infty)/I.
\]
A single regular element
\[
 f=\alpha w_0+\beta w_1w_2+\gamma ABw_2\in\mathfrak k
\]
would prove that this quotient is zero.

Two global modular colon attempts were made:

- \(p=32003\), \(f=w_0+w_1w_2+ABw_2\), `slimgb`;
- \(p=101\), the same \(f\), `std`.
- \(p=101\), the simpler candidate \(f=w_0\), `slimgb`.

The first two exceeded five minutes and the last exceeded four minutes
before producing a quotient basis.  A timeout is not evidence either for or
against saturation.

Exact untruncated characteristic-zero fiber calculations do give two useful
negative counterexample probes:

\[
(e_1,e_2,t,w_0,w_1,w_2)=(1,2,3,0,5,0)
\]
on the odd-core divisor \(w_0=w_2=0\) has an eight-element basis and
reduces the synchronization defect to zero.  The generic \(A=0\) rational
probe
\[
(1,1,3/5,0,0,1)
\]
has a nine-element basis and also reduces the defect to zero.  These fiber
calculations test the distinguished defect, not the whole saturation module,
and therefore cannot replace the global colon.

### Finite-jet embedded class

The divisor-stratified calculation has now produced an embedded class.  Put
\(\mathfrak m=(u,v)\) after the three exact unit pivots and
\[
 C_n=R/(I+\mathfrak m^n),\qquad
 T_n=((I+\mathfrak m^n):\mathfrak k^\infty)/(I+\mathfrak m^n).
\]
Over \(\mathbf F_{101}\), with all six base variables retained,
\[
 T_5=0,\qquad T_6\ne0.
\]
The standard-basis sizes are \(59/59\) at \(n=5\) and \(199/220\) at
\(n=6\).  Generator-by-generator reduction is essential here: Singular's
bulk ideal reduction can incorrectly hide the displayed quotient witness.

The first extracted class \(c_6\) lies in the fifth normal layer.  Its
cyclic annihilator has five generators and
\[
 \sqrt{\operatorname{Ann}(c_6)}
 =(u,v,w_0,w_1,w_2).
\]
After eliminating \(u,v\), the smallest base prime is therefore
\[
 (w_0,w_1,w_2),
\]
the deepest sevenfold monomial vertex.  The earlier class \(w_0^2v\) in
\((I+\mathfrak m^6):w_0\) is not relevant: its base annihilator is only
\((w_0)\), so it is not \(\mathfrak k\)-torsion.

This phenomenon is not confined to positive characteristic.  At the exact
characteristic-zero specialization
\[
 (e_1,e_2,t)=(1,2,3),\qquad A=-12,\quad B=1359,
\]
the sixth normal jet again has a strict saturation.  Both standard bases
have size \(39\), but the first nonzero generator normal form occurs at
index \(29\); its cyclic annihilator again has radical
\[
 (u,v,w_0,w_1,w_2).
\]
Thus there is an exact characteristic-zero **truncated embedded-prime
counterexample**.  The corresponding formal object is the nonzero
\((u,v,w_0,w_1,w_2)\)-primary cyclic submodule generated by \(c_6\) in the
fifth normal layer of \(C_6\).  The driver prints the full polynomial
representative; its size makes the normal-layer description more intrinsic
than copying its coefficients here.

The all-order conclusion is not yet available.  At the same specialization
modulo \(101\), an order-seven torsion class maps onto \(c_6\):
\[
 T_7\longrightarrow T_6,\qquad c_7\longmapsto c_6.
\]
The checked residual/saturation basis sizes are \(39/39\) at \(n=6\) and
\(44/43\) at \(n=7\), and the direct lift-membership test returns one.  The exact
characteristic-zero transition did not finish in fifteen minutes.  Hence
the current result is stronger than a moving-socle observation but is not
yet a proof that
\[
 (I:\mathfrak k^\infty)/I\ne0
\]
in the completed or polynomial residual algebra.

An explicit modular lift \(c_7\) was also extracted.  Direct reduction
certifies both \(c_7\ne0\) in \(C_7\) and
\[
 c_7\equiv c_6\pmod{I+\mathfrak m^6}.
\]
This is the input for rational reconstruction; it is not merely an
existence conclusion from the saturation basis sizes.

The reproducible discovery driver is
[`scripts/search_degree42_support_saturation_jets.py`](scripts/search_degree42_support_saturation_jets.py).
The principal commands are

```bash
.venv/bin/python scripts/search_degree42_support_saturation_jets.py \
  --prime 101 --normal-power 6 --analyze-witness
.venv/bin/python scripts/search_degree42_support_saturation_jets.py \
  --prime 0 --base-values 1,2,3 --normal-power 6 --analyze-witness
.venv/bin/python scripts/search_degree42_support_saturation_jets.py \
  --prime 101 --base-values 1,2,3 --normal-power 6 \
  --transition --print-lift
```

The next decisive computation is rational reconstruction of the modular
order-seven lift followed by direct exact membership, which avoids another
full characteristic-zero saturation.

## C. Three-boundary Keller suspensions

The existing three-factor Cox ledger already realizes the determinant
mechanism:
\[
 \widehat\mu:(Y\times\mathbb A^1_z)\longrightarrow
 (T\times\mathbb A^1_Z),\qquad Z=z/r_{12},
\]
with constant Jacobian and boundary lattice
\[
\ker B=\mathbb Z(1,1,1,-2).
\]
It is a finite étale torsor-Keller morphism, not a polynomial self-map of
affine space.  The obstruction is now precise: the primitive row of the
unimodular Cox ledger is realized by division by the nonconstant unit
\(r_{12}\).  Affine space has no such unit.

The smallest plausible algebraization is therefore an affine modification,
not another reciprocal coordinate.  Introduce a modification relation
\[
 z=r_{12}Z
\]
and ask whether the resulting smooth affine fourfold is actually
\(\mathbb A^4\), while the ordered-factor map and residue form extend across
\(r_{12}=0\).  A Danielewski or flexible modification can make \(Z\)
regular, but it must pass three independent tests:

1. factoriality and trivial Makar--Limanov/Derksen obstructions compatible
   with affine space;
2. no new divisor in the Jacobian ledger over \(r_{12}=0\);
3. a polynomial coordinate system, not merely stable or flexible
   equivalence.

This turns the three-boundary proposal into a concrete affine-modification
recognition problem.

## D. Non-rational critical normalization

There is an immediate obstruction for every one-parameter marked-line
chart.  If the normalized critical curve receives a nonconstant morphism
from \(\mathbb A^1\), then it cannot have genus one: the map extends to
\(\mathbb P^1\) and Riemann--Hurwitz excludes a nonconstant map to an
elliptic curve.  It also cannot land in
\(\mathbb P^1\setminus\{0,1,\infty\}\), because both \(f\) and \(1-f\)
would have to be units in \(k[t]\), hence constant.  The same unit argument
excludes a nonconstant map from \(\mathbb G_m\) to the three-punctured line.

Therefore a successful genus-one or three-puncture example must already
place nontrivial geometry in the source critical divisor; it cannot be
created by reparametrizing the existing \(\mathbb A^1\) or
\(\mathbb G_m\) critical charts.

For singular rational critical divisors the situation is different.
Polynomial normalization by \(\mathbb A^1\) is compatible with a nontrivial
conductor (cusps and nodes), so the first bounded ansatz should retain an
\(\mathbb A^1\) normalization and prescribe a conductor pair in the
incidence algebra.  The determinant identity does not see this conductor;
the finite normalization/reconstruction algebra does.  This is the most
accessible of the three nonstandard normalization searches.

## E. Support-saturated BCW minimization

Two exact audits were rerun on the current essential 21-variable cubic
witness:

1. its 65 Jacobian coefficient matrices have only the known proper row
   module, constant on the collision; the quotient action is the full
   \(M_{20}(\mathbb Q)\);
2. the 2,484 coefficient equations for an affine translation symmetry have
   full rank \(441\), so there is no nonzero affine vector-field symmetry.

Consequently no further collision-preserving **linear** quotient or affine
translation/LND slice exists.  Any improvement by support saturation must
be genuinely nonlinear.

A usable homological search object is a bounded-degree intertwining pair of
derivations
\[
 D_{\rm src},D_{\rm tgt},\qquad
 D_{\rm src}(F_i)=D_{\rm tgt}(y_i)\big|_{y=F}.
\]
If both are locally nilpotent and possess polynomial slices, their
contractible orbit directions can be quotiented.  Without the LND and slice
conditions, a local-cohomology torsion class only removes a formal or
localized variable; it does not produce a polynomial map of a smaller
affine space.

The next finite search should therefore solve the intertwining equations for
quadratic source and target vector fields, then test local nilpotence and
slices.  The completed linear and affine audits prove that degree zero and
one contain no new quotient, so quadratic degree is the first meaningful
case.
