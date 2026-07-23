# Cubic critical-boundary marking extraction

This note advances the geometric-degree-three case of the
[minimal-boundary gateway and classification conjecture](MINIMAL_BOUNDARY_CLASSIFICATION.md).
It proves that the primitive marking on a one-place cubic plane core is
automatic, constructs the exact two-place defect atlas showing why the same
argument cannot work from puncture count alone, and excludes the most direct
diagonal reciprocal lift of those defects.

The results concern coordinate-preserving plane cores.  They do not extract
such a core from an arbitrary Zariski--Main package and do not straighten the
positive threefold chart.

Work over an algebraically closed field `k` of characteristic zero.

## 1. The one-place cubic marking is automatic

### Theorem 1.1 -- unmarked one-place cubic core rigidity

Let

\[
\chi:\mathbb A^2_{w,q}\longrightarrow\mathbb A^2_{q,t},
\qquad
\chi(w,q)=(q,T(w,q))
\tag{1.1}
\]

be a coordinate-preserving polynomial plane map of generic degree three.
Assume:

1. its critical divisor is one reduced irreducible curve `E`;
2. `E` is isomorphic to `A^1`;
3. the restriction `E -> A^2_(w,q)` is the given closed embedding.

Then a triangular source automorphism preserving `q`, followed by a target
shear preserving `q` and a nonzero scaling of `t`, puts `chi` in the tangent
normal form

\[
\boxed{
(w,q)\longmapsto(q,wq-H(w)),
\qquad \deg H=3.
}
\tag{1.2}
\]

No primitive coordinate on `E` needs to be supplied separately.

#### Proof

Generic degree three gives

\[
\deg_w T=3.
\tag{1.3}
\]

Since

\[
\det D\chi=-T_w,
\]

the unique reduced critical equation `D` satisfies

\[
T_w=cD,\qquad c\in k^*,
\tag{1.4}
\]

and `deg_w D=2`.  The irreducible equation `D(w,q)=0` is the minimal
polynomial of `w` over `k(q)`.  Consequently the finite polynomial
`q|_E:E=A^1 -> A^1` has degree two.

Choose a parameter `tau` on `E`.  Write the embedded coordinate functions as

\[
q=q(\tau),\qquad w=w(\tau),
\qquad k[q(\tau),w(\tau)]=k[\tau],
\tag{1.5}
\]

with `deg q=2`.  The degree theorem in the Abhyankar--Moh epimorphism theorem
says that one of `deg q` and `deg w` divides the other.  If `deg w=1`, it is
already a coordinate.  Otherwise `deg w=2d`.  A scalar multiple of `q^d`
has the same leading term as `w`; replacing

\[
w\longmapsto w-P(q)
\tag{1.6}
\]

for a suitable polynomial `P` lowers `deg w` and preserves both `q` and the
equality in (1.5).  Reapply the degree theorem.  Every intermediate degree
greater than one is even, so iteration terminates at degree one; degree zero
would contradict (1.5).

Thus, after a triangular source automorphism preserving `q`, the restriction
of `w` is an affine coordinate on `E`.  Its equation is

\[
q=h(w),\qquad \deg h=2.
\tag{1.7}
\]

Equations (1.4) and (1.7) give

\[
T_w=c\{q-h(w)\}.
\]

If `H'=h`, integration yields

\[
T=c\{wq-H(w)\}+g(q).
\]

A target shear removes `g`, and a target scaling removes `c`.  This proves
(1.2).  QED

The theorem closes the primitive **plane-core** marking gap on the
one-place cubic branch.  It does not prove that an arbitrary minimal-boundary
Keller map supplies the coordinate-preserving presentation (1.1), nor that
its positive polynomial threefold chart is a weighted chart.

### Corollary 1.2 -- uniqueness of the cubic one-place core

Every core in Theorem 1.1 is polynomially left--right equivalent to the
plane tangent core of

\[
H_0(W)=W^2(1-W).
\tag{1.8}
\]

Indeed, translation of `w` removes the quadratic term of an arbitrary cubic
`H`; its linear and constant terms are removed by target shears, and nonzero
source and target scalings normalize the cubic coefficient.  Equivalently,
one may place the unique critical point of `H'` at the marked origin and a
second affine point at `1`, recovering (1.8).

## 2. Cubic degree does not force the two-place marking

The two-place case behaves differently even before attempting a threefold
lift.

### Proposition 2.1 -- infinite toric cubic defect atlas

For every positive odd integer `b`, put

\[
D_b(s,q)=1-s^2q^b
\tag{2.1}
\]

and define

\[
\chi_b:\mathbb A^2_{s,q}\longrightarrow\mathbb A^2_{q,t},
\qquad
\chi_b(s,q)=
\left(q,s-\frac13q^bs^3\right).
\tag{2.2}
\]

Then:

1. `chi_b` has generic degree three;
2. its Jacobian is `-D_b`, so its critical divisor is reduced;
3. the critical normalization is `G_m`, with
   \[
   s=Y^{-b},\qquad q=Y^2;
   \tag{2.3}
   \]
4. the restriction of `chi_b` to the critical curve is birational onto its
   image;
5. the two displayed unit vectors are
   \[
   v(s)=(-b,b),\qquad v(q)=(2,-2).
   \tag{2.4}
   \]

#### Proof

Differentiation gives

\[
\det D\chi_b
=-\frac{\partial}{\partial s}
 \left(s-\frac13q^bs^3\right)
=-(1-s^2q^b).
\tag{2.5}
\]

The second target coordinate is cubic in `s`, proving the generic degree.
Since `gcd(2,b)=1`, the binomial `s^2q^b-1` is irreducible.  Equations (2.3)
parameterize it.  Conversely,

\[
Y=sq^{(b+1)/2}
\tag{2.6}
\]

in its coordinate ring.  Both `s` and `q` are units there, so (2.6) and its
inverse identify that ring with `k[Y,Y^{-1}]`.

On `D_b=0`, the second target coordinate is `2s/3`.  Hence the critical
image is parameterized by

\[
\left(q,t\right)
=\left(Y^2,\frac23Y^{-b}\right).
\]

The coprimality of `2` and `b` makes this parameterization birational.
Finally (2.4) follows at the two ends of the `Y`-line.  QED

For `b>=3`, neither displayed function is a primitive generator of the unit
lattice.  The two vectors generate a saturated rank-one lattice jointly
because `gcd(2,b)=1`, but no distinguished primitive affine-linear lift has
been selected.  Thus:

\[
\boxed{
\text{cubic degree}+\text{reduced }G_m\text{ critical curve}
\not\Longrightarrow
\text{primitive cancellation marking}.
}
\tag{2.7}
\]

The maps (2.2) are plane cores, not Keller maps of `A^3`.  They do not
disprove the gateway and classification conjecture.  Instead, they show
that `PC` cannot be proved from degree, ramification, and puncture rank
alone: its witness must use the threefold boundary package and its ambient
embedding.

## 3. Reciprocal polynomiality excludes the toric defects

The atlas (2.2) cannot occur in any saturated reciprocal lift with polynomial
core outputs.  This is valuation-theoretic and does not depend on a diagonal
chart.

### Proposition 3.1 -- intrinsic reciprocal valuation obstruction

Let `V` be the boundary DVR of a reciprocal link and suppose

\[
D_b=1-s^2q^b=A^{-1},
\qquad b\ge1,
\tag{3.1}
\]

where `q` and

\[
T=s-\frac13q^bs^3
\tag{3.2}
\]

are regular core outputs.  Then no such link exists.

#### Proof

Put

\[
u=v_A(s),\qquad v=v_A(q)\ge0.
\]

Since `v_A(D_b)=-1`, the nonconstant monomial in (3.1) has valuation `-1`:

\[
2u+bv=-1.
\tag{3.3}
\]

In particular `u<0`.  The two terms of `T` have valuations

\[
v_A(s)=u,\qquad
v_A(q^bs^3)=bv+3u=u-1.
\tag{3.4}
\]

The second is uniquely lower, so it cannot cancel.  Therefore

\[
v_A(T)=u-1<0,
\]

contradicting regularity of the polynomial core output `T`.  QED

Thus the infinite plane-core atlas proves that abstract puncture data do not
select the marking, but Proposition 3.1 removes the entire atlas as soon as
the reciprocal height-one link and polynomial target core are extracted.
Other nonmonomial Laurent embeddings can have several lowest valuation terms
and are not excluded by this one-monomial argument.

The following diagonal calculation remains useful as a chart-level
regression and generalizes to arbitrary coprime monomial exponents.

### Proposition 3.2 -- diagonal toric chart obstruction

Let `a,b>=1`, put

\[
A=1+x^ay^b,
\tag{3.5}
\]

and choose integers `p,v` with

\[
ap+bv=-1.
\tag{3.6}
\]

Set

\[
s=xA^p,\qquad q=yA^v.
\tag{3.7}
\]

Then

\[
1-s^aq^b=A^{-1}.
\tag{3.8}
\]

For the primitive

\[
T(s,q)=s-\frac1{a+1}q^bs^{a+1},
\tag{3.9}
\]

the source pullbacks are

\[
q=yA^v,
\qquad
T=\frac{xA^{p-1}(aA+1)}{a+1}.
\tag{3.10}
\]

The first is polynomial exactly when `v>=0`, while the second is polynomial
exactly when `p>=1`.  These inequalities contradict (3.6).  Therefore no
diagonal reciprocal monomial link makes both preserved core coordinates
polynomial.

#### Proof

Equation (3.8) follows from

\[
s^aq^b=x^ay^bA^{ap+bv}=(A-1)A^{-1}.
\]

Substitution into (3.9) gives (3.10).  The prime `A` divides neither `x` nor
`y`, and `aA+1` is nonzero modulo `A`; the two polynomiality criteria are
therefore exact.  But `p>=1` and `v>=0` imply `ap+bv>=a>0`, contradicting
(3.6).  QED

For the cubic atlas take `a=2` and odd `b`.  Proposition 3.1 is stronger:
it eliminates every reciprocal lift with regular `q,T`, while Proposition
3.2 records the explicit monomial chart failure.

There is also a positive extraction result: for an arbitrary quadratic
cubic-core equation, polynomiality determines all three boundary valuations.

### Theorem 3.3 -- cubic reciprocal coefficient extraction

Let `V` be the boundary DVR of a reciprocal cubic suspension.  Suppose its
reduced critical equation and core primitive have been normalized to

\[
D=1+d_1s+d_2s^2,
\qquad
T=s+\frac12d_1s^2+\frac13d_2s^3,
\tag{3.11}
\]

where `d_1,d_2 in V`, `v(D)=-1`, and `T in V`.  Put

\[
v(s)=-n,\qquad n\ge1.
\]

Then

\[
\boxed{
v(d_1)=n-1,\qquad v(d_2)=2n-1.
}
\tag{3.12}
\]

If the quadratic coefficient is transversely saturated,

\[
d_2\in\mathfrak m_V\setminus\mathfrak m_V^2,
\tag{3.13}
\]

then `n=1`, `d_1` is a unit, and, with

\[
Q=-d_1,\qquad P=d_2,\qquad Y=Q-Ps,
\tag{3.14}
\]

one has

\[
\boxed{
D=1-s(Q-Ps),\qquad Y|_{D=0}=s^{-1}.
}
\tag{3.15}
\]

Thus the primitive `G_m` unit and its affine-linear ambient lift are
automatic in the transversely saturated cubic reciprocal branch.

#### Proof

Write

\[
a_i=v(d_i)\ge0.
\]

The three summands of `T` have valuations

\[
\tau_0=-n,\qquad
\tau_1=a_1-2n,\qquad
\tau_2=a_2-3n.
\tag{3.16}
\]

Since `T` is regular while `tau_0<0`, its least valuation must occur at
least twice.  It cannot be `tau_0=tau_1`: that equality gives `a_1=n`, so
the first two terms of `D` are regular; the pole `v(D)=-1` would then have
to come from `d_2s^2`, making

\[
\tau_2\le-1-n<\tau_0,
\]

a contradiction.  The case `tau_0=tau_2` is symmetric, and a three-way tie
makes every term of `D` regular.  Hence

\[
\tau_1=\tau_2<\tau_0.
\tag{3.17}
\]

Put

\[
L=a_1-n=a_2-2n.
\tag{3.18}
\]

This is the common valuation of `d_1s` and `d_2s^2`.  It must be negative,
because `D` has a pole.  If `L<-1`, those two leading terms must cancel in
`D` to raise its valuation to `-1`.  If their nonzero initial residues are
`r_1,r_2`, this says

\[
r_1+r_2=0.
\]

Regularity of `T` at the still lower valuation `L-n` says

\[
\frac12r_1+\frac13r_2=0.
\]

The two equations force `r_1=r_2=0`, a contradiction.  Therefore `L=-1`,
which is exactly (3.12).  At the leading pole of `T`, regularity imposes

\[
\frac12r_1+\frac13r_2=0;
\]

then `r_1+r_2` is nonzero, consistently giving `v(D)=-1`.

Condition (3.13) and (3.12) give `2n-1=1`, hence `n=1`; the first equation
of (3.12) makes `d_1` a unit.  Substitution of (3.14) proves the first
identity in (3.15), and restriction to `D=0` gives `sY=1`.  Since
`v(s)=-1`, this is the primitive two-place unit.  QED

Theorem 3.3 replaces two formerly separate assumptions--the primitive unit
and the affine-linear lift `Y=Q-sP`--by one intrinsic conormal condition on
the quadratic coefficient.  After (3.15), the completed reciprocal
cancellation theorem applies with `(m,r)=(1,1)`.

## 4. An exact positive cubic chart criterion

The positive threefold problem can be compressed to two conormal quotient
classes.

### Proposition 4.1 -- weighted quotient-tower criterion

Let `R=k[x,y,z]`, let `a in k`, `b in k^*`, and put

\[
\gamma=1+axy+bx^2z,\qquad
C=x\gamma,\qquad
W=(1+xy)\gamma.
\tag{4.1}
\]

Then

\[
\det\frac{\partial(W,\gamma,C)}{\partial(x,y,z)}
=b x^3\gamma^2.
\tag{4.2}
\]

Conversely, suppose a positive polynomial chart supplies elements
`W,gamma,C` and boundary equations `x,gamma` such that

\[
x=\frac C\gamma,\qquad
y=\frac{W-\gamma}{C},\qquad
z=\frac{\gamma-1-axy}{b x^2}
\tag{4.3}
\]

belong to its source ring and form a polynomial coordinate system.  Then the
chart is exactly (4.1).  In particular, it is a weighted vertical chart.

#### Proof

Substitution of (4.3) gives successively

\[
C=x\gamma,\qquad
W=\gamma+yC=(1+xy)\gamma,
\]

and

\[
\gamma=1+axy+bx^2z.
\]

This is (4.1).  Direct differentiation gives (4.2).  QED

For the cubic tangent primitive `H=W^2(1-W)`, weighted admissibility fixes

\[
a=-\frac32,
\tag{4.4}
\]

and a scaling of `z` removes `b`.  Hence (4.3) is already a uniqueness
certificate for the foundational positive chart.

The first quotient in (4.3) is forced once the intrinsic UFD divisor ledger
has identified

\[
\operatorname{div}(C)=\operatorname{div}(x)
 +\operatorname{div}(\gamma):
\]

unique factorization gives `C=cxgamma`, and a scalar normalization removes
`c`.  The unresolved positive-chart theorem is therefore narrower than an
arbitrary affine-modification classification.  It must prove that the
intrinsic boundary package selects the two congruences

\[
W\equiv\gamma\pmod C,
\qquad
\gamma\equiv1+axy\pmod{x^2},
\tag{4.5}
\]

with primitive quotients, and that those quotients complete to polynomial
coordinates.  Proposition 4.1 then performs all remaining straightening.

The coordinate hypothesis in Proposition 4.1 can itself be replaced by the
same primitive-LND and one-component-boundary conditions that complete the
reciprocal branch.

### Theorem 4.2 -- positive cubic chart straightening under intrinsic saturation

Let `R=k^[3]` and suppose a normalized cubic tangent-suspension square
supplies `x,y,gamma in R` with:

1. `x,y` algebraically independent, `x` prime, and
   \[
   R[x^{-1}]=k[x,y,\gamma,x^{-1}];
   \tag{4.6}
   \]
2. the Jacobian derivation
   \[
   \partial=\operatorname{Jac}(x,y,-)
   \tag{4.7}
   \]
   is primitive along `x=0`;
3. the boundary `V_R(x)` is not contracted by `(x,y)`, and its general
   fiber over `V(x) in A^2_(x,y)` is irreducible;
4. the source determinant ledger is
   \[
   \partial(\gamma)=b x^2,\qquad b\in k^*;
   \tag{4.8}
   \]
5. with
   \[
   C=x\gamma,\qquad W=(1+xy)\gamma,
   \tag{4.9}
   \]
   the normalized cubic tangent target quotients
   \[
   B=\frac{1+2u-3u^2\gamma}{x},
   \qquad
   A=\frac{u+u^2-2u^3\gamma}{x^2},
   \qquad u=1+xy,
   \tag{4.10}
   \]
   are polynomial.

Then there is a polynomial coordinate `z` on `R` such that, after a source
shear preserving `x,y`,

\[
\boxed{
\gamma=1-\frac32xy+bx^2z,\quad
W=(1+xy)\gamma,\quad C=x\gamma.
}
\tag{4.11}
\]

Thus the positive vertical chart is the foundational weighted chart.

#### Proof

On the localization in (4.6), equation (4.8) gives

\[
\partial=bx^2\frac{\partial}{\partial\gamma}.
\tag{4.12}
\]

Every element of `R` becomes a polynomial in `gamma` after this
localization, so a sufficiently high power of `partial` kills it.  Injectivity
of localization shows that `partial` is locally nilpotent on `R`.

The kernel argument of Propositions 8.4--8.5 in
[`LOG_GEOMETRY_OF_SUSPENSIONS.md`](LOG_GEOMETRY_OF_SUSPENSIONS.md) applies
verbatim with boundary equation `x`: noncontraction identifies
`ker(partial)=k[x,y]`, primitivity removes divisorial content, and
irreducibility of the general boundary fiber makes the plinth ideal a unit.
Hence `partial` has a slice `z` and

\[
R=k[x,y,z],\qquad \partial=\frac{\partial}{\partial z}
\tag{4.13}
\]

after scaling `z`.  Equation (4.8) integrates to

\[
\gamma=bx^2z+g(x,y).
\tag{4.14}
\]

Polynomiality of `A` in (4.10) now forces its numerator to vanish to order
two at `x=0`.  The constant term is

\[
2-2g(0,y),
\]

so `g(0,y)=1`.  Write

\[
g(x,y)=1+xg_1(y)+x^2h(x,y).
\]

The coefficient of `x` in the numerator of `A` is

\[
-2g_1(y)-3y.
\]

It vanishes exactly when `g_1(y)=-3y/2`.  Finally the source shear

\[
z\longmapsto z+h(x,y)/b
\]

absorbs the remaining `x^2h` term and gives (4.11).  The numerator defining
`B` already has zero constant term, so it introduces no further condition.
QED

Theorem 4.2 is a straightening theorem, not merely a recognition criterion.
Its hypotheses are finite intrinsic candidates: localized chart equality,
primitive conormal Fitting content, noncontraction, and boundary Stein degree
one.  The remaining positive extraction problem is to derive precisely these
labels and the target quotients (4.10) from the bare Zariski--Main package.

## 5. Updated cubic frontier

The degree-three classification now has three logically separate steps:

1. **Suspension extraction (open):** recover a coordinate-preserving
   divisor-minimal plane core from the intrinsic Zariski--Main package.
2. **Core marking:**
   - one-place branch: proved by Theorem 1.1;
   - two-place branch: Theorem 3.3 makes the primitive unit and its
     affine-linear lift automatic once the quadratic coefficient is the
     primitive conormal class; Proposition 2.1 shows why that transverse
     threefold saturation cannot be omitted.
3. **Threefold straightening:**
   - reciprocal branch: already complete once its primitive height-one link
     is extracted;
   - positive branch: Theorem 4.2 proves straightening once the intrinsic
     package supplies the primitive quotient/Stein labels in its hypotheses.

Thus the next useful theorem should use the ambient boundary algebra to rule
out the vectors `(-b,b),(2,-2)` for `b>=3`, or prove that such a defect cannot
occur in a polynomial Keller suspension.  Searching only the abstract
critical curve cannot succeed.

There is also a branch-free alternative.  The
[cubic normalization frontend](CUBIC_NORMALIZATION_FRONTEND.md) starts from
the rank-three canonical normalization algebra.  It proves flatness away
from a zero-dimensional Fitting defect and proves that intrinsic
curvilinearity of the closed collision fibers eliminates that defect.  The
equivalent test is cyclicity of the fiberwise relative cotangent module, a
unit first Fitting ideal already visible in the scheme package.
Equivalently, the primitive conormal class must generate each collision
nilradical; this is the closed-point extension of the generic conormal
marking already used branchwise.  Flatness then extracts the
Deligne--Faddeev binary cubic.  If its coefficient
morphism is affine-linear of full rank modulo polynomial Tschirnhausen
gauge and there is no second unramified boundary divisor, the
hyperplane-orbit theorem gives the foundational map directly.  This
replaces the branch markings in that subcase, but leaves cotangent
cyclicity, coefficient straightening, and exclusion of the second divisor
to be used to construct and verify the `PC`, `CS`, and `LC` parts of the
finite-normalization witness.

The Hartogs extension criterion narrows the first item further: it is enough
that the pure two-dimensional scheme-theoretic ramification support satisfy
`S_2`, its rank-one full-support relative cotangent module satisfy `S_1`,
and the primitive conormal class generate in codimension one.  Depth then
extends the generator through the closed collision points.  The two
remaining closed-point failures are exactly
`Ext_A^2(T,A)` and `Ext_A^3(Omega_{B/A},A)`.

## 6. Reproduction and external input

Run

```bash
.venv/bin/python scripts/verify_cubic_marking_frontier.py
```

The checker verifies the toric atlas, its normalization and valuation
vectors, and the diagonal lift obstruction in bounded exact ranges.  The
one-place proof uses the degree-divisibility part of:

- S. S. Abhyankar and T.-T. Moh, *Embeddings of the line in the plane*,
  J. reine angew. Math. 276 (1975), 148--166,
  [doi:10.1515/crll.1975.276.148](https://doi.org/10.1515/crll.1975.276.148).

The checker is regression evidence for the explicit identities; it does not
replace that theorem.
