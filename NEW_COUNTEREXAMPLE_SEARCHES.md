# New counterexample searches: first exact pass

This note records the first active pass on the five proposed searches.  It
separates exact conclusions from computational evidence and from new search
targets.

## Universal moment-search preflight

Before enumerating coefficients for any Gaussian, Image-Mathieu, or weighted
constant-term candidate, test whether its algebra and functional admit an
exact intertwining embedding into

\[
 \Gamma_r(u^nz^\mu)=n!\mathbf1_{\mu=0}.
\]

Christopher D. Long's
[factorially weighted multitorus theorem](extended-geometry/FACTORIALLY_WEIGHTED_MULTITORUS_THEOREM.md)
proves that every such kernel is Mathieu--Zhao, for arbitrary angular torus
rank.  This removes the entire one-factorial-radial architecture from
counterexample searches:

1. if the origin is outside the torus-weight convex hull, support separation
   forces eventual mixed vanishing;
2. if the origin is inside, one fixed power has nonzero pure moments at every
   sufficiently large prime dilation.

Finite zero prefixes cannot evade the second branch.  A retained search must
exhibit several independent radial directions, a non-prime-separating radial
functional, or a proved failure of functional/multiplicative intertwining.
Do not apply the sieve to ordinary multi-pair Gaussian or SIC contraction
solely because their angular torus has low rank: their product-factorial
radial functionals are outside its hypotheses.

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

This theorem is now part of the `OP-RITT` architecture rather than an
isolated counterexample search.  The
[Hessian--Ritt deformation complex](extended-geometry/HESSIAN_RITT_DEFORMATION_COMPLEX.md#31-the-frobenius-cell-module)
uses \(V_q(H,R)\) as a separately labelled Frobenius cell module.  Its
characteristic-\(p\) edge labels distinguish the full \(H'=0\) kernel, the
separable but Hessian-invisible \(H'=a\) row, and the ordinary zero-defect
row.  The same note records the additional denominator, pivot, tame-degree,
and no-Frobenius conditions needed before a characteristic-zero
synchronization certificate may be claimed to survive reduction.

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

The next decisive computation was initially expected to be rational
reconstruction of the modular order-seven lift.  The subsequent search
replaced it by the simpler single-colon test below.

The reconstruction search sharpens this substantially.  The fixed integer
characteristic-zero witness reduces to an order-seven class represented by
the same normal-form polynomial in characteristic \(101\), and
\[
 c_6\ne0\pmod {I+\mathfrak m^7},\qquad
 w_0c_6=0\pmod {I+\mathfrak m^7}.
\]
Both the original `slimgb` basis and an independently re-standardized basis
give the same annihilation result.  Equivalently, the transition problem is
the single-colon membership
\[
 c_6\in (I+\mathfrak m^7:w_0)+(I+\mathfrak m^6).
\]
This modular calculation has basis sizes \(39\) at order six and \(44\) at
order seven; the \(w_0\)-colon has size \(24\), and all four direct checks
(lift, nonvanishing, transition, and annihilation) return one.

This also explains why naive coefficientwise CRT failed.  Different primes
select different representatives in the nontrivial transition kernel; even
the pure sixth-layer corrections are not canonical.  The correct exact
target is therefore the displayed single-colon identity, not reconstruction
of arbitrary saturation representatives.

Two characteristic-zero direct attempts were made: first with the redundant
order-six basis and then with only the order-seven basis and the fixed
integer polynomial \(c_6\).  Both hit a twenty-minute cap.  Thus order seven
is still modular rather than an exact characteristic-zero theorem.  The
remaining exact computation is the fixed direct identity
\(w_0c_6\in I+\mathfrak m^7\) and no longer requires support saturation or
an unknown lift.  The direct-identity regression is
[`scripts/verify_degree42_order7_known_witness.py`](scripts/verify_degree42_order7_known_witness.py).
The colon-transition and graded-correction diagnostics are
[`scripts/verify_degree42_order7_colon_lift.py`](scripts/verify_degree42_order7_colon_lift.py)
and
[`scripts/search_degree42_order7_graded_lift.py`](scripts/search_degree42_order7_graded_lift.py).

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

## F. LND ideal-image Mathieu search

This search is attached to the Image/Mathieu chain, but it is not a
compression of the known two-pair counterexample.  For a global slice
\(D=\partial_s\), the exact primitive criterion closes three large classes:

1. all zero-dimensional carrier-free ideals, including nonreduced schemes;
2. ideals \(I=qJ\) with \(q\) monic of positive \(s\)-degree and \(B/J\)
   finite-dimensional, including nonreduced \(J\);
3. \(q=x^c(xs-1)^d\) with the finite residual scheme supported at its
   degree-drop point.

The proofs and exact replays are in
[`LND_MATHIEU_SLICE_CONDUCTOR_FRONTIER.md`](extended-geometry/LND_MATHIEU_SLICE_CONDUCTOR_FRONTIER.md).
Thus the old broad instruction “search small \(qJ\)” is too wasteful.  The
first exact pass should use only:

\[
 \deg_s q>0,\ q\text{ nonmonic},
\]

excluding the normalized repeated carrier \(x^c(xs-1)^d\).  Support \(J\)
over a multiple or reducible vanishing divisor of the leading
\(s\)-coefficient.  This is where roots escape the monic slice and where
valuation/conductor saturation can still change the image.

The bounded compiler should output, for each primary chart:

1. the primitive-membership matrix modulo \(J\);
2. the leading-coefficient and carrier valuations;
3. the normalization branches and conductor quotient;
4. the exposed support faces of \(f\);
5. pure and mixed sequences in recurrence-ready form.

Reject a chart immediately if it has a strict support separator, only
reduced interval moments, or a unique valuation face.  Promote nothing
without an all-order primitive certificate and an infinite mixed
nonmembership certificate.

The first exact bounded compiler
`scripts/search_lnd_nonmonic_degree_drop.py` audits the normalized carriers
\(x^c(xs-1)^d\) for \((c,d)=(0,1),(0,2),(1,1)\), three small primary
residual schemes, and 256 sparse seeds.  Its primitive-carrier assertion is
a regression for the support-weight theorem; every finite search result is
only a replay or candidate audit.

The next compiler, `scripts/search_lnd_plinth_ideal_images.py`, crosses the
plinth divisor of \(D=x\partial_y+y\partial_z\) by constructing each
homogeneous space \(D(I_n)\) exactly.  Five ideals and 45 sparse seeds leave
thirteen pure-prefix survivors, all divisible by the visible plinth
coordinate \(x\), and no bounded mixed-tail obstruction.

The exact reducible-plinth compiler
`scripts/search_lnd_reducible_plinth.py` then treats

\[
 D=x(x-1)\partial_y+y\partial_z,
\]

using the weight grading \((0,1,2)\) and Singular standard bases over
\(\mathbf Q[x]\).  Five ideals coupling \(x=0\) and \(x=1\) and 205
branch-aware seeds leave fourteen pure-prefix survivors.  They are exactly
the displayed support-cone forms \(x(x-1)y^r\), \(1\leq r\leq3\), and none
has a bounded mixed-tail obstruction.  The same compiler now also treats

\[
 D=x(x-1)\partial_y+y^2\partial_z,
\]

with grading \((0,1,3)\).  Five coupled ideals and 210 seeds leave
seventeen pure-prefix survivors, all \(x(x-1)y^r\) for
\(1\leq r\leq4\), and again no bounded mixed-tail obstruction.  The first
branch-asymmetric target

\[
 D=x(x-1)\partial_y+(y^2+x)\partial_z,
\]

admits an exact normalized-primitive compiler despite breaking the grading.
Across eight ideals, including three mixed-weight jet charts, and 350
mixed-weight seeds, it leaves 48 chartwise survivors representing only
seven distinct forms.  All seven are divisible by \(x(x-1)\), and none
has a bounded mixed-tail obstruction.

Two disjoint points on the invariant \(x\)-line therefore do not provide
enough conductor interaction.  The crossing-plinth LND

\[
 D=uv\partial_y+(y^2+u)\partial_z
\]

on \(k[u,v,y,z]\) makes the plinth components \(u=0\) and \(v=0\)
intersect.  An exact normalized-primitive search over \(\mathbf Q[u,v]\)
tests five zero-dimensional crossing ideals and 956 seeds.  Its 177
chartwise survivors represent 36 distinct forms, all divisible by \(uv\);
no mixed-tail obstruction survives, and every non-\(uv\) seed already
fails by its square.  This removes the simplest crossing-divisor route.

The nonprincipal plinth ideal begins with

\[
 D=u\partial_x+v\partial_y,\qquad \ker D=k[u,v,uy-vx].
\]

Here the local slices \(x/u\) and \(y/v\) have a genuine overlap cocycle
along the codimension-two plinth locus \((u,v)\), rather than a single
product denominator.  The exact finite-residual audit tests five ideals
and 1,055 seeds.  It leaves 494 chartwise survivors representing 100
distinct forms, all in \((u,v)A\), and no mixed obstruction.  Every seed
outside that support cone fails at its first power.

The positive-dimensional continuation begins with \(I=(u,v,x)\).  Here
there is an exact identity

\[
 D(I)=\operatorname {Im}D\cap
 \ker\bigl(h\mapsto[v^1]h(0,v,0,y)\bigr).
\]

More generally, with \(\mathfrak p=(u,v)\), the base-degree grading gives

\[
 \operatorname {Im}D\cap\mathfrak p^N
 =D(\mathfrak p^{N-1}).
\]

The known rank-two linear-derivation theorem for \(\operatorname {Im}D\)
therefore proves all-order that \(D(I)\) is Mathieu--Zhao for every ideal
containing a power of \(\mathfrak p\).  This includes
\((u,v,\ell(x,y))\) with \(\ell\ne0\) and every nilpotent jet below.

Four positive-dimensional nilpotent jets were then searched exactly
without truncating their free \(y\)-direction.  Across all five charts,
1,055 seeds leave 289 chartwise survivors representing 75 forms, all in
\((u,v)A\), and no mixed obstruction.  The theorem independently proves
all five charts safe.  The live counterexample target is now \(I=(x)\):
modulo \(x\), the invariant image is \(k[u,v,uy]\), has positive module
rank, and is not controlled by plinth-power saturation.

For this ideal, primitive membership is the exact support inequality
\(a\ge c\) on every residue monomial \(u^av^by^c\).  The first
valuation-face search tests 1,055 seeds.  Forty-eight survive six pure
powers, all divisible by \(u\), with no mixed-tail obstruction; every
failed seed dies by its square.  Finite Neumann-series inversions prove
all-order safety when \(\nu_u(f)>\deg_yf\), or when
\(\nu_v(f)>\deg_xf\) and the \(u\)-minus-\(y\) support weight is positive.
The local-slice argument generalizes further: every
\(f=u^rF\), \(r\ge1\), with \(F\in\ker(D)[x]\), has the mixed-power
conclusion.  All 48 survivors lie in \(u\ker(D)[x]\), including the three
affine forms left by the strict cones, so the entire survivor list is
proved safe to all orders.

The suggested square gate is false after an algebraic scalar extension.
For
\[
 p(t)=2t-1+\frac{\sqrt{-15}}3(6t^2-6t+1),
\]
the first two moments on \([0,1]\) vanish but the third does not.
Homogenizing and substituting \(z=uy\), \(w=uy-vx\) gives an exact
\(f\) with \(f,f^2\in D((x))\) and \(f^3\notin D((x))\).  This falsifies
the finite shortcut, not LNED.

The full eventual-power hypothesis does close the ideal.  In the local
slice \(A\subset k[u,v,uy-vx]_u[x]\), a nonpositive lowest \(u\)-face
would make all sufficiently large interval moments of one nonzero
polynomial vanish.  Passing to a fixed sufficiently high power and using
the one-variable polynomial moment lemma gives a contradiction.  The
lowest \(u\)-order is therefore positive, putting \(f\) in
\(u\ker(D)[x]\), where the existing local-slice theorem supplies every
mixed tail.  Hence \(D((x))\), and by simultaneous linear change
\(D((\ell(x,y)))\) for every nonzero linear form \(\ell\), is
Mathieu--Zhao.  Divisibility bootstrapping then proves the same for every
thickening \(D((\ell^d))\), \(d\ge1\).

The first reducible target \(I=(xy)\) also closes all-order, for a
different reason.  Over the invariant fraction field, \(xy\) has the two
distinct local-slice roots \(0\) and \(-(uy-vx)/v\).  A primitive in
\((xy)\) vanishes at both, so eventual pure-power membership gives all
large interval moments between those roots.  The one-variable moment
lemma forces the seed to be zero.  More generally, every principal carrier
with two distinct roots on the generic additive-group orbit has zero
Mathieu radical.  This also closes \(I=(x(x-1))\).

The first invariant-content controls also close.  Since
\(D(uxA)=uD(xA)\), the local-slice proof for \(D((x))\) tolerates the
single fixed multiplier denominator \(u^{-1}\).  In the crossed case,
pure membership in \(uD(yA)\) forces
\(f\in uA\cap v\ker(D)[y]=uv\ker(D)[y]\); its \(u^mv^m\) content clears
both fixed denominator losses.  By symmetry,
\[
 D((a\ell))\text{ is Mathieu--Zhao}
 \quad(a\in\{u,v\},\ \ell\in\{x,y\}).
\]

The first genuinely rational ladder also closes.  For
\(q_n=u^nx+(uy-vx)\), \(n\ge1\), normalization of \(A/(q_n)\) gives the
exact support cone \(u^av^bt^c:a\ge(n+1)c\).  The apparent second-face
obstruction \(y=v(q_n-w)/u^{n+1}+w/u\) invalidates divisibility of an
arbitrary coefficient, but the actual lowest \(u\)-face retains a
\((q_n-w)\)-factor.  Eventual pure powers therefore force
\(f\in u\ker(D)[q_n]\) by the one-variable moment lemma, and this content
clears every fixed mixed-tail denominator.  Hence every
\(D((q_n))\) is Mathieu--Zhao.  For \(q_1\), the exact bounded census has
seventeen survivors on powers \(4,5,6\), all in the proved safe cone.

The first tied family closes as well.  For
\(q_{r,s}=u^rv^sx+(uy-vx)\), \(r,s\ge1\), quotient normalization gives
the simultaneous exact cone
\[
 a\ge(r+1)c,\qquad b\ge sc.
\]
Independent lowest-face moment arguments force positive \(u\)- and
\(v\)-content, hence \(f\in uv\ker(D)[q_{r,s}]\); the resulting
\((uv)^m\) clears every fixed mixed denominator.  For \(q_{1,1}\), eight
of 1,055 seeds survive powers \(4,5,6\), all in the proved safe cone.

The intercept can in fact be arbitrary: replacing \(w\) by any
\(b\in\ker(D)\) merely replaces the moment endpoint \(w\) by \(b\).
The zero-intercept case reduces to invariant content times the
moving-linear theorem.  Hence every
\[
 q=u^rv^sx+b,\qquad r\ge1,\ s\ge0,\ b\in\ker(D),
\]
is safe.

Generic nonmonomial slopes close prime by prime.  For coprime
\(q=ax+b\), every prime factor of \(a\) forces growing content, and the
\(u\)-face does the same unless
\[
 v(b\bmod u)=(a\bmod u)w.
\]
Off this locus, \(D((q))\) is Mathieu--Zhao.  On it, \(q\in uA\) despite
\(\gcd(a,b)=1\) in \(B[x]\): this is a genuine conductor cancellation.

The next principal target is therefore the aligned conductor-descent
locus, plus removal of nontrivial common invariant content.
