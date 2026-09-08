# Universal RR genus-two family: descent-object preflight

## Decisive correction: all smooth RR members are already rationally pointed

**New deduction, independently verified application.** For this entire net,
not only the retrospectively selected first-witness pencil,

\[
\boxed{\mathcal C_{u,v}(\mathbf Q)\ne\varnothing
\quad\text{for every smooth geometrically integral rational member}.}
\]

The point is inherited from a generic section. This rules out **raw curve
solubility**, everywhere-local solubility, or nonemptiness of the curve's
fake2-Selmer set as a discriminator between302 and control RR members.
It does not rule out a *specified residual class* or an independent-point
criterion. No numerical Selmer panel is claimed.

**Subsequent verified application.** The
[first-unlock Jacobian certificate](DET1092_FIRST_UNLOCK_JACOBIAN_CLASS_2026-09-08.md)
now proves that the rational first witness, based at this inherited point,
is independent of the full rank17 NS-restriction image on its retrospective
RR member. Their combined Jacobian subgroup has rank18. This closes that
specific non-generic rational-class test, not the blind control comparison.

**Subsequent control result.** The
[generic-point specificity panel](DET1092_RR_GENERIC_POINT_SPECIFICITY_2026-09-08.md)
now gives an explicit generic-input family indexed by the original parameter,
`C_t: s^2=c*q(T;r0(t),0)`. All eight V4 controls and a generic-point302 control
have the same non-generic rational Jacobian-class property, while their
marked elliptic points are in MW17. This supplies a precise obstruction to
using that property alone as an elliptic-extra-point certificate.

**Proof from the existing divisor data.** On the smooth K3 surface,

\[
D=2O+5F+\phi(w)=O+P_w,\qquad
D\cdot O=2(-2)+5=1.
\]

A smooth integral member does not contain `O`. Its intersection with `O`
is an effective divisor of degree1 defined over Q, hence a rational point.
This uses only the generic net, not the exceptional302 point. In fact the
44 distinct generic sections in the already certified22 pairs all have
`D.P_x=1`: they give44 marked rational intersections, distinct generically
but potentially coincident at special members. No enumeration is repeated.

**Explicit uniform construction.** Retain the exact `A_i,B_i,h,c,q` below.
Polynomial division, independently checked, gives

\[
A_2=\alpha h,\qquad B_2=(\beta_0+\beta_1T)h,
\quad \alpha\ne0,\quad\deg h=3.
\]

Set

\[
u_*=-\beta_0/\alpha,\quad v_*=-\beta_1/\alpha,
\qquad U=u-u_*,\quad V=v-v_*.
\]

There is an explicitly recorded nonzero rational constant `g` such that

\[
B_1+(u_*+v_*T)A_1=g h^2.
\]

The member at `(u_*,v_*)` is exactly the already certified reducible
`O+P_w`; its vertical-chord coefficients agree with the old certificate.
Writing `h_hom(T,Z)=Z^3 h(T/Z)`, the uniform weighted-projective point is

\[
\boxed{[T:Z:s]=[-U:V:\ \pm g^2h_{\rm hom}(-U,V)]}
\quad\text{in weights }(1,1,3).
\]

The only base parameter of this formula is `(U,V)=(0,0)`, the excluded
reducible member. In particular it also handles `V=0` and points above
`T=infinity`; no unwritten finite-denominator assumption is needed.

**Exact identity and replay.** Put
`t0=-(beta0+alpha*u)`, `z0=beta1+alpha*v`,
`H=h_hom(t0,z0)`. Homogenizing `f1` to degree6 gives

\[
f_{1,\rm hom}(t_0,z_0)=gH^2,\qquad
cq_{\rm hom}(t_0,z_0;u,v)=g^4H^2.
\]

All coefficients, the constant `g`, the unique excluded parameter, and
an exact affine point on the formerly source-selected member `u=v=0`
are in the
[uniform-point certificate](../../artifacts/generated-results/elliptic-curves/det1092_rr_inherited_point_v1.json).
The [independent replay](../../artifacts/generated-results/elliptic-curves/det1092_rr_inherited_point_replay_v1.json)
reconstructs the centre from the generic17 sections, checks the complete
sextic identity, checks the homogeneous point coefficientwise, and identifies
the excluded member with the earlier vertical-chord construction. It does
not import the constructor or use any exceptional point coordinates.
Both bounded jobs finished in about one second under25-second caps.

**Established descent interpretation.** A rational point `P0` embeds the
curve in its Jacobian by `P -> [P-P0]`. Pulling back multiplication by2
gives a two-cover containing the rational lift `(P0,0)`. Thus the curve's
two-Selmer set, and its fake two-Selmer set, are nonempty. This proves
nonemptiness without computing either set; see
[Bruin--Stoll](https://arxiv.org/abs/0803.2052) and the
[two-cover definition](https://magma.maths.usyd.edu.au/magma/handbook/text/1619).
The full Jacobian2-Selmer group, other cover classes and Cassels--Tate data
are **not computed** by this argument.

**Evidence correction.** The previous double-plane note left rational and
everywhere-local solubility of `C_(0,0)` unknown. Both are now proved,
uniformly for the net. Earlier immutable certificates remain historical
evidence of what those computations checked; their negative statements
about elliptic quotients are unchanged. An elliptic quotient was never a
necessary prerequisite for testing an empirical arithmetic correlation.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_rr_inherited_point.sage
```

## What is constructed, and what is not yet a valid blinded panel

**Verified application.** The universal genus-two family from the calibrated
historical RR net is now explicit:

\[
\mathcal C_{u,v}:s^2=c\,q(T;u,v),\qquad\deg_Tq=6.
\]

All65 nonzero coefficients of the integral polynomial `q`, the rational
squareclass factor `c`, and exact maps are supplied in the
[construction certificate](../../artifacts/generated-results/elliptic-curves/det1092_universal_rr_descent_preflight_v1.json).
The [independent replay](../../artifacts/generated-results/elliptic-curves/det1092_universal_rr_descent_preflight_replay_v1.json)
checks the full identity, not an interpolation at the first point.
The equation construction reads only generic parent/RR data; the centre
itself is calibrated from the historical run. Sealed point information is
read only afterward for the labelled retrospective class diagnostic.

**Updated comparison boundary.** The full net is indexed by RR parameters
`(u,v)`, not automatically by the original elliptic-fibre parameter. The new
generic-input rule `(u(tau),v(tau))=(r0(tau),0)` does define a control family,
and its local/Jacobian/rational-Kummer comparison is complete. It marks a
known generic elliptic point, so it is not an extra-point predictor. Full
residual2-Selmer groups and a source-selected class sufficient to encode an
independent point on `E_tau` remain unconstructed. These are distinct from
the now completed specificity control.
For the raw solubility predicate the theorem above is stronger than such a
panel: every smooth rational assignment would give the same positive answer.

## Three different objects that must not share one label

**Verified application.** For clarity, use `tau` for a fixed original fibre,
`T` for the variable along a multisection, and `(u,v)` for its RR parameters.

| Object | Dimension/genus | What its rational points mean |
|---|---|---|
| `w^2=F_tau(lambda)` | genus1 | Points of the already pointed elliptic curve `E_tau` |
| `s^2=c*q(T;u,v)` for fixed `(u,v)` | genus2 generically | Points on the surface, usually on different original fibres `E_T` |
| `s^2=c*q(tau;u,v)` for fixed `tau,u,v` | zero-dimensional degree2 fibre | Splitting of this particular multisection above `tau` |

Freezing a genus-two member leaves one Jacobian, not a varying Jacobian
`J_tau`. Varying `(u,v)` varies genus-two Jacobians, but must be specified
without the exceptional witness if the comparison is to be blind.
Existence of a point on a multisection somewhere also does not imply a point
above the specific marked parameter `tau`.

**New deduction: the marked square condition is a different test.** Fix a
smooth member and a finite nonbranch inherited point `(T0,s0)`. In the
sextic etale algebra `A=Q[theta]/(q(theta))`, for finite `tau` with
`q(tau)!=0`, set

\[
\beta_\tau=\frac{\tau-\theta}{T_0-\theta},\qquad
N_{A/\mathbf Q}(\beta_\tau)
=\frac{q(\tau)}{q(T_0)}=\frac{c q(\tau)}{s_0^2}.
\]

Hence square norm here is exactly the already explicit splitting condition
`c*q(tau)` is a rational square. For this *zero-dimensional* marked fibre,
solubility in every completion is equivalent to rational solubility: a
nonzero rational number square everywhere has even prime valuations and
positive sign, hence is a rational square. This does not give a Hasse
principle for genus-two curves or for their two-covers. Nor does norm square
prove that the resulting elliptic point is independent of the inherited
MW17 subgroup. The sextic algebra is not automatically the cubic algebra
of the elliptic fibre or the earlier MW16 ideal-class construction.

## Exact universal equation and maps

**Verified application.** With the generic norm10 centre `C=P_w`, short
coordinates `(c_x,c_y)`, short coefficient `a=-c4/48`, and the saved RR
basis `A,B`, put

\[
f_i(T;u,v)=B_i(T)+(u+vT)A_i(T),\qquad
m=-f_1/f_2+a_1/2,
\]

and let `h^2=den(x_C)`. The complete checked identity is

\[
f_2^4\left(m^4-6c_xm^2-8c_ym-3c_x^2-4a\right)
=h^6c\,q(T;u,v).
\]

The specialization `(u,v)=(0,0)` has a squarefree sextic, proving generic
genus2 on a nonempty open parameter set. The scalar `c` is retained: dropping
it could change rational solubility even though geometric genus is unchanged.

Off the displayed denominators, `s^2=c*q` gives

\[
W=h^3s/f_2^2,\quad X=(m^2-c_x+W)/2,
\quad Y=m(X-c_x)-c_y.
\]

The inverse short-model transport gives a point of `E_T`. The varying `T`
in that last statement is essential. This is not a map to a fixed `E_tau`.

## Why pointed-quartic solubility cannot be the discriminator

**Verified application and established interpretation.** The universal
fibrewise quartic is

\[
F_t(\lambda)=\lambda^4-6c_x(t)\lambda^2-8c_y(t)\lambda
-3c_x(t)^2-4a(t).
\]

Its discriminant is `256*Delta_E`. Its leading coefficient is1, giving the
two rational weighted-projective points `(1:0:+/-1)` at infinity. Consequently
every smooth specialization is rationally, and hence everywhere locally,
soluble, on302 and on every smooth determinant1092 control. This is an
exact uniform statement; no finite-prime test is being promoted to an
everywhere-local theorem.

Its binary-quartic invariants satisfy `I=-48a`, `J=-1728b`. Thus the usual
Jacobian coefficients `(-27I,-27J)` are `(6^4 a,6^6 b)`, and the Jacobian is
Q-isomorphic to `E_t`. These identities are independently checked.

The degree-two marking is still meaningful. Under the identification with
`E_t`, it is the divisor `O+C`, and the associated2-covering map has the form
`P -> 2P-C`. Its Kummer class is the already inherited `delta(C)`; its
underlying torsor class in `H1(Q,E_t)` is zero, and its class modulo
`delta(M17)` is zero. The **new point** can give a new Kummer class, but that
is not the same as the already known covering class of its search chart.
See the repository's
[pointed-quartic proof](POINTED_QUARTIC_CENTRE_CLASSES_2026-09-07.md) and
[Fisher's discussion of genus-one models with degree-n divisor classes](https://www.dpmms.cam.ac.uk/~taf1000/papers/minbqtc.pdf).

## The first witness is not strict modulo MW17

**Retrospective verified application.** The exact first historical point
has public31-basis word

```text
-1,2,0,0,-1,0,-1,-2,0,0,0,1,0,0,0,1,0,-1,0,0,1,0,-2,1,-1,0,1,0,0,0,0.
```

Let `K` be the previously certified kernel of all20 bad-place and real
Kummer characters. It has dimension10; the generic image `M17` has
dimension17 and intersects `K` trivially. The checker verifies an explicit
linear functional annihilating `M17+K` but taking value1 on the first word.
Therefore

\[
[P_{\rm first}]\notin M_{17}+K.
\]

No translation by a generic MW17 point makes this first direction strict
at the stated places. This does not contradict the existence of the later
strict directions in the completed recovery; it rules out identifying this
particular first direction with one of them under the same local definition.

The explicit elliptic Kummer element is also recorded. For literal302
coordinates `(x,y)`, let `X=4x`, `Y=8y+4x+4`, and write
`X=n/d^2`, `Y=b/d^3`. In the retained cubic algebra

\[
K_3=\mathbf Q[\theta]/\bigl(\theta^3+5\theta^2+
(16a_4+8)\theta+64a_6+16\bigr),
\]

the element `beta=n-d^2 theta` has norm `b^2`. The checker recomputes that
resultant exactly. This is a concrete Kummer class, but not an ideal-class
identification with a separate MW16 reference. Such a comparison needs an
exact algebra map and compatible local conditions.

## What descent can establish

**Established literature.** Genus-two curve2-cover descent uses the etale
algebra of its sextic (and a Selmer *set* of the curve), not automatically
the cubic algebra of an elliptic fibre. An empty fake2-Selmer set proves no
rational point; a nonempty set does not prove a rational point. See
[Bruin--Stoll, Two-cover descent on hyperelliptic curves](https://arxiv.org/abs/0803.2052)
and the [documented implementation](https://magma.maths.usyd.edu.au/magma/handbook/text/1619).

Similarly, a Jacobian2-Selmer group is not a certificate that its curve has
a rational point. Cassels--Tate pairings can obstruct specified classes;
surviving classes can still be insoluble. The existing
[soluble-versus-Sha panel](EXCEPTIONAL_SOLUBLE_VS_SHA_PANEL_2026-09-05.md)
already proves this distinction and is not rerun here.

**Unverified hypothesis.** A source-constructed residual arithmetic class
could discriminate rank incidence. Testing that requires frozen class/cover
selection outside the inherited subgroup, an explicit map to the claimed
elliptic point, and independent rank certification. Control search misses
must remain **unknown arithmetic rank**, not labels asserting no extra point.
Otherwise the proposed comparison still measures bounded visibility.

## A fixed-fibre repair: an exact branch-divisor obstruction

**New deduction from established double-cover theory, verified on the parent.**
One concrete way to obtain a genus-two curve mapping to the *marked* elliptic
fibre would be to use the RR intersection pair as its branch divisor. This
particular repair is now decided, rather than left as an unspecified selector.
It is not a classification of all possible maps from genus-two curves to the
elliptic fibre.

Write `C=P_w(t)` for the generic centre. A reduced RR intersection pair has
the form

\[
D=P+(C-P),\qquad [D]=[O+C]\quad\hbox{in }\operatorname{Pic}^2(E_t).
\]

For any field `k` of characteristic zero containing these data,

\[
\boxed{\begin{aligned}
&\text{a smooth degree-two genus-two cover of }E_t
  \text{ branched exactly at }D\text{ exists over }k\\
&\hspace{35mm}\Longleftrightarrow C\in 2E_t(k).
\end{aligned}}
\]

**Established literature and proof.** For a double cover, the trace-zero
summand of the pushforward of its structure sheaf is a line bundle `L^-1`,
and multiplication gives `L^2 = O(D)`. Conversely this datum constructs the
cover, which is smooth for reduced `D` and geometrically connected because
it ramifies. These are the degree-two building data of
[Pardini's abelian-cover construction](https://eudml.org/doc/153330).
Here `deg L=1`. Since the elliptic curve has the rational origin `O`, every
rational degree-one line-bundle class is represented by a rational point
`Q`. The square-root condition is exactly `2Q=C`. Riemann--Hurwitz gives
`2g-2=2`, hence genus two. A constant twist cannot repair the absence of this
rational line bundle.

**Verified application.** For the historical centre at302, the half quartic
is not merely rootless modulo17: it is irreducible there,

\[
H_0(m)\bmod17=m^4+2m^2+8m+14.
\]

Independent Frobenius checks give
`gcd(H_0,m^(17^2)-m)=1` and `m^(17^4)=m mod H_0`. Its rational coefficients
are17-integral and it is monic, so it is irreducible over `Q` as well.
Every half of this centre therefore has degree four. The proposed branched
genus-two cover does not exist over `Q`, or over any extension of degree
less than four. This conclusion holds for every reduced RR pair in this
centre's degree-two linear system, independent of the chosen member.

### Explicit construction over the half field

**New deduction, with exact maps independently checked.** Adjoin a root
`mu` of `H_t`, and set

\[
q_x=(\mu^2-c_x)/2,\qquad q_y=\mu(q_x-c_x)-c_y,
\qquad Q=(q_x,q_y),\quad 2Q=C.
\]

Let `n` be the short-model slope of the selected RR chord through `-C`.
Put

\[
e=q_x+\frac{2q_y}{n-\mu},\qquad f(X)=X^3+aX+b.
\]

Over the half field, a genus-two cover with the required branch pair is

\[
\mathcal B_{t,n,\mu}:y^2=f(z^2+e)
=z^6+3ez^4+(3e^2+a)z^2+f(e),
\qquad (z,y)\longmapsto (z^2+e,y)+Q\in E_t.
\]

The chord identity after translation by `Q` is

\[
n(R+Q)=\mu+\frac{2q_y}{x(R)-q_x}.
\]

Thus the two branch points `z=0` map to exactly the desired RR pair. The
coefficient replay proves

\[
(n-\mu)^4f(e)=q_y^2H_t(n),\qquad
\operatorname{disc}_z f(z^2+e)
=-64f(e)(-4a^3-27b^2)^2.
\]

The displayed affine construction excludes singular elliptic fibres,
`q_y=0`, `n=mu`, `H_t(n)=0`, and poles of these charts. In particular it
does not count a repeated RR intersection as a reduced branch divisor.

**New deduction.** Its Jacobian is isogenous over the half field to `E_t`
times the explicit second elliptic curve

\[
E': V^2=U^3+(3e^2+a)U^2+3ef(e)U+f(e)^2.
\]

The second quotient map is `U=f(e)/z^2`, `V=f(e)y/z^3`. Together with the
first quotient, the pulled-back invariant differentials are multiples of
`dz/y` and `z dz/y`, a basis of the genus-two differentials. Hence the
induced map of abelian surfaces is an isogeny. Both quotient equations are
checked exactly; no Jacobian rank or Selmer group is inferred from this.

**Verified application and limitation.** The untwisted sextic is monic and
has two rational infinity points over the half field. Both map to `Q`,
whose double is the inherited centre. They do not give an independent
Mordell--Weil direction, and their existence does not explain the first
extra rational point at302. This genus-two construction is different from
the RR multisection `C_(u,v)` at the top of this note.

### Frozen local comparison:302 and ten unchanged V1 controls

**Verified application, calibration-informed design with outcome-blind
arithmetic execution.** The new
[protocol](../../artifacts/generated-results/elliptic-curves/det1092_rr_fibre_cover_gate_protocol_v1.json)
was emitted before its arithmetic run. It takes exactly the eleven entries
of the completed V1 roster, transports their parameters to the original
parent, uses the one historical generic centre, and tests all eighteen
previously fixed primes through197. No control is replaced and the prime
list is not extended. The arithmetic worker uses hashed case labels and
does not read ranks, search results, or exceptional point coordinates.
The designer already knows302 and the previous aggregate control outcomes;
this is not a claim of fully unseen experimental design.

| Fixed population | Exact outcome for this branch-cover construction |
|---|---|
|302 | Obstructed over `Q_17`, among other retained places |
|All ten V1 determinant1092 controls | Each has a recorded local obstruction |

The [full exposure certificate](../../artifacts/generated-results/elliptic-curves/det1092_rr_fibre_cover_gate_v1.json)
retains all198 prime/case trials, including coefficient-denominator
exclusions and positive residue roots. Positive residue roots are labelled
`ROOT_MOD_P_ONLY`, not everywhere-local solubility. The
[independent replay](../../artifacts/generated-results/elliptic-curves/det1092_rr_fibre_cover_gate_replay_v1.json)
reconstructs the centre from the generic seventeen sections, verifies every
control's parameter transport and elliptic `j`, and checks every root count
using Frobenius gcds instead of the constructor's residue enumeration.

**New deduction and scope.** This proposed fixed-fibre genus-two repair
fails its required positive case,302, as well as all controls. It therefore
cannot provide the desired arithmetic discriminator. This is not a failure
of genus-two descent: the specified cover does not exist over `Q` to begin
with. Nor is it a non-Selmer certificate for `delta(C)`: that known Kummer
class comes from a rational point and already lies in the elliptic Selmer
group. The obstruction says that `C` has no rational half, not that its
pointed genus-one covering is insoluble. Other branch divisors, higher-degree
maps and other correspondences remain open.

## Computation and next gate

The subsequent [double-plane and Jacobian gate](DET1092_RR_PLANE_AND_JACOBIAN_GATE_2026-09-08.md)
identifies the elliptic fibres as lines through the branch node and the
genus-two members as lines missing it. The source-selected member `B` has
an absolutely simple Jacobian, witnessed by exact reduction at17. So does
the geometric generic RR member: neither has an elliptic quotient of any
degree. Special split-Jacobian loci and other correspondences remain open;
this is not a genus-two Selmer or rational-solubility certificate.

One universal RR sextic and one symbolic fixed-fibre bielliptic construction
were checked. The bounded branch-divisor supplement ran198 small-prime
trials on the unchanged V1 roster. No Selmer group, global unit/class group,
Cassels--Tate matrix, new control point search, or pilot mutation was run.
A first independent implementation hit its25-second cap while reducing
multivariate rational functions; clearing the known denominators avoided
that unnecessary operation, and the corrected coefficient replay passed.
The branch-divisor replay initially had the wrong sign in the multiplier
of a cleared chord identity; that assertion failed before any replay
certificate was written. Correcting the multiplier to `-y+qy-m*(x-qx)`
gave the exact identity and a passing replay. All supplement jobs finished
in under one second each. No background job is running.

The next necessary gate is an oracle-free specification connecting the
genus-two family and its marked rational points to each original control
fibre. The universal equation above does not itself supply that selector or
correspondence, and the one-RR-pair branch construction above cannot supply
it over `Q` on the fixed panel. Until a suitable construction exists, there
is no honest blinded panel called
“the universal genus-two first-unlock curve C_t”. The updated goal remains
open, not completed by these preflight results.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_universal_rr_descent_preflight.sage
sage -python research/elliptic-curves/cas/verify_det1092_rr_fibre_cover_gate.sage
```
