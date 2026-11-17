# Nonlinear sweep monodromy and the remaining construction gates

Work over algebraically closed constants of characteristic zero.
This note extends the [ruled-sweep analysis](SWEEP_MONODROMY.md).
The source framework is Gao,
[arXiv:2608.00222v2](https://arxiv.org/html/2608.00222v2),
Sections 4.3, 4.4.2 and 4.6, retained in the
[source manifest](sources/sweep-monodromy-manifest.json).
The monodromy arguments below are written separately from that construction.

**KMFI26.** Every nontrivial four-dimensional curve-type sweep in
Gao's normalized family, including every nonlinear transverse
potential, has full natural symmetric monodromy.
More precisely, for any P(w,z), phi(w) and c!=0, put

\[
 T=-P_w+(2P+\phi')P_z/c,\quad
 X=(T,\ P+wT,\ 2wP+cz+\phi+w^2T),\quad
 S=X+\gamma(1,w,w^2).
 \tag{1}
\]

Then det JS=c gamma. If its generic degree N>1, the group is S_N.
Any Gao polynomial affine twist has the same group.
An explicit nonlinear, rational-coefficient A4 Keller map of
degree four and group S4 is constructed in Section 4.
Here A4 denotes affine four-space, not the alternating group.

**KMFI27.** For a finite normalized sweep S:A^(n-1)->A^(n-1)
with Jacobian c gamma^k, gamma a source coordinate and k>=1,
its unique branch divisor has inertia

\[
 (k+1)^\delta 1^{N-(k+1)\delta}, \tag{2}
\]

where delta is the generic degree of X=S|_(gamma=0) onto its image.
These inertia elements normally generate the geometric group.
For an actual Gao affine Keller twist, N>(k+1)delta.
For PSL2(11) this leaves exactly the rows in Section 5:
the natural degree-twelve action first requires dimension at least six.

**KMFI28.** In the finite tautological three-variable sweep with
Delta=(1,w,z), suppose its nonzero rank-one nilpotent normalization
matrix has polynomial entries

\[
 P_0=\lambda\begin{pmatrix}-h(w)&1\\-h(w)^2&h(w)\end{pmatrix},
 \qquad h\in k[w]\text{ nonconstant},\quad\lambda\in k[w,z]\setminus0.
 \tag{3}
\]

Then the generic degree delta of X onto its branch surface is
at most two. If delta=1 the group is A_N.
If delta=2 the sweep has a quadratic source extension over
a proper intermediate field, hence has imprimitive monodromy.
Consequently every primitive nontrivial finite sweep in this
whole one-parameter polynomial-multiplier class has group A_N.
In particular none realizes either PSL2(11) action.

These are construction-class results, not a classification of
all polynomial Keller maps. Rational multipliers, genuinely
two-parameter polynomial multipliers, other direction fields,
nonfinite sweeps and the mixed-fiber Shabat completion remain separate.
The abstract PSL2(11) occurrence questions remain open.

## 1. The nonlinear collision curve

The polynomial change z -> z+phi(w)/c replaces P by
P(w,z-phi(w)/c) and removes phi from (1).
We may therefore take phi=0. Let e=deg_z P.
For e<=1, X is ruled and KMFI24 applies.
Assume e>=2, and write P=sum a_j(w)z^j, a_e!=0.

The inverse equation for a target Y=(Y1,Y2,Y3) is

\[
 z_Y(w)=(Y_3-2wY_2+w^2Y_1)/c,\quad
 R_Y(w)=P(w,z_Y(w))+wY_1-Y_2. \tag{4}
\]

It determines z and gamma=Y1-T uniquely for each finite root.
Moreover

\[
 R_Y'=\gamma+2P_zR_Y/c. \tag{5}
\]

For two distinct parameter values w,u, put d=w-u.
The ordered-pair incidence is isomorphic, on d!=0, to

\[
 G=c(z-t)+d\{P(w,z)+P(u,t)\}=0, \tag{6}
\]

by recovering
Y1=(P(u,t)-P(w,z))/d and then Y2,Y3 from (4).
In particular no roots with the same w are being discarded.

The generic projective curve (6) over k(w,u) is smooth and
geometrically irreducible. Here is an all-degree degeneration proof.
Set d=epsilon^(e-1), u=w-epsilon^(e-1), z=Z/epsilon,
t=Tvar/epsilon, and multiply G by epsilon.
Its degree-e projective limit is

\[
 F_0=a_e(w)(Z^e+Tvar^e)+c(Z-Tvar)H^{e-1}. \tag{7}
\]

It is smooth at H=0, since its Z,Tvar derivatives cannot both
vanish at a projective point there.
On H=1, with a=a_e(w), F=a(Z^e+Tvar^e)+c(Z-Tvar),
write J=(Z^(e-1)-Tvar^(e-1))/(Z-Tvar). The identity

\[
 c(e-1)(F_Z-F_T)-eaJ(eF-ZF_Z-Tvar F_T)
   =2c^2(e-1) \tag{8}
\]

proves that F,F_Z,F_T have no common zero.
The projective plane curve (7) is therefore smooth; it is also
geometrically irreducible, since distinct positive-degree plane
components would meet and make a singularity.

Smoothness of a projective family is open: the projected singular
locus is closed by properness.
The substituted generic fiber is smooth, and this descends across
the finite extension k(w,u) -> k(w,epsilon).
Thus the generic curve (6) is smooth and geometrically irreducible.

The polynomial G has no content in k[w,u].
A common factor of its coefficients would divide both
d a_e(w) and d a_e(u), hence d; the coefficients reducing to
c and -c on d=0 exclude this.
Gauss's lemma makes the total ordered-pair incidence irreducible.
It dominates Y-space because the generic root degree is at least four.
The group is therefore doubly transitive.

## 2. A generic branch value has only one double root

The restriction X(w,z) is generically finite onto a surface D.
Indeed T has leading z-term (2e a_e^2/c)z^(2e-1).
The two-by-two Jacobian of (T,P+wT) has leading term
-(2e-1)(2e a_e^2/c)^2 z^(4e-3), which is nonzero.
The simple Jacobian divisor gamma=0 gives local ramification
index two at its general point.

We show that X has generic degree one.
At a root, the target normal obtained from (4) is

\[
 \nu(w,z)=\bigl(w+a w^2,\ -1-2aw,\ a\bigr),
 \qquad a=P_z(w,z)/c. \tag{9}
\]

Suppose a general Y in D had two different ramified points
(w,z),(u,t). Their normals to the same smooth branch surface
are proportional. Pairing them with (1,s,s^2), whose roots
are w,u, gives

\[
 c+dP_z(w,z)=0,\qquad -c+dP_z(u,t)=0. \tag{10}
\]

Thus G=G_z=G_t=0.
The ordered ramified-pair variety would have dimension two
and would dominate each copy of the ramification surface.
Consequently w,u are general, a_e(w)a_e(u)!=0.
For fixed such w,u, the last two equations in (10) are nonzero
polynomials of degree e-1 separately in z and t, so have finite
solution sets. The pair variety must therefore dominate the
(w,u)-plane. That contradicts the generic smoothness of (6).
Hence X is birational onto D.

No missing infinite root can add inertia at general D.
Give w,z weights one and two. The degree N of R is the weighted
degree of P, and its leading coefficient is P_top(1,Y1/c),
a nonzero polynomial in Y1 alone.
The degree-drop locus is a finite union of planes Y1=constant.
But Y1=T(w,z) is nonconstant on D, so D is none of these planes.
After inverting that leading coefficient, (4) makes the sweep
finite and all roots remain finite.

It follows that generic inertia at D is exactly one transposition.
Combined with double transitivity, this proves S_N and KMFI26.
The proof allows arbitrary nonlinear P and does not assume
global finiteness of S.

## 3. Finite normalized sweep inertia

For KMFI27, finiteness makes the smooth source the entire
normalization of the target in its function field.
Its only ramification divisor is gamma=0.
In characteristic zero its different exponent equals e_ram-1;
Jacobian c gamma^k gives e_ram=k+1.
There are delta geometric points of this divisor over the
generic branch point. All other source branches are unramified.
This proves (2).

The normal subgroup generated by finite affine inertia is the full
group: after quotienting by it, purity makes the finite cover of
affine space etale, hence trivial. This is the same argument as KMFI1.
All branch meridians here are conjugate.

A Gao target twist is an isomorphism on C!=0:
Yi=C^(e_i) times the polynomial target output.
The branch divisor of S persists on that affine torus chart,
with one independent C coordinate.
An actual polynomial Keller map must have fixed-sheet inertia
there by KMFI1. Thus (2) must have at least one fixed sheet.
This conclusion concerns this target twist architecture; it is
not a claim about arbitrary birational contractions.

## 4. An actual nonlinear S4 Keller control

Take c=1, phi=0 and

\[
 P=z^2+w^3-w^2,\quad
 T=4z^3+4w^3z-4w^2z-3w^2+2w,\quad
 X=(T,P+wT,2wP+z+w^2T).
\]

The sweep has Jacobian gamma and generic degree four.
At Y=(4,1,1), its root polynomial is

\[
 R_Y=w^2(16w^2-15w+11),
\]

with exactly one double root and a squarefree nonzero residual.

Here is a full affine lift, not only this branched sweep.
On source coordinates (x,y,z_src,t_src) put

\[
 v_1=xy,\quad v_2=xz_{\rm src},\quad v_3=x^2t_{\rm src},
\]
\[
 \gamma=1+3v_1+v_3,\quad u_1=1-2v_1,\quad u_2=v_2,
 \quad w=\gamma u_1,\quad z=\gamma u_2,\quad C=\gamma x.
\]

Define

\[
 F=(C,\ S_1/C,\ S_2/C^2,\ S_3/C). \tag{11}
\]

All components are polynomials.
After substituting w=gamma u1,z=gamma u2, the S_i are divisible
by gamma^(1,2,1), respectively.
The resulting quotients vanish at (gamma,u1,u2)=(1,1,0).
For the middle quotient, its gradient there is (-2,-3,0).
The two stage columns of weight one are (3,-2,0) and (0,0,1),
both in this kernel. This is exactly the extra first-jet condition
needed for divisibility by x^2; the other outputs need only x.

The affine stage determinant is -2.
Its determinant chain is x^4*(-2)*gamma^2*gamma^2/C^4=-2.
Thus (11) is a polynomial Keller self-map of ordinary A4.
The source inverse recovers x=C/gamma, u1=w/gamma,u2=z/gamma,
v1=(1-u1)/2,v2=u2,v3=gamma-1-3v1, and then y,z_src,t_src.
The target inverse is Y=(C F2,C^2 F3,C F4).
The field pair is the sweep plus independent C, so its geometric
and arithmetic groups are S4.

This is a concrete control of the nonlinear theorem, without
a novelty-priority claim.

## 5. The exact PSL2(11) finite-sweep gates

For the natural degree-twelve action the only nontrivial
fixed-sheet uniform-cycle types are 5^2 1^2 and 11 1.
For either exceptional degree-eleven action they are
2^4 1^3, 3^3 1^2 and 5^2 1.
The full cycle dictionaries are retained in KMFI3 and KMFI10
and the canonical PSL2(11) benchmark.

In Gao's n-dimensional architecture, the branch index k
satisfies 1<=k<=n-2. Equation (2) therefore gives:

| action | k | delta | fixed sheets | minimum affine dimension |
|---|---:|---:|---:|---:|
| natural degree 12 | 4 | 2 | 2 | 6 |
| natural degree 12 | 10 | 1 | 1 | 12 |
| exceptional degree 11 | 1 | 4 | 3 | 3 |
| exceptional degree 11 | 2 | 3 | 2 | 4 |
| exceptional degree 11 | 4 | 2 | 1 | 6 |

These are necessary gates, not constructed maps.
In particular finite normalized sweeps in dimensions three,
four or five cannot give the natural degree-twelve action.
The curve-type four-dimensional route is further excluded by KMFI26.
The one-parameter multiplier route to the k=2 row is excluded next.

## 6. The Gauss degree of a one-parameter multiplier

Assume (3) and finite S=X+gamma(1,w,z).
Finiteness gives a branch surface D and generic degree delta
for X onto D. Its normal is

\[
 (\sigma,-h,1),\qquad \sigma=h(w)w-z.
\]

At smooth points this normal is a rational function of X.
Consequently

\[
 k(h,\sigma)\subset k(X)\subset k(w,\sigma)=k(w,z). \tag{12}
\]

Write d=deg h. The outer extension has degree d, so delta divides d.
All intermediate fields descend from k(w)/k(h) after the independent
variable sigma is adjoined: take its normal closure and use the
unchanged finite subgroup correspondence.
Luroth and the unique infinity place of the polynomial h then give

\[
 k(X)=k(q(w),\sigma),\quad h=a(q),\quad
 q,a\text{ polynomials},\quad\deg q=\delta. \tag{13}
\]

Every X_i is polynomial in q,sigma.
Indeed k[w,sigma] is integral over the normal ring k[q,sigma],
so its intersection with k(q,sigma) is exactly k[q,sigma].

Put F=X1 in this ring, and use D0=partial_w+h partial_z.
The matrix (3) gives

\[
 D_0X_1=h'\lambda,\quad (X_1)_z=-D_0\lambda,\quad
 D_0X_2=h'\lambda w,\quad D_0X_3=h'\lambda z.
 \tag{14}
\]

For delta>=3, 1,w,w^2 are independent over k(q,sigma).
Writing A=F_q/a', B=F_sigma gives lambda=A+wB.
Equations (14) force

\[
 A_q=0,\quad a'A_\sigma+B_q=0,\quad B_\sigma=0.
\]

Since F_(q sigma)=a'A_sigma=B_q, both mixed terms vanish,
and A,B are constants alpha,beta.
Thus F=alpha h+beta sigma+constant, lambda=alpha+beta w.
The X2 equation, compared in 1,w,w^2, forces beta=0 and
(X2)_q=0,(X2)_sigma=alpha.
The X3 equation forces

\[
 (X3)_q=-\alpha a'\sigma,\quad (X3)_\sigma=\alpha h.
\]

Equality of mixed partials gives -alpha a'=alpha a', hence alpha=0.
This contradicts lambda!=0. Therefore delta<=2.

If delta=1, generic branch inertia is a single three-cycle.
Finiteness and Section 3 make the whole group generated by its
conjugate three-cycles. A transitive group generated by three-cycles
is A_N, by the elementary support-hypergraph proof in KMFI12.

## 7. Delta two gives an imprimitive action

When delta=2, an affine change in w and a linear target change
normalize q=w^2; the corresponding normal coordinate becomes
sigma=h(w)w-z, with h=a(q).
Repeat (14), now using w^2=q. They give

\[
 F=f(q)+g(\sigma),\quad
 A'= -a'q g'',\quad
 (X2)_q=a'q g',\quad (X2)_\sigma=A .
\]

Mixed partials of X2 force g''=0 and A constant.
Mixed partials of X3 then force that constant to vanish.
After harmless target translations and a common nonzero scaling,
the surface is exactly

\[
 X=(\sigma,\ b(q),\ d(q)-\sigma^2/2),\qquad
 b'=q a',\quad d'=a b',\qquad \lambda=w. \tag{15}
\]

For its sweep the target satisfies

\[
 Y_1=\sigma+\gamma,\quad
 Y_2=b(q)+\gamma w,\quad
 Y_3+Y_1^2/2=d(q)+a(q)\gamma w+\gamma^2/2.
\]

The polynomial involution

\[
 (w,\sigma,\gamma)\longmapsto(-w,\sigma+2\gamma,-\gamma)
\]

fixes Y and q. Conversely, from Y,q,w recover gamma=(Y2-b(q))/w
and sigma=Y1-gamma.
Thus the full source field is quadratic over K(q).
Since inertia has two three-cycles, N>=6, so this is a proper
intermediate field and the action is imprimitive.
This completes KMFI28.

For a concrete finite delta-two control, choose h=w^2,lambda=w:
X=(w^3-z,w^4/2,w^6/3-(w^3-z)^2/2).
Its degree-eight eliminant is

\[
 w^8-12Y_2w^4+12(Y_1^2+2Y_3)w^2-12Y_2^2.
\]

The source is integral over the target (the other equation is
monic quadratic in z), and it has the exhibited quadratic block.
This is a finite-sweep control, not an asserted affine Keller lift.

## 8. What remains to construct

KMFI26 closes all nonlinear transverse potentials for the
four-dimensional curve-type field, extending the ruled exclusion.
KMFI28 closes primitive finite tautological sweeps with every
polynomial multiplier depending only on w.
[KMFI32--34](SWEEP_ASYMPTOTICS.md) now close the entire finite
tautological k=2,delta=3 route, without a multiplier assumption:
its asymptotic-direction bound is delta<=2.
They also show every polynomial multiplier is independent of z and
classify the complete homogeneous rational multiplier class.
A different direction field, another architecture, or a separately
analyzed nonfinite sweep is required for this exceptional action.
The k=1 and k=4 gates need other data satisfying the finite-sweep
and fixed-sheet constraints. No gate is sufficient for completion.

The historical marked candidate is retained in
[sweep_frontier_next_candidate.json](sweep_frontier_next_candidate.json).
In (15), take a(q)=alpha q^2+q with
16alpha^2+30alpha+15=0.
Then b=2alpha q^3/3+q^2/2 and
d=2alpha^2 q^5/5+3alpha q^4/4+q^3/3.
At w=1, z=alpha/3+1/2, gamma=-2alpha/3-1/2,
the sweep reaches the origin with gamma!=0.
Scaling weights (1,3) and output exponents (1,2,4) satisfy
the determinant balance and gamma divisibility.
This checkpoint originally left the affine stage and generic group open.
[KMFI29](SIGNED_ALTERNATING.md) now completes both and embeds it in an
infinite family of actual signed-alternating Keller maps. The retained
checkpoint flags describe its original evidence, not current status.

[Checker](verify_sweep_frontier.py) · [model](sweep_frontier_model.json) ·
[receipt](sweep_frontier_receipt.json).
Written proofs plus exact controls are not formal verification
or independent whole-proof review. The full goal remains open.
