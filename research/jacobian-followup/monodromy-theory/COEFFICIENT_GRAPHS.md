# Polynomial coefficient graphs: affine completion and the root-group gate

The preceding [marked-root family](MARKED_WREATH_ACTIONS.md) has
target polynomial
\[
 P(T)=qT^n+rT^{n-1}+\sum_{j=2}^{n-2}\lambda_jT^j+T-b,\qquad m=n-2.
                                                               \tag{1}
\]
Its construction is due to Harish; Fonollosa proves coordinate descent
for constant lower coefficients. Here the lower coefficients may vary
polynomially with b. The polynomial completion lemma used below is
Fonollosa's Lemma3.1, with new first jets and a different initial step.
The prior [primary pin](sources/marked-wreath-manifest.json) is retained.

[Constructor](coefficient_graph.py) · [checker](verify_coefficient_graph.py) ·
[model](coefficient_graph_model.json) · [receipt](coefficient_graph_receipt.json) ·
[decomposition-group input](sources/coefficient-graph-manifest.json).

## Statements

**KMFI73 (polynomial b-graph completion).** For every n>=3 and arbitrary
polynomials h_j(b), 2<=j<=n-2, the full inverse image of
lambda_j=h_j(b) under the original marked-root map is isomorphic to
ordinary affine3. There is an explicit recursive polynomial source
chart, fixing x, in which all
\[
                  v_j=C_j-h_j(A_0)
\]
are coordinate functions (C2 means A2). Restricting v_j=0 gives a
polynomial Keller map
\[
                  E_{n,h}:\mathbb A^3\longrightarrow
                          \mathbb A^3_{r,q,b}
\]
with full determinant (-1)^(n-2)*(n-2)/2 in the specified ordering.
Its generic degree is n(n-2), not merely the degree of the root polynomial.
All h_j have arbitrary degree; no truncation is imposed.

**KMFI74 (full monodromy under a weighted graph gate).** More generally
allow h_j(q,r,b) polynomial. Set r=0. Suppose every nonzero monomial
q^a b^c in h_j(q,0,b) satisfies
\[
                         c+j-1>a(n-1).                       \tag{2}
\]
Then P has root group S_n over algebraically closed characteristic-zero
constants, and a divisor of simple double-root collision with q a unit.
Consequently the field pair K(t,a), a^m=P'(t), has exactly the group
G_n, arithmetic group and cyclotomic constants of KMFI69--70.
In particular (2) holds for every q-independent graph h_j(r,b),
and therefore for every actual map in KMFI73.
No such graph supplies a proper root group, including exceptional
degree11 PSL2(11). This is a constructor/field-family exclusion,
not a new absolute exclusion of that action.
For graphs outside KMFI73 the affine source is not asserted.

**KMFI75 (atomic cyclic deck witnesses in dimension3).** Every map
E_(n,h) in KMFI73 is absolutely and stably atomic. Its complete proper
block sizes are d|m, all within one root block, and its geometric generic
deck group is C_m. Its generic deck group over k is mu_m(k);
its regular polynomial fixed-target source deck group is trivial.
Thus every finite cyclic geometric generic deck group has a fully
constructive atomic witness on ordinary affine3. Independent products
give every finite abelian generic deck group on affine3s, where s is
the number of nontrivial cyclic factors. No minimum-dimension claim
or abelian-monodromy realization is made.

**KMFI76 (an exact nonlinear-graph countercontrol).** For every n>=4,
take
\[
 h_2(q,r,b)=-qb^{n-2}-rb^{n-3},\qquad h_j=0\quad(j>2).         \tag{3}
\]
Then T=b is a rational root and its marked component has field
\[
 E=K(a),\qquad
 a^m=D=1+m q b^{m+1}+(m-1)r b^m.                              \tag{4}
\]
Over algebraically closed constants E/K is cyclic of degree m>=2,
and both E and K are rational. This exact pair has no Keller model
in any dimension, even after common independent stabilization.
The corresponding actual component of the ambient inverse-image
graph is smooth affine rational3, but has a nonconstant unit and is
not stably affine space. Thus arbitrary polynomial graph restriction
does not preserve the ordinary-affine source.
Condition(2) fails at equality, so replacing its strict inequality by
a weak one would be false.

These claims do not finish OP-KMON or decide exceptional PSL2(11).

## 1. The elementary polynomial-coordinate lemma

Let R be a polynomial coefficient ring, c a nonzero constant, and
\[
 f=f_0(x,u)+x^2z,\qquad
 f_0=cu+\beta+xH(u)\pmod{x^2},\quad \beta\in R.
\]
The map
\[
 \psi(v)=(v-\beta)/c-(x/c)H((v-\beta)/c),\qquad
 (u,z)\longmapsto\left(v=f,\ w=(u-\psi(f))/x^2\right)           \tag{5}
\]
is polynomial. Its inverse is
\[
 u=\psi(v)+x^2w,\qquad z=(v-f_0(x,\psi(v)+x^2w))/x^2.          \tag{6}
\]
Both numerators are divisible by x^2 since psi and f0 invert each
other modulo x^2. Substitution proves both inverse identities.
The determinant in (u,z) is -1. This is the cited completion lemma,
and it applies with arbitrary polynomial auxiliary parameters.

## 2. A new first step and the all-degree induction

Use the original coordinates x,y,z3,...,zn, p=1+xy.
The coefficient formulas are
\[
 A_2=-\tfrac m2y-\tfrac{n(n+1)}2xy^2-\tfrac{3n}2x^2y^3
           -\tfrac32xp^2z_3+x^2z_4,
 \quad
 b=A_0=x-\tfrac n2x^2y-\tfrac12x^3z_3.                       \tag{7}
\]
Terms z4 are omitted when n=3. For n>=4,
f=A2-h2(b) has exactly x^2z4 as its z4 dependence.
Over R=k[z3], its first jet is
\[
 f=-\tfrac m2y-h_2(0)
       +x\left(-\tfrac{n(n+1)}2y^2-\tfrac32z_3-h'_2(0)\right)
          \pmod{x^2}.
\]
Apply (5) with u=y,z=z4,c2=-m/2. This produces coordinates
(x,v2,z3,s4,z5,...), and old z4 is, modulo x, a linear function
of s4 of slope -c2, plus a polynomial in v2,z3.
Unlike a rational elimination, (6) is regular on x=0.

For k>=3 the original coefficient Ck has the form
\[
 C_k=p^{e_k}z_k-2pxz_{k+1}+x^2z_{k+2}+g_k(x,y),
 \qquad e_3=3,\ e_k=2\ (k>=4).                               \tag{8}
\]
After the earlier changes, original y,z3,...,z_(k+1) involve only
the already used coordinates. In particular b involves no unused
z_(k+2), so Ck-hk(b) still has coefficient exactly x^2 in this
highest unused coordinate. Also every change fixes x, and
\[
 b=x\pmod{x^2}.
\]
Thus subtracting hk(b) changes its zeroth and first jets only by
the constants -hk(0) and -h'_k(0).

Inductively, in the current source slots (u=s_k,z=z_(k+2)),
the zeroth jet of Ck-hk(b) is c_k*u plus a polynomial in the
already chosen v's; c3=1, and each application of (5) gives
c_(k+2)=-c_k. All slopes are nonzero constants.
The offset for the next current slot depends only on earlier v's
and that slot's predecessor: substitute u=(v_k-offset)/c_k modulo x.
Consequently (5)--(6) apply through k=n-2, with no degree bound
on any graph polynomial.

The resulting full source coordinates are
(v2,...,v_(n-2),x,s_(n-1),s_n). There are n-3 pair exchanges of
determinant -1. Moving x past the n-3 parameters cancels that sign,
so this full coordinate change has determinant1 in this ordering.

On the target, subtracting h_j(b) is a polynomial automorphism.
Moving A2 before C3,...,Cn introduces sign (-1)^(n-2).
The ambient determinant m/2 from KMFI69 therefore becomes
(-1)^(n-2)m/2. In the new source and target coordinates the first
n-3 outputs are the identical parameters v_j, so the full
three-variable fiber Jacobian is the same constant.
The full inverse image of v=0 is precisely affine3, and the
explicit substitutions are polynomial at every source point.

On its dense xp!=0 open, t=x/p,a=1/p have exactly the old two-sided
root-and-marking reconstruction. Hence its generic source field is
K(t,a). The next section proves its degree is nm.

## 3. A uniform root-group proof, allowing q-dependent graphs

Under (2), substitute r=0 and
\[
 q=\epsilon^{-(n-1)},\quad b=\epsilon B,\quad T=\epsilon Z,
 \qquad F_\epsilon(Z,B)=\epsilon^{-1}P(\epsilon Z).
\]
A monomial q^a b^c T^j becomes
\[
                  \epsilon^{c+j-1-a(n-1)}B^cZ^j.
\]
Thus F_epsilon is monic and polynomial in epsilon, and
\[
                         F_0=Z^n+Z-B.                       \tag{9}
\]
This identity is an exact coefficient statement, not a numerical limit.

The polynomial Z^n+Z has n-1 distinct critical points alpha,
n*alpha^(n-1)=-1, and distinct critical values (n-1)*alpha/n.
The connected degree-n cover of the B-line is therefore generated
by transpositions and has group S_n.

Over the epsilon-adic DVR with residue field k(B), (9) is separable.
The splitting extension is unramified there: its roots and their
distinct residues lift over the unramified extension of the henselian
DVR. Its decomposition group is therefore the residue splitting
group S_n. The decomposition/inertia facts used here are
[Stacks, Lemma15.114.4](https://stacks.math.columbia.edu/tag/09E3)
and the elementary simple-root Hensel lifting.
Hence the generic group of F_epsilon contains S_n.
Base change q=epsilon^(-(n-1)) and unramified specialization r=0
cannot enlarge a root group; the original group is S_n as well.
This argument works after any algebraically closed characteristic-zero
constant extension.

The discriminant of (9) has n-1 simple zeros in B.
They lift as simple zeros of disc(F_epsilon) over k[[epsilon]].
At each such branch value the polynomial has exactly one double
root and the other n-2 roots simple. The simplicity and transverse
B derivative persist by Hensel/implicit lifting.
Consequently the original target discriminant has a multiplicity-one
prime divisor, away from q=0, with this generic root pattern:
a repeated global factor could not restrict to a simple zero.
In the root splitting field the colliding difference has valuation1,
and every other root difference and q is a unit.

These are precisely the pair valuations in the Kummer argument of
KMFI69. They force r_i+r_j=0 modm for every putative derivative-root
relation, leaving only the single even relation.
The identical derivative-product identity then gives G_n, and
KMFI70's arithmetic size bounds give constants k(mu_m) and its
exact cyclotomic parity twist.
The full marked extension is transitive of degree nm.

In particular all q-independent h_j(r,b) satisfy (2), since
c+j-1>=1 for every term surviving r=0. To obtain a proper root
group from a polynomial coefficient graph, at least one term with
r-degree zero must have positive q degree and violate the strict
weight inequality. This is necessary for that route, not sufficient.

## 4. Blocks, generic deck groups and source poles

The action and marking field are the same as in KMFI69.
The group-theoretic arguments in KMFI71 give exactly the fields
K(t,a^d), d|m, and cyclic upper extensions of degrees d.
KMFI1--2 then exclude every proposed nontrivial polynomial inner
factor; independent stabilization preserves this argument.
This proves absolute and stable atomicity on the actual affine3
models in Section2.

The geometric generic deck rotations are a->zeta*a.
Their old source functions satisfy
\[
                    x\mapsto x/\zeta,\qquad
                    y\mapsto y+(1-\zeta)/x.
\]
The new source chart fixes x, and original y is polynomial in its
three coordinates. Hence the latter function still has a genuine
pole on x=0 whenever zeta!=1. No nonidentity deck rotation extends
to a regular source endomorphism. Over k only mu_m(k) rotations
are defined. Independent products give the claimed finite abelian
deck groups, without identifying deck symmetry with monodromy.

## 5. A nonlinear graph that actually destroys the affine source

For (3), direct evaluation gives P(b)=0, and
\[
 P'(b)=1+(n-2)q b^{n-1}+(n-3)r b^{n-2}=D.
\]
Over k(b,r), D has a simple linear zero as a polynomial in q.
The corresponding valuation has value1 on D, so X^m-D is
Eisenstein and has degree m. With algebraically closed constants,
adjoining a is a cyclic Galois extension.
Moreover q=(a^m-1-(m-1)rb^m)/(m b^(m+1)), so E=k(b,r,a)
is purely transcendental, as is K=k(b,r,q).
KMFI1--2 forbid this nontrivial Galois pair as a Keller field pair,
in every stabilized or birationally changed affine presentation.

There is an explicit stronger check on the actual source component.
It is a component of the smooth affine base change of the ambient
etale marked-root map to the target graph. Generically t=b,
so x=b*p. This identity holds everywhere on that component, and
\[
                p(1-by)-1=y(x-bp)=0.                        \tag{10}
\]
Thus p is a unit with inverse1-by. It is nonconstant:
generically p=1/a and a^m=D is nonconstant.
The component is rational by the calculation of E, but its
nonconstant unit persists under every polynomial stabilization.
It is therefore not isomorphic, even stably, to ordinary affine space.

Finally h2(q,0,b)=-q b^(n-2) gives exponent0 in(2), and the scaled
polynomial is Z^n-B^(n-2)Z^2+Z-B, with rational root Z=B.
This proves the equality-boundary countercontrol for every n>=4.
Failure of (2) alone does not prove a proper group or a bad source;
(3)--(10) are the exact evidence for this particular family.

The finite checker records explicit polynomial coordinate divisions,
first jets, graph identities and whole three-variable Jacobians,
plus exact scaling and boundary countercontrols. The induction and
decomposition/Kummer proofs supply the unbounded results.
Formal verification, independent whole-proof review and a full
classification of Keller monodromy remain unclaimed.
