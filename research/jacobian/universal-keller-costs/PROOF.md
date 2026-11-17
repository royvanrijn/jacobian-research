# The separate costs of universality, canonicity and absolute atomicity

14 November 2026. These results are independent of the quadratic-gauge formulas. They give a quantitative obstruction to canonical linear constructions, a stronger obstruction when a commuting torus is retained, a necessary nonproperness rule, and sharp distinctions between the Frobenius constructions.

The canonical prime-power construction from the previous note simplifies from a plane bundle to a line bundle. In composite prime-power degree this simpler construction loses absolute atomicity. A fixed universal map and a canonical map varying functorially with its input are different requirements.

## 1. The ambient cost of a natural symmetric fibre

**UACANCOST1.** Let k be any field, N>=2, and let F:V->W be an S_N-equivariant polynomial Keller map between linear representations of common dimension d. Suppose a full-group-fixed target y has complete fibre the natural N-point S_N set, with its specified identification. Then
\[
 \boxed{d-\dim V^{S_N}\ge N-1,\qquad
        \dim V^{S_N}\ge2,\qquad d\ge N+1.}           \tag{1}
\]
There is no restriction on the characteristic. If d=N+1 in characteristic zero, the fixed-space restriction is a nonsurjective plane Keller map; it would disprove JC(2). Conditional on JC(2), d>=N+2.

These bounds concern a natural marked construction on linear representation bundles, not an unmarked field-versal map. No existence at the lower bound is asserted.

### Proof from the actual point orbit

It suffices to work after algebraically closing k. Write the fibre points x1,...,xN in V and let M=k^N be the natural permutation module. There is an equivariant linear map
\[
 \varphi:M\longrightarrow V,\qquad e_i\longmapsto x_i.
\]
Its kernel K lies in the scalar line k*1. Indeed, if a vector in K has unequal coordinates i,j, its transposition difference is a nonzero multiple of e_i-e_j, forcing xi=xj. Thus K is zero or the scalar line; this argument needs no representation classification or averaging.

Put P=im(phi), T=V^(S_N), and t=dim T. If K=0, P has dimension N and P^(S_N) has dimension one, hence
\[
 d\ge\dim(P+T)=N+t-1.
\]
If K=k*1, P=M/(k*1) has dimension N-1 and has no invariant vector. For N>=3, an invariant class represented by v has every transposition difference in the scalar line; a third coordinate is zero in that difference, so the difference must be zero. All coordinates of v are equal and its class is zero. For N=2 in characteristic other than two the quotient is the nontrivial sign line. In characteristic two this kernel is impossible, since e1-e2=1 would identify the two fibre points. Thus in the second possible case
\[
 d\ge\dim(P+T)=N-1+t.
\]
This proves the first inequality in (1).

UAFIXALL1 gives a Keller restriction F^(S_N):V^(S_N)->W^(S_N) with empty fibre over y. In dimensions zero and one a Keller self-map is surjective over an algebraic closure in every characteristic. Therefore t>=2, proving (1).

When equality d=N+1 holds, t=2. The two possible source representation types are
\[
 M\oplus k,\qquad (M/(k*1))\oplus k^2,               \tag{2}
\]
with the second absent at N=2 in characteristic two. This follows directly from P+T=V: split off the extra invariant vector in the first case, or use P intersect T=0 in the second. In characteristic zero both types are
\[
 V\simeq H\oplus k^2,\qquad H=\{\sum u_i=0\},
 \quad\dim H=N-1.                                  \tag{3}
\]
The fixed restriction is then the stated empty-fibre plane map.

For characteristic-zero affine-linear source actions the same bound holds: averaging an orbit gives a fixed origin and converts the affine action into a linear representation. A source affine bundle without an equivariant origin in wild characteristic is outside this assertion.

### Doubly transitive groups

The same bound holds for any doubly transitive finite group G on N labels when |G| is invertible in k. Over an algebraic closure the permutation module splits as 1 plus an irreducible augmentation module H of dimension N-1. One elementary proof is that its commuting endomorphism matrices are constant on the two orbits of ordered label pairs, so End_G(M) has dimension two; semisimplicity and the single trivial summand force H to be irreducible. Distinct fibre points make phi restricted to H nonzero, hence injective. Its invariant space is zero, and the fixed-space argument again adds at least two dimensions. This extension is not asserted for all modular doubly transitive groups.

### A uniform stabilization limit

A marked linear construction in m original coordinates with r added linear coordinates must have
\[
 r\ge N+1-m.                                       \tag{4}
\]
Thus a fixed number of stabilizing coordinates cannot make three-dimensional constructions canonical in this sense across all ranks. Atomicity under identity stabilization and canonical marked descent remain separate.

At rank seven the existing unmarked universal minimum is four over the specified characteristic-zero bases. Every marked linear construction instead needs at least eight coordinates, and at least nine conditional on JC(2). These are lower bounds for the specified canonical architecture, not lower bounds for ktdim.

## 2. Retaining a commuting torus costs another coordinate

**UATORUSCOST1.** In characteristic zero, assume the hypotheses of UACANCOST1 and suppose F is also equivariant for a nontrivial linear G_m action on source and target commuting with S_N. Then
\[
 \boxed{d\ge N+2.}                                 \tag{5}
\]
No assumption of JC(2) is needed. The torus need not fix the marked target y.

The same conclusion holds for a doubly transitive G in characteristic zero. In particular a candidate in dimension N+1 must break every nontrivial commuting linear multiplicative-group symmetry.

### The planar input is established literature

Shaska's 2026 Theorem3.4 in arXiv:2607.20210v2 proves that a planar Keller map equivariant for a nontrivial G_m action is an automorphism. For the linear actions needed here there is a short algebraic proof over every characteristic-zero field.

The derivative at the origin matches the source and target weight multisets. After a linear target identification they agree. For a nonzero two-dimensional weight vector:

- For two weights of the same sign, reverse signs if needed and order 0<a<=b. The map is linear if a=b; otherwise its components are cx and dy+e*x^(b/a), with the last term possible only if a divides b. This is triangular and invertible.
- For weights (0,a), a!=0, the components are P(x) and y*Q(x). Constancy of P'*Q forces P affine and Q constant nonzero.
- For opposite primitive weights (a,-b), a,b>0, put z=x^b*y^a. The components are x*f(z), y*g(z), and
\[
 \det DF=fg+bz f'g+az fg'.                          \tag{6}
\]
If r=deg f and s=deg g, its top coefficient has nonzero factor 1+b*r+a*s. A constant nonzero determinant therefore forces r=s=0.

The zero weight vector is excluded: its action is trivial and it would say nothing about JC(2). The all-signature planar theorem is not claimed as a new result here.

### Applying the planar theorem to the canonical minimum

Suppose d=N+1. By (3), V=H plus a two-dimensional fixed summand T. The augmentation representation H is absolutely irreducible, so a commuting G_m acts on H by one character.

If the action on T is nontrivial, the fixed-plane Keller restriction is G_m-equivariant and hence an automorphism. It cannot have the required empty fibre.

If the action on T is trivial, nontriviality of the full action forces a nonzero scalar weight on H. Weighted equivariance then makes the H component of F linear in the H variables, with arbitrary polynomial coefficients in T: a monomial of H-degree e has weight e times that nonzero weight, so e must equal one. UANORMAL1 excludes precisely this affine-linear moving component, irrespective of nonlinear invariant outputs. This contradiction proves (5).

A nontrivial commuting split torus also yields (5), by choosing a nontrivial one-parameter subgroup. The proof does not exclude nonlinear group actions or a receiver which discards the torus.

## 3. A finite universal map cannot carry varying tame resolvents

Let ell be invertible in k. Kummer theory gives, for any k-scheme S,
\[
 0\to O(S)^*/O(S)^{*\ell}
  \to H^1(S,\mu_\ell)\to\operatorname{Pic}(S)[\ell]\to0.
 \tag{7}
\]
Consequently, if O(S)^*=k* and Pic(S)[ell]=0, every mu_ell torsor is constant. For affine space this always applies.

**UAFINITE1.** A finite Keller self-map of affine space in characteristic other than two has one constant discriminant square class for all its fibre algebras. It cannot be field-versal for rank-N etale algebras, N>=2. Thus every fixed universal S_N Keller map in odd characteristic is necessarily nonfinite and nonproper, in every ambient dimension.

This is a restriction on every construction, not on the boundary of one compiler.

Indeed a finite Keller map is finite etale of constant rank N. The sign of its full frame S_N torsor is a mu2 torsor on the whole affine target, and (7) makes it constant. Over k(s), the algebra k(s)^(N-2) times k(s)(sqrt(s)) has nonconstant discriminant class s, so it cannot be a fibre.

The same argument applies to any finite etale rank-N family over a base S with only constant units and no two-torsion in Pic(S). On a stack, nontrivial Picard torsion is an actual way around this particular obstruction; the statement is not asserted for arbitrary stack bases.

### Framed groups and the exact scope of the extension

For a finite group G with faithful N-point action, a G-framed Keller realization here means a fixed Keller self-map on A^d, a nonempty open U of its affine target on which the complete cover is finite etale of rank N, and a G-torsor P/U identifying that cover with P times^G the N labels. Field-versality means that P realizes every G-torsor over every extension field. At G=S_N the frame is intrinsic.

If the Keller map is finite everywhere, any such G structure extends from U to the entire target. To see this, take its global finite etale S_N frame scheme Q. A reduction to G is a section of the finite etale scheme Q/G over U. The closure of that section is a finite etale component birational to the normal affine target, hence is the whole target. The reduction therefore extends.

If G has a nontrivial cyclic quotient C_ell of order prime to the characteristic, the associated quotient torsor is constant after adjoining the roots of unity. Such a family cannot be field-versal. No split quotient hypothesis is needed: after algebraically closing constants, use the generic G-torsor k(X_g)/k(X_g)^G from the faithful regular permutation action. Its quotient by ker(G->C_ell) is a nontrivial C_ell field extension, whereas every quotient from the family is trivial.

The framing is part of this statement. A family of unstructured fibres known only to admit some subgroup reduction is not silently substituted for a specified torsor family on U.

### Why a finite chart needs a divisor

If U is an open subset of affine space with complement of codimension at least two, then O(U)^*=k* and Pic(U)=0: normal Hartogs gives the units and factoriality gives the divisor class group. Hence (7) again prevents a varying prime-to-characteristic cyclic resolvent. An S_N field-versal finite chart in odd characteristic therefore cannot have only codimension-two omissions.

This does not rule out completing a source as a nonfinite Keller map. It says that filling the finite chart while preserving a globally finite universal cover is impossible. It also explains the alternatives in (7): a nonconstant unit on an open target, or nontrivial Picard torsion on another receiver.

## 4. Sharp Frobenius costs and a simpler canonical construction

Let q=p^e and k contain F_q; write E_q=(F_q,+) and H_q=AGL1(F_q).

**UAFROBCOST1.** The canonical marked bundle construction over BH_q simplifies to the finite etale map
\[
 \boxed{\Phi_q:\ A^1_{\rm aff}\longrightarrow A^1_\chi,
                  \qquad x\longmapsto x-x^q.}       \tag{8}
\]
Its relative Jacobian is one and its zero fibre is the universal natural q-point algebra. A canonical marked affine-line construction exists for every q, and it is absolutely stably atomic when q is prime. In composite prime-power degree, no one-dimensional canonical marked affine-line construction can be absolutely atomic; the plane construction UAFROB1 is dimension-minimal for that stronger requirement.

Separately, the minimum dimension of one fixed field-versal H_q-framed Keller map, in the preceding open-chart sense, is one for q=2 and two for q>2. For E_q torsors it is one for every q. The fixed-map and canonical-family requirements must not be identified.

### Rank-one polynomial descent

For (a,b) in H_q, act by
\[
 x\longmapsto ax+b\quad\text{on the source},\qquad
 t\longmapsto at\quad\text{on the target}.
\]
Since a^q=a and b^q=b,
\[
 \Phi_q(ax+b)=a\Phi_q(x).                            \tag{9}
\]
For an H_q torsor P/S, the source descends to an affine-line bundle E_P under the character line L_chi, and the target to L_chi itself. The derivative-one identity glues because both tangent lines have character a. The zero fibre is P times^(H_q) F_q, or equivalently the universal cover B(F_q^*)->BH_q.

No primitive generator line or affine origin has been chosen. The source has no equivariant origin, since translations include x->x+1, and taking its product with vector bundles cannot create one. The character line has a canonical (q-1)-st-power trivialization; the tame quotient can live in that Picard torsion. Thus this finite family over BH_q does not contradict UAFINITE1 for one fixed map over an affine-space target.

At q=3, H_3=S3. Every cubic etale algebra in characteristic three therefore has this canonical marked one-variable Keller realization on bundles. After choosing fieldwise trivializations it is a univariate polynomial whose coefficients depend on the algebra. UA1DISC1 still excludes a single fixed universal univariate polynomial.

### Absolute atomicity forces dimension two in composite degree

The split map x-x^q has geometric monodromy the regular translation group E_q. For e>1 its factorization is already explicit:
\[
 x-x^{p^e}=(T_e(x))-(T_e(x))^p,\qquad
 T_e(x)=x+x^p+\cdots+x^{p^{e-1}}.                    \tag{10}
\]
Both factors have nonzero constant derivative and degrees greater than one.

The same obstruction holds for every one-dimensional canonical marked realization over BH_q, not just (8). On the split object the source is an affine line with H_q action faithful on its q marked points. Its elementary abelian p subgroup E_q acts by a faithful translation group W of size q: a p group has trivial multiplicative linear part in characteristic p. The target has a marked fixed point, so after translating that point its action is linear and E_q acts trivially. The univariate map f is W-invariant. Since a univariate Keller map is finite, its degree equals the marked fibre rank q. Thus f-f(0) is a nonzero multiple of
\[
 P_W(x)=\prod_{w\in W}(x-w).
\]
Its geometric monodromy is regular W.

To exhibit a polynomial factorization without appealing to a converse of primitive-monodromy atomicity, choose W0=F_p*a inside W. The polynomial
\[
 h(x)=x^p-a^{p-1}x
\]
is additive, has kernel W0 and constant nonzero derivative. Its image h(W) has size q/p, and
\[
 P_W(x)=\prod_{v\in h(W)}(h(x)-v)
\]
is a nontrivial decomposition when e>1. The outer polynomial has constant nonzero derivative by the chain rule and dominance of h. This disproves absolute atomicity, and the decomposition persists under every identity stabilization.

When e=1, the geometric degree is prime, so the characteristic-free primitive-monodromy argument gives absolute stable atomicity. In composite degree the two-dimensional canonical construction K_q from UAFROB1 has primitive geometric H_q monodromy and supplies the matching upper bound two.

### One fixed universal map: the prime-to-p quotient changes the minimum

Every one-dimensional Keller map is finite everywhere. Its framed family extends to A1 by the argument above. A connected geometric finite etale Galois cover of A1 can only have a quasi-p group: quotient by the normal subgroup generated by the p subgroups gives a cover of prime-to-p degree. A smooth compactification is tamely branched only at infinity, so Riemann--Hurwitz would give
\[
 2g-2=-d-r<-2\quad(d>1),
\]
where r>=1 counts the points above infinity. This is impossible. This is the elementary necessary direction of the classical Abhyankar theory; no sufficiency theorem for Keller maps is claimed.

A field-versal framed family must have full geometric group: a regular generic G-torsor over an extension field cannot come from a family with image inside a proper geometric subgroup. Thus one-dimensional field-versal framed groups must be quasi-p.

For H_q with q>2 the quotient H_q/E_q=F_q^* is nontrivial and prime to p. Dimension one is impossible. UAFROB1 gives the matching two-dimensional fixed map
\[
 K_q(u,x)=(u,x+(u-1)x^q).
\]
Its H_q torsor frame is defined on u!=1, and it realizes all H_q torsors over every extension field. For q=2, H_2=E_2 and (8) is already one-dimensional and field-versal.

For E_q alone, additive Hilbert90 proves that x-x^q is one fixed field-versal map in every q: an E_q-valued cocycle becomes gamma(mu)-mu, and the constant mu-mu^q gives the target. Its minimum dimension is one. Composite q makes this one-variable universal map nonatomic, as in (10).

The nonproper locus of K_q is exactly u=1. It is finite etale of degree q off that line, and has one point over every target on it; finiteness near such a target would force the same degree q, a contradiction. The multiplicative resolvent on the finite chart is
\[
 \lambda^{q-1}=-1/(u-1).
\]
For q>2 its valuation prevents extension across u=1. The boundary is necessary for this varying tame quotient, even though K_q is surjective.

### Arithmetic atomicity is weaker than absolute atomicity

Over k=F_q(a), with a transcendental constant, the univariate map
\[
 f_a(x)=x+a x^q
\]
has arithmetic monodromy H_q over k(v), but geometric monodromy E_q after closing the constant field. Its arithmetic action is primitive, so it is stably atomic over k. For e>1 it becomes decomposable after extending constants and is not absolutely atomic. Promoting a to a moving coordinate gives the absolutely stably atomic plane family. This is coefficient promotion, not identity stabilization.

### The three distinct minima

Here canonical means a marked affine-bundle construction over BH_q; fixed means one field-versal H_q-framed Keller map as defined above.

| q | Canonical marked affine-bundle dimension | Canonical dimension with absolute stable atomicity | Fixed framed field-versal dimension |
|---|---:|---:|---:|
| 2 | 1 | 1 | 1 |
| prime p>2 | 1 | 1 | 2 |
| p^e with e>1 | 1 | 2 | 2 |

These are exact results for the stated structured groups and bundle categories. They are not characteristic-zero S_N target-dimension values.

## 5. A modular control for the fixed-plane gate

The characteristic-zero qualification in the JC(2) consequence is essential. For any prime p, in characteristic p the plane map
\[
 f=x+(1+x)^p y^p,\qquad
 g=y-1-(1+x)^{p-1}y^p                              \tag{11}
\]
has Jacobian one and omits the origin. Put a=1+x and w=ay. Then
\[
 ag=w-1-f.
\]
On f=0, a cannot vanish, a=1-w^p=(1-w)^p, and
\[
 g=(w-1)/a=-(1-w)^{1-p},
\]
which never vanishes. This is a written elementary control, not a novelty claim or a claimed canonical realization at the ambient lower bound.

## 6. What this changes

Canonical linear marked construction now has a compiler-independent cost of at least N+1, and retaining a nontrivial commuting torus in characteristic zero raises it to N+2 without assuming JC(2). A general canonical construction at lower dimension must change the source-action category, not just replace the coefficient section.

Every fixed universal S_N map in odd characteristic must retain nonproperness. The affine-bundle source and the stack's character-line torsion explain how a canonical finite family can instead exist in the modular case.

The Frobenius line construction simplifies canonical descent, while its exact atomicity limitation explains why the plane construction remains useful. The characteristic-zero nonlinear canonical receiver, arbitrary-base dimension minima and general affine completions remain open.

The [replay](verify.py) supplies exact controls, not formalization or independent mathematical review of the all-degree arguments. The [literature map](LITERATURE_AND_GAPS.md) distinguishes the published planar theorem, classical Kummer and covering theory, and their explicit applications here.
