# Coordinate-level obstructions and a canonical Frobenius construction

14 November 2026. This closes the previously unclassified B=constant level for the independently promoted coefficient families. It also supplies a constructive alternative in positive characteristic: canonical marked Keller maps on affine bundles for affine finite-field lines, including every cubic etale algebra in characteristic three.

The classical ingredients are finite-field spreading, Hilbert 90, additive polynomials and associated bundles. The claims concern their explicit application and exact limits in this programme; no literature-priority assertion is made.

## 1. The actual B levels in every degree

Let k have characteristic other than two. Retain an arbitrary finite subset J of seed indices j>=4, put r=|J|, and promote each a_j independently as an affine coordinate. All other seeds are zero. For generic degree N>=4, retain N=max J; J empty is the cubic map. Use
\[
 t=1+xy,\quad q=t^2z+y^2(1+3t),\quad \Pi=tq,
\]
\[
 B=y+3xq+\sum_{j\in J}j a_jt^2x^{j-2}q^j,\qquad
 C=x(5-3t)-x^3z-\sum_{j\in J}(j-2)a_j(xq)^j.          \tag{1}
\]
The normalized target is (a_j;Pi,-B/2,C). For beta in k set
\[
 S_\beta=\{B=\beta\}\subset A^{r+3}.
\]
It is smooth of dimension r+2, since the promoted compiler is etale and B is a nonzero scalar multiple of one target coordinate.

**UABLEVEL1.** The following classification holds over every such field, including after any constant-field extension:
\[
 \boxed{S_\beta\text{ is stably affine space }
 \quad\Longleftrightarrow\quad
 \operatorname{char}k=3\text{ and }3\mid j
          \text{ for every }j\in J.}                 \tag{2}
\]
The condition is vacuous for J empty. In the exceptional case B=y and S_beta is literally A^(r+2). In every other case S_beta times any affine space is not an affine space. Thus B is a coordinate, or even a stable coordinate, exactly in the exceptional case.

This statement concerns independent promoted coefficients. Relations such as a4=a5=w, nonlinear monomial coefficient sections, fixed nonzero seed values, and arbitrary rational profiles define other source levels and require their own proofs.

### A four-stratum computation

Assume first that some j0 in J is nonzero in k, so r>=1. Write L=[A1] in the Grothendieck ring of varieties. The four disjoint strata are:

| Stratum | Description | Class |
|---|---|---:|
| x=0 | y=beta, seeds and z free | L^(r+1) |
| t=0 | y=-beta/2, x=2/beta; empty if beta=0 | L^(r+1) if beta!=0 |
| x*t*q!=0 | use x,t,q as torus coordinates and solve a_(j0) | (L-1)^3 L^(r-1) |
| x*t!=0, q=0 | y=beta; x!=0 and 1+beta*x!=0; seeds free | (L-2)L^r if beta!=0, otherwise (L-1)L^r |

For the third row, t is a unit so q replaces z, and the coefficient j0*t^2*x^(j0-2)*q^j0 is a unit. The remaining seeds are free. No coefficient-size or degree calculation is involved.

Consequently
\[
 [S_\beta]=
 \begin{cases}
 L^{r+2}+L^r-L^{r-1},&\beta\ne0,\\
 L^{r+2}-L^{r+1}+2L^r-L^{r-1},&\beta=0.
 \end{cases}                                         \tag{3}
\]
Over any finite field F_Q in which j0 is nonzero this is also an exact point count, with L replaced by Q. The differences from Q^(r+2) are
\[
 Q^{r-1}(Q-1),\qquad -Q^{r-1}(Q-1)^2,                 \tag{4}
\]
respectively. Both are nonzero. For active coefficients the nonzero-level polynomial has value one at L=1, so the ordinary characteristic-zero Euler characteristic alone would miss the obstruction.

At beta=0 there is also an elementary unit certificate. Modulo t one has xy=-1, q=y^2 and xB=2. Thus t is invertible on B=0 because two is a unit. Outside the affine exception, t is nonconstant on the displayed torus stratum (or on the cubic Gm^2 chart in the inactive case), so stabilization cannot make that level affine space.

If every retained j vanishes in k and char k!=3, the seed terms in B vanish and B=y+3xq. On x*t!=0 one solves uniquely for z, giving Gm^2 times A^r. The x=0 and t=0 rows are unchanged. Thus
\[
 [S_\beta]=
 \begin{cases}
 L^r(L^2+1),&\beta\ne0,\\
 L^r(L^2-L+1),&\beta=0.
 \end{cases}                                         \tag{5}
\]
These counts also differ from affine space over every finite residue field of characteristic other than two and three. In characteristic zero this case only occurs for J empty. If char k=3 and every j vanishes, (1) simply has B=y, as asserted.

### The nonzero levels pass several simpler affine-space tests

In the active-coefficient case and beta!=0, S_beta is geometrically integral and rational, factorial, and has only constant units. Indeed its dense open x*t*q!=0 is Gm^3 times A^(r-1). Every other stratum has smaller dimension, while the smooth hypersurface has pure dimension r+2, so no additional component exists, even geometrically.

The three boundary divisors are prime and reduced: x=0 and t=0 are affine spaces of dimension r+1, while
\[
 q=0\simeq\operatorname{Spec}k[x,1/(1+\beta x)]
                     \times A^r.
\]
The last description uses y=beta and solves z since t is a unit. These divisors are respectively the principal divisors of x,t,q, each of multiplicity one. Any Weil divisor becomes principal on the factorial torus-times-affine open; its difference from that principal divisor is supported on these three already principal divisors. Hence every divisor class is zero, and the smooth affine coordinate ring is factorial.

A unit restricts on that open to a constant times x^a*t^b*q^c. The valuations at the three boundary divisors are respectively a,b,c, since the other two coordinates are generically nonzero there. A unit has all three valuations zero, so it is constant. Thus smoothness, rationality, factoriality, constant units and Euler characteristic one can all hold while (4) still excludes stable affine space.

### Why finite counts exclude stable affine space over arbitrary fields

The count is used with a spreading argument, not as an inference from a bounded list of finite fields or as an assertion that L is algebraically independent in every Grothendieck ring.

All nonzero beta levels are isomorphic to S_1: the retained scaling sends B to lambda^(-1)*B, so take lambda=beta. Thus only beta=0,1 need consideration.

In the active-coefficient case define the formulas over R=Z[1/(2*j0)] in characteristic zero, or over the prime field in positive characteristic. In the inactive noncubic case use Z[1/6] for J empty in characteristic zero, or the prime field of the stated characteristic. A hypothetical isomorphism
\[
 S_\beta\times A^m\simeq A^{r+2+m}
\]
over k, together with its inverse, uses finitely many coefficients. It descends, after enlarging that coefficient ring, to a nonzero finitely generated R-algebra D. The defining inverse identities also descend, since the coordinate rings are finitely presented. D has a closed point with finite residue field F_Q, and the inverted integers stay units there. Specialization would give an isomorphism over F_Q. But multiplying any count discrepancy in (4) or (5) by Q^m cannot make it zero. This contradiction proves (2).

The finite-residue-field and descent facts are standard; primary references are given in the literature map. This is an obstruction over every field extension, not just over C.

### Completing the three-coordinate normalization test

For these same independent promoted families, no level Pi=pi0 or C=c0 is stably affine space, in any characteristic other than two.

If pi0!=0, t is a unit on Pi=pi0. It is nonconstant: set all seeds zero, x=1,y=u-1,t=u, and
\[
 z=\bigl(\pi_0/u-(u-1)^2(1+3u)\bigr)/u^2,\qquad u\in G_m.
\]
Thus every affine stabilization retains a nonconstant unit. If pi0=0, the comaximal nonempty t=0 and q=0 pieces are disjoint open-and-closed components.

Likewise C=x*K, where K modulo x is two. At c0!=0, x is a unit and is nonconstant on the curve
\[
 a_j=0,\quad x=u,\quad y=0,\quad z=(2u-c_0)/u^3.
\]
At c0=0, x=0 and K=0 are nonempty disjoint open-and-closed components, because (x,K)=1. For example the latter contains a_j=0,x=1,y=0,z=2.

Therefore in characteristic zero **none of the three moving target-coordinate levels produces an affine-space source, even after stabilization**. This closes those direct coordinate-normalization routes for the named atlas, not every coefficient section or another affine completion.

## 2. A universal lower bound from univariate discriminants

Use the programme's definition of field-versality: one fixed Keller self-map with a family of complete finite-etale rank-N fibres realizing every rank-N algebra over every field extension of k. We use the identical definition in positive characteristic.

**UA1DISC1.** For every field k of characteristic other than two and N>=2, no one-dimensional Keller map is field-versal for rank-N finite etale algebras. Thus the corresponding Keller target dimension, whenever finite, is at least two.

In characteristic zero a univariate Keller polynomial is linear, so the assertion is immediate. In characteristic p, a univariate Keller polynomial f has derivative a1 in k*, and its degree is either one or divisible by p. It is finite, and every complete geometric fibre has that same degree and is etale. Thus a rank-N family requires deg f=N.

Write the leading coefficient aN. The discriminant of f(X)-v is independent of v:
\[
 \operatorname{Disc}(f-v)
   =(-1)^{N(N-1)/2}a_N^{\,N-2}a_1^N.                \tag{6}
\]
One can verify this without degree conventions for resultants: multiply f'(r_i)=a1 over all roots, use the root-difference product, and the usual leading-coefficient normalization. The discriminant of the monic polynomial differs by a square, so every fibre algebra has the same discriminant square class delta in k*/k*^2.

Over L=k(s), the algebra
\[
 L^{\,N-2}\times L(\sqrt{s})
\]
has discriminant class s. Its quotient by the constant delta has odd s-adic valuation and cannot be a square. It is not a fibre. This disproves field-versality, even if the original k has few quadratic extensions.

This is a lower bound for every one-dimensional Keller construction, not a boundary invariant of the quadratic-gauge compiler. It gives a matching minimum in the next rank-three case.

## 3. Prime-power Frobenius maps: atomicity, versality and canonical descent

Let q=p^e be any prime power and let k contain F_q. On two-dimensional affine space define
\[
 \boxed{\mathcal K_q(u,x)=
       \bigl(u,\ x+(u-1)x^q\bigr).}                  \tag{7}
\]
Its Jacobian is one. On u!=1 it is finite etale of rank q; on u=1 it has rank one.

**UAFROB1.** The map (7) has arithmetic and geometric monodromy
\[
 H_q=\operatorname{AGL}_1(F_q)=F_q\rtimes F_q^*
\]
in its natural q-point action. It is absolutely and stably atomic, including at composite prime-power degrees. Its rank-q fibres are exactly the finite etale algebras whose geometric q-element sets admit a Galois-invariant affine F_q-line structure.

Moreover it has a canonical marked descent over BH_q as a polynomial Keller morphism **from an affine-plane bundle to a vector bundle**, with complete marked fibre the universal q-point algebra. At q=3 this gives a canonical marked construction for every cubic etale algebra in characteristic three, without retaining a generator line. It does not choose trivializations of these bundles as fixed affine spaces.

### Exact monodromy

Put z=u-1 and use target coordinate v. Choose lambda with
\[
 \lambda^{q-1}=-1/z,\qquad x=\lambda y.
\]
Then the inverse equation becomes
\[
 y^q-y=-v/\lambda.                                  \tag{8}
\]
The Kummer extension has degree q-1 by the valuation of z at zero. Over k(lambda,v), (8) has degree q: v=lambda*(y-y^q) is a rational function of degree q in y, and its derivative is nonzero. Its q roots differ by F_q, so it is Galois with additive group F_q. The full splitting field is k(lambda,y), a rational field, and has degree q(q-1) over k(z,v). Transformations of the roots are exactly x_i -> x_(a*i+b). Thus the group is H_q, over k and its algebraic closure.

The affine group is doubly transitive: there is a unique affine transformation carrying one ordered pair of distinct field elements to another. Hence it is primitive even when q is composite. The characteristic-free primitive-monodromy argument recorded in UAALLATOM1 proves absolute stable atomicity.

### The exact field-realization criterion

A fibre with u!=1 has roots forming an affine F_q-line, since root differences satisfy z*D^q+D=0. Thus it has the stated structure.

Conversely let a finite etale algebra A/k of rank q have such a structure. Its finite Galois action can be written
\[
 i\longmapsto a_\gamma i+b_\gamma,
 \qquad a_\gamma\in F_q^*,\quad b_\gamma\in F_q.
\]
Multiplicative Hilbert 90 supplies lambda!=0 with gamma(lambda)=a_gamma*lambda. The elements b_gamma*lambda form an additive cocycle, so additive Hilbert 90 supplies mu with
\[
 \gamma(\mu)-\mu=b_\gamma\lambda.
\]
The distinct elements mu+lambda*i transform by the required q-point action. Their polynomial is
\[
 P(T)=T^q-\lambda^{q-1}T+
                  \lambda^{q-1}\mu-\mu^q=T^q+aT+b
 \in k[T],\quad a\ne0.
\]
At z=1/a, v=-b/a, or u=1+1/a, the full fibre of (7) is k[T]/P, canonically identified with A once the chosen affine labeling data is used. This proves the criterion for disconnected algebras and finite fields as well as field extensions.

The Hilbert 90 ingredients and generic additive-polynomial framework are established literature. The point here is the complete-fibre interpretation and its compatible canonical bundle form.

### An explicit canonical marked bundle construction

For g=(a,b) in H_q use these actions:
\[
 g_{\rm src}(u,x)=(u,ax+b),\qquad
 g_{\rm tgt}(u,v)=(u,av+bu).                         \tag{9}
\]
The first is affine linear; the second is linear. Since a^q=a and b^q=b,
\[
 \mathcal K_q\circ g_{\rm src}
       =g_{\rm tgt}\circ\mathcal K_q                 \tag{10}
\]
as polynomial identities on the entire planes. The origin in the target is H_q-fixed, and its complete fibre is
\[
 u=0,\quad x^q-x=0,
\]
the natural H_q set F_q.

For an H_q torsor P/S, form the associated affine-plane bundle E_P and rank-two vector bundle W_P using (9). Equivariance descends (7) to a polynomial morphism E_P->W_P. The determinants of the two linear parts are both a, so their relative determinant bundles identify; the Jacobian-one identity descends. Its fibre over the zero section is exactly P times^(H_q) F_q. All of this commutes with base change. Equivalently the equivariant map gives a representable etale morphism [A2_src/H_q]->[A2_tgt/H_q] over BH_q. Its pullback along the target zero section is [F_q/H_q], naturally equivalent to B(F_q^*) via the stabilizer of zero. Thus the universal degree-q root cover B(F_q^*)->BH_q is exactly the marked Keller fibre; for q=3 this is BS2->BS3.

Over a field both bundles admit trivializations, by the linear and affine forms of Hilbert 90; these trivializations are choices. Over BH_q the source has **no equivariant origin**: the subgroup of translations includes x->x+1, which has no fixed point. It is an affine bundle, not a vector bundle supplied with an origin. The lack of an equivariant origin persists after taking its product with any linear H_q representation: a fixed point in that product would project to a fixed source point. Thus vector-bundle stabilization cannot turn this canonical source into a vector bundle over BH_q.

This is the precise escape from the earlier linear-representation obstruction. Inertia on the marked fibre now contains the full H_q. The construction forgets the good generator line and changes the source-bundle type. It does not contradict UAINERT1, UANORMAL1 or the fixed-space theorem below.

### Universal cubics and an exact minimum in characteristic three

At q=3, H_3=S3 as permutation groups. Every three-element Galois set has the required affine structure, so (7) realizes every cubic finite etale algebra over **every** field of characteristic three, finite or infinite, and over every extension of that field.

There is also a short intrinsic generator proof. For A of rank three, the perfect trace pairing makes ker Tr two-dimensional; it contains the scalar line k*1 because Tr(1)=3=0. Choose theta in ker Tr outside that scalar line. Its characteristic polynomial is
\[
 T^3+d_1T+d_0.
\]
If d1=0, theta would be purely inseparable over k; an element of an etale algebra is separable, so it would be scalar, a contradiction. Thus d1!=0 and the characteristic polynomial is separable, making theta primitive. Taking z=1/d1 and v=-d0/d1 gives the whole algebra as a fibre.

With the positive-characteristic extension of the programme's field-versal definition, UA1DISC1 gives the matching lower bound:
\[
 \boxed{\operatorname{ktdim}_k(S3)=2
               \quad\text{for every }\operatorname{char}k=3.} \tag{11}
\]
The generic monodromy is S3 and the map is absolutely stably atomic. This is not a characteristic-zero plane Keller construction.

At q=4, H_4 is A4 in its natural action. Over fields containing F4, the same canonical affine-bundle construction covers quartic algebras equipped with an A4 reduction. For q>=4 it is not a universal S_q construction; the affine-structure condition must be retained.

### Relation to the coordinate-level exception

In characteristic three, if all retained indices are divisible by three then B=y, and restricting the promoted compiler to y=beta gives an actual Keller map on A^(r+2). This is the exceptional source in (2), not an omitted divisor completion. The generic inverse still has its retained top degree. It does not automatically realize every algebra in that degree; the retained coefficient support remains a separate arithmetic restriction.

For the cubic at beta=0, t=1,q=z and the restricted map is (z,2x-z*x^3), which becomes (z,x+z*x^3) after negating the second target coordinate in characteristic three. Thus the low-rank universal plane map is the explicit modular coordinate-level exception, up to elementary coordinate changes.

## 4. Fixed-space restrictions without a characteristic assumption

**UAFIXALL1.** For any field, any finite group G, and any G-equivariant polynomial Keller map F:V->W between linear representations, every subgroup-fixed restriction
\[
 F^H:V^H\longrightarrow W^H
\]
is Keller between equal-dimensional affine spaces. At a G-fixed target with finite etale fibre X, its complete fibre is exactly X^H.

The differential at the source origin is an equivariant isomorphism, so the dimensions agree. At every H-fixed point, its differential and inverse preserve invariant tangent vectors, giving an isomorphism on the fixed tangent spaces. The restricted determinant has no geometric zeros, hence is a polynomial unit. Fixed-point equalizers commute with the fibre product. These arguments require no averaging and remain valid when the characteristic divides |G|.

If m is the geometric number of H-fixed labels, the low-dimensional rules are:
\[
 \begin{array}{c|c}
 \dim V^H&\text{necessary fixed-label count}\\ \hline
 0&m=1\\
 1,\ \operatorname{char}k=0&m=1\\
 1,\ \operatorname{char}k=p&m=1\text{ or a positive multiple of }p.
 \end{array}                                        \tag{12}
\]
Indeed a univariate Keller map is a nonconstant finite polynomial, surjective over an algebraic closure, of degree one or a positive multiple of p; in characteristic zero its degree is one. In particular **m=0 always forces dim V^H>=2**, in every characteristic.

For G=S_N and its natural N-point fibre, N>=2, the full-group fixed set is empty. Thus the permutation representation bundle alone is excluded in every characteristic, including wild characteristic and characteristic two. In characteristic zero, dimension two still entails the earlier JC(2) counterexample gate; that conclusion is not imported into positive characteristic.

The affine source action in (9) has no full-group-fixed point and is not a linear representation with an origin. It therefore falls outside this theorem exactly as stated. For a finite affine-linear action in characteristic zero, averaging any orbit produces a fixed point and turns the affine space into a linear representation after translating that point. The modular source-without-origin construction does not provide this escape in characteristic zero.

## 5. Remaining gates and assurance

The B-level problem is closed for the independent promoted atlas, and together with the Pi and C arguments every moving target-coordinate normalization has a source obstruction in characteristic zero. Nonlinear coefficient sections and genuinely different completions remain open.

The canonical affine-bundle construction is effective over BH_q, and over BS3 for cubic algebras in characteristic three. A characteristic-zero presentation-free construction, or a marked construction using fixed linear representation bundles, still faces the existing inertia, fixed-space and nonlinear-moving requirements.

The [replay](verify.py) gives exact representative counts, whole cubic fibre algebras, covariance identities and finite-field group controls. The all-degree spreading, Hilbert 90 descent, canonical associated-bundle and atomicity arguments are written proofs. Formalization, independent whole-proof replay, external review and novelty priority are not claimed.
