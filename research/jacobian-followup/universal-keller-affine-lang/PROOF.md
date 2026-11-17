# Affine canonical bounds, automatic atomicity and Lang realizations

14 November 2026. This extends the [symmetry classification](../universal-keller-symmetry/PROOF.md) and the retained [affine inertia theorem](../monodromy-theory/THEORY.md). All maps and fibres below are whole polynomial maps and complete scheme-theoretic fibres. The category and degree hypotheses are part of each assertion.

## 1. Affine sources cannot make the high-rank canonical cost small

**UAAFFCAN1.** Let N>=5. Suppose an S_N-equivariant polynomial Keller map between d-dimensional affine representations has a fixed target whose complete marked fibre is the natural N-point S_N set. Then
\[
 \boxed{d\ge N+1\quad\text{if the characteristic does not divide }N;
 \qquad d\ge N\quad\text{otherwise}.}              \tag{1}
\]
In the divisible-characteristic case, d>=N+1 still holds if the source has an equivariant origin. No assumption on the generic degree is made. These are marked affine-bundle bounds, not unmarked field-versal dimension bounds.

### Barycentres and the actual affine orbit

Work geometrically. If N is invertible, the barycentre of the N fibre points is a fixed source point, even if the characteristic divides N!. Translate it and the marked target to origins and apply UACANCOST1. The same argument applies whenever a source origin is already fixed.

Now suppose p divides N. Extend the affine source action on V to the linear action on k plus V, and map the permutation module M=k^N by
\[
 e_i\longmapsto(1,x_i).
\]
The kernel is zero or k*1: a vector with unequal coordinates in the kernel would, by its transposition difference, identify two fibre points. Thus the affine span P of the orbit is one of
\[
 \{\textstyle\sum z_i=1\}\subset M,\quad\dim P=N-1,
 \qquad
 \{\textstyle\sum z_i=1\}/(k*1),\quad\dim P=N-2.    \tag{2}
\]
The second quotient is defined because N=0 in k.

Choose 2<=r<=N-2 with p not dividing r: use r=2 for odd p and r=3 for p=2 (then N>=6). Let H=S_r times S_(N-r). It fixes no natural label. In the first model P^H is an affine line: the two constant block coordinates c,d satisfy r(c-d)=1. In the second model it is one point, after quotienting by simultaneous constant translation. There are no extra fixed classes in the quotient: every block transposition difference has zero coordinates outside its two entries, so cannot be a nonzero scalar vector.

The ambient fixed locus is nonempty because P^H is nonempty. Its dimension t satisfies
\[
 t\le\dim P^H+d-\dim P=d-N+2.                      \tag{3}
\]
At a fixed source point the differential identifies source and target H-fixed tangent spaces. Their complete fixed loci are affine spaces, and the restriction is Keller by UAFIXGEOM1. Its marked fibre is empty, so t>=2 in every characteristic. Equation(3) gives d>=N.

This proof uses one two-block stabilizer uniformly; no catalogue of small modular representations is required. The rank-three characteristic-three and rank-four characteristic-two examples below explain why the N>=5 boundary matters. At high rank, losing an origin saves at most one coordinate from the linear lower bound, and is possible only when the characteristic divides the rank.

## 2. Missing sheets bound affine divisorial inertia

**UAINERTIACOST1.** Let F:A^d->A^d be Keller over an algebraically closed field, with geometric generic degree n and faithful monodromy on n sheets. For an affine target prime divisor D let r_D be the number of points in its complete geometric generic fibre. Its inertia group I_D satisfies
\[
 \boxed{1\le r_D\le n,\qquad
        I_D\text{ embeds in }S_{n-r_D}.}            \tag{4}
\]
The bounds and fibre ranks persist under polynomial left-right changes and identity stabilization. In particular:

- A loss of at most one sheet permits no nontrivial inertia.
- A loss of exactly two sheets, if ramified, forces a tame transposition when the characteristic is not two.
- Nontrivial wild inertia in characteristic p requires loss of at least p sheets.
- In characteristic zero a noninvertible map has some ramified affine divisor and loses at least two sheets there. If its monodromy is contained in A_n it loses at least three at some such divisor.
- In odd characteristic, a nonsquare geometric discriminant forces a divisor with odd sign inertia and r_D<=n-2.
- If the generic group is F_q semidirect C_m in its natural affine action, m>1 dividing q-1, at least one ramified divisor has r_D=1. In fact every ramified affine divisor has rank one and cyclic inertia in a conjugate of C_m.

### Proof and the difference between ramification and omitted branches

The divisorial argument in KMFI1 works in every characteristic: for D=V(h), h(F) is nonconstant, and every source prime factor gives a divisor dominating D. Its local map is etale. This proves r_D>=1. After strict henselization each actual source point gives a separate unramified sheet fixed by inertia. Hence I_D fixes at least r_D labels and acts faithfully on the remaining n-r_D, proving(4).

There may be omitted unramified branches in the finite normalization, so equality between r_D and the number of inertia-fixed sheets is not asserted. If a p subgroup acts nontrivially on a finite set it moves at least p points. The remaining small-support assertions follow immediately. In characteristic zero KMFI1 says affine inertia normally generates the geometric group, giving a nontrivial I_D unless the map is birational and hence an automorphism. Nonidentity even permutations move at least three labels.

For the odd-characteristic discriminant assertion, factor its square class in the rational function field of affine space. If all divisorial valuations are even, it is a square over algebraically closed constants. A nonsquare class therefore has odd valuation at some D, giving a nontrivial sign inertia character.

For the affine group with m>1, its nontrivial C_m quotient likewise has some divisorial Kummer ramification. Every affine inertia group fixes a label, so it lies in that label's stabilizer C_m. A nontrivial multiplier fixes exactly one label. Thus each nontrivial I_D is tame cyclic and r_D=1. A constant cyclic quotient could evade the existence argument before closing constants; the stated group is geometric.

A useful control is the translation-family map (u,x)->(u,x-u^(q-1)*x^q). Its normalization is
\[
 z-z^q=uv,
\]
finite etale over the whole target. The original affine source is the open retaining z=0 over u=0 and omitting the other q-1 unramified branches; z=ux. Thus the actual boundary fibre has rank one while affine inertia is trivial. A rank drop by itself is not a ramification certificate.

### A characteristic-free Kummer consequence

For all integers a>=2, the same divisor pullback argument gives
\[
 b(F)\in k(X)^{*a}\ \Longrightarrow\ b\in k(Y)^{*a}.              \tag{5}
\]
Indeed target prime divisors pull back with nonempty disjoint squarefree supports. Compare valuations at one source prime over each, and then use algebraically closed constants. This is the retained KMFI1 power-injectivity argument, extended to every characteristic, not a new priority claim for it. Kummer theory implies that the generic source field has no nontrivial cyclic prime-to-characteristic Galois intermediate field. In particular a nontrivial regular C_ell action with ell prime to the characteristic is excluded.

## 3. Full-degree canonical symmetry usually supplies atomicity automatically

**UACANATOM1.** Suppose a polynomial Keller map has generic degree N>=2, a complete N-point natural marked S_N fibre at a fixed target, and compatible S_N actions by polynomial automorphisms on source and target. Then its geometric monodromy K is a transitive normal subgroup of S_N. It is primitive, and the map absolutely stably atomic, except for the possible case
\[
 \boxed{N=4,\quad\operatorname{char}k=2,\quad K=V_4.}             \tag{6}
\]
In characteristic zero the possible groups are S3 at rank three and A_N or S_N for N>=4; rank two is impossible. This does not require affine or linear group actions. The full generic degree hypothesis is essential: a marked rank-N boundary fibre of a higher-degree map is not covered.

### Why the marked target lies on a finite chart

Use Zariski's Main Theorem to put the source into the finite normalization of the affine target in its function field. Over the strict henselization at the marked target, the N distinct etale source points supply N rank-one factors. They consume the entire generic degree N. A further normalization component would have positive generic rank, impossible. Thus no normalization point over the target is omitted. Finiteness holds on a neighbourhood; intersect its finitely many S_N translates to obtain an invariant finite-etale open containing the fixed target.

With that target as basepoint, functoriality of the etale fibre functor gives
\[
 \rho(g_*\gamma)=g\rho(\gamma)g^{-1}.
\]
Thus the S_N label action normalizes K. Since it is the full symmetric group on those same N sheets, K is normal in S_N. Connectedness of the source gives transitivity.

For N>=5 the only possibilities are A_N and S_N, both primitive. At rank three C3 and S3 are primitive. At rank two C2 is primitive. At rank four the transitive normal possibilities are V4,A4,S4, and only V4 is imprimitive. If the characteristic is not two, a regular V4 field would contain a cyclic quadratic Galois intermediate field, excluded by(5). The same argument excludes regular C2 in odd characteristic and C3 outside characteristic three; KMFI1 excludes all nontrivial regular actions in characteristic zero. Apply the characteristic-free primitive-monodromy argument retained in UAALLATOM1.

The exception is real: in characteristic two, (x,y)->(x-x^2,y-y^2) has natural S4 equivariance through AGL2(F2), complete degree four, geometric group V4, and two quadratic polynomial Keller factors. A generic imprimitive group alone is still not a factorization theorem in other examples.

## 4. Lang maps give canonical affine-space local systems in every rank

Let q=p^e and k contain F_q. Write G=AGL_n(F_q) on N=q^n labels.

**UALANG1.** The map
\[
 \Phi_{n,q}(X)=X-X^{[q]}:A^n_{\rm aff}\longrightarrow A^n
                                                               \tag{7}
\]
is finite etale of degree q^n with Jacobian one. Its source affine G action and target linear action satisfy Phi(AX+b)=A Phi(X). Twisting gives a canonical marked affine-bundle realization of every affine F_q-space local system, without a primitive generator or source origin. It is decomposable when en>1, and absolutely stably atomic when en=1.

A fixed field-versal G-framed Keller map with absolute stable atomicity is
\[
 \boxed{\mathcal L_{n,q}(U,X)=(U,X-U X^{[q]})}                    \tag{8}
\]
on ordinary affine space of dimension n^2+n. Its arithmetic and geometric groups are G. It is finite etale of rank q^n on det(U)!=0; its source and target are whole affine spaces globally. On a base F_q, the finite-field source budget gives dimension at least n+1 for fixed framed universality whenever GL_n(F_q) is nontrivial. No matching minimum is asserted in general.

The canonical minimum for(7)'s group is exactly n over every base containing F_q when q>2, and over F2 when q=2. For q>2 choose a prime ell dividing q-1. The label stabilizer contains (C_ell)^n. At its fixed labelled source point its faithful affine action becomes a faithful linear tangent action; diagonal characters require at least n dimensions. Over F2 the q^n split points require at least n coordinates. For q=2 over larger bases, the general minimum is not deduced here. The case n=2 has exact minimum two over every characteristic-two base by UA1CANCLASS1.

### A whole-cover frame and every-field versality

On det(U)!=0, the equations X-U X^[q]=V are equivalent to monic equations X_i^q=(U^-1(X-V))_i. Their coprime leading monomials give a free basis with 0<=exponent_i<q. The Jacobian is the identity. Lang's theorem gives an invertible matrix B with U B^[q]=B; the kernel is B*F_q^n. Its affine frames are pairs(B,mu) with
\[
 UB^{[q]}=B,\quad\det B\ne0,\quad\mu-U\mu^{[q]}=V.              \tag{9}
\]
The frame changes (B,mu)->(BA,mu+Bb) make this a G torsor, not merely a splitting-field description.

For a G cocycle(A_gamma,b_gamma) over any L/k, matrix Hilbert90 gives gamma(B)=B A_gamma. Additive Hilbert90 gives gamma(mu)-mu=B b_gamma. The quantities
\[
 U=B(B^{[q]})^{-1},\quad V=\mu-U\mu^{[q]}                         \tag{10}
\]
lie in L and give the prescribed frame. This includes disconnected torsors and finite fields. Applying this after closing constants to a generic connected G torsor forces the geometric monodromy to be full G; otherwise every fibre would reduce to one proper generic subgroup. The affine action is doubly transitive, hence primitive, so(8) is absolutely stably atomic.

For canonical descent of(8), write T=U-I. The actions are
\[
 (T,X)\mapsto(ATA^{-1},AX+b),\quad
 (T,V)\mapsto(ATA^{-1},AV-ATA^{-1}b).                              \tag{11}
\]
The source action is affine and the target action linear. Both determinant characters are det A. The fixed target(T,V)=(0,0) has the natural F_q^n fibre, and translations prevent any source origin, even after vector-bundle stabilization.

The elementary canonical map(7) has regular translation monodromy. For n>1 its coordinate product factors; for n=1,e>1 use the trace factorization of X-X^q. Coefficient promotion to(8) changes geometric monodromy; identity stabilization does not.

### The basis-free receiver, including the recent literature extension

More generally take a rank-n vector bundle E with a specified isomorphism Frob_q^*E->E, and an affine E torsor Z with compatible specified Frobenius-periodicity isomorphism theta:Frob_q^*Z->Z. Then
\[
 z\longmapsto z-\theta(z^{[q]})\in E                              \tag{12}
\]
is a globally defined finite etale degree-q^n map with derivative the identity. Local polynomial expressions glue without selecting an origin. Its zero scheme is the corresponding affine F_q local system. The periodicity isomorphisms are data, not a canonical choice supplied merely by knowing that a bundle is periodic.

This is the affine-group instance of the classical Lang correspondence and of Lei Zhang's September2026 classifying-stack formulation (Theorem I). That paper permits arbitrary fibered categories and imperfect ground fields. Over a characteristic-p base not containing F_q, the finite etale group scheme need not be constant; it must not be silently identified with the constant group AGL_n(F_q). The direct formulas above concern constant groups over bases containing F_q.

## 5. A globally finite universal quartic map in characteristic two

**UA4CHAR2.** Over every characteristic-two field, the one fixed map
\[
 \boxed{\mathcal Q(s,t,x,y)=
 (s,t,\ x+(1+st)x^2+s y^2,\ y+t x^2+y^2)}                         \tag{13}
\]
is globally finite etale of rank four, has Jacobian one, arithmetic and geometric S4 monodromy, realizes every quartic etale algebra over every extension field, and is absolutely stably atomic. Thus
\[
 \operatorname{ktdim}_{F2}(S4)\in\{3,4\}.                        \tag{14}
\]
The lower bound is UAFINCOUNT1, not a claimed matching minimum. Over any characteristic-two field the fixed-map dimension is between two and four.

Separately, the canonical marked affine-bundle dimension for quartic algebras in characteristic two is exactly two, by(7) with(n,q)=(2,2) and S4=AGL2(F2). Every two-dimensional canonical marked affine-bundle realization is decomposable; its absolute-atomicity minimum therefore lies between three and six, the upper bound supplied by(8). Neither the four-dimensional fixed map nor its basis normalization is asserted to be a canonical bundle construction.

### A two-parameter normal form for every rank-two Frobenius module

Set
\[
 A(s,t)=\begin{pmatrix}1+st&s\\t&1\end{pmatrix},\qquad\det A=1.
                                                               \tag{15}
\]
Every A0 in GL2(L), for every characteristic-two field L, is Frobenius-equivalent to one of(15): A0~C^-1 A0 C^[2]. Put D=det A0 and first use C0=diag(D^-1,1); the new determinant is one.

Write the resulting matrix as [[a,b],[c,d]]. If b!=0, use
\[
 C_1=\begin{pmatrix}1&0\\(d-1)/b&1\end{pmatrix}.
\]
The bottom-right entry of C1^-1*A0*C1^[2] is one, and its determinant stays one. If b=0, d!=0 and C1=diag(d,d^-1) gives the same conclusion. Every determinant-one matrix with bottom-right entry one has form(15). These two rational charts use no root extraction and cover all fields, including F2.

Start with the affine S4 cocycle construction(10) for n=2,q=2 and replace B,mu by C^-1 B,C^-1 mu. Then its U has form(15), and its V lies in L^2. This realizes the original S4 torsor as a whole fibre of(13). Since the natural affine permutations of F2^2 are all24 permutations, this covers every quartic etale algebra, without an additional structure choice.

Because A^-1 is polynomial and det A=1, the monic equations for x^2 and y^2 give the free basis 1,x,y,xy over the entire target ring. Thus the map is globally finite etale of rank four. The generic frame is S4 by every-field versality as above; primitivity gives absolute stable atomicity.

Over F2 all16 complete fibres have the following splitting-type distribution:
\[
 1^4:1,\quad 2+2:3,\quad2+1+1:4,\quad3+1:4,\quad4:4.            \tag{16}
\]
The replay verifies these as complete quotient algebras using Frobenius point counts, not by treating a list of rational preimages as the whole fibre.

A fixed univariate quartic Keller polynomial in characteristic two has form aX^4+bX^2+cX+d, c!=0. Over closed constants its geometric group is V4; the three nonzero differences are constant. Its finite-etale S4 frame therefore cannot realize a connected generic S4 torsor over every extension. This proves the general fixed-map lower bound two without applying the odd-characteristic discriminant argument.

### Why the two-dimensional canonical map cannot be absolutely atomic

Let an affine S4 action on a plane have a natural four-point orbit. Its linear part kills V4: if the commuting unipotent V4 acted nontrivially in dimension two, its common fixed line would be S4-stable, making the whole representation upper triangular. Its diagonal characters are trivial in characteristic two because S4_ab=C2, so its image would be abelian and kill A4, a contradiction. The affine action is faithful on labels; its translation kernel is therefore exactly V4, since a larger normal subgroup cannot inject into the additive plane.

Shift one marked point to zero. The other points are its V4 translation vectors W. They span the plane: otherwise the faithful S3 stabilizer action on W would be scalar, impossible. Choosing an F2 basis of W identifies it with F2^2. The target affine action has a fixed marked point, hence is linear after translation, and also kills V4. Any equivariant polynomial map is thus invariant under the two independent translations and factors through
\[
 (x,y)\longmapsto(x^2+x,y^2+y).
\]
That map has two nontrivial quadratic polynomial Keller factors. The remaining factor is Keller by the chain rule. Hence every such canonical plane realization is decomposable, regardless of a possible larger generic degree. This does not classify arbitrary nonlinear plane actions or the unknown canonical dimensions three through five.

### Why a finite quartic map does not contradict the tame obstruction

The odd-characteristic UAFINITE1 obstruction does not apply in characteristic two. S4's only nontrivial abelian quotient is C2, now wild. Conversely A4 has a tame C3 quotient, so a globally finite A4-framed field-versal family is excluded by UAFINITE1. Map(13) still realizes every individual quartic algebra induced from an A4 torsor. It has no global A4 reduction of its generic frame. This is an actual control demonstrating why the framed-group qualification cannot be dropped.

## 6. What remains open

The affine canonical bound and normality lemma remove two independent uncertainties: wild origins cannot give low-dimensional natural S_N models at high rank, and full-degree canonical symmetry normally supplies atomicity automatically. They do not construct the characteristic-zero receiver.

The Frobenius-periodic receiver gives effective canonical descent for every affine finite-field local system, including all quartics in characteristic two. The fixed universal quartic map is finite and atomic in four coordinates, but the three-dimensional finite-base minimum and canonical absolute-atomicity dimensions three through five remain open.

General unmarked characteristic-zero dimension minima, more general non-affine or nonlinear receivers, and affine completions outside the retained families remain open. The [literature map](LITERATURE_AND_GAPS.md) records classical inputs and the recent stack formulation. The [replay](verify.py) supplies exact controls, not formalization, independent whole-proof review or novelty certification.
