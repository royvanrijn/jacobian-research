# Winding-one phases through sphere degree six

Let K have characteristic zero and put
\[
 A_K=K[x,y,t]/(t^2+2xy-1),\qquad
 \Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^1
 f(x,(1-t^2)/(2x),t)\,dt .
\]
The integral is the rational polynomial rule. Pure means all positive power moments vanish, eventually pure means all sufficiently large power moments vanish, and full-marker safe means every fixed Q in A_K has eventual mixed vanishing. Degree is the sphere filtration degree.

## 1. The exact classification

Consider
\[
 f=U(t)+xS(t)+yR(t),\qquad
 \deg U\le6,\quad \deg S,\deg R\le5.                    \tag{1}
\]
The zero polynomials are allowed. This is the full angular-support interval [-1,1] inside sphere degree at most six, a nineteen-dimensional vector space.

**Theorem 1.** For a phase in (1), purity, eventual purity and full-marker safety are equivalent. Their exact common locus consists of:

1. U=0 and at least one of S,R is zero;
2. f=a y+b t+c x with b^2+2ac=0.

The overlap and zero phase are included. Every such safe phase obeys
\[
                    \Lambda(Qf^m)=0\qquad(m>\deg Q) .   \tag{2}
\]
In particular this angular-support sector contains no pure phase unsafe for the full sphere functional. No reduction of arbitrary quintics to (1) is asserted.

The proof combines an elementary exposed-face obstruction with ordinary reduction of cyclic covers of genus at most five, whose Jacobians decompose into factors of dimension at most two. It uses only the abelian-surface ordinary theorem among the deep ordinary inputs; no full-endomorphism/simple-signature hypothesis is needed.

## 2. Two exposed monomials or one zero-charge vertex

Assume first that coefficients are in a number field, and S,R are nonzero. Write
\[
 r=\deg R,\quad s=\deg S,\quad
 a=[t^r]R\ne0,\quad c=[t^s]S\ne0,\quad
 H=(r+s+2)/2 .                                        \tag{3}
\]
Let F(x,t) be the literal Laurent polynomial after eliminating y. For an angular exponent i and height j, assign weight j+rho*i, where
\[
                        \rho=(r+2-s)/2 .
\]
The positive and negative ends attain exactly the common height H:
\[
            c x t^s,\qquad -(a/2)x^{-1}t^{r+2}.        \tag{4}
\]
Every other term of xS+yR has smaller weight. This includes the lower term a*x^{-1}t^r/2 from the factor 1-t^2.

Suppose deg U<H, including U=0. For m=2k, angular extraction has t-degree at most k(r+s+2), and its coefficient there is exactly
\[
 \binom{2k}{k}(-ac/2)^k .                              \tag{5}
\]
Choose arbitrarily large good rational primes p=1 modulo L=r+s+2 and k=(p-1)/L. Since L>=2, 2k<=p-1; hence (5) is a unit at a prime ideal above p. In the normalized interval integral of the angular constant term, the only denominator divisible by p is p itself, at t^{p-1}. Thus p*Lambda(f^{2k}) reduces to (5), and that moment is nonzero. The orders are unbounded.

If q=deg U>H, choose a rational rho strictly between r+2-q and q-s. Then the unique maximal-weight monomial of F is u_q*t^q, at angular exponent zero. The t-degree of CT F^m is at most qm and its leading coefficient is u_q^m. At good primes p=1 modulo q, take m=(p-1)/q. The same unique-denominator argument gives infinitely many nonzero moments. Here q>=2.

It remains only the tie q=H, which requires the integer H in {1,2,3,4,5,6}. If r+s is odd no tie is possible.

The required infinitude of primes1 modulo any fixed positive integer uses the elementary cyclotomic argument, as in the angular-face classification. All coefficient denominators and specified nonzero coefficients exclude only finitely many primes.

## 3. The tied face and its exact cyclic-cover scalar

For q=H define b=[t^H]U. The maximal face is
\[
              F_\infty=t^H L(z),\qquad
              L(z)=-a/(2z)+b+cz,\qquad z=x/t^\rho .     \tag{6}
\]
Here rho is an integer. For every m, all angular-constant terms have height at most Hm, and
\[
               [x^0t^{Hm}]F^m=\operatorname{CT}_zL(z)^m .
                                                               \tag{7}
\]
Thus at good odd p=1 modulo H and m=(p-1)/H,
\[
                 p\Lambda(f^m)\equiv
                         \operatorname{CT}L^m \pmod{\mathfrak p}.
                                                               \tag{8}
\]
All lower coefficients are unrestricted.

For H>=2 put
\[
 h_H(u,v)=u^{H-1}v^{H-1}
                 \bigl(cu^2+buv-a v^2/2\bigr).          \tag{9}
\]
Its degree is 2H and z^{-H}h_H(z,1)=L(z). Since a,c are nonzero, the quadratic has no root at zero or infinity.

If the quadratic has a double root alpha, necessarily alpha!=0. Then L(z)=c*z^{-1}(z-alpha)^2, so CT L^m=c^m*binom(2m,m)*(-alpha)^m. For m=(p-1)/H and H>=2 one has2m<=p-1, and this is a unit at all sufficiently good primes. No ordinary theorem is used in this repeated-root branch.

When the quadratic has distinct roots, let C_H be the smooth projective normalization of
\[
                              w^H=h_H(u,v),             \tag{10}
\]
where w has homogeneous weight two. It is connected and fully ramified at four distinct base points with multiplicities H-1,H-1,1,1. Its genus is H-1. All branch multiplicities are smaller than H, so the character-one eigensheaf of the normalized covering algebra is O(-2), with one-dimensional H^1. At an absolute degree-one good prime p=1 modulo H, its Frobenius scalar on the Cech generator w/z is
\[
 [z^{p-1}]h_H(z,1)^{(p-1)/H}
                  =\operatorname{CT}L(z)^{(p-1)/H}.     \tag{11}
\]
This follows by raising w/z to its p-th power and projecting to H^1(O(-2)). Ordinary reduction of the Jacobian makes Frobenius bijective, so (11) is nonzero. This is a sufficient ordinary condition for the exact scalar needed in (8), not an assertion that every cyclic-cover character is one-dimensional.

## 4. Every tied cover has ordinary primes after finite extension

The only deep ordinary input is Sawin, arXiv:1506.04784v2: every abelian surface has density-one ordinary primes after a finite extension. It also covers an elliptic curve by applying it to its square. Products and isogenies preserve ordinary reduction away from finitely many primes.

For H=2 the genus is one; for H=3 it is two. Sawin's input therefore applies immediately.

For H=4 the genus is three. The square of the deck generator is an involution, and its quotient is the normalization of
\[
             v^2=z^3(cz^2+bz-a/2),
\]
equivalently, after dividing v by z,
\[
                         v_1^2=z(cz^2+bz-a/2).
\]
The latter is a smooth genus-one curve because its three finite roots are distinct. Let iota be this involution on the genus-three curve. On its Jacobian, the images of1+iota and1-iota have dimensions one and two: the invariant H^1 is the pullback from the genus-one quotient. Their sum is the entire Jacobian up to isogeny since their sum as operators is2, and their intersection is finite since iota acts with opposite signs there. Thus J(C_4) is isogenous to an elliptic curve times an abelian surface after finite extension. Sawin's input again supplies ordinary primes.

For H=5 the branch exponents are 4,4,1,1. Pair each exponent4 with an exponent1. Every four-point configuration has a projective involution interchanging the two chosen pairs. It negates the branch exponents modulo5 and lifts, after finite extension, to an automorphism tau inverting the deck generator g. Since tau^2 is a deck transformation and conjugation by tau both fixes its own square and inverts that odd-order group, tau^2=1. The two base fixed points are unramified, and the reflection on each five-point fibre fixes exactly one point. Riemann--Hurwitz gives genus two for C_5/<tau> and C_5/<g*tau>. Their pullback H^1 spaces each have dimension four, and their intersection is deck-invariant and zero. Their Jacobians thus have product isogenous to J(C_5). This is the complementary-exponent four-root argument, with no generic endomorphism assumption.

For H=6 the cover has genus five. Let g be its deck generator. The quotient by the order-two subgroup <g^3> is a normalized cyclic-cubic cover with branch exponents2,2,1,1, of genus two. The quotient by the order-three subgroup <g^2> is a normalized double cover branched at four points, of genus one. On J(C_6), the rational algebraic idempotents
\[
 e_2=(1+g^3)/2,\qquad e_3=(1+g^2+g^4)/3
\]
have these two quotient images. Their product is the full deck average, zero because the deck quotient is P1. The remainder1-e_2-e_3 is an algebraic rational idempotent. For each deck character i=1,...,5, the eigensheaf formula gives
\[
 -1+2\{i/6\}+2\{5i/6\}=1
\]
for its dimension on H^1(O), so each character has Betti dimension two. The remaining characters are1,5, with total Betti dimension four; after clearing the idempotent denominators their image is an abelian surface. The three images have product isogenous to J(C_6). Sawin's input applies to all three factors.

Enlarge the number field to contain the covering data, the needed isogenies and zeta_H, and to have density-one ordinary primes for all the stated factors. Absolute degree-one primes have density one; their rational characteristics are1 modulo H. Excluding finitely many bad primes leaves infinitely many ordinary primes usable in (8)--(11). Consequently every tied case H>=2 has infinitely many nonzero moments.

## 5. The linear exception, one-sided branches and all coefficient fields

If the tie has H=1, then r=s=0 and U=b t+d. Thus
\[
                         f=d+a y+b t+c x .
\]
For L=a y+b t+c x, formal Gaussian contraction or orthogonal trace invariance gives
\[
 \Lambda(L^{2n+1})=0,\qquad
 \Lambda(L^{2n})=(b^2+2ac)^n/(2n+1).                    \tag{12}
\]
If b^2+2ac!=0, choose its square root over a finite extension. The moments of d+L are exactly those of the nonzero affine polynomial d+sqrt(b^2+2ac)*t under the interval rule. The one-variable interval-moment theorem of Mueger--Tuset, arXiv:1910.12148, Theorem2.1, gives infinitely many nonzero moments. If b^2+2ac=0, all pure moments are d^m, so eventual purity forces d=0. A nonzero isotropic L is carried to a scalar multiple of x by a linear sphere isometry over an extension. The trace is invariant, and moving the isometry onto a fixed marker proves (2), since a degree-d marker has angular weights in [-d,d].

If S=0 or R=0, the angular constant term of f^m is U(t)^m. Mueger--Tuset gives eventual purity exactly when U=0. The remaining one-sided phases have angular weight exactly m or -m in their m-th powers, so (2) follows directly.

For an arbitrary characteristic-zero field, assume eventual purity and generate a finitely generated Q-algebra from all coefficients. In each alleged excluded branch, invert the finitely many leading coefficients, the nonzero discriminant when used, and the coefficients needed to preserve the indicated degrees and inequalities. ONE number-field specialization preserves that branch and all eventual moment identities. The preceding contradiction then applies. For the repeated-root and isotropic branches the defining equalities are polynomial identities and survive specialization. In the allowed branches safety can be proved directly after an algebraic extension, and descends by the faithful scalar extension of A_K and Lambda. This proves Theorem1 over every K.

## 6. The complete linear-negative part of the stripped quintic chart

Inside
\[
 S_K=K[y,t,q],\qquad q=xt^2,
\]
consider the full class
\[
 f=U(t)+yR(t)+qQ(t),\qquad
 \deg U\le5,\quad\deg R\le4,\quad\deg Q\le2.             \tag{13}
\]
It has dimension fourteen. Taking S=t^2Q places it in (1), in its degree-at-most-five subspace (which has dimension sixteen). A nonzero S then has degree at least two, so the mixed isotropic-linear exception cannot occur. Hence purity, eventual purity and full-A marker safety are equivalent to
\[
                U=0,\qquad R=0\ \text{or}\ Q=0 .       \tag{14}
\]
All safe phases obey (2). Its even-t subspace has dimension eight and obeys the same exact classification.

This is the ENTIRE linear-in-y part of the stripped chart, allowing arbitrary height coefficients and both parities. It does not contain nonlinear negative winding y^2R_2(t)+..., and does not settle the even14-dimensional/full24-dimensional stripped chart or the global primitive Image quintic question.

## 7. Source comparison and proof limits

The exposed-face method and full-ambient marker cutoff extend the written univariate angular-face argument. The tied coefficient uses the same normalized-cover Cech calculation as the spherical leading-form proof. The H=5 complementary-exponent splitting and H=6 rational-idempotent decomposition are also written in full above. Ordinary reduction of arbitrary abelian surfaces and Mueger--Tuset's one-variable interval theorem are the only external mathematical premises. Neither Cantoral's simple-signature theorem nor a CM-subfield claim is used. The exact primary hypotheses, including the ordinary-after-finite-extension quantifier, must be retained in any application.

Every exclusion is an all-order infinite-nonvanishing argument. No bounded moment search, generic-position assumption, or formalization claim is substituted for the proof.
