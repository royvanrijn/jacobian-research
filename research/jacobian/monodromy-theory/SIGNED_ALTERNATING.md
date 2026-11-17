# Actual signed-alternating Keller actions in affine four-space

Work geometrically over an algebraically closed characteristic-zero field.
For an even integer n>=6, put

\[
 V_n=\{(\epsilon_1,\ldots,\epsilon_n)\in{\bf F}_2^n:
             \sum\epsilon_i=0\},\quad G_n=V_n\rtimes A_n.
\]

Its action is on the 2n signed letters (+1,-1,...,+n,-n):
A_n permutes the pairs, and V_n flips an even number of pairs.
The action is faithful; it is transitive because A_n moves between
pairs and a two-pair flip exchanges the signs in any chosen pair.

**KMFI29.** Every one of these faithful transitive actions is realized
by a polynomial Keller self-map of ordinary affine four-space.
The maps below have generic degree 2n and an explicit nonzero
constant Jacobian. Thus occurrence is no longer confined to the
previous symmetric and alternating natural actions and product examples.

Write n=2m+2, m>=2. The coefficient field is
K_m=Q(sqrt(-3(2m+1))).
The arithmetic group is the inverse image, in C2 wr S_n, of

\[
 \operatorname{Gal}\bigl(
 K_m(\sqrt{(-1)^m(2m+1)},\sqrt{-(2m+1)})/K_m\bigr)
 \subset C_2\times C_2,
 \tag{1}
\]

under the characters (sign of the block permutation, parity of flips).
This displayed field is exactly the constant field of the normal closure.

In particular, n=6 gives degree twelve, geometric group 2^5:A6
of order 11520, arithmetic group C2 wr S6 of order 46080,
and constants Q(sqrt(-15))(i,sqrt(5)) over its coefficient field.

**KMFI30.** These maps are absolutely and stably atomic under polynomial
Keller composition, although their monodromy action is imprimitive.
Their only proper nontrivial block system is the system of n pairs.
Its intermediate field has source degree two and cannot be
algebraized as an intermediate affine-space Keller factor.

**KMFI31.** For the signed action, the lcm C of all permutation
cycle-length gcds is n when n is a power of two, and 2n otherwise.
Consequently the actual maps with n=2^a, a>=3, admit no independently
stabilized curve-cover presentation. For other n this invariant
passes and gives no such exclusion.
This uses the curve-only KMFI9--10 obstruction.

The earlier delta-two checkpoint is now completed by this construction.
No assertion of novelty priority, independent whole-proof review,
formal verification or a complete Keller action classification is made.
The specified PSL2(11) cases remain open.

## 1. Explicit coefficients, surface and marked point

Choose either root alpha of

\[
 A_m\alpha^2+B_m\alpha+1/24=0,\qquad
 A_m={m^2\over2(2m+1)(m+1)^2},\quad
 B_m={m\over2(m+1)(m+2)}.
 \tag{2}
\]

Its discriminant is

\[
 -{m^2(m-1)^2\over12(m+1)^2(m+2)^2(2m+1)}.
\]

Thus alpha is nonreal and generates K_m. One explicit choice is

\[
 \alpha=-{(2m+1)(m+1)\over2m(m+2)}
 +{(m-1)(m+1)\over6m(m+2)}\sqrt{-3(2m+1)}.
\]

Set

\[
 a(q)=\alpha q^m+q,\quad
 b(q)={\alpha m\over m+1}q^{m+1}+q^2/2,
\]
\[
 d(q)={\alpha^2m\over2m+1}q^{2m+1}
       +{\alpha(m+1)\over m+2}q^{m+2}+q^3/3.
 \tag{3}
\]

These satisfy b'=q a' and d'=a b'.
On source parameters w,z put q=w^2,sigma=a(q)w-z and

\[
 X=(\sigma,b(q),d(q)-\sigma^2/2),\quad
 S(\gamma,w,z)=X+\gamma(1,w,z).
 \tag{4}
\]

Its Jacobian is gamma^2.
The origin has the marked preimage

\[
 w_0=1,\quad z_0=a(1)-b(1),\quad\gamma_0=-b(1)\ne0.
 \tag{5}
\]

Indeed the remaining equation is
d(1)-a(1)b(1)+b(1)^2/2=0, exactly (2).
The nonvanishing follows since alpha is nonreal.

## 2. A completely specified polynomial affine lift

Use stage variables theta=(gamma,u,v) and define polynomial quotients

\[
 E_i(\theta)=S_i(\gamma,\gamma u,\gamma^3v)/\gamma^{e_i},
 \qquad(e_1,e_2,e_3)=(1,2,4).
 \tag{6}
\]

These are polynomials: the source weights of w,z are 1,3;
the lowest weights of X1,X2,X3 are respectively at least 3,4,6,
and those of gamma(1,w,z) are 1,2,4.

Put
theta0=(gamma0,1/gamma0,z0/gamma0^3).
Every E_i vanishes there.
Let J be the three-by-three matrix DE(theta0), L=J^{-1},
and c_i the columns of L.
The determinant chain at the marked point gives
det J=gamma0^{-1}, so all these data are defined over K_m.

Write directional derivatives of E3 at theta0 as
H_ij=D_(c_i)D_(c_j)E3 and T111=D_(c1)^3 E3.
Define the three explicit constants

\[
 f_{20}=-H_{11}/2,\quad f_{11}=-H_{12},\quad
 f_{30}=-T_{111}/6-f_{20}H_{13}.
\]

The polynomial stage is

\[
 \theta=\theta_0+
 L\begin{pmatrix}v_1\\v_2\\v_3+
                  f_{20}v_1^2+f_{11}v_1v_2+f_{30}v_1^3\end{pmatrix}.
 \tag{7}
\]

It is an affine isomorphism composed with a triangular shear;
its inverse and constant Jacobian det L are explicit.
These derivative formulas are a finite algebraic recipe, not an
unsolved coefficient system.

On ordinary source coordinates (x,y,z_src,t_src) take
v1=xy,v2=x^2 z_src,v3=x^4 t_src, C=x gamma, and put

\[
 F_m=(C,\ E_1(\theta)/x,\ E_2(\theta)/x^2,\ E_3(\theta)/x^4).
 \tag{8}
\]

Here is a full polynomiality proof.
Give v1,v2,v3 weights 1,2,4.
The linear terms of E(theta0+Lv) are v_i.
For E1 only its zero constant term is needed.
For E2 the terms of weight below two are its constant and v1 terms,
both zero. The shear cannot change them.
For E3 the entire list of weights below four is
1,v1,v2,v1^2,v1*v2,v1^3.
The first three vanish by the marked point and JL=I.
The last three vanish by the displayed f20,f11,f30 formulas:
the only correction to the cubic coefficient is f20 H13.
Taylor terms of order four or more have weight at least four.
Thus every quotient in (8) is polynomial everywhere.

The full Jacobian factors are

\[
 x^7\cdot\det L\cdot\gamma^4\cdot\gamma^3/C^7
       =\det L=\gamma_0\ne0.
\]

This proves the Keller identity on the nonzero chart, then globally
because all coordinates are polynomials.

The field inverse recovers x=C/gamma,u=w/gamma,v=z/gamma^3.
Apply L^{-1} to theta-theta0, undo the triangular shear in (7),
and divide v1,v2,v3 by x,x^2,x^4.
The target inverse is
Y1=C F2,Y2=C^2 F3,Y3=C^4 F4.
Hence (8) has exactly the sweep field pair with independent C added.

For m=2, (2) is 16alpha^2+30alpha+15=0 and the data simplify to

\[
 \theta_0=(-2\alpha/3-1/2,\ 16\alpha+18,\ -24),
\]
\[
 L=\begin{pmatrix}
0&1/8-\alpha/4&17\alpha/72+13/48\\
0&-32\alpha-29&-8\alpha/3-1\\
-96\alpha-84&336\alpha+414&104\alpha+94
\end{pmatrix},
\quad(f_{20},f_{11},f_{30})=(48\alpha+42,0,0).
\]

Equations (3)--(8) are therefore a compact explicit degree-twelve map.

## 3. Finite sweep and the alternating block group

Let A=Y2 and B=Y3+Y1^2/2.
The source equation reduces to

\[
 B=d(q)+a(q)(A-b(q))+(A-b(q))^2/(2q).
\]

Define

\[
 Q(q)=(A-b)^2+2qa(A-b)+2qd-2qB. \tag{9}
\]

Its degree is n=2m+2, constant term A^2, and leading coefficient

\[
 L_n=-{\alpha^2m^2\over(2m+1)(m+1)^2}\ne0.
\]

The full eliminant is Q(w^2) of degree 2n.
It makes w integral over the target.
After eliminating gamma, the third sweep equation is monic in z
up to the constant 1/2:

\[
 z^2/2+Y_1z+d(w^2)-w^2a(w^2)^2/2-Y_3=0.
\]

Thus z, then gamma are integral as well: S is finite.
When A!=0 every root w of Q(w^2) is nonzero and uniquely determines
z and gamma through the second and first sweep equations.
Its generic degree is 2n.

For fixed generic A, differentiate the rational function of q:

\[
 {dB\over dq}=-(A-b(q))^2/(2q^2).
\]

There are m+1 finite critical points with local degree three.
Their values are distinct: at such a point (A,B)=(b(q),d(q)),
and the map q -> (b,d) is birational since its degree divides
both m+1 and 2m+1, which are coprime.
Infinity has pole partitions (n-1,1), giving the remaining
Riemann--Hurwitz contribution.

Finite inertia generates the curve group and consists of single
three-cycles. The support-hypergraph lemma of KMFI12 gives A_n.
There is no enlargement from the spectator parameter A:
the identity

\[
 qQ'-Q=-(A-b(q))^2 \tag{10}
\]

makes the discriminant a square up to a constant.
Thus over algebraically closed constants the block group is A_n.

## 4. The full signed kernel

Choose square roots w_i of the n roots q_i.
The action embeds in C2 wr S_n.
The product of q_i is

\[
 \prod q_i=-(2m+1)\left({(m+1)A\over\alpha m}\right)^2.
 \tag{11}
\]

Over geometric constants, the product of the chosen w_i is
therefore in the base field. Only even sign flips occur.
By Section 3 the block action is A_n, so G<=V_n:A_n.

Let W be its sign kernel.
At B=infinity the w-curve has pole orders 2n-2 and 2.
Raising that inertia element to n-1, which is odd, gives
j=(1,...,1), the simultaneous flip of all n pairs.
Thus 0!=<j><=W.

Any A_n-invariant subspace of V_n containing a nonconstant vector
is all of V_n. Take a three-cycle moving one support point and
two points outside the support; adding its translate gives a
weight-two vector. Its A_n orbit spans V_n.
Consequently W is either <j> or V_n.

The first alternative is impossible for n>=6.
If W=<j>, every base three-cycle has a unique lift of order three:
of the two lifts, exactly one has cube one.
It fixes every signed letter outside its three pairs.
Take the lifts of (i,n-1,n), i=1,...,n-2.
They satisfy x_i^3=1 and (x_i x_j)^2=1:
each relation lies in <j> and fixes the pairs outside four labels,
whereas j moves those pairs.
Carmichael's classical presentation therefore gives a section A_n->G.

The finite sweep's only ramification divisor is gamma=0.
Its map to the branch image has degree two:
k(X)=k(q,sigma), because gcd(deg b,deg d)=1, whereas q=w^2.
Hence finite inertia is two three-cycles, and these conjugates
normally generate G by KMFI27.
Under the assumed central kernel, each such conjugate is the
unique order-three lift and lies in the section.
They cannot generate j, a contradiction.

Therefore W=V_n, and the full geometric group is exactly
V_n:A_n in its signed 2n-letter action.
The presentation input is equation (3.1) of
[Guralnick--Kantor--Kassabov--Lubotzky, JEMS 2011](https://ems.press/content/serial-article-files/31768),
attributed in the [source manifest](sources/signed-alternating-manifest.json).
No finite-field sample or group-order guess replaces this proof.

## 5. Exact arithmetic constants

From (10), resultant identities give

\[
 \operatorname{disc}Q
 = {(-1)^{n(n-1)/2}\over L_n}\,
     {\operatorname{Res}(Q,A-b)^2\over A^2}.
\]

Its square class is (-1)^m(2m+1).
Equation (11) supplies the flip-parity square class -(2m+1).
They describe the two characters of
(C2 wr S_n)/(V_n:A_n), a group C2 times C2.
Their square roots lie in the constant field.

Since the geometric group is already the full common kernel
of those two characters, there are no further constants.
This proves (1), including cases where one class is square
or the two classes coincide in K_m.

For m=2 the classes 5 and -5 are independent over Q(sqrt(-15)):
neither they nor their product -1 is a square there.
Thus the arithmetic group is the full C2 wr S6 and the constant
extension is K_2(sqrt(5),sqrt(-5))=K_2(i,sqrt(5)).

## 6. Atomic maps with imprimitive monodromy

Fix the sheet +1 and let H be its stabilizer in G_n.
Let K be the stabilizer of the pair {+1,-1}.
Then [K:H]=2, and its block system is the n pairs.
This is the only proper intermediate subgroup H<J<G_n.

Indeed H intersect V_n is the (n-2)-dimensional subspace with
first coordinate zero.
If J is not contained in K, its projection contains A_(n-1)
and another element, hence all A_n, since the natural action is primitive.
Then J intersect V_n is A_n-invariant and contains H intersect V_n,
which is not contained in <j>.
Section 4's elementary submodule argument makes it V_n;
therefore J=G_n.
Inside K no subgroup lies strictly between H and K.

A nontrivial polynomial factorization of F_m through ordinary
affine four-space would make both factors Keller:
the product of their polynomial Jacobian determinants is a unit.
A birational Keller factor is an automorphism, so two noninvertible
factors give a proper intermediate field, as in the retained atomicity theorem.
The unique intermediate field would force the first factor
to have degree two. Its function-field extension is then Galois,
contradicting the Galois-case Keller obstruction in KMFI1.
Thus F_m is atomic. The argument persists after algebraically closed constant extension,
every identity stabilization and polynomial left-right equivalence,
since the action and its subgroup interval are unchanged and the
degree-two obstruction holds in every affine dimension.
This is an actual failure of the converse
"atomic implies primitive monodromy", including its stable version.

There is also a distinction between generic and regular source symmetry.
The normalizer N_G(H) is K: H is normal in K, while a larger
normalizer would be G and make the faithful nontrivial point
stabilizer normal. Thus the generic deck group is K/H=C2.
Its involution is
(w,sigma,gamma)->(-w,sigma+2gamma,-gamma), with C fixed.
In the actual source coordinates it sends x=C/gamma to -x.
It cannot extend to a regular source map.
Indeed gamma is polynomial and restricts to the constant gamma0
on x=0. An extension preserving x=0 would pull it back to gamma0
there, whereas the generic identity sends it to -gamma0.
Since gamma0!=0, this is impossible.
Thus the polynomial source deck group is trivial even though
the generic deck group is C2.

## 7. Exact cycle-gcd invariant

For a signed block cycle of length ell, the action has either
two cycles of length ell (no net flip) or one cycle of length
2ell (net flip).

Suppose n=2^b is a power of two.
A cycle gcd equal to the degree 2n would require a single 2n-cycle.
Its block permutation would be an n-cycle, which is odd and
does not belong to A_n. Every possible cycle gcd divides 2n,
so its lcm is at most n. The value n is attained by two block
cycles of length n/2, both with net flip.
Thus C=n in this case.

Otherwise write n=2^b s with s>=3 odd.
Partition the blocks into cycles of lengths 2^b and n-2^b,
and flip both. Their block permutation is even and their flip
parity is even. On signed letters the two cycle lengths are
2^(b+1) and 2^(b+1)(s-1), with gcd 2^(b+1).
For every odd prime power p^a dividing n, partition the blocks
into p^a-cycles with no flips. This lies in G_n and has cycle
gcd p^a. These witnesses attain every prime-power divisor of 2n,
so C=2n.

KMFI9--10 require C=degree for independently stabilized
curve-cover Keller pairs. The power-two subfamily fails this
necessary condition, beginning with an actual degree-sixteen map
with group 2^7:A8. Other even n pass this invariant; no curve
presentation follows merely from passing it.
This proves KMFI31.

## 8. Status and evidence

The [formula module](signed_alternating.py) implements every
coefficient and finite stage operation in Sections 1--2.
The [checker](verify_signed_alternating.py),
[model](signed_alternating_model.json) and
[receipt](signed_alternating_receipt.json) retain exact controls.
The universal degree, group, atomicity and curve-obstruction
arguments are written above; bounded controls are not their substitute.

This completes the earlier marked candidate and supplies an infinite
family of actual faithful transitive actions. It does not settle
PSL2(11), natural A4, all imprimitive groups, or the requested
necessary-and-sufficient classification. The full goal stays open.
