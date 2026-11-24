# All polynomial fibres of a fixed odd-height endpoint coordinate

Work over a characteristic-zero field \(K\). Fix an odd integer \(a\geq1\)
and a rational function

\[
u\in K(x),\qquad u(0)=0,\quad u'(0)=1.
\tag{F1}
\]

The fibre considered here consists of all normalized polynomial pairs
\(C,D\in K[x]\), \(C(0)=-D(0)=1\), for which the regular sphere phase

\[
p=t^a\frac{C(x)+t^2D(x)}x,
\qquad t^2+2xy=1,
\tag{F2}
\]

has every positive pure moment zero and its selected endpoint coordinate
\(x(-D/C)^{a/2}/C\) equals the **fixed** function \(u\). The square root
has value 1 at zero. This is an unrestricted statement within the exact
two-height profile (F2), with no coefficient-degree bound. The coordinate
\(u\) is held fixed throughout; classifying all possible \(u\), arbitrary
sphere phases, and unrestricted primitive degrees 5 and 7 remain separate
problems.

## 1. The classification

Put

\[
b=\frac{u}{xu'},\qquad b(0)=1.
\tag{F3}
\]

Let \(j(x)=-x+O(x^2)\) be the formal branch defined by
\(u(j(x))=-u(x)\). Let \(\mathcal C_u\) be the irreducible component of
\(u(x)+u(y)=0\) through \((0,0)\), and use its smooth projective model
when discussing orders at points. This component is uniquely determined:
the equation has a nonzero \(y\)-derivative at the marked point. In
particular, choosing an unrelated component is not permitted. The
function \(j\) need not be rational.

A normalized **reciprocal unit** means a rational function \(\eta\) with

\[
\eta(0)=1,\qquad \eta(x)\eta(y)=1
\quad\hbox{on }\mathcal C_u.
\tag{F4}
\]

This terminology describes an identity in the function field of the
marked component. It does not assert an involution on \(K(x)\).

Define a monic polynomial \(J_{a,u}\) by its geometric root orders:

\[
\begin{aligned}
\operatorname{ord}_\alpha J_{a,u}
 &=\left\lceil\frac{s}{a}\right\rceil-1
 &&\text{if }\alpha\ne0,\quad \operatorname{ord}_\alpha u=s>0,\\
 &=\operatorname{ord}_\alpha u'
 &&\text{if }u(\alpha)\notin\{0,\infty\},\quad u'(\alpha)=0,\\
 &=0&&\text{otherwise.}
\end{aligned}
\tag{F5}
\]

The data are Galois invariant, so \(J_{a,u}\in K[x]\). Its roots avoid
zero and the poles of \(u\). The value \(J_{a,u}=1\) is allowed.

**Theorem FixedCoordinateFibreClassification.** The following statements
hold over every characteristic-zero \(K\).

1. The fibre is in bijection with the reciprocal units \(\eta\in K(x)\)
   satisfying these two finite local conditions:

   * all finite poles of \(\eta\) occur at poles of \(u\); at a pole of
     \(u\) of order \(m\), their order is at most
     \(1+\lfloor m/(a+2)\rfloor\);
   * \(\eta\equiv-1\pmod {J_{a,u}}\), interpreted in the local rings at
     the roots of \(J_{a,u}\).

   The pair is recovered exactly by

   \[
   X=\frac b2(1+\eta),\qquad
   C=\frac{xX^a}{u},\qquad D=-\frac{xX^{a+2}}u.
   \tag{F6}
   \]

2. Before imposing the congruence, those reciprocal units form either a
   finite set, or a disjoint union of finitely many sets

   \[
   \{\eta_i R^k:k\geq0\},
   \tag{F7}
   \]

   where \(R\in K[x]\) is a unique nonconstant polynomial with
   \(R(0)=1\). More precisely, the polynomial reciprocal units are
   either just \(1\), or exactly \(R^k\), \(k\geq0\). In the first case
   the bounded-pole set is finite. In the second case (F7) applies.

3. In each set (F7), the surviving indices are empty, a singleton, or
   one congruence class modulo a positive integer, intersected with
   \(\mathbb Z_{\geq0}\). If there are infinite sets of survivors, they
   have the **same modulus** \(d\): the multiplicative order of
   \([R]\) in \((K[x]/J_{a,u})^*\). When \(J_{a,u}=1\), the congruence
   is vacuous and \(d=1\).

4. Whenever \(R\) exists, there is a rational \(U\in K(z)\) such that

   \[
   u=U(R),\qquad U(1/z)=-U(z),\qquad
   \deg R\mid\deg(u).
   \tag{F8}
   \]

   Every rational reciprocal unit belongs to \(K(R)\). Thus an infinite
   fibre has a canonical polynomial quotient on which return is
   \(z\mapsto1/z\), although return on the original \(x\)-line can be
   nonrational.

5. At fixed \(u,a\), only finitely many pairs have least sphere degree
   at most a specified bound. If the fibre is infinite, its set of exact
   least degrees is, up to finitely many exceptions, a finite union of
   arithmetic progressions. Their common step is

   \[
   (a+2)d\deg R.
   \tag{F9}
   \]

The finite family in (F7) is an existence conclusion. The proof also gives
an exact finite-support integer-lattice description and a finite local
test. It does **not** supply a general numerical bound for the number or
degrees of the \(\eta_i\), or claim a completed algorithm that enumerates
them over an arbitrary presented field. In particular, the use of
Dickson's lemma below is not an effective bound.

## 2. Endpoint purity becomes a reciprocal-unit equation

Write \(X=\sqrt{-D/C}\), with \(X(0)=1\). Rationality of this square root
in the present fixed-coordinate setting is elementary: since \(a\) is
odd and \(u=xX^a/C\) is rational,

\[
X=\frac{uC/x}{(X^2)^{(a-1)/2}}\in K(x).
\tag{F10}
\]

Let \(v\) be the formal inverse of \(u\). The retained all-height endpoint
identity gives

\[
\Lambda(x^q p^m)=\kappa_{a,m}[s^m]v(s)^q\Phi(s),\quad
\Phi(s)=T(v(s)),\quad T=\frac{X^{a-1}}{Cu'}=\frac bX
\tag{F11}
\]

for even \(m\), where
\(\kappa_{a,m}=\int_0^1 t^{am}(1-t^2)^m\,dt\ne0\).
Odd powers have zero moment by parity. Consequently purity is exactly
\(T(x)+T(j(x))=2\).

Set \(\eta=2/T-1=2X/b-1\). Since \(T(0)=1\), \(\eta(0)=1\). Directly,

\[
\frac2{1+\eta(x)}+\frac2{1+\eta(y)}=2
\quad\Longleftrightarrow\quad \eta(x)\eta(y)=1.
\tag{F12}
\]

All denominators here are units at \((0,0)\). An identity on the formal
branch holds on its irreducible component, and conversely. This proves
the purity part of (F6), in both directions. The choices \(\eta=-1\) or
a norm constant \(-1\) are excluded by the marked normalization.

The only retained moment input is the exact endpoint identity, not a
finite-prefix inference. The older all-order rationality theorem is
needed only for the optional observation that **every** pure odd-height
unit-center phase has some rational \(u\), and hence belongs to one of
the fibres classified here.

## 3. Exact polynomiality and the finite support

At a finite point \(\alpha\), let \(s=\operatorname{ord}_\alpha u\) and
\(\delta=1\) if \(\alpha=0\), otherwise \(0\). Formula (F6) makes both
\(C,D\) regular exactly when

\[
\operatorname{ord}_\alpha X\geq
L_\alpha:=\max\left\{
\left\lceil\frac{s-\delta}{a}\right\rceil,
\left\lceil\frac{s-\delta}{a+2}\right\rceil\right\}.
\tag{F13}
\]

At zero this is automatic, and (F6) gives \(X(0)=C(0)=-D(0)=1\).
At every other finite zero or pole of \(u\), \(b\) has a simple zero.
At a point with finite nonzero \(u\), its pole order is precisely the
zero order of \(u'\). These facts follow by differentiating a local
leading term; characteristic zero is essential here.

At a pole of \(u\) of order \(m\), (F13) becomes
\(\operatorname{ord}X\geq-\lfloor m/(a+2)\rfloor\). Since
\(\operatorname{ord}b=1\), it is equivalent to

\[
\operatorname{ord}\eta\geq-1-\lfloor m/(a+2)\rfloor.
\tag{F14}
\]

There is no possible cancellation issue in this equivalence: the lower
bound in (F14) is negative, and a negative order of \(1+\eta\) is exactly
the order of \(\eta\). At a zero of \(u\) of order \(s>0\), (F13) is
\(\operatorname{ord}(1+\eta)\geq\lceil s/a\rceil-1\). At a critical
point with finite nonzero \(u\), it is
\(\operatorname{ord}(1+\eta)\geq\operatorname{ord}u'\). At all other
finite points it only says that \(\eta\) is regular. This proves exactly
the local conditions in statement 1. A rational function regular at every
finite point is a polynomial; nothing is required at infinity.

The possible zeros of \(\eta\) also lie in a fixed finite set. Put
\(c=u(\infty)\in\mathbb P^1(K)\), with \(-\infty=\infty\), and let

\[
S=u^{-1}(\infty)\ \cup\ u^{-1}(-c)\ \cup\ \{\infty\}.
\tag{F15}
\]

If \(\eta\) has a zero at \(\alpha\), choose any point of the smooth
model of \(\mathcal C_u\) above \(x=\alpha\). The identity (F4) forces
\(y\) to be a pole of \(\eta\). Such a pole is either a pole of \(u\)
or infinity. Therefore \(u(\alpha)=\infty\) or \(-c\). Both projections
of the marked curve are nonconstant finite maps, so the chosen point
exists. Thus every zero and pole of \(\eta\) lies in (F15). If \(c=0\),
this set can contain zero, but (F4) prescribes order zero there.

Over an algebraic closure write the finite nonzero points of \(S\) as
\(\alpha_1,\ldots,\alpha_s\). Every possible normalized unit is uniquely

\[
\eta=\prod_{i=1}^s(1-x/\alpha_i)^{n_i},\qquad n_i\in\mathbb Z.
\tag{F16}
\]

The reciprocal condition is an integer lattice condition. On the marked
curve form the divisors of
\((1-x/\alpha_i)(1-y/\alpha_i)\). They have finite support; their orders
give a finite integer matrix \(M\). Then (F4) is precisely \(Mn=0\).
Indeed a product with zero divisor on the smooth projective curve is
constant, and its value at the marked point is 1. The pole restrictions
are just coordinate lower bounds \(n_i\geq\ell_i\): use (F14) at poles
of \(u\), and \(\ell_i=0\) elsewhere.

This describes the pre-congruence set as

\[
\{n\in\ker_{\mathbb Z}M:n_i\geq\ell_i\}.
\tag{F17}
\]

It is not permissible to omit multiplicities in the pullback divisors,
or replace the marked component by the full reducible correspondence.
For \(K\)-rational units one may either impose Galois invariance on
(F16), or use normalized irreducible factors in \(K[x]\) from the start.
Both give the same subgroup and coordinate bounds.

## 4. There is at most one polynomial generator

**Lemma PolynomialReciprocalUnits.** The normalized polynomial solutions
of (F4) are either \(1\), or the nonnegative powers of one uniquely
determined polynomial \(R\). If \(R\) exists, then (F8) holds and every
rational solution of (F4) is in \(K(R)\).

First suppose \(P\) is a nonconstant polynomial solution. Embed
\(E=K(u,P)\subset K(x)\) in \(K((s))\) using \(x=v(s)\). Replacement
\(s\mapsto-s\) preserves \(E\): it sends \(u\mapsto-u\) and
\(P\mapsto P^{-1}\). Thus it induces a genuine involution \(\sigma\) of
this intermediate field. No extension to \(K(x)\) is being asserted.

By Lüroth, \(E=K(w)\). The generator can be chosen polynomial. To see
this, the polynomial \(P\) has only one pole on the \(x\)-line. Its
descended rational function on the \(w\)-line therefore has a unique
pole \(Q\), and the entire fibre above \(Q\) is the point \(x=\infty\).
The place \(Q\) is \(K\)-rational, since it is the image of that rational
point. Choose a coordinate \(w\) with its only pole at \(Q\), and
translate so that \(w(0)=0\). Its pullback has no finite poles, and is
therefore a polynomial. The descended \(P=\bar P(w)\) is polynomial too.

Since \(\sigma(P)=P^{-1}\), \(\bar P\) has just one zero, at
\(\alpha=\sigma(Q)\), and that zero is finite. It follows that
\(\bar P=c(w-\alpha)^e\). The involution is Möbius and exchanges
infinity with \(\alpha\). It fixes the marked point \(w=0\), because
formal replacement fixes the place \(s=0\). Thus

\[
S=\frac{w-\alpha}{-\alpha}\in K[x],\qquad
S(0)=1,\qquad \sigma(S)=S^{-1},\qquad P=S^e.
\tag{F18}
\]

Here \(\alpha\ne0\), since \(P(0)=1\). These equalities also hold on
the original marked component by its formal embedding.

Choose a polynomial solution \(R\) of least positive degree. For any
other polynomial solution \(P\), repeat this argument with
\(E=K(u,R,P)\). It gives a polynomial \(S\) with \(R=S^e\). Minimality
forces \(e=1\), so \(E=K(R)\). The descended polynomial \(P\) has its
only zero at \(R=0\), hence \(P=R^k\); normalization removes its scalar.
This proves existence, cyclicity and uniqueness.

The same argument applies to \(E=K(u,R,\eta)\) for **any rational**
reciprocal unit \(\eta\): the field is stable under
\((u,R,\eta)\mapsto(-u,R^{-1},\eta^{-1})\). Polynomiality of \(R\)
again gives (F18), and its minimal degree again forces \(E=K(R)\).
Therefore \(u=U(R)\), \(\eta\in K(R)\), and the involution is inversion
in this coordinate. Rational-map degrees multiply, giving (F8).

The lemma works over the original field \(K\). It also shows that a
geometric generator over \(\overline K\), if one exists, is defined over
\(K\): uniqueness and normalization make it Galois invariant. No choice
of roots or a split original return is added to the hypotheses.

## 5. Finite translates and the exact integer filters

Apply the elementary form of Dickson's lemma: any infinite subset of
\(\mathbb N^s\) contains two distinct coordinatewise comparable vectors.
For completeness, an infinite sequence of nonnegative integers has an
infinite nondecreasing subsequence, by taking an infinitely repeated
value or successively larger values. Applying this to each of finitely
many coordinates proves the assertion by induction.

If (F17) were infinite and there were no nonconstant polynomial unit,
two comparable vectors would have a nonzero nonnegative difference in
\(\ker M\). Their quotient would be a nonconstant polynomial reciprocal
unit. This is a contradiction. Hence that case is finite.

Now let \(r\) be the finite valuation vector of the generator \(R\).
It is nonzero and coordinatewise nonnegative. Call a vector \(n\) in
(F17) minimal if \(n-r\) violates one of its lower bounds. Every vector
can be reduced to a minimal one by repeatedly subtracting \(r\): some
positive coordinate of \(r\) ensures termination. There are only
finitely many minimal vectors. Otherwise two would be comparable;
their nonnegative difference would be the vector of a polynomial unit,
hence \(kr\), \(k\geq1\), contradicting minimality of the larger vector.
This proves (F7). Its minimal representatives have disjoint forward
orbits, by the same subtraction argument.

It remains to filter one orbit. All finite poles of \(\eta_i\) are
poles of \(u\), so its denominator is invertible modulo \(J=J_{a,u}\).
If \(\eta_i\) is not a unit modulo \(J\), no index works. If \(R\) is
not a unit modulo \(J\), no positive index works, leaving only \(k=0\)
to test. In every other case the complete condition is

\[
[R]^k=-[\eta_i]^{-1}\quad\hbox{in }(K[x]/J)^*.
\tag{F19}
\]

The map \(\mathbb Z\to(K[x]/J)^*\), \(k\mapsto[R]^k\), is a group
homomorphism. A nonempty fibre is a coset of its kernel, and a subgroup
of \(\mathbb Z\) is either \(0\) or \(d\mathbb Z\), \(d\geq1\).
This proves statement 3, including its common-modulus claim. No
Skolem--Mahler--Lech theorem or zero-set theorem for general recurrence
sequences is needed: (F19) is one equation in a cyclic subgroup.

The order test has a useful exact geometric form. In characteristic
zero, \([R]\) has finite order modulo \(J\) exactly when, at every
geometric root \(\alpha\) of \(J\) of order \(q\),

\[
R(\alpha)\text{ is a root of unity},\qquad
\operatorname{ord}_\alpha(R-R(\alpha))\geq q.
\tag{F20}
\]

In that case its order is the least common multiple of the orders of
the values \(R(\alpha)\). Indeed the first nonzero term of
\(R^d-R(\alpha)^d\) has the same order as \(R-R(\alpha)\), because its
linear coefficient is \(dR(\alpha)^{d-1}\ne0\). Thus checking only the
root values, and ignoring jets at multiple roots of \(J\), would be
incorrect.

## 6. Exact degree growth and extension fields

For any polynomial pair in (F2), its least sphere degree is

\[
\deg_S p=\max\{\deg C+a-1,\deg D+a+1\}.
\tag{F21}
\]

One way to verify exactness is to reduce \(t^2=1-2xy\). The unique
normal form with \(t\)-degree at most one is
\(t(1-2xy)^{(a-1)/2}((C+D)/x-2yD)\). Its highest powers of \(y\)
distinguish the top \(D\)-term from the \(C\)-term. If their total
degrees coincide, those powers still prevent cancellation. Reduction
does not increase total degree, so this is the least degree, not only
an upper bound.

In one infinite orbit write \(g=\deg R\),
\(s_\infty=\operatorname{ord}_\infty u\), and

\[
\delta_i=-\operatorname{ord}_\infty b
          -\operatorname{ord}_\infty\eta_i.
\]

For all sufficiently large \(k\), \(\eta_iR^k\) has a pole at infinity,
so its leading term cannot cancel with 1. Formula (F6) gives
\(-\operatorname{ord}_\infty X=\delta_i+kg\). Hence

\[
\begin{aligned}
\deg C&=1+a(\delta_i+kg)+s_\infty,\\
\deg D&=1+(a+2)(\delta_i+kg)+s_\infty,\\
\deg_Sp&=(a+2)(\delta_i+kg)+s_\infty+a+2.
\end{aligned}
\tag{F22}
\]

Restricting \(k\) to a surviving congruence class gives precisely the
step in (F9). All leading-cancellation exceptions have bounded \(k\).
Together with the finite case, this proves statement 5.

There is also no continuous coefficient freedom hidden in a fixed
coordinate. Over any extension of \(K\), (F15)--(F16) place all zeros
and poles of \(\eta\) at points algebraic over \(K\), and normalization
fixes its remaining scalar. Thus every coefficient of every pair in
the geometric fibre is algebraic over \(K\). The finite-degree fibre
is finite even over an algebraic closure. Its \(K\)-pairs are exactly
the Galois-fixed pairs. This does not claim that coefficients of
different endpoint coordinates are discrete.

## 7. Complete fibres that illustrate the distinctions

### 7.1 Every degree-one endpoint coordinate

**Corollary DegreeOneEndpointFibres.** If \(u=x\), its fibre is just
\(C=1,D=-1\). If

\[
u=\frac{x}{1+\lambda x},\qquad\lambda\in K^*,
\]

its entire fibre, at every odd \(a\geq1\), is

\[
\begin{aligned}
R&=1+2\lambda x,& B&=1+\lambda x,& L_k&=\frac{1+R^k}{2},\\
X_k&=BL_k,& C_k&=B^{a+1}L_k^a,&
D_k&=-B^{a+3}L_k^{a+2},\qquad k\geq0.
\end{aligned}
\tag{F23}
\]

The exact least degrees are

\[
\deg_Sp_k=(a+2)(k+2).
\tag{F24}
\]

Proof. For \(u=x\), finite poles of a candidate \(\eta\) are forbidden,
and \(\eta(x)\eta(-x)=1\) makes a polynomial \(\eta\) constant.
For \(\lambda\ne0\), return is
\(j=-x/(1+2\lambda x)\), or \(R\mapsto R^{-1}\). The support (F15)
in the \(R\)-coordinate is \(\{0,-1,\infty\}\). At the fixed point
\(R=-1\), (F4) forces order zero. Thus \(\eta=R^k\). A negative \(k\)
has a finite pole at \(R=0\), which is not a pole of \(u\), and is
forbidden. There are no jet conditions: \(J=1\). This proves (F23),
and \(\deg X=k+1\), (F21) give (F24), including \(k=0\).

These phases are unsafe by the retained all-marker dichotomy, except
the separate \(u=x\) safe phase. The corollary does not make a new
claim about minimal unsafe degree among arbitrary phases.

In fact **no nonlinear polynomial \(u\)** has a nonempty fibre at any
odd height. For polynomial \(u\), (F15) is just \(\{\infty\}\), so the
only candidate is \(\eta=1\). If \(\deg u=n>1\), its simple zero at zero
implies that it has at least two distinct geometric zeros. If it has
\(h\) distinct zeros, the zeros of \(u'\) at those points have total
order \(n-h\), whereas \(\deg u'=n-1\). Thus \(u'\) has a zero with
nonzero \(u\)-value. The resulting nonconstant factor of \(J\) makes
\(1\equiv-1\pmod J\) impossible. This proof does not bound the
polynomial degree or the number of roots.

For a finite nonempty pre-congruence set that is eliminated by the jet
condition, take \(u=x/(1-x^2)\). Its marked return is \(y=-x\), even
though the full correspondence has the other component \(xy=1\).
The only possible finite zeros and poles of \(\eta\) are \(\pm1\).
Reciprocity makes their orders opposite, and (F14) bounds them by 1
for every \(a\). Hence the complete candidate list is
\(\eta=((1+x)/(1-x))^k\), \(k=-1,0,1\). Here \(J=1+x^2\), and none
of these three units is \(-1\) modulo \(J\). The whole fibre is empty.

### 7.2 Recovery of the full hyperbolic fibre

Take \(R=1+x\) and \(u=(R-R^{-1})/2\). Then

\[
b=\frac{R(R+1)}{R^2+1},\qquad J=R^2+1.
\]

The support is \(\{R=0,\infty\}\), so every candidate is \(\eta=R^k\),
with \(k\geq-1\). Modulo \(J\), \(R\) has order 4. The exact filter
\(R^k\equiv-1\pmod J\) gives

\[
k=2+4r,\qquad r\geq0.
\tag{F25}
\]

This recovers the previously classified \(N=1\) hyperbolic coordinate
fibre. It is a specialization of the general argument, not a new claim
for that already retained family.

### 7.3 A quotient of degree two with nonrational original return

Put

\[
R=1+x+x^2/8,
\qquad u=\frac{2(R-1)}{R+1}
       =\frac{2x(x+8)}{(x+4)^2}.
\tag{F26}
\]

Here \(b=(x+4)(x+8)/32\) and \(J=1\). The entire fibre is

\[
\eta=R^k,\quad X=\frac{(x+4)(x+8)(1+R^k)}{64},
\quad C=\frac{(x+4)^2X^a}{2(x+8)},\quad D=-CX^2,
\qquad k\geq0.
\tag{F27}
\]

Every displayed division is polynomial. The exact degree is
\((a+2)(2k+3)\).

To check exhaustiveness without assuming a rational return, \(R\) is
a polynomial reciprocal unit and the generator of the lemma has degree
dividing \(\deg u=2\). A linear generator would give a Möbius return
on \(x\). Such a return must fix the marked zero 0, the other zero
\(-8\), and the unique pole \(-4\); hence it would be the identity,
which cannot send \(u\) to \(-u\). Thus \(R\) itself is the minimal
generator, up to the already fixed normalization. Every unit descends
to \(K(R)\). The support in that line is again \(\{0,-1,\infty\}\),
so the same order argument as in (F23) gives exactly \(R^k\), \(k\geq0\).
This example checks a case in which a rational original return really
does not exist.

### 7.4 The finite jets can obstruct lowering the starting height

Use the same quadratic \(R\), but now take

\[
u=(R-R^{-1})/2.
\tag{F28}
\]

For every \(k=2+4r\), \(r\geq0\), define \(X=b(1+R^k)/2\) and use
(F6). At the two roots of \(R\), the only poles of \(u\), the local
conditions hold. At the roots of \(R^2+1\), the congruence in (F25)
cancels the simple critical poles of \(b\). At \(x=-4\), however,
\(u\) has a zero of order two and \(R=-1\). For \(a=1\), (F5) requires
\(\eta(-4)=-1\), whereas \(R(-4)^k=1\). Thus these candidates fail
polynomiality at height one. For every odd \(a\geq3\), that extra
condition is absent and the displayed candidates are polynomial pure
phases. Their exact degrees are

\[
2(a+2)k+a.
\tag{F29}
\]

Only the displayed candidates, not the entire fibre of (F28), are being
classified in this paragraph. In particular its height-one rejection
is not a claim that the complete height-one fibre is empty. The point
is that the finite local orders in (F5) retain a genuine higher-height
phenomenon that a blanket polynomial reduction to height one would lose.

## 8. Sources, proof scope and exact controls

The indispensable endpoint source is
[ALL_START_HEIGHT_MARKER_SAFETY.md](ALL_START_HEIGHT_MARKER_SAFETY.md),
Section 1, equations (S2)--(S3). Its all-marker safety classification is
used only in the last paragraph of Section 7.1. The optional coverage
statement at the end of Section 2 uses OddStartRationalCenter in
[ODD_START_HEIGHT_RATIONALITY_AND_QUINTIC_SAFETY.md](ODD_START_HEIGHT_RATIONALITY_AND_QUINTIC_SAFETY.md).

The field import is J.S. Milne, *Fields and Galois Theory*, version 5.10,
Theorem 9.19 (printed page 117), Lemma 9.24 (page 118), and Example 9.25
(page 119), in the retained
[primary PDF](sources/milne_fields_galois_v510.pdf). The polynomial
generator and reciprocal-unit deductions from Lüroth are proved here.

For the smooth projective model and finite projections, the precise
standard facts used are the Stacks Project,
[Section 53.2, tag 0BXX](https://stacks.math.columbia.edu/tag/0BXX),
Theorem 53.2.6, Lemmas 53.2.2, 53.2.4, 53.2.8, and 53.2.9(3).
The actual statements were inspected and the
[primary page is retained](sources/stacks_curves_function_fields_0bxx.html).
They apply because the marked component has a smooth rational point
and characteristic zero; it is geometrically integral. In particular,
a zero divisor means a constant rational function on this projective
curve, with that constant fixed by the marked value.

Two primary comparisons guided the argument but are **not theorem
imports**. F. Pakovich, *On rational functions orthogonal to all powers
of a given rational function on a curve*, arXiv:0910.2105v1,
[Lemma 4.1 and Theorem 4.2, pages 18--19](https://arxiv.org/pdf/0910.2105),
connect identities among inverse branches with common compositional
factors and impose explicit monodromy hypotheses for their moment
criterion. Those hypotheses are not substituted for the marked
two-branch condition here. Barroero, Capuano and Zannier,
[*Betti maps, Pell equations in polynomials and almost-Belyi maps*](https://doi.org/10.1017/fms.2022.77),
Section 2 and Theorem 2.3, describe the ramification of an actual
polynomial Pell equation and the associated prescribed-branch covering
correspondence. No polynomial Pell equation or Belyi property is assumed
for \(u\) in this note. The finite support, cyclic unit and jet arguments
above give the required conclusions directly. No novelty claim is made.

The [checker](verify_odd_start_fixed_coordinate_fibres.py) uses exact
rational arithmetic for the displayed fibres, endpoint reconstruction,
polynomiality and degrees, selected marked identities, and raw sphere
moments. It includes invalid exponents, a multiple-root jet with a
root-of-unity value but infinite order, a wrong correspondence, and a
corrupted coefficient. Those bounded controls do not prove the
universal classification, the finite-translate assertion, or the
infinite progression statement; the written proofs do.

The [source ledger](odd_start_fixed_coordinate_fibres_source_ledger.json)
and [receipt](odd_start_fixed_coordinate_fibres_verification.json) pin the
proof, checker and exact imported sources. Replay from the repository root:

```bash
python3 research/foundations/extended-geometry/verify_odd_start_fixed_coordinate_fibres.py \
  --repo-root . --output /tmp/odd_start_fixed_coordinate_fibres_review.json
```
