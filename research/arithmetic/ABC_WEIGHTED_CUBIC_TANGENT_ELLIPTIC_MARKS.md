# Elliptic marked fibres on a weighted cubic tangent line

Every smooth weighted cubic tangent line gives a
genus-one double cover with three marked fibres.
Its elliptic curve is a cubic twist of \(y^2=x^3+1\).
The marked points are two translates of its rational
order-three subgroup. The translation is torsion
exactly on the symmetric profile. An asymmetric
profile has two genuine prime triangles in the same
rank-one quadratic twist, so the unweighted
primitivity obstruction requires more than the
three marked fibres.

## The double cover and its three marks

Use \(s=(y,z,x)\), the positive normalized weights
of the [exact reconstruction](ABC_COHERENT_NORM_PENCILS.md#exact-reconstruction-from-three-signed-coefficients),
and the nonzero rational tangent normal \(n\) of
[TL1](ABC_CUBIC_TANGENT_RING_CLASS_LIFTS.md#fixed-profile-tangent-squares-and-the-cofactor-root).
Put
\[
 a_i=n_i^3/s_i^2,\qquad
 X_i=p_i s_i^2/n_i^2,\qquad
 k=a_2/a_1,\qquad w=k(k+1),\qquad m=2k+1 .
\]
The conditions on \(n\) and the prime line give
\(\sum_i a_i=\sum_i a_iX_i=0\).
Thus \((X_1,X_2,X_3)=(u,u-(k+1)h,u-kh)\)
for rational \(u,h\). On a genuine fibre \(h\ne0\).
The same notation below is used with \(u,h\)
as homogeneous variables. Put \(t=u/h\).

**Theorem WM1.** The weighted tangent field is
\(\mathbb Q(\sqrt{f_k(t)})\), where
\[
 f_k(t)=-3t^4+4mt^3-6wt^2+w^2 .
\]
The double cover \(U^2=f_k(t)\) is
\[
 E_w: Y^2=X^3+16w^2,\qquad
 T=(4w,-4wm),\qquad D=(0,4w)
\]
with involution \(P\mapsto T-P\).
Its three fibres at \(t=0,k,k+1\) are respectively
\(\{O,T\}\), \(\{-D,T+D\}\), and \(\{D,T-D\}\).
The point \(T\) is torsion if and only if
\(k\in\{1,-2,-1/2\}\), the profiles for which two
of \(a_1,a_2,a_3\) agree.

In the Heron form, the weighted inputs
\(s_i^2p_i^3\) equal \(a_i^2X_i^3\).
With \(A,B,C\) denoting these inputs,
expanding gives
\[
 A^2+B^2+C^2-2(AB+BC+CA)=a_1^4 h^6 w^2 f_k(t).
\]
This has the opposite sign from the Heron
numerator. All displayed prefactors are squares,
as is the denominator in the exact field
reconstruction. The condition \(h=0\) is the
tangent contact and gives Heron value zero.
The discriminant of \(f_k\) is \(-6912w^4\).
In particular the cover is smooth as a
projective curve for all \(k\ne0,-1\).

Write \(b=(U+w)/t^2\) and
\(v=(b^2+3)t-2m\).
The quartic gives
\[
 (b^2+3)t^2-4mt-2wb+6w=0,\qquad
 v^2=2w(b-1)^3+4 .
\]
Consequently \((X,Y)=(2w(b-1),2wv)\), with
inverse relation
\[
 t=\frac{2w(X-4w)}{Y-4wm},\qquad U=bt^2-w .
\]
At \(t=0\) the two points specialize to
\(O,T\), so the deck involution is \(T-P\).
At \(k\), the sign \(U=-k\) gives \(-D\).
At \(k+1\), the sign \(U=k+1\) gives \(D\).
The other points are their images under
the deck involution.

The rational torsion on \(E_w\) has order
dividing \(6\). One way to see this in the
present model is to use good primes
\(q\equiv2\bmod3\), at which the cube
map is a bijection and \(\#E_w(\mathbb F_q)=q+1\).
Dirichlet's theorem allows such primes
outside the finite set of bad primes
with \(q+1\) prime to any fixed
\(\ell\ge5\), with \(v_2(q+1)=1\), and
with \(v_3(q+1)=1\), respectively.
The reduction injections on prime-to-\(q\)
torsion give the bound. The subgroup
of order three is \(\{O,D,-D\}\).
Since \(X(T)=4w\ne0\), torsion for
\(T\) requires rational two-torsion and
thus \(4w\) to be a rational cube.
Scaling \(E_w\) to \(y^2=x^3+1\)
then sends \(X(T)\) to \(\sqrt[3]{4w}\).
By the torsion bound and its six
displayed rational torsion points,
the only nonzero possibilities are
\(-1\) and \(2\). Thus
\(4w=-1\) or \(8\), giving the
three stated profiles.

## A degree-three map merges the three fibres

**Proposition WM2.** The quotient by
\(\langle D\rangle\) gives a degree-three
rational map on the parameters. In coordinates
\[
 Q(t)=(w+1)t^2-wmt+w^2,\qquad
 r(t)=\frac{t(t-k)(t-k-1)}{2Q(t)}
\]
it sends all three marks to zero. Its other
double cover is the quotient of
\(E'_w:Y^2=X^3-432w^2\) by \(T'-P\),
where \(T'=(4(w+1),-4m(w-2))\).
Its field polynomial is
\[
 g_k(r)=1-24(w+1)r^2+32m(w-2)r^3-48(w+1)^2r^4 .
\]
One has the exact identity
\(g_k(r(t))=4r'(t)^2 f_k(t)\).

Use the \(a=0,D=1\) case of
[Cohen–Pazuki, Section 1](https://www.impan.pl/shop/publication/transaction/download/product/82564).
Writing \(B=4w\), the isogeny is
\[
 X'=\frac{X^3+4B^2}{X^2},\qquad
 Y'=Y\frac{X^3-8B^2}{X^3}.
\]
It sends \(T\) to \(T'\). The quotient coordinate
\((X'-X(T'))/(Y'+Y(T'))\) has zero fibre
\(O,T'\). Substitute the inverse relation
for \(t\) from WM1 in the isogeny. Eliminating
\(Y\) by \(Y=4wm+2w(X-4w)/t\) and reducing
in \(X\) gives the displayed \(r(t)\).
The quadratic \(Q\) has discriminant
\(-3w^2\) and none of the marked roots is
a pole. In homogeneous coordinates the
poles are infinity and the two roots of \(Q\).

Writing \(X(T')=X_0,Y(T')=Y_0\), the equation
\(Y'=(X'-X_0)/r-Y_0\) cuts out a quadratic
in \(X'\) after removing the fixed point
\(-T'\). Its discriminant is
\(1-6X_0r^2-8Y_0r^3-3X_0^2r^4\).
This gives \(g_k\).
The isogeny pulls back the regular differential
on the second double cover to a constant times
the regular differential on the first.
Thus \(g_k(r(t))=c^2r'(t)^2f_k(t)\).
At \(t=0\), \(g_k(0)=1\),
\(f_k(0)=w^2\) and \(r'(0)=1/(2w)\).
It follows that \(c^2=4\).
The polynomial identity also verifies this
at ramified and exceptional parameter values
by clearing denominators. There the field
comparison is read on the smooth covers.

## The cyclotomic orbit and a necessary reduction test

**Proposition WM3.** Over
\(K=\mathbb Q(\sqrt{-3})\) put
\(\zeta=(-1+\sqrt{-3})/2\) and
\(J=[\zeta]-[\zeta^2]\) in the endomorphism
ring of \(E_w\). For each odd integer \(n\)
define
\[
 P_n=-[\zeta]T+\left[\frac{n+1}{2}\right] JT .
\]
Its trace is \(T\), and \(2P_n-T=[n]JT\).
For a nonsymmetric profile this gives
infinitely many finite rational parameters
with field \(K\), including infinitely
many with all three homogeneous markings
positive. If \(q>3\) is a good marking
prime and \(P_n\) reduces to the \(i\)-th
marked fibre, then
\[
 [3n^2+1]T=[2]D_i\bmod q,\qquad
 (D_1,D_2,D_3)=(O,D,-D).
\]
Here the indices correspond to
\(t=0,k+1,k\), the zeros of
\(p_1,p_2,p_3\), respectively.
A good marking prime excludes the
denominators and the primes at which
the three marked fibres are not
distinct and étale.

The action is \([\zeta](X,Y)=(\zeta X,Y)\).
It satisfies \(1+[\zeta]+[\zeta^2]=0\),
\(J^2=[-3]\), and \(JD=O\).
The point \(-[\zeta]T\) is one
point at infinity in the \(t\)-cover.
Its trace is \(T\) and
\(2(-[\zeta]T)-T=-JT\).
The anti-invariant point \(JT\)
is nontorsion for a nonsymmetric
profile by WM1. The trace coset
therefore contains infinitely many
distinct points with rational \(t\).
A multiple of \(JT\) lies in the
identity circle of the compact
complex anti-invariant subgroup.
Its multiples are dense in that
circle. At the point at infinity
the map \(t\) has a simple pole,
since the quartic has nonzero
leading coefficient. Thus there
are infinitely many finite parameters
\(t>\max(0,k,k+1)\) by approaching
that point along the appropriate side
of the trace circle. These give
positive rational linear markings,
which can be scaled to coprime
integers.

At a good marking prime use the two
possible reductions \(P_n=D_i\) and
\(P_n=T-D_i\). They give
\([n]JT=\pm(2D_i-T)\). Apply \([n]J\).
Using \(JD_i=0\) and \(J^2=[-3]\)
gives \([-3n^2]T=T-2D_i\),
independent of the sign.
This is a necessary condition at
each of the three incident primes.
It concerns the rational reduction
of \(T\), so the resulting equality
does not depend on the prime above \(q\).

## Two genuine triangles on a rank-one twist

**Proposition WM4.** The normalized
weights \((x,y,z)=(10,1,9)\) and the
prime line \(p_1+9p_2-10p_3=0\) give
two genuine triangles over
\(\mathbb Q(\sqrt{-3})\), with
\[
\begin{array}{c|r|r|r|r}
 n & p_1 & p_2 & p_3 & h\\ \hline
 3&56611&6451&11467&5016\\
 5&8826476898331&1502124583051&2234559814579&732435231528
\end{array}
\]
Their anti-invariant points \(2P_n-T\)
are \(3JT\) and \(5JT\). The companion
quadratic twist has rank one. Hence the
rank-one uniqueness conclusion of
[FC2](ABC_CUBIC_TANGENT_FIELD_COLLISIONS.md#the-three-marked-fibres-force-primitivity)
for unweighted progressions does not
extend to arbitrary fixed weighted
tangent lines.

For every primitive integer triple
\((a_1,a_2,a_3)\) of nonzero sum zero,
one can take \(n_i=a_i\) and \(s_i=|a_i|\).
These \(s_i\) are pairwise coprime and
the prime line is \(\sum a_i p_i=0\).
For \(i,j,l\) distinct, reduction modulo
\(a_i\) gives \(a_j=-a_l\) and then
\(p_j=p_l\bmod a_i\). So
\(a_j^2p_j^3=a_l^2p_l^3\bmod a_i\),
giving all full trace congruences.
In particular this applies to \((1,9,-10)\).
Both rows have six certified primes in
total, all coprime to \(2\cdot3\cdot5\).
The exact Heron scale \(\delta\) in
\(E=3\delta^2\) is \(732435231528\)
for the first row and
\(5486777363970309410097697565862059160\)
for the second. In each row \(E\) is
coprime to all incident primes.
This is enough for the exact
coherent ideal reconstruction.
Evaluating \(t(P_3),t(P_5)\) for
\(k=9\) gives the two ratios \(p_1/h\).

For \(k=9\), \(E_w\) is isomorphic over
\(\mathbb Q\) to \(y^2=x^3+2025\).
The exact two-descent bounds from
[PARI/GP ellrank](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html)
on its minimal model both equal one.
The \(-3\) quadratic twist is isogenous
to \(E_w\) by WM2, so it also has
rank one. The list of primes, both
reconstructions, orbit relations,
reduction tests and rank bounds are
reproducible with
[the Sage checker](verify_weighted_cubic_tangent_elliptic_marks.py).
The rank statement uses the unconditional
upper and lower descent bounds of
that implementation.

The mechanism distinguishes the torsion
and nontorsion shifts of the involution.
In the symmetric family the marked
points are torsion, and good reduction
preserves their prime divisibility
under odd multiplication. In the
asymmetric family the shift is
nontorsion and the reduction test
depends on \(3n^2+1\).
An upper count for the indices with
three prime markings on this orbit
would require additional control of
those varying division conditions.
