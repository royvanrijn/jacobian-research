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

## Prime markings force prime CM index norm

**Theorem WM5.** Fix a nonsymmetric weighted
profile and the orbit \(P_n\) of WM3 for odd
\(n\). For all sufficiently large \(|n|\),
a single prime marking in the primitive
integer direction of \(t(P_n)\) requires
\[
 \operatorname{N}\!\left(\frac{1+n\sqrt{-3}}2\right)
       =\frac{3n^2+1}{4}
\]
to be a rational prime. Up to a finite
set of rational primes fixed by the
profile, the three markings are the
norms of the three CM elliptic
denominator ideals for the respective
base points \(T,T+D,T-D\) and the same
index \((1+n\sqrt{-3})/2\).

Let \(\lambda=(1+n\sqrt{-3})/2\) in
\(\mathcal O_K=\mathbb Z[\zeta]\).
The explicit orbit formula is
\(P_n=[\lambda]T\).
Since \(\lambda+\bar\lambda=1\),
the two index ideals are coprime.
On the rational group generated by \(D\),
\(\zeta\) acts trivially, so
\(\lambda\) and \(\bar\lambda\) both
act as \(-1\bmod3\). In particular
\[
 P_n-D_i=[\lambda](T+D_i),\qquad
 \overline P_n-D_i=[\bar\lambda](T+D_i).
\]
Choose a fixed set \(S\) containing
the bad primes, the coefficient
denominators, the primes of collision
of the marked fibres, and the primes
where a base point \(T+D_i\) reduces
to the origin. All three base points
are nontorsion. Fix also \(S\) to
contain primes dividing a nonzero
maximal minor of the map from the
two homogeneous parameters to the
three integer markings, after
clearing denominators. Normalizing
the markings to a primitive
direction then divides by a unit
outside \(S\).

For \(q\notin S\), a vanishing marking
places \(t\) in the corresponding
finite étale fibre. The two points
in this fibre are rational and
distinct, so the field \(K\) splits
at any such \(q\). The valuation of
the marking is the valuation of a local
parameter of \(P_n\) at \(D_i\) or
\(T-D_i\). At the prime of \(K\)
where it is \(D_i\), this is the
elliptic denominator valuation of
\([\lambda](T+D_i)\). At the
conjugate prime the reduction of
\(P_n\) is \(T-D_i\), and
\([\lambda](T+D_i)\) does not reduce
to the origin there. More generally,
if both primes were denominators
then \(\lambda\) and \(\bar\lambda\)
would both annihilate \(T+D_i\)
at one of the primes, which is
impossible outside \(S\).
Thus the absolute value of the
marking with its \(S\)-part removed
equals the ideal norm of the
denominator with its \(S\)-part
removed. This uses denominator
exponents \(v_{\mathfrak q}(z)\)
for local parameter \(z=-X/Y\)
at the origin, equivalently half
the negative valuation of \(X\).

Apply [Streng's main theorem](sources/streng-2008-cm-elliptic-divisibility.pdf)
to \(E_w/K\) and one of the
three fixed base points \(T+D_i\).
Its endomorphism ring over \(K\)
is the full \(\mathbb Z[\zeta]\),
so all nonzero ideals are invertible
and principal. For all ideals of
sufficiently large norm there is a
primitive denominator prime. Each
fixed prime is primitive at at
most one index by strong divisibility
(Streng, Lemma 3.3), so the primitive
prime can be required outside \(S\)
for sufficiently large index norms.

If \((\lambda)\) is composite of large
norm, it has a proper divisor ideal
\(\mathfrak a\) of norm at least
\(\sqrt{\operatorname{N}(\lambda)}\),
by dropping a smallest prime ideal
factor. The two indices
\(\mathfrak a\) and \((\lambda)\)
then have distinct primitive
denominator primes. Both divide
the denominator at \((\lambda)\).
They cannot be the two primes
over one rational prime outside
\(S\), by the conjugation argument
and \(\lambda+\bar\lambda=1\).
Hence the marking has at least
two distinct prime divisors.
The asserted prime marking
forces \((\lambda)\) to be prime.
It cannot lie over an inert or
ramified prime of \(K\), since then
it would equal its conjugate.
Thus its ideal norm is a rational
prime.

**Corollary WM6.** For each fixed
nonsymmetric tangent profile, the
number of odd indices \(|n|\le N\)
on the orbit of WM3 with all three
markings prime is \(O_k(N/\log N)\)
as \(N\) tends to infinity. The
number of those triangles with all
three vertices at most \(Z\) is
\(O_k(\sqrt{\log Z}/\log\log Z)\).

With \(n=2j+1\), the required rational
prime is \(3j^2+3j+1\). Its number
of roots modulo \(q\) is two for
\(q\equiv1\bmod3\) and zero at all
other primes. The local root weights
therefore have sieve dimension one
by the prime number theorem in
the fixed residue classes modulo
three. Apply the upper bound sieve
in [Richert, Theorem 11.1](sources/richert-sieve-methods-1976.pdf)
at level \(z=N^{1/8}\). For
squarefree \(d\) the remainder
of the divisibility count is
\(O(2^{\omega(d)})\).
The remainder sum is
\(O_\epsilon(z^{2+\epsilon})\)
after changing \(\epsilon\),
and the local product is of
order \(1/\log z\). The prime
values at most \(z\) contribute
\(O(\sqrt z)\). This proves the
index bound.

For the height bound, the function
\(X(2P-T)\) on \(E_w\) is invariant
under \(T-P\). It is therefore a
degree-four rational function of
\(t\), since \([2]\) has degree four
and \(X\) and \(t\) have degree two.
Its height is \(4h(t)+O_k(1)\).
The canonical height of \(2P_n-T\)
equals \(3n^2\widehat h(T)\) because
\(J\) has degree three. Comparing to
half the \(X\)-height gives
\(h(t(P_n))=\tfrac32 n^2\widehat h(T)+O_k(1)\).
The positive height of \(T\) is
fixed. Each parameter from primes
at most \(Z\) has height at most
\(\log Z+O_k(1)\), by its fixed
rational linear reconstruction.
Thus \(|n|\ll_k\sqrt{\log Z}\),
which yields the last bound.
The argument is confined to the
specified orbit and does not
estimate a sum over varying
tangent profiles.

## The entire cyclotomic trace coset in rank one

**Theorem WM7.** Fix a nonsymmetric weighted
tangent profile for which \(E_w(\mathbb Q)\)
has rank one. In the entire field
\(K=\mathbb Q(\sqrt{-3})\), the number
of genuine prime triangles on the line with
all vertices at most \(Z\) is
\[
 O_k\!\left(\frac{\sqrt{\log Z}}{\log\log Z}\right).
\]
After choosing an Eisenstein generator
of \(E_w(K)\) modulo torsion, all trace
parameters have indices on one integral
affine line in \(\mathbb Z[\zeta]\).
For all sufficiently large indices,
prime triples require prime index norm.
The three marks together allow a fixed
choice whose associated denominator
has no torsion shift.

First \(E_w(K)_{\rm tors}=\langle D\rangle\).
At good primes \(q\equiv2\bmod3\),
the residue field over \(K\) is
\(\mathbb F_{q^2}\) and the reduction
order is \((q+1)^2\). Varying \(q\)
as in WM1 excludes all primary torsion
except at \(2,3\).
A root over \(K\) of \(X^3+(4w)^2\)
would imply a rational root, since
a cubic with a degree-two root also
has a rational root. This requires \(4w\) to be a
rational cube. Since \(m^2=4w+1\),
that would yield a rational point
on \(y^2=x^3+1\) with \(x^3=4w\).
This last curve has rank zero.
Indeed its translated two-isogeny
models have coefficients
\((a,b)=(-3,3)\) and \((6,-3)\).
The first squareclass image is
\(\{1,3\}\), by positivity after
removing \(x\). The second contains
\(-3\), and class \(-1\) would give
\[
 n^2=-u^4+6u^2v^2+3v^4,\qquad (u,v)=1.
\]
Modulo \(3\) this forces \(3\mid u\),
after which the right side has
valuation one. Thus the second image
is \(\{1,-3\}\) and the same
two-isogeny formula as in
[FC1](ABC_CUBIC_TANGENT_FIELD_COLLISIONS.md#the-fixed-elliptic-double-cover)
proves rank zero. Its six torsion
points give only \(4w=-1,0,8\).
These are the symmetric or excluded
profiles, so there is no two-torsion.

The polynomial for three-torsion on
\(E_w\) is \(3X(X^3+4(4w)^2)\).
An additional root over \(K\) would
force \(8w\) to be a rational cube,
and hence \(w=k(k+1)\) to be a
rational cube. Write \(k=a/b\)
in lowest terms. Then \(ab(a+b)\)
would be an integer cube, and its
three pairwise coprime nonzero
factors would each be cubes.
[Fermat's exponent-three theorem](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#witness-selection)
excludes this. Thus the
three-primary group over \(K\)
is cyclic and its subgroup of
order three is rational.
Conjugation acts trivially on
the entire cyclic group, since
an involution acting trivially
modulo three cannot act as minus
one. Higher three-power torsion
would then be rational, which
WM1 excludes. This proves the
torsion claim.

The quadratic twist by \(-3\)
is isogenous to \(E_w\) by WM2.
The free subgroup over \(K\)
therefore has rank two over
\(\mathbb Z\), and is a free
rank-one module over the
principal ideal domain
\(\mathbb Z[\zeta]\).
Indeed its endomorphism action
has no nonzero torsion, because
a nonzero endomorphism has finite
kernel. Choose a generator \(G\)
modulo torsion. On the free module
conjugation sends \(G\) to a
unit multiple of \(G\). By
changing the generator by a unit,
the multiplier can be made
\(\epsilon=1\) or \(-1\).
The remaining torsion difference
in \(\sigma G-\epsilon G\)
can be taken zero as well.
For \(\epsilon=1\) conjugating
twice already kills that difference,
since it lies in a rational
group of order three. For
\(\epsilon=-1\) a translate by
an appropriate multiple of \(D\)
kills it.

Now write \(T=[\mu]G+jD\),
where \(\mu\ne0\). The rationality
of \(T\) gives \(\mu=\epsilon\bar\mu\).
Every \(P=[\lambda]G+bD\)
of trace \(T\) satisfies
\[
 \lambda+\epsilon\bar\lambda=\mu,\qquad
 2b=j\bmod3 .
\]
The torsion component \(bD\) is
fixed across the coset. It is one
of \(D_1,D_2,D_3\) in WM3, so in
that marking the denominator is
the unshifted denominator of
\([\lambda]G\). Repeat the
étale valuation comparison of
WM5 at this one mark.
The conjugate point on the
unshifted sequence is indexed
by \(\epsilon\bar\lambda\).
A denominator at both primes
over \(q\) would force the
reduction of the fixed
nontorsion point \([\mu]G\)
to be the origin. The primes
where this happens are added
to the finite excluded set.
Streng's theorem at
\((\lambda)\) and a large proper
divisor therefore forces
\((\lambda)\) to be prime
when the corresponding
marking is prime and the
index is large. An inert or
ramified prime ideal would
equal its conjugate and
divide \((\mu)\), so all
large indices under
consideration have prime
rational norm.

If \(\epsilon=1\), write
\(\mu=c\) and
\(\lambda=(c+ns)/2\), where
\(s=\sqrt{-3}\), \(c\ne0\)
and \(n\equiv c\bmod2\).
If \(\epsilon=-1\), write
\(\mu=cs\) and
\(\lambda=(n+cs)/2\)
with the same parity constraint.
The norm is respectively
\((c^2+3n^2)/4\) or
\((n^2+3c^2)/4\).
On each progression of \(n\)
with the required parity
this is an integral quadratic
with discriminant \(-3c^2\).
Outside a finite set of primes
the weights for the sieve are
the same as in WM6. A fixed
prime divisor would itself
give only finitely many prime
values. Finally \(2P-T\)
has free index \(2\lambda-\mu\),
so its canonical height grows
as a positive constant times
\(n^2\). The same degree-four
height comparison as in WM6
bounds \(|n|\) by \(O_k(\sqrt{\log Z})\).
This proves the claimed
count on the whole trace coset,
including the profile of WM4.
