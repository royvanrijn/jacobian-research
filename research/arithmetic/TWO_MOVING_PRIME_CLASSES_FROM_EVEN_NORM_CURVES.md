# Two moving prime classes from an even norm curve

A middle term in the even norm curve makes independent
prime Frobenius prescriptions possible at both ends of a
genuine ratio. Two fixed auxiliary primes ramify to remove
their valuation obstructions. The resulting mirror layer
has an everywhere unramified \(C_m^2\) quotient with
complete splitting at three. In the imaginary base fields,
the two actual moving prime ideals have independent
\(3\)-primary parts of order \(m\), and their even parts
satisfy an exact relation with the fixed ramified primes.

## Coefficients and local neighborhoods

Fix \(m=3^n\), \(n\ge1\). Choose two auxiliary primes
\(p,q\), and put
\[
 A=\frac{p-q}{4},\qquad B=-\frac{p+q}{2},\qquad
 z=u^m,\qquad R^2=F(u)=A^2u^{2m}+B u^m+1 .
\]
For finite nonzero \(u\) in a quadratic specialization,
define
\[
 f_\pm=A z\pm1+R,\qquad
 \gamma_\pm=\frac{f_\pm}{\bar f_\pm},
 \qquad N(f_+)=p z,\quad N(f_-)=q z.
 \tag{1}
\]
Let \(L=2n+3\) and \(r_0=1-3^{2n+2}\).
Take
\[
 \begin{split}
 p&\equiv1\pmod8,\qquad q\equiv5\pmod8,\\
 p,q&\equiv r_0\pmod{3^L},\qquad
 q\equiv4\pmod p,\qquad p/4\in\mathbb F_q^{\times m}.
 \end{split}
 \tag{2}
\]
These conditions have solutions for every fixed \(n\).
The prime \(p\) is chosen by Dirichlet. For \(q\), use
complete splitting in
\(\mathbb Q(\zeta_m,\sqrt[m]{p/4})\).
The cyclic radical subgroup has order \(m\) by the
valuation at \(p\), and complex conjugation acts on
it by inversion. Its maximal abelian subfield over
\(\mathbb Q\) is therefore \(\mathbb Q(\zeta_m)\).
The simultaneous rational progression in (2) asks
only for a Frobenius agreeing with this complete
splitting on the intersection, since \(r_0\equiv1\bmod m\).
Chebotarev in the compositum supplies \(q\ne p\).

**Proposition TP1.** The two independent order-\(m\)
boundary classes of \(R^2=F(u)\) have norm-one
representatives \(\gamma_+,\gamma_-\). There is a
nonempty local slice for \(u\), with \(u\) a unit at
all primes dividing \(6Apq\), giving all of the
following properties on quadratic specializations.

- The field \(K=\mathbb Q(\sqrt{F(u)})\) splits at
  \(2,3\) and ramifies at \(p,q\).
- Both \(\gamma_\pm\) are \(m\)-virtual and are local
  \(m\)-th powers above \(3\).
- The slice supports either real sign.

The four rational boundary points at \(u=0,\infty\)
give the two independent torsion classes of
[Gillibert–Levin, Lemma 3.3](sources/gillibert-levin-torsion-line-bundles.pdf),
after rescaling \(R\) by \(A\). This applies since
\(A\ne0\) and \(B^2-4A^2=pq\ne0\).
Dividing the functions \(f_\pm\) by their conjugates
doubles the corresponding divisor classes, so
retains independence for odd \(m\).

At \(2\), \(A,B\) are odd. Choose \(u\equiv -B^{-1}
\bmod8\); the odd power \(z\) then gives \(F(u)\equiv1\bmod8\).
At \(3\) impose \(u\equiv1\bmod3^L\). Then
\(A\equiv0\), \(B\equiv-r_0\bmod3^L\), and
\(F(u)\equiv3^{2n+2}\bmod3^L\). It is a square
locally and the two roots have valuation \(n+1\).
Both ratios in (1) are units congruent to one
modulo \(3^{n+1}\), hence local \(m\)-th powers.

At \(p\) take \(z\equiv-1/A\bmod p\), which is possible
with \(u\equiv1\bmod p\) by \(q\equiv4\bmod p\).
The identity \(F=(A z+1)^2-pz\) gives valuation one.
At \(q\), \(z\equiv1/A\bmod q\) is an \(m\)-th power
since \(1/A\equiv4/p\bmod q\). Now
\(F=(A z-1)^2-qz\) gives valuation one.
All norm-one values are units in these nonsplit
completions. At the other primes dividing \(2A\),
unit \(u\) gives integral \(f_\pm\) and unit norms
by (1).

At any remaining prime, if \(u\) is integral then
\(f_\pm\) are integral. For unit \(u\) their norms
are units. For positive \(u\)-valuation, \(R\) is
locally close to either \(1\) or \(-1\); one
conjugate of each \(f_\pm\) is a unit, and the
norm has valuation divisible by \(m\). For a pole
of \(u\) the field splits by the square leading
coefficient \(A^2\). On the branch with
\(R/(A z)\) close to one, \(f_\pm\) have valuation
\(m v(u)\), and the other branches are units by
the norms. This proves virtuality at every prime.
All conditions can be imposed simultaneously by
weak approximation. The quadratic polynomial in \(z\)
has two distinct positive roots, since \(B<0\),
\(B^2-4A^2=pq>0\) and its constant term is one.
Positive real \(u\) realizes either sign.

## Independent moving Frobenius and field counts

For \(m>1\), let \(V_m(K)\) denote the elements of
\(K^\times/K^{\times m}\) with all finite valuations
divisible by \(m\). At a split place with residue
\(1\bmod m\), remove the local \(m\)-power of a
uniformizer from \(\gamma_\pm\), and take the
two \(m\)-power residue characters.

**Theorem TP2.** With \(p,q\) fixed as in (2), there
are positive-density prime sets giving genuine
distinct \(P,Q\equiv1\bmod m\), \(u=P/Q\), on the
slice of TP1. At \(P\) and \(Q\) the two-dimensional
Frobenius vectors can be prescribed independently,
up to a unit scalar on each vector. In particular
one can make the two vectors a basis modulo \(m\).
Then \(\gamma_+,\gamma_-\) generate a copy of
\(C_m^2\) in \(V_m(K)\), and the mirror layer
\[
                 M=\mathbb Q(\sqrt{-3D}),\qquad
                 E=M\mathbb Q(\zeta_m)^+
                 \quad(D=\operatorname{disc}K)
\]
has a cyclic-factor quotient \(C_m^2\) of its ordinary
class group, completely split above three. Both
moving rational primes split completely in \(E\),
where their oriented Frobenius images give the
prescribed two vectors. For either base sign
the distinct fields satisfy
\[
 \#\{|D|\le X\}\gg_{m,p,q}
            \frac{X^{1/m}}{(\log X)^{C_m}},\qquad
 \#\{|\operatorname{disc}E|\le Y\}\gg_{m,p,q}
            \frac{Y^{3/m^2}}{(\log Y)^{C_m}} .
 \tag{3}
\]

At \(P\), use the branch \(R\) near \(1\).
The leading terms of \((\gamma_+,\gamma_-)\) give
the residue vector represented by
\[
              \left(\frac4p,\frac q4\right)
              \quad\text{after removing }
                     (z^{-1},z).
 \tag{4}
\]
At \(Q\), use \(R/(Az)\) near one. Here the residues
are represented by
\[
              \left(\frac{4A^2}p,\frac{4A^2}q\right)
              \quad\text{after removing }(z,z).
 \tag{5}
\]
Each pair gives independent Kummer directions over
\(\mathbb Q(\zeta_m)\) by valuations at the two fixed
primes \(p,q\), which are unramified in this
cyclotomic field and do not divide \(2A\).
Its governing field has group
\(C_m^2\rtimes(\mathbb Z/m)^\times\) over \(\mathbb Q\).
In each case the maximal abelian subfield is
\(\mathbb Q(\zeta_m)\) because inversion generates all
translations in the commutator. Thus the Frobenius
choices in these governing fields are compatible
with any fixed rational progression away from three,
and with \(P,Q\equiv1\bmod3^L\). Use such progressions
to impose all unit ratio conditions in TP1.
Chebotarev supplies \(\gg U^2/(\log U)^2\) pairs in
fixed rectangles \(P,Q\asymp U\) inside either sign
interval. Unit-scalar ambiguity is the cyclotomic
conjugation action on a vector; it preserves the
basis condition.

Independence modulo cubes follows already from the
two residue vectors modulo three. A noncube
combination over \(K\) stays noncubic in
\(K(\zeta_m)\) by the nonnormal cubic argument in
[SD2](CYCLOTOMIC_MIRROR_TORSION_FROM_SMALL_DICKSON_TRACE.md);
\(K\ne\mathbb Q(\zeta_3)\) because \(3\) splits in it.
Kummer theory then gives independent order-\(m\)
translation directions over the compositum. The
simultaneous involution inverts both radicals and
\(\zeta_m\), so is central. Its fixed field is a
\(C_m^2\) extension over \(E\), unramified by
virtuality away from three and completely split at
three by TP1. The odd degree also handles real
places of \(E\). As \(P,Q\equiv1\bmod m\), they
split in \(K(\zeta_m)\) by the indicated branches,
so the relative Frobenius vectors over \(E\) are
the residues just prescribed.

The squarefree binary form of \(F\) has degree
\(2m\), because its quadratic polynomial in \(z\)
has distinct nonzero roots. The prime-factor
truncation and the fixed twist-Jacobian multiplicity
argument of [PP2–PP3](ABC_COHERENT_PRIME_PROGRESSION_NORMS.md)
give at most a power of \(\log U\) retained pairs
per field, as in SD3. Here the genus is \(m-1\ge2\),
so uniform Mordell–Lang applies. The size of \(D\)
is \(O_{m,p,q}(U^{2m})\).
The unique quadratic subfield of \(E\) is \(M\).
As \(3\) splits in \(K\), the discriminant is
\[
                 |\operatorname{disc}E|
                       =3^{(2n-1)m/3}|D|^{m/3},
\]
by the local full cyclotomic completion at three
and the unramified real cyclotomic factor away
from three. This proves (3).
The quotient at \(m=3\) over the same \(M\) also
shows that every marked mirror has class \(3\)-rank
at least two.

## The exact absolute relations in the imaginary fields

**Corollary TP3.** On the imaginary base sign in TP2,
for \(|D|>4pq\), the actual two oriented moving prime
classes generate
\[
                         C_{4m}\times C_{2m}.
\]
Each moving prime class has full order \(4m\), and the
two \(3\)-primary parts have independent exact order
\(m\). Neither orientation makes the two prime
classes belong to a common coset of the \(m\)-torsion
subgroup. These assertions and the field count
refer to the same genuine prime ratios.

Choose the branch orientations as in (4) and (5),
and write \(\mathfrak p,\mathfrak q\) for the fixed
ramified primes over \(p,q\). The full principal
ideal formulas are
\[
      (f_+)=\mathfrak p(\bar{\mathfrak P}/\mathfrak Q)^m,
      \qquad
      (f_-)=\mathfrak q(\mathfrak P/\mathfrak Q)^m .
      \tag{6}
\]
Indeed at \(P\) one plus conjugate and the opposite
minus conjugate have valuation \(m\). At \(Q\)
the two selected conjugates have valuation \(-m\).
At \(p,q\) the norm in (1) gives the unique ramified
valuation one for the indicated function. At all
other primes \(u\) and the norms are units and
the functions are integral. In additive class
notation (6) gives
\[
 m([\mathfrak P]+[\mathfrak Q])=[\mathfrak p],
 \qquad
 m([\mathfrak P]-[\mathfrak Q])=[\mathfrak q].
 \tag{7}
\]
The sign on the second ramified class is immaterial,
since it has order dividing two.

In imaginary \(K\) all units are cubes, since
\(3\) splits, so the independent virtual classes
inject into the \(m\)-torsion of the class group.
Their divided ideals give the two differences and
sums of the moving classes multiplied by two.
As \(2\) is invertible modulo \(m\), this gives the
claimed independent \(3\)-primary parts.

The fixed ramified classes \([\mathfrak p]\) and
\([\mathfrak q]\) are independent in the imaginary
fields under consideration. A nonempty product
would be an integral ideal of norm \(p,q\) or \(pq\).
A nonrational algebraic integer in an imaginary
quadratic field has absolute norm at least
\(|D|/4\), and a rational generator has square norm.
Thus no such product is principal. Equation (7)
now gives each moving class exact even part of
order four, with shared double \([\mathfrak p]+[\mathfrak q]\)
and a difference of order two outside the cyclic
group generated by either even part. This proves
the subgroup and rules out the odd torsion coset
for every orientation, since each signed sum
has nonzero \(m\)-multiple by (7).
In the real base sign, a cyclic unit kernel may remove
one direction from the ordinary class group, while
both directions survive in \(W/E\) of TP2.

## Scope

The Jacobian subgroup is the classical one of
Gillibert–Levin, Lemma 3.3. Their Corollary 3.2
requires a rational Weierstrass point on the
chosen model, whereas this fixed even-degree
curve has no rational branch point: at the
auxiliary prime \(p\) the only possible root
would reduce to the double-root center and have
valuation one. The present construction uses
two explicit ordinary local norm conditions
and the genuine prime residues to make the
reflection and ideal classes simultaneous.
Yamamoto’s Part I, Theorem 2 gives quadratic
torsion with prescribed rational prime types.
The fixed auxiliary relation (7), the two moving
Frobenius prescriptions and the higher cyclotomic
mirror quotient are the additional conclusions here.
[Yamamoto, Part I, Theorem 2](https://ir.library.osaka-u.ac.jp/repo/ouka/all/10165/ojm07_01_04.pdf).

The constants and auxiliary primes vary with the fixed
exponent. No uniform summation over \(n,p,q\) is asserted.
The [exact checker](verify_two_moving_prime_classes_even_norm_curves.py)
replays norm identities, local coefficient fixtures, the two moving
residue frames, dyadic chart identities, genuine prime pairs,
and elliptic and ray class fixtures.

## Dyadic compact-type reduction and exact generic torsion

**Theorem TP4.** Let \(m>1\) be odd, \(d=v_2(A)\ge0\),
and \(v_2(B)=0\), with \(A\ne0\) and
\(\gcd(m,1+d)=1\). For the smooth projective curve
\(C:R^2=A^2u^{2m}+Bu^m+1\), the odd rational Jacobian
torsion is exactly \(C_m^2\), generated by the two
boundary directions of TP1. When \(m=3^n\) and
\(d=0\), this is the full rational torsion. For all
the specified coefficients, every rational norm-one
function with \(m\)-divisible divisor is generated
modulo \(m\)-th powers by \(\gamma_+,\gamma_-\).

Over two put \(p=2A-B\), a local unit, and
\(w=R+Au^m+1\). The curve is a tame cyclic cover of
the \(w\)-line:
\[
                   u^m=\frac{w(2-w)}{p-2Aw}.
                   \tag{8}
\]
Its branch points are \(0,2,p/(2A),\infty\).
A semistable model of the marked line has two
affine component charts with \(x=w/2\) and
\(z=p/(2Aw)\), and a node chart with \(v=1/x\),
\(s=1/z\), \(vs=4A/p\). The branch pairs reduce to
\(0,1\) on each component. The model is proper
as the closure of \((4A/p)xz=1\) in
\(\mathbb P^1\times\mathbb P^1\).

Make the tame extension adjoining \(\zeta_m\) and
\(\pi\) with \(\pi^m=2\), and put \(A_0=A/2^d\).
On the first component chart put \(U=u/\pi^2\).
On the second put \(V=1/(\pi^{2+2d}u)\).
Normalizing (8) gives the two equations near
their finite special points
\[
              U^m=\frac{x(1-x)}{p-4Ax},\qquad
              V^m=\frac{A_0^2 z(z-1)}{4Az-p}.
                    \tag{9}
\]
Both denominators are units on the corresponding
finite special charts. The special components
have equations \(U^m=x(1-x)\) and \(V^m=z(z-1)\).
They are smooth on those charts since the
derivative of the right hand side is nonzero in
characteristic two. Their smooth compactifications
have genus \((m-1)/2\), and each has one point at
infinity, since the orders at that point are
coprime to \(m\).

In the completed node chart, (8) is
\[
     u^m=\frac{s}{v}\frac{v-1}{A(1-s)},
                    \qquad vs=4A/p .
\]
Multiplying \(u\) by \(\pi^d\) and absorbing a unit
\(m\)-th root given by Hensel's lemma makes
\(u_{\mathrm{new}}^m=s/v\). After an unramified base
change take \(\delta^m=4A/p\). The normalized node
chart is
\[
       ab=\delta,\qquad v=a^m,\quad s=b^m,\quad
       u_{\mathrm{new}}=b/a .
\]
Indeed the ring \(ab=\delta\) is normal, and the
invariants under \(a\mapsto\zeta_m a,\, b\mapsto
\zeta_m^{-1}b\) are the base node ring. The
ratio \(b/a\) is a primitive eigenfunction,
so gives the required cyclic extension of degree
\(m\). The two components of (9) therefore meet
at exactly one node. The model has compact type.
Resolving the node of thickness \(2+d\) inserts a
chain of rational components. The relative Picard
theorem identifies the good special fiber of the
Jacobian with the product of the two component
Jacobians.
[Bosch–Lütkebohmert–Raynaud, §9.2 and Theorem 9.5/4](https://www.math.stonybrook.edu/~kamenova/homepage_files/Bosch_Raynaud_Neron_Model_tc.pdf).

The root of unity is unramified at two. A tame
inertia generator \(\pi\mapsto\zeta_m\pi\) fixes
the base coordinates and acts on \(U,V\) in (9)
by \(\zeta_m^{-2},\zeta_m^{-2-2d}\), up to choosing
the inverse generator. Both scalars are primitive
by the hypothesis on \(d\). Prime-to-two torsion
specializes injectively and equivariantly to
the product. On each factor the fixed subgroup
of this cyclic automorphism has order \(m\),
by the invariant-divisor argument of
[SD4](CYCLOTOMIC_MIRROR_TORSION_FROM_SMALL_DICKSON_TRACE.md)
applied to the three branch points with orders
\(1,1,-2\). Therefore odd rational torsion has
order at most \(m^2\). The two rational boundary
directions from Gillibert–Levin attain it.

When \(m=3^n\) and \(d=0\), the branch polynomial in \(u\)
reduces at two to
\[
                      u^{2m}+u^m+1=\Phi_{3^{n+1}}(u).
\]
It is irreducible over \(\mathbb F_2\) because
\(2\) has order \(2\cdot3^n\) modulo \(3^{n+1}\).
The \(2m\) branch points form a single Galois
orbit. An even subset representing a rational
\(2\)-torsion class cannot be exchanged with its
complement, as this would give the odd size \(m\).
It is therefore invariant as a subset and gives
the zero class. This eliminates the remaining
\(2\)-primary torsion.
For a norm-one function with \(m\)-divisible
divisor, its divided divisor class is in the
two-generator rational group just proved.
After removing these two classes, the residual
function is a rational constant times an
\(m\)-th power. The norm-one condition makes
that constant an \(m\)-th power, as in SD4.

The base change and normalization are the tame
cyclic-cover mechanism used in
[Helminck, Theorem 1.2 and §3](https://link.springer.com/article/10.1007/s00209-021-02959-5).
The explicit charts determine the tame inertia
characters as well as the reduction graph, which
are both needed for the rational torsion bound.

## Two nonisogenous elliptic twists at the same conductor

**Theorem TP5.** For every fixed \(m=3^n\), there
are two fixed elliptic curves over \(\mathbb Q\)
which are not \(\mathbb Q\)-isogenous and a fixed
degree-\(2m\) norm curve with full rational
Jacobian torsion \(C_m^2\), with the following
properties. On both signs of the genuine prime
ratios there are \(\gg_m X^{1/m}/(\log X)^{C_m}\)
distinct fundamental discriminants \(|D|\le X\)
for which the two fixed curves have positive-rank
\(D\)-twists of common conductor \(36D^2\). The
independent moving Frobenius prescriptions and
mirror quotient of TP2 hold, and the exact class
relations of TP3 hold on the imaginary side for
\(|D|>4pq\). The constants may depend on the
coefficient primes chosen for the fixed \(m\).

Use \(A=a^m\) with \(a>0\) even and divisible by
\(3^L\), \(p-q=4a^m\), and \(p,q\equiv r_0\bmod3^L\).
Such primes exist by Tao–Ziegler's relative
polynomial progression theorem with \(0\) and
\(4(2\cdot3^L r)^m\) as the two polynomials.
The residue class \(r_0\) has positive relative
prime density, and both polynomials have constant
term zero.
[Tao–Ziegler, Theorem 5](sources/tao-ziegler-narrow-progressions-primes.pdf).
Keep one choice for each \(m\).

For this version of TP1 impose \(u\equiv-1/a\bmod p\)
and \(u\equiv1/a\bmod q\). The two norm identities
give ramification at these primes as before.
At two, \(A\) is even and \(B\) is odd. The function
\(z\mapsto A^2z^2+Bz+1\) takes the odd residues
modulo \(32\) bijectively onto the even residues:
the difference of its values is \(z-z'\) times
the odd factor \(A^2(z+z')+B\). The odd power map
on the units modulo \(32\) is bijective. Choose
\(u\) odd with \(F(u)\equiv4\bmod32\), which gives
a square in \(\mathbb Q_2\). The norms in (1) are
odd. The slice at three, the remaining unit
conditions and the two sign intervals are those
of TP1. In this alternative slice the arguments
for virtuality, the residue frames (4)–(5), the
prime ratio counts and the oriented ideals (6)
all apply. In particular, \((p,q)\) stay coprime
to \(2A\) and the two valuation matrices stay
unimodular.

The odd torsion is \(C_m^2\) by TP4 since
\(m\mid v_2(A)\). To determine the two-primary
part, over \(\mathbb Q(\sqrt{pq})\) put
\(\beta=(-B+\sqrt{pq})/(2A)\) of norm one and
\(c=a^{m/3}\), so that \(A=c^3\).
If \(\theta^3=\beta\) in this quadratic field,
then \(N(\theta)=1\). The number \(c\theta\) is
integral because its cube \(A\beta\) is integral.
Thus \(T=c(\theta+\bar\theta)\) is an integer
with \(T^3-3c^2T=-B\). This would give
\[
           q=-B-2c^3=(T+c)^2(T-2c).
\]
Since \(q>0\) and \(c\ge2\), its displayed square
factor exceeds one, contradicting primality.
Now \(x^m-\beta\) is irreducible over the quadratic
field. Induct on \(n\): the norm of any root of
the irreducible \(x^{3^j}-\beta\) to the quadratic
field is \(\beta\), so the root is itself a
noncube in that intermediate field. A root of
\(F(u)\) contains the quadratic field and
satisfies \((au)^m=\beta\) or its conjugate.
Hence \(F(u)\) is irreducible over \(\mathbb Q\),
and the even-subset argument of TP4 eliminates
the rational two-primary torsion.

Write \(w=c u^{m/3}\) and \(\lambda=B/A\).
The intermediate genus-two equation is
\(R^2=w^6+\lambda w^3+1\). It has involutions
\((w,R)\mapsto(1/w,\pm R/w^3)\). For
\(s=w+w^{-1}\), the invariant functions
\(v_\pm=R(w\pm1)/w^2\) satisfy
\(v_\pm^2=(s\pm2)(s^3-3s+\lambda)\).
Putting \(h=3c\) gives the elliptic targets
\[
      E_+: Y^2=X^3+(hX+p)^2,\qquad
      E_-: Y^2=X^3+(hX-q)^2.
                    \tag{10}
\]
Their points over \(K=\mathbb Q(\sqrt D)\) are
\[
 \begin{split}
 X_+&=\frac{-pw}{c(w+1)^2},&Y_+&=\frac{-pR}{(w+1)^3},\\
 X_-&=\frac{-qw}{c(w-1)^2},&Y_-&=\frac{-qR}{(w-1)^3}.
 \end{split}
\]
Here \(w\ne\pm1\) on the chosen slice since \(w\)
has positive valuation at two. The points are
anti-invariant under the involution of \(K\).
With \(L_+=hX_++p\) and \(L_-=hX_--q\), direct
substitution gives
\[
       \frac{(Y_++L_+)^2}{X_+^3}=-\gamma_+^{-1},
       \qquad
       \frac{(Y_-+L_-)^2}{X_-^3}=-\gamma_- .
                    \tag{11}
\]
The left hand side is the class of \(Y-L\) modulo
cubes. The usual three-isogeny descent is a
homomorphism into \(K^\times/K^{\times3}\); by
(11) and a basis prescription at the moving
primes, both images are nonzero. A nonzero
anti-invariant point of order three would have
rational \(X\) and \(Y^2=-(hX+3b)^2/3\), with
\(b=p\) or \(-q\), by the division polynomial.
The point at \(X=0\) has rational nonzero \(Y\)
and is invariant. Since \(3\) splits in \(K\),
the quadratic field \(K\) is different from
\(\mathbb Q(\sqrt{-3})\). Thus there is no
anti-invariant point of order three. Any
anti-invariant torsion point of order divisible
by three would give one by multiplication.
Both points in (11) therefore have infinite
order and give positive rank on their twists.
[Cohen–Pazuki, Proposition 1.4 and Lemma 2.1](https://www.impan.pl/shop/publication/transaction/download/product/82564).

For \(Y^2=X^3+(hX+b)^2\) the invariants are
\(\Delta=16b^3(4h^3-27b)\) and
\(c_4=16h(h^3-6b)\). The discriminants of (10)
are \(-432p^3q\) and \(-432q^3p\). Both curves
are multiplicative at \(p,q\) since \(c_4\) is
a unit there. To see the exact conductor at
two, set \(y'=Y-(hX+b)\) and, over \(\pi^3=2\),
put \(X=\pi^2 x,\ y'=\pi^3y\).
The resulting good reduction is \(y^2+y=x^3\),
since \(h\) is even and \(b\) is odd. Tame inertia
acts by \(x\mapsto\zeta_3^{-2}x\) and fixes
\(y\). Its nonidentity action on the Tate
module has determinant one and no fixed
direction. The exponent at two is two.

At three, \(h\equiv0\bmod9\) and \(b^2\equiv1\bmod9\).
Over \(\pi^4=3\), put \(X=-1+\pi^2 x,\ Y=\pi^3y\).
The constant term after division by \(\pi^6\) is
\(((b-h)^2-1)/\pi^6\), which reduces to zero.
The linear coefficient reduces to one and the
quadratic coefficient to zero. The reduction
\(y^2=x^3+x\) is good. Tame inertia acts by
\((x,y)\mapsto(\zeta_4^{-2}x,\zeta_4^{-3}y)\),
whose square is negation. Its exponent is also
two. Both curves have good reduction outside
\(2,3,p,q\). Since \(D\) splits at two and three,
these conductor exponents are unchanged on
twisting. Each odd prime dividing \(D\)
contributes exponent two, including \(p,q\)
where the ramified twist of a multiplicative
curve has that exponent. This gives \(36D^2\).

The three-isogenous partner of \(E_+\) is
\(Y^2=X^3-3(hX+3q)^2\), the \((-3)\)-twist of
\(E_-\) by the explicit isogeny in
Cohen–Pazuki. A multiplicative prime gives
\(E_+\) a nonintegral \(j\)-invariant, so it has
no complex multiplication. A nonzero rational
isogeny of \(E_+\) to its \((-3)\)-twist would
give a nonzero endomorphism of \(E_+\) taken
to its negative by the quadratic involution.
This is impossible in its geometric endomorphism
ring \(\mathbb Z\). Thus \(E_+\) and \(E_-\) are
not rationally isogenous, and common twisting
preserves that distinction. On each sign the
distinct common conductors up to \(Z\) number
\(\gg_m Z^{1/(2m)}/(\log Z)^{C_m}\).

## The two isogeny Selmer groups and an inert dyadic slice

**Theorem TP6.** Keep the curves of TP5. Let \(D\)
be a nontrivial fundamental discriminant of either
sign, with \(2\) unramified, \(3\) split, and \(p,q\)
ramified in \(K=\mathbb Q(\sqrt D)\). Put
\(M=\mathbb Q(\sqrt{-3D})\). For
\(E=E_\pm^{(D)}\), let
\(\phi:E\to\widehat E\simeq E_\mp^{(-3D)}\)
be the three-isogeny marked by \(X=0\).
Both choices of \(E\) have exactly the same two
isogeny Selmer groups under their respective
Kummer identifications. Write
\(S_K=\operatorname{Sel}_{\widehat\phi}(\widehat E)\)
and \(S_M=\operatorname{Sel}_{\phi}(E)\).
Each group is the cubic character group of the
group in the corresponding entry:
\[
\begin{array}{c|cc}
 &S_K&S_M\\ \hline
 D\equiv1\pmod8&
 \mathrm{Cl}_{(6)}(M)&
 \mathrm{Cl}(K)/\langle[\mathfrak t_K]\rangle\\
 D\equiv5\pmod8&
 \mathrm{Cl}_{(3)}(M)/\langle[\mathfrak t_M]\rangle&
 \mathrm{Cl}_{(2)}(K).
\end{array}
\tag{12}
\]
The ray class groups are ordinary and \((b)\)
denotes the integer modulus \(b\) in that field.
The ideal \(\mathfrak t_F\) denotes either prime
above two in the indicated split field, with its
class in the appropriate ordinary or ray group.
Their cube-quotient dimensions satisfy
\[
             \dim S_K-\dim S_M
                  =\mathbf1_{D>0}+
                     \begin{cases}
                       1 & D\equiv1\bmod8,\\
                      -1 & D\equiv5\bmod8 .
                     \end{cases}
\tag{13}
\]
In particular the finite \(3\)-Selmer dimensions
of both curves and both isogenous partners are
odd for negative \(D\) and even for positive
\(D\).

There is a genuine prime-ratio slice of TP5 with
\(D\equiv5\bmod8\). It has the same field and
common-conductor counts, higher mirror quotient
and moving prime relations as in TP5. On its
imaginary side all four positive-rank curves
\(E_\pm^{(D)},E_\pm^{(-3D)}\) have conductor \(36D^2\)
and finite \(3\)-Selmer dimension at least three.
On the real side all four dimensions are at
least two. The four curves have distinct
\(\mathbb Q\)-isomorphism classes.

To construct this slice impose \(F(u)\equiv20\)
modulo \(32\) instead of \(4\). The same bijection on
odd residues proves existence. Then \(K_2\) is the
unramified quadratic extension of \(\mathbb Q_2\).
The ratios \(\gamma_\pm\) are units reducing to
one there: \(R\) has valuation one and \(Au^m\)
is even. They are local cubes, since cube
exponentiation is bijective on principal units
at two. In the split slice they are local
cubes as well. The valuation, class and field
arguments of TP2–TP3 use virtuality at two and
complete splitting at \(3\). The two special
residue frames at \(P,Q\) stay the same. An
unramified quadratic twist at two preserves the
conductor exponent in TP5. This proves all the
previous conclusions on the new slice.

For (12), use the Kummer models from
Cohen–Pazuki. At every odd prime \(\ell\ge5\)
dividing \(D\), both local quadratic algebras
are ramified and their anti-invariant cube
quotients vanish: conjugation fixes the
valuation and the residue field, and the
principal units are cubes. This includes
the multiplicative coefficient primes \(p,q\).
Elsewhere at \(\ell\ge5\) both elliptic models
have good reduction, so the prime-to-\(\ell\)
local images are the unramified groups.
At infinity both local cohomology groups
vanish.

At two, when the marking algebra is split,
the model is \(Y^2=X^3+(hX+b)^2\) with \(h\)
even and \(b\) odd, up to local square
coordinate changes. The marked point has
Kummer class \(2b\) or its inverse. It
generates the full valuation line of
\(\mathbb Q_2^\times/\mathbb Q_2^{\times3}\).
In the nonsplit algebra the local model is
\(Y^2=X^3-3(hX+3b)^2\) for an odd \(b\).
An even integral \(X\) would give \(Y^2\equiv5\)
modulo \(8\), and odd integral \(X\) gives
even \(Y\). The Kummer argument
\(Y-\sqrt{-3}(hX+3b)\) in that case is a unit
with residue one. For nonintegral \(X\), the
ordinate dominates the irrational term and
its valuation is a multiple of three, so
the normalized residue is again rational.
This gives the zero image on the
unramified quadratic side.

At three the split model has Kummer image
the unit line. For integral \(X,Y\), the
factors \(Y\pm(hX+b)\) of \(X^3\) cannot both
be divisible by three. For nonintegral \(X\),
the ordinate dominates \(hX+b\) and its
valuation is divisible by three. At the
marked point \(2b\) is a noncube unit since
\(b\equiv\pm1\bmod9\).

For the ramified model at three, put
\(\pi=\sqrt{-3}\) and \(L_X=hX+3b\), where
\(b\equiv\pm1\bmod9\). The Kummer class has
norm a cube, so its square is represented
by \((Y-\pi L_X)/(Y+\pi L_X)\).
For \(v_3(X)\le0\), the valuation of
\(\pi L_X/Y\) in the \(\pi\)-scale is at
least three, with equality for unit \(X\).
Valuation \(v_3(X)\ge2\) is impossible,
since \(Y^2\) would have valuation three.
If \(X=3t\) with \(t\) a unit, then
\[
         Y^2=27\bigl(t^3-(b+ht)^2\bigr).
\]
The residue modulo three requires \(t\equiv1\)
modulo three. The difference is then zero
modulo nine, so \(v_3(Y)\ge3\), allowing
\(Y=0\). Now \(Y/(\pi L_X)\) has valuation
at least three, and one may discard the
rational sign in the ratio. Thus the image
is contained in the line of \(1+\pi^3\mathcal O\).
This line is unramified and nonzero. Indeed,
taking a cube root of \(1+\pi^3 c_0\) reduces
after the substitution \(1+\pi z\) to the
separable equation \(z^3-z-\bar c_0=0\).
For \(c_0\) a unit this cubic is irreducible
over \(\mathbb F_3\), and if \(c_0\) is in the
maximal ideal Hensel's lemma supplies a cube
root. The unit \(X=1\) supplies a point with
valuation exactly three in the ratio, since
\(1-3(h+3b)^2\) is a square in \(\mathbb Q_3\).
This proves equality with the unramified line.

Reflect the two Kummer spaces through
\(K(\zeta_3)=M(\zeta_3)\) as in
[CD3](ABC_CUBIC_TANGENT_COMPANION_DESCENT.md).
The split unit line at three reflects to the
anti-invariant character line of conductor
dividing \((3)\) in \(M\), and the unramified
line in \(M_3\) reflects to unramified
characters in \(K\). At two, the split
valuation line reflects to the tame
anti-invariant characters of conductor
dividing \((2)\) in the unramified quadratic
field. The zero line in that quadratic
algebra reflects to complete splitting
in the split marking field. At all other
finite primes the conditions are unramified.
The cubic characters of the listed ray
groups are anti-invariant since a rational
norm prime to \(6\) can be signed to be
\(1\) modulo the rational modulus. This
gives the full groups (12), with no
additional conditions at real places.

The size ratio for \(\phi\) at three is one:
its local Kummer image has order three and
its kernel has three rational points.
At two these sizes are \(1,3\) on the split
\(K\) side and \(3,1\) on the inert \(K\)
side. All other finite ratios are one by
the local-image description. At infinity
the ratio is \(1/3\) for real \(K\) and one
for imaginary \(K\). Neither isogeny kernel
has a nonzero global rational point, and
neither curve has rational three-torsion,
by the division polynomial argument in
TP5 applied to both nonsquare marking
fields. Cassels' ratio formula gives
(13), and the three-isogeny parity theorem
gives the finite Selmer parities for the
curves and their partners.
[Bhargava–Klagsbrun–Lemke Oliver–Shnidman, §9](sources/bhargava-klagsbrun-lemke-oliver-shnidman-three-isogeny.pdf).

Both independent \(\gamma\) classes are in
\(S_K\) on the genuine slices: they have
unramified Kummer images away from two and
three and are local cubes at two and three.
For imaginary \(K\) with \(2\) inert, (13)
then gives \(\dim S_M\ge3\), which injects
into \(\operatorname{Sel}_3(E)\). On the
partners, \(S_K\) injects into their
finite three-Selmer groups, so dimension
at least two and odd parity give dimension
at least three. For real \(K\) with \(2\)
inert, both isogeny dimensions are at least
two and give the two finite Selmer bounds.
Positive ranks and conductor are preserved
by isogeny. Within either isogeny class,
the two curves are nonisomorphic since a
degree-three self-isogeny of a non-CM
elliptic curve is impossible. Across the two
isogeny classes, use nonisogeny in TP5.

## The exact rational elliptic part of the Jacobian

**Theorem TP7.** Keep \(p=2A-B\) and \(q=-2A-B\)
distinct primes at least five, with the dyadic
hypotheses of TP4. For \(m=3^n\), every rational
elliptic quotient of the Jacobian factors up
to isogeny through the genus-two quotient
\[
             R^2=A^2 v^6+Bv^3+1,\qquad
             v=u^{m/3}.
\tag{14}
\]
If \(A\) is a rational cube, this quotient
accounts for exactly the two rational elliptic
isogeny factors \(E_+,E_-\) of (10), each with
multiplicity one. If \(A\) is a noncube, the
genus-two Jacobian is simple over
\(\mathbb Q(\zeta_3)\), and there is no rational
elliptic quotient. In particular, passing to
higher \(n\) in TP5 supplies no further
elliptic direction from a rational morphism of
the fixed norm curve.

Use the two components of the compact-type
model in TP4. On each, the subgroup generated
by the third power of the inertia generator
has order \(m/3\). Its quotient curve is
\(U_1^3=x(1-x)\), or its counterpart on
the other component, and has genus one.
On the rational Tate module for an auxiliary
prime ℓ different from two, the invariants of
the subgroup are exactly
the Tate module of the quotient Jacobian:
pullback and norm compose to multiplication
by the subgroup order on the quotient,
and to the sum of the automorphisms upstairs.
Consequently the inertia subgroup of order
\(m/3\) has a fixed space of dimension exactly
four in the Tate module of the whole Jacobian.
The genus-two quotient (14) supplies four
such dimensions, since its own dyadic
inertia has order three by the same charts.
The complementary rational isogeny factor
therefore has no such invariants.
[Milne, Proposition 10.1](https://jmilne.org/math/CourseNotes/AV.pdf).

For an elliptic isogeny factor, dyadic
inertia acts on its potentially good
reduction through automorphisms preserving
the origin, of order dividing \(m\). Good
reduction after the same tame extension follows
from the Tate-module criterion.
[Serre–Tate, Theorem 1](https://annals.math.princeton.edu/1968/88-3/p05).
Such an
automorphism \(\tau\) satisfies a monic
integral quadratic: from the endomorphism
degree formula,
\(\tau+\tau^{-1}=[2-\deg(1-\tau)]\).
The formula also holds when \(\tau=1\).
Thus an automorphism of odd three-power
order has order dividing three. Its Tate
module lies in the four-dimensional space
already supplied by (14). This proves
the first assertion via the rational
isogeny decomposition.

For the genus-two analysis put
\(C^0:R^2=w^6+(B/A)w^3+1\).
Its two elliptic quotients are defined over
\(\mathbb Q\), with equations (10) after
replacing \(h,p,q\) by \(3,p/A,q/A\).
The two reciprocal involutions act by
opposite signs on each of the two
differential eigenspaces. In the basis
\(dw/R,w\,dw/R\), the first involution
interchanges the basis elements and
negates both, and the second acts as
its negative. Thus the Jacobian of
\(C^0\) is isogenous over \(\mathbb Q\)
to the product of the two quotients.
At \(p\) and \(q\) both have multiplicative
reduction by the elementary \(c_4,\Delta\)
calculation in TP5. The isogeny computation
of TP5 gives two non-CM curves which are
nonisogenous over \(\mathbb Q\) and
isogenous over \(F=\mathbb Q(\zeta_3)\).
Over \(F\) this Jacobian is isogenous to the
square of a non-CM elliptic curve. All elements
of its geometric endomorphism algebra
\(M_2(\mathbb Q)\) are therefore defined over \(F\).

The isomorphism from (14) to \(C^0\)
takes \(w=cv\) with \(c^3=A\).
For \(A\in\mathbb Q^{\times3}\) it is rational
and gives the stated two factors. Otherwise
\(F(c)/F\) has degree three, since adjoining a
root of a rational noncube first has degree
three. A generator acts on the isomorphism
by the automorphism \(\rho:w\mapsto\zeta_3w\)
of \(C^0\). Hence its action on geometric
endomorphisms of (14) is conjugation by
the order-three matrix of \(\rho\).
This matrix is nonidentity, since
\(C^0/\rho\) is a conic and its invariant
Tate module is zero. Its polynomial
over \(\mathbb Q\) is \(t^2+t+1\) and its
centralizer in \(M_2(\mathbb Q)\) is the
field \(\mathbb Q(\zeta_3)\), which has no
nontrivial idempotent. The genus-two
Jacobian of (14) is therefore simple over
\(F\). This excludes all rational
elliptic quotients in the noncube case.
