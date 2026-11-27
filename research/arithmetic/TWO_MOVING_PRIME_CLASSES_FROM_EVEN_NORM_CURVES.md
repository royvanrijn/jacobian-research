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
residue frames, dyadic chart identities and the elliptic formulas.

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
