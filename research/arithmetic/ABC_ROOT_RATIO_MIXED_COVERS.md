# Mixed valuation grids from the ratio of two power roots

Write a primitive positive triple in one of its three signed orientations as
\[
 U=A x^d,\qquad V=B y^d,\qquad
 T=\epsilon U+\eta V>0,\qquad \epsilon,\eta\in\{1,-1\},
\]
where \(d\ge2\) and \(A,B,x,y\) are positive integers. Fix a prime \(p\) and
orders \(m_0,m_\infty,m_T\) that are powers of \(p\), allowing order \(1\).
For a squarefree integer \(S\), impose
\[
 m_0\mid v_q(x),\qquad m_\infty\mid v_q(y),\qquad
 m_T\mid v_q(T)\qquad(q\nmid S).
 \tag{RR1}
\]
The relevant hyperbolicity condition is
\[
          \frac1{m_0}+\frac1{m_\infty}+\frac d{m_T}<d .
 \tag{RR2}
\]
Fix \(d,p\) before the height limit. For fixed \(d,p\) the constants can be
uniform in the three orders. Indeed, the triples of their exponents that
satisfy (RR2) form an upward-closed subset of \(\mathbf Z_{\ge0}^3\). It
has finitely many minimal elements by Dickson's lemma, so one can always
use one of finitely many covers.

## The height-free support bound

For fixed \(A,B\) and either orientation satisfying (RR2), the number of
primitive triples satisfying (RR1), with no height bound, is at most
\[
                         \exp O(1+\omega(S)).
 \tag{RR3}
\]
The constant may depend on the fixed data. More generally, take \(A,B\) to
be the canonical \(d\)-free kernels. Let \(d_0\) be the part of \(d\) prime
to \(p\). Impose that the class of \(-\eta B/(\epsilon A)\) in
\(\mathbf Q^\times/\mathbf Q^{\times d_0}\) lies in a prescribed **fixed finite
set**. For each squarefree \(H\), the number of solutions with
\(\operatorname{rad}(ABS)\mid H\) is then at most \(\exp O(1+\omega(H))\),
even with varying \(A,B,S\) and with no height bound. If \(d_0=1\) the class
condition is automatic. In this case the signature is already covered by
the prime-primary bound on the original line.

Summing over the support gives at most \(O_\xi(Z^{1+\xi})\) primitive
triples with combined support radical at most \(Z\), for every \(\xi>0\).
For fixed kernels one may instead take \(Z\) to bound just the defect
radical \(S\). In a small-budget class with \(c\le X\) and
\(Z\le X^{f\lambda}\), the bound is \(O_{\lambda,\varepsilon}
(X^{f\lambda+\varepsilon})\), with the fixed data suppressed. Thus
\(f\le 1/3-\delta\) yields a strict saving with that fixed margin on
this sector. On the least-repair formulation the same count holds for
primitive quotients and for their uniquely chosen least repairs. If
\(Z\le Y^{f/M}\) in repaired height \(Y\), the resulting exponent is
\(f/M+\varepsilon\). The count has no multiplier factor.
Under the balanced repair hypotheses, if \(Z\le c^{f\lambda}\),
\(\operatorname{rad}(c)\ge c^{(1/3-\tau)\lambda}\) and \(H_M(c)\le Y\),
then \(Z\le Y^{f/[M(1/3-\tau)]}\). This gives a coefficient below
\(1/M\) for \(f<1/3-\tau\).
[Least target repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

For example, take two fixed weighted squares (\(d=2\)). With odd \(p\ge5\)
one can impose \(p\)-divisibility on the first square root and on the
third letter, leaving the second square root free. This realizes the
mixed signature \((2p,2,p)\) away from \(S\). With \(p=3\), taking
\((m_0,m_\infty,m_T)=(3,1,9)\) or \((3,3,3)\) covers \((6,2,9)\)
or \((6,6,3)\) away from \(S\). These conclusions allow coefficients
in the high root powers and the third letter to vary on \(S\).
None of these signatures has a hyperbolic prime-primary signature on
the original triple ratio. With fixed weighted cubes one can use
\((m_0,m_\infty,m_T)=(1,1,5)\), for example, yielding \((3,3,5)\).
The sign choices include both same-side and opposite-side power letters.

## An auxiliary cover of the root ratio

Fix \(\rho\) with \(\rho^d=-\eta B/(\epsilon A)\), and put
\[
 z=x/(\rho y),\qquad T=\eta B y^d(1-z^d).
 \tag{RR4}
\]
For \(A,B\) fixed, \(z\) determines the primitive triple since
\(U/V=(A/B)(\rho z)^d\) determines the two primitive integral letters,
and the sign orientation determines \(T\).

On the projective \(z\)-line assign index \(m_0\) at zero,
\(m_\infty\) at infinity and \(m_T\) at each of the \(d\) roots of
unity. There exists a geometrically connected regular cover
\(\pi:C\to\mathbf P^1\) of genus at least two with a finite \(p\)-group
of deck transformations, whose inertia orders divide these indices.
The curve and map depend only on the fixed orders and \(d\).
Riemann–Hurwitz gives the necessity of (RR2).

For sufficiency discard sections of index \(1\), and let \(r\) be the
number left. Condition (RR2) forces \(r\ge3\). For \(r=3\), use the
three-point construction, after a fixed Möbius change of variable.
[Three-point prime-group covers](ABC_MISMATCHED_GRID_NUMBER_FIELDS.md#an-ell-group-covering-with-controlled-inertia)
If \(p\) is odd and \(r\ge4\), choose four sections and use a cyclic
degree-\(p\) cover with local Kummer exponents \(1,-1,1,-1\).
Its genus is \(p-1\). For \(p=2,r\ge6\), a double cover branched at any six
sections has genus two. For \(p=2,r=5\), use two independent quadratic
equations with inertia vectors
\((1,0),(1,0),(1,0),(0,1),(1,1)\) at the five sections. The vectors sum
to zero, their coordinates prescribe the two squarefree branch products,
and their span gives degree four. Riemann–Hurwitz gives genus two. For
\(p=2,r=4\), one index must be at least four. Double-cover the line
branched at the four sections, take an unramified degree-two isogeny
onto the resulting genus-one curve, and double-cover the isogeny source
branched at the two lifts of the selected higher-index point. This gives
genus two and indices at most four at that section and two elsewhere.
The last double cover exists by divisibility of the degree-zero Picard
group over \(\overline{\mathbf Q}\). Take the Galois closure of any
tower used here. It is still a \(p\)-cover, and its local indices are
the least common multiples of the initial local indices. All data
descend to a fixed number field \(K_0\), which we enlarge to contain
the branch coordinates, all deck transformations, \(\mu_p\),
the coordinates of a basepoint on \(C\), and \(J[p]\) for
\(J=\operatorname{Jac}(C)\).

## Specialization and the moving kernel cost

Put \(F=K_0(\rho)\). Include the primes over \(ABS\) and a fixed set
containing \(dp\), bad reduction and the exceptional primes of the
specialization theorem in \(\mathcal S_F\). Away from this set, \(\rho\)
is a unit and the \(d\) roots of unity have distinct reductions. At a
prime where \(x\) or \(y\) is divisible, primitivity makes the other a
unit. The contact order of \(z\) with zero or infinity is thus the local
valuation of \(x\) or \(y\). When \(z\) meets a root of unity, both
\(x,y\) are units and exactly one of those roots can be met. Its contact
order equals the local valuation of \(T\) by (RR4). Extending a valuation
from \(\mathbf Q\) to \(F\) multiplies these orders by the ramification
index, preserving the needed divisibility. The specialization theorem
gives a fibre field \(E/F\) unramified outside \(\mathcal S_F\).
The fibre is a torsor under the constant deck group, so \(E/F\) is
Galois of bounded \(p\)-power degree.
[Conrad, Theorem 0.1 citing Beckmann, Proposition 4.2](https://math.stanford.edu/~conrad/papers/inertiapaper.pdf)

The following bound makes the varying-base cost explicit:
\[
 \#\{\text{primitive triples for }A,B,S\}
 \ \le\
 \exp O_{d,p,m_0,m_\infty,m_T}\bigl(
      1+\omega(ABS)+\dim_{\mathbf F_p}\operatorname{Cl}(F)[p]\bigr).
 \tag{RR5}
\]
Indeed, let \(G\) be the Galois group of the maximal pro-\(p\) extension
of \(F\) unramified outside \(\mathcal S_F\). Its generator number is
the dimension of \(\operatorname{Hom}(G,\mathbf F_p)\) by the pro-\(p\)
Frattini theorem. By Kummer theory the latter is bounded by the unit
rank, \(\#\mathcal S_F\) and the \(p\)-rank of \(\operatorname{Cl}(F)\).
Here \([F:\mathbf Q]\) is bounded and \(F\) contains \(\mu_p\).
Each specialization is a homomorphism from \(G\) to a fixed finite
deck group after choosing a geometric fibre point. Therefore the
number of fields \(E\) is at most the group order raised to the
generator bound. The open subgroup over \(E\) has generator number
at most the index times that bound plus one, by Schreier's
inequality. This subgroup is the full Galois group of the maximal
pro-\(p\) extension of \(E\) unramified away from the primes over
\(\mathcal S_F\). In fact the normal closure over \(F\) of a finite
Galois \(p\)-extension of \(E\) with that ramification is still a
\(p\)-extension with that ramification, since \(E/F\) is Galois of
\(p\)-power degree.

The Kummer injection for \([p]\) on \(J(E)\) is unramified outside
these places, by finite étaleness of \([p]\) on the good reduction
abelian scheme. Since \(J[p]\) is constant, the Kummer target is
bounded in dimension by \(2g(C)\) times the generator number of the
open subgroup. Thus \(\operatorname{rank}J(E)\) is bounded by a
constant times the quantity in (RR5). Uniform Mordell–Lang bounds
\(\# C(E)\) exponentially in the rank. At least one point of \(C(E)\)
lies over each occurring \(z\), which proves (RR5).
[Gao–Ge–Kühne, Theorem 1.1′](https://arxiv.org/pdf/2105.15085v4)

For fixed \(A,B\), (RR5) immediately gives (RR3). For the finite
class variant write \(d=p^j d_0\), and write
\(-\eta B/(\epsilon A)=\kappa w^{d_0}\) with \(\kappa\) fixed in the
chosen class and \(w\in\mathbf Q^\times\). One may take
\(\rho=\kappa^{1/d} w^{1/p^j}\). After enlarging the fixed \(K_0\)
to contain \(\kappa^{1/d}\) and \(\mu_d\), \(F/K_0\) is Galois of
bounded \(p\)-power degree, unramified outside fixed primes and those
dividing \(AB\). The same open-subgroup argument, starting over
this fixed base and using Kummer there, bounds
\(\dim\operatorname{Cl}(F)[p]\) by \(O(1+\omega(AB))\).
At each prime of \(H\), the canonical \(d\)-free kernel exponents
and inclusion in \(S\) give a bounded number of choices. This proves
the assertion for \(H\) and hence the support-budget count, using
\(\exp O(\omega(H))\ll_\xi H^\xi\).

For fully varying kernels, (RR5) reduces the additional cost to a
weighted class-group moment over the Kummer fields \(K_0(\rho)\).
Already for \(d=2\) and odd \(p\) it involves the \(p\)-torsion of
moving quadratic extensions of the fixed base. To obtain the
defect-budget exponent by (RR5), it suffices to bound
\[
 \sum_{\substack{A,B\text{ canonical }d\text{-free}\\
                  (A,B)=1,\ \operatorname{rad}(AB)\le Z}}
       \frac{|\operatorname{Cl}(K_0(\rho))[p]|^K}{\operatorname{rad}(AB)}
       \ \ll_{\xi,K}\ Z^\xi
 \qquad(\xi>0)
 \tag{RR6}
\]
for the fixed moment order \(K\) needed in (RR5). In fact for
\(r=\operatorname{rad}(AB)\) at most \(Z/r\) possible total supports
\(H\le Z\) are divisible by \(r\). One can apply (RR5) with \(S=H\),
absorbing \(\exp O(\omega(H))\) in an arbitrarily small power of \(Z\).
No estimate for (RR6) is assumed in the unconditional statements.

For the height-free support bound in signature \((6,6,3)\), the
moving fields in this criterion can be taken to be biquadratic.
The bounded-height defect-budget bound for fully moving kernels
already follows from the Euclidean \((3,3,3)\) signature, since
each order is a multiple of three.
[Euclidean grid point count (EP6–EP7)](ABC_EUCLIDEAN_GRID_POINT_COUNTS.md#primitive-quotients-on-two-euclidean-signatures)

Take the cyclic cubic cover
\[
                         w^3=z/(1-z^2)
\]
branched at \(0,\infty,1,-1\), of genus two. The cover and its deck
group are defined over \(K_0=\mathbf Q(\zeta_3)\). A ramified point
gives a rational basepoint, and \(J[1-\zeta_3]\) is constant over
\(K_0\) because it is generated by differences of ramification
points. One may apply the branch-divisor description after a
fixed Möbius change moving infinity away from the branch locus.
Use the Kummer map for \(1-\zeta_3\) in place of \([3]\). The map
is finite étale away from \(3\) at good primes, and
\([3]\) is a unit times \((1-\zeta_3)^2\), so the same rank argument
and (RR5) hold over \(F=\mathbf Q(\zeta_3,\rho)\).
[Poonen–Schaefer, Proposition 6.2](https://math.mit.edu/~poonen/papers/descent.pdf)

Here \(\rho^2\) has squareclass \(\kappa=-\eta AB/\epsilon\).
If the two quadratic extensions are distinct and nontrivial,
the odd class torsion of \(F\) splits into the three character parts
from its quadratic subfields. Indeed extension and norm identify
the invariants of each order-two subgroup with the torsion of the
fixed field, since two is invertible on the \(3\)-torsion. The full
invariant part is zero by the same argument over \(\mathbf Q\).
The only variable factors are thus the \(3\)-torsion of
\(\mathbf Q(\sqrt{\kappa})\) and \(\mathbf Q(\sqrt{-3\kappa})\).
If \(F\) is already quadratic its contribution is fixed. For
each \(AB=r\), the number of kernel allocations is at most
\(2^{\omega(r)}\). Cauchy–Schwarz and dyadic summation show that a
subpower normalized moment bound for the \(2K\)-th powers of
quadratic \(3\)-torsion would suffice for (RR6) in this signature.
Each squarefree parameter \(\kappa\) or the squarefree part of
\(-3\kappa\) has absolute discriminant comparable to \(r\). The
first-moment theorem for quadratic \(3\)-torsion alone does not
supply this higher-moment bound.
[Davenport–Heilbronn, Theorem 3](https://www.math.utoronto.ca/~ila/DH2.pdf)

## Cyclic isogeny descent for \((2p,2,p^a)\)

Fix an odd prime \(p\) and an integer \(a\ge1\) with
\(\varphi(p^a)\ge4\), and put \(n=p^a\).
For signature \((2p,2,n)\), take \(d=2\) in (RR1) and use the cover
\[
           w^n=z^{n/p}(z-1)(z+1)^{n-n/p-1}.
\]
The branch indices at \(0,1,-1\) are \(p,n,n\), respectively.
The point at infinity is unramified because the polynomial has
degree \(n\). The exponent at \(1\) proves geometric connectedness.
Riemann–Hurwitz gives genus \(g=\varphi(n)/2\).
The cover is regular over \(K_0=\mathbf Q(\zeta_n)\), with deck
generator \(\sigma:w\mapsto\zeta_n w\).

The quotient by the subgroup of order \(p\) has equation
\[
           (w^p)^{n/p}=z^{n/p}(z-1)(z+1)^{n-n/p-1}.
\]
The only possible branch points are \(1,-1\), with index \(n/p\)
at each, so the quotient has genus zero. The norm followed by
pullback on the Jacobian is
therefore the zero endomorphism
\[
              1+\sigma^{n/p}+\cdots+\sigma^{(p-1)n/p}
                 = \Phi_n(\sigma).
\]
Consequently \(\mathbf Z[\zeta_n]\) acts on the Jacobian \(J\), with
\(\zeta_n\) acting as \(\sigma\). In this ring \(p\) is a unit times
\((1-\zeta_n)^{\varphi(n)}\). Thus \(\phi=1-\sigma\) is an
isogeny and \([p]=u\phi^{\varphi(n)}\) for an automorphism \(u\).
Taking degrees gives \(\deg\phi=p\), since
\(\deg[p]=p^{2g}\). The two unique points \(P_+\), \(P_-\) over
\(1,-1\) are rational over \(K_0\) and fixed by \(\sigma\).
Their difference defines a nonzero class because a degree-one
function would make the curve rational. It generates \(J[\phi]\),
so the kernel is constant over \(K_0\).

With this isogeny the specialization argument for (RR5) uses
\(F=K_0(\sqrt{\kappa})\), where \(\kappa=-\eta AB/\epsilon\)
up to rational squares. At good places away from \(p\), the
\(\phi\)-Kummer image is unramified. The identity for \([p]\) bounds
the dimension of \(J(E)/pJ(E)\) by \(\varphi(n)\) times the dimension
of \(J(E)/\phi J(E)\). The remaining field count and rank bound in
(RR5) therefore apply with this \(F\), without adjoining \(J[p]\).

At \(p=3\), the field \(\mathbf Q(\zeta_{3^a},\sqrt{\kappa})\)
is a bounded cyclic \(3\)-extension of the field
\(\mathbf Q(\zeta_3,\sqrt{\kappa})\), unramified outside \(3\).
The open-subgroup argument applied over that base with the places
above \(3\) allowed bounds the upper \(3\)-class rank by a constant
times one plus the lower \(3\)-class rank. The biquadratic
decomposition and the dyadic summation above give the same
quadratic \(2K\)-th moment criterion for each fixed signature
\((6,2,3^a)\), with the required moment order depending on \(a\).
Scholz reflection allows both quadratic factors to be bounded by
the \(3\)-torsion of the imaginary member of the reflected pair.
[Scholz reflection in Ellenberg–Venkatesh, Section 3.3](https://people.math.wisc.edu/~ellenberg/Scholz-submit-july17.pdf)

For \((6,2,3^a)\) this moment would give a moving-kernel defect-budget
bound that is also new at bounded height. The signature has no
Euclidean or hyperbolic prime-primary signature dividing its three
orders. Joint signed dictionaries and the first moment already give
an unconditional small-budget saving when the odd-kernel heights
are small or their degrees concentrate in a bounded number of windows.
[Unconditional moving-kernel sector](ABC_JOINT_FREE_LETTER_DICTIONARIES.md#an-unconditional-sector-with-fully-moving-square-kernels)

The earlier cyclic cubic cover for \((6,6,3)\) strengthens the
height-free statement only.

## Why splitting the remaining branch is insufficient by itself

For \(d=2\) with \(m_0=m_\infty=1\), (RR2) always fails. This is
reflected in actual primitive solutions. For any fixed odd \(q\ge3\)
and every odd integer \(s>1\), the integers
\[
 x=(s^q-1)/2,\qquad y=(s^q+1)/2
\]
give \(y^2-x^2=s^q\). The roots are consecutive and
\(\gcd(s,xy)=1\), so the triples are primitive. Both square kernels
are one and the third letter is an exact \(q\)-th power. The empty
defect support admits infinitely many triples with no height bound.
At least one further divisibility condition or another arithmetic
input is required for a height-free fibre bound on this signature.

Darmon and Granville prove qualitative finiteness for fixed
coefficients of hyperbolic generalized Fermat equations and for
suitable binary forms. Their specialization theorem controls
ramification. Here the counting estimate comes from a cover with
a single prime in the deck group and a bound exponential in the
number of support primes. The additional step from the three-point
prime-primary bound is to pass to the ratio of two roots before
building the cover.
[Darmon–Granville, Theorems 1 and 2](https://www.math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf) It supplies a
mixed-order support estimate and an explicit varying-base cost.
The unrestricted small-budget problem still requires a way to
find such grids cheaply, or a joint bound for the remaining
families with diffuse valuation spectra.

## Bounded-degree composita and class records

For each integer \(D\ge1\) there is a constant \(C_D\) such that
for every extension \(L/K\) with \([L:\mathbb Q]\le D\) and
every \(m\ge1\),
\[
             |\mathrm{Cl}(L)[m]|\ \ge\
                   C_D^{-1}|\mathrm{Cl}(K)[m]|.              \tag{RR7}
\]
Thus a bounded-degree compositum preserves the size of the torsion
subgroup used in a descent up to a bounded factor.

To prove (RR7), pass to the normal closure \(E/\mathbb Q\)
of \(L\), of degree at most \(D!\). The kernel of extension of
classes to \(L\) is contained in the kernel for \(E\). Put
\(G=\operatorname{Gal}(E/K)\), \(n=|G|\) and \(U=\mathcal O_E^\times\).
If a fractional ideal \(\mathfrak a\) of \(K\) becomes \((\alpha)\)
in \(E\), the unit cocycle \(s\mapsto s(\alpha)/\alpha\) defines
a class in \(H^1(G,U)\), independent of the generator and the ideal
representative. If it is a coboundary, \(\alpha\) can be divided
by a unit to lie in \(K^\times\). Comparison of the extended
ideal valuations then makes \(\mathfrak a\) principal in \(K\).
This injects the capitulation kernel into \(H^1(G,U)\).
For a cocycle \(u_s\), the product \(v=\prod_{s\in G}u_s\)
satisfies \(u_s^n=v/s(v)\), so \(n\) kills \(H^1(G,U)\).
The cocycle group is a subgroup of \(U^n\). By Dirichlet's unit
theorem it has rank at most \(n([E:\mathbb Q]-1)\), and its
torsion is a subgroup of a product of \(n\) cyclic groups.
Its quotient by \(n\)-th powers therefore has size at most
\(n^{n[E:\mathbb Q]}\), bounding the capitulation kernel in terms
of \(D\) alone. Restrict extension of ideals to \(\mathrm{Cl}(K)[m]\)
to obtain (RR7).
[Mazur–Ullom, Section 4](https://www.numdam.org/item/10.5802/jtnb.621.pdf#page=12)

There is a separate obstruction to restricting odd-order class
records by requiring a common rational prime. Let \(K_1,\ldots,K_r\)
be distinct quadratic fields, \(L\) their compositum and
\(c_i\in\mathrm{Cl}(K_i)\) any classes of odd order. Infinitely many
rational primes \(\ell\) split completely in \(L\) and admit a prime
\(\mathfrak P\mid\ell\) of \(L\) whose contractions to \(K_i\)
have classes \(c_i\) simultaneously.

Write the class groups additively and let \(j_i\) and \(N_i\) be
extension and norm between \(K_i\) and \(L\). For \(i\ne j\),
take an automorphism of \(L\) fixing \(K_j\) that is nontrivial on
\(K_i\). Quadratic conjugation negates ideal classes, since the
product of an ideal and its conjugate is principal. Hence
\(N_jj_i(c_i)=-N_jj_i(c_i)\), which makes this odd-order class zero.
Also \(N_i j_i=[L:K_i]\), a power of two. Choose an odd-order class
\(c_i'\) with \([L:K_i]c_i'=c_i\). Then \(b=\sum_i j_i(c_i')\)
has \(N_i(b)=c_i\) for each \(i\). This norm calculation is the
odd class group decomposition in the biquadratic case.
[Lemmermeyer, Section 1](https://arxiv.org/pdf/2607.27034#page=2)

The Hilbert class field of \(L\) is Galois over \(\mathbb Q\),
since \(L/\mathbb Q\) is Galois and its Hilbert class field is
preserved by every automorphism over \(\mathbb Q\).
Apply Chebotarev over \(\mathbb Q\) to the conjugacy class of
the Artin image of \(b\). For each resulting rational prime choose
a prime upstairs whose Frobenius is that element. Its restriction
to \(L\) is the identity, so the rational prime splits completely
there, and the restricted prime \(\mathfrak P\) in \(L\) has class
\(b\). The relative norms \(N_i(\mathfrak P)\) are the contractions
with exponent one, which proves simultaneous realization.
[Stevenhagen, Hilbert class fields](https://websites.math.leidenuniv.nl/algebra/Stevenhagen-CFT.pdf#page=6)
· [Lenstra, Chebotarev density theorem, Theorem 3.1](https://websites.math.leidenuniv.nl/algebra/Lenstra-Chebotarev.pdf#page=3)

For a concrete instance, \(41\) splits in both
\(\mathbb Q(\sqrt{-23})\) and \(\mathbb Q(\sqrt{-31})\).
Writing \(\omega_d=(1+\sqrt{-d})/2\), the elements
\(181+64\omega_{23}\) and \(99+80\omega_{31}\) each have norm
\(41^3\) and each is divisible by just one of the two primes over
\(41\) in its field. The norm forms \(a^2+ab+6b^2\) and
\(a^2+ab+8b^2\) do not represent \(41\), as completing the square
reduces the check to \(|b|\le2\). The resulting prime ideals have
exact order three, with their orientations compatible in the
biquadratic compositum. A further restriction on the class records
in (RR5) must therefore use more than coincidence of the rational
prime supports. The split-prime construction imposes no additive
identity or global height constraint on the norm generators.
