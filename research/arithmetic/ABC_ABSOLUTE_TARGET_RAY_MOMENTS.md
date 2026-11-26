# Absolute target residues and half-price norm counts

A deep divisor of a quadratic norm omits the prime ideal where the
norm vanishes. The remaining prime does not differentiate a rational
scalar in a projective ray group. An absolute residue of one input
restores an ordinary ray condition there. A list of such residues
enters the count at the square root of its size. The same argument
counts algebraic integers over any field of fixed degree.

## Algebraic integers with a ray residue dictionary

Fix a degree \(d\ge1\), constants \(B,\nu>0\), and \(g,\ell\ge0\).
For every \(\eta>0\) and sufficiently small fixed \(L>0\), let \(K\)
be a number field of degree \(d\) with
\[
                    |\operatorname{Disc}(K)|\le X^{gL}.
\]
The field can vary with \(X\). Fix an integral ideal \(\mathfrak m\)
with \(X^\nu\le N\mathfrak m\le X^B\), a positive integer \(P\le
X^B\), and a set \(\mathcal U\subseteq
(\mathcal O_K/\mathfrak m)^\times\) of size \(R\ge1\). Count the
nonzero \(\alpha\in\mathcal O_K\) satisfying
\[
 |\sigma(\alpha)|\le X^B\quad(\sigma:K\hookrightarrow\mathbb C),
 \qquad \alpha\bmod\mathfrak m\in\mathcal U,\qquad
 \prod_{\substack{\mathfrak p\mid(\alpha)\\
                    \mathfrak p\cap\mathbb Z=(p),\ p\nmid P}}
             N\mathfrak p\le X^{\ell L}.
\]
The count is
\[
                   \ll_{d,B,\nu,g,\ell,\eta,L,\varepsilon}
                       R^{1/2}X^{(\ell/2+\eta)L+\varepsilon}.
                                                               \tag{AT1}
\]
The small-budget threshold is independent of the field, modulus
and residue set. One can replace the discriminant bound by
\(\ll_d X^{gL}\) with the same conclusion.

### Ideals and the compact ray group

Write \(V\) for the real vector space of embedding logs with
conjugate coordinates equal and coordinate sum zero. Write
\(\mathcal L(a)_\sigma=\log|\sigma(a)|-d^{-1}\log|N a|\),
counting each complex coordinate twice in the sum. Use
\[
 G=(I(\mathfrak m)\times V)/
   \{((a),-\mathcal L(a)):a\equiv1\bmod\mathfrak m\},
\]
where \(I(\mathfrak m)\) is the discrete group of fractional ideals
prime to \(\mathfrak m\). The congruence is local at every prime
dividing \(\mathfrak m\), and there is no positivity condition.
Write \(g(\mathfrak a)=[(\mathfrak a,0)]\) and \(h(t)=[(1,t)]\).
The denominator subgroup is closed because its intersection with
each discrete ideal fibre is empty or a translate of the log
lattice of units congruent to one. These units have finite index in
the full unit group. The image \(h(V)\) is therefore a torus and
the quotient by it is the ordinary ray class group with finite
modulus \(\mathfrak m\). This is the Arakelov ray class group in
centered logarithmic coordinates, using the negative of the
archimedean logarithm convention in Bartel–Pagano. Their definition
applies to arbitrary finite ideal moduli. [Bartel–Pagano, §2.1](https://arxiv.org/html/2509.20185v1#S2.SS1) · [Projective variant (PR3)](ABC_PROJECTIVE_RAY_CLASS_NORM_MOMENTS.md#the-compact-group-and-the-diagonal-moment)

For every fixed positive integer \(e\),
\[
                           |G[e]|\le |\mathrm{Cl}(K)[e]|X^{o_X(1)}.
                                                               \tag{AT2}
\]
Indeed the torus contributes at most \(e^d\), and the finite ray
group contributes at most the product of the \(e\)-torsion sizes of
the ordinary class group and the residue unit group. The latter
bound follows from the ray exact sequence and the fact that the
size of the \(e\)-torsion in a finite abelian group equals the size
of its quotient by \(e\)-multiples. For a local residue field of
characteristic \(p\nmid e\), its multiplicative group contributes
at most \(e\) roots and the principal unit group has no \(e\)-torsion.
For \(p\mid e\), there is also a bound depending only on \(d,e\),
uniform in the modulus exponent. In the normalized local
valuation \(v_{\mathfrak p}\), on the subgroup
\(1+\mathfrak p^{d+1}\) binomial expansion gives
\[
             v_{\mathfrak p}(a^e-1)
                    =v_{\mathfrak p}(a-1)+v_{\mathfrak p}(e)
\]
for \(a\ne1\). Apply the expansion first to the prime-to-\(p\)
part of \(e\), then successively to multiplication by \(p\).
The initial depth exceeds \(v_{\mathfrak p}(p)/(p-1)\), which
makes the first term the unique term of minimum valuation at
each step. The quotient before depth \(d+1\) has size at most
\(p^{d(d+1)}\) and \(v_{\mathfrak p}(e)\le d v_p(e)\). These
bounds also cover modulus depths below \(d+1\). Multiplying over
the primes gives a height subpower since each \(e\) is fixed in
the height limit and \(N\mathfrak m\le X^B\).

There are \(X^{o_X(1)}\) possibilities for the full \(P\)-supported
part of \((\alpha)\), uniformly in \(P\) and \(K\). To see this,
enumerate its norm on the fixed rational support and then its ideal
factorization. For each integer \(n\), there are at most
\(\tau(n)^d\) integral ideals of norm \(n\). The norm and all
conjugates of \(\alpha\) have polynomial height uniformly, with
\(|\sigma(\alpha)|\ge X^{-(d-1)B}\) from its nonzero integral norm.
Also fix the full remaining factor supported on valuations above
an integer cutoff \(J\). Its squarefree ideal norm is at most
\(X^{dB/(J+1)}\). Choose \(J\) at fixed \(L\) so this exponent is at
most \(\eta L/8\). Each such radical ideal gives a rational support
of polynomial height, so its complete supported part has only
height-subpower many possibilities by the same enumeration.
[Uniform fixed-support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

We have
\[
         (\alpha)=\mathfrak d\prod_{e=1}^{J}\mathfrak a_e^e,\qquad
         N\mathfrak a_e\le X^{L\theta_e},\qquad
         H_0:=\sum_e\theta_e\le\ell+\eta/4 .
\]
Here \(\mathfrak d\) is fixed and all ideals are prime to
\(\mathfrak m\). The \(\mathfrak a_e\) are squarefree ideals.
The caps come from a mesh of total width at most \(\eta/4\), fixed
at the budget before the height limit. The number of ideals of
norm at most \(X^{L\theta_e}\) is at most
\(X^{L\theta_e+o_X(1)}\) by the bound on ideals of any given norm.
Slots of zero cap contribute only the unit ideal.

### The two diagonal comparisons

For \(r\in\mathcal U\) choose an integral lift \(\beta_r\) and set
\[
                         c_r=[((\beta_r),-\mathcal L(\beta_r))].
\]
This is independent of the lift since the quotient of two lifts
is congruent to one. If \(\alpha\equiv r\bmod\mathfrak m\), then
\(g((\alpha))=c_r+h(\mathcal L(\alpha))\).
Since \(\mathcal L(\alpha)=O_{d,B}(\log X)\), cover its range by
a bounded number of translates of the cube \(C_{W/2}\) in \(V\)
of sup norm radius \(W/2\), with \(W=\nu\log X/(32d)\).
For each translate, \(\sum_e e\,g(\mathfrak a_e)\) must lie in
one of the corresponding translates of \(c_r+h(C_{W/2})\).

Let \(K_0\) be the autocorrelation of the indicator of \(h(C_W)\)
divided by its Haar measure in \(G\). It is supported on
\(h(C_{2W})\), has Fourier coefficients \(\widehat K_0\ge0\)
summing to one, and is bounded below by a positive constant
depending on \(d\) on \(h(C_{W/2})\). Indeed the intersection
when translating by a point of \(h(C_{W/2})\) contains
\(h(C_{W/2})\), and \(C_W\) can be covered by a bounded number of
translates of \(C_{W/2}\). These facts apply to the zero-dimensional
torus as well. The Fourier series is absolutely convergent by
Parseval for the indicator.

With \(S(\chi)=\sum_{r\in\mathcal U}\chi(c_r)\), we have
\[
                         \sum_\chi\widehat K_0(\chi)|S(\chi)|^2
                                  \le R X^{o_X(1)} .          \tag{AT3}
\]
For a pair contributing on the left we have \(c_r-c_s=h(t)\)
with \(\|t\|_\infty\le2W\). By the definition of \(G\), some global
unit \(\epsilon\) then satisfies
\(\epsilon\bmod\mathfrak m=s/r\) and
\(\mathcal L(\epsilon)=t\). For each \(r\) the number of possible
\(s\) is at most the number of units in that log box, which is
\(O_d((1+\log X)^d)\) uniformly in \(K\). For this uniformity,
a unit all of whose conjugates are bounded by two belongs to a
fixed finite set of algebraic numbers of degree at most \(d\)
by the bounds on minimal polynomial coefficients. The nonzero
log vectors therefore have uniform positive separation, by
applying this observation to quotients of units. The zero log
vectors give a uniformly bounded group by the same argument.
Packing gives the claimed unit bound and (AT3).

For each moving ideal slot \(e\), let \(\mathcal A_e\) be its
dictionary and put \(F_e(\chi)=\sum_{\mathfrak a\in\mathcal A_e}
\chi(e\,g(\mathfrak a))\). Choose a constant \(T\ge1\) in terms
of \(g,\eta\), and set \(m_e=\lceil TH_0/\theta_e\rceil\) when
\(\theta_e>0\). Then \(\sum_e 1/m_e\le1/T\) and
\(m_e\theta_e\le(T+1)H_0\). As in (PR5), a \(2m_e\)-th moment
groups products \(\mathfrak B\) corresponding to a fixed first
product into at most \(|G[e]|\) classes. In each class two such
products have
\[
 \mathfrak B_1/\mathfrak B_2=(\gamma),\qquad
 \gamma\equiv1\bmod\mathfrak m,\qquad
 \|\mathcal L(\gamma)\|_\infty\le4W/e .
\]
In fact the equation \(e\,g(\mathfrak B)-e\,g(\mathfrak A)=h(t)\)
with \(\|t\|_\infty\le2W\) gives
\(g(\mathfrak B)-g(\mathfrak A)\in G[e]+h(C_{2W/e})\),
so differences in one class give the displayed condition.

Put \(D_2=N\mathfrak B_2\). The integer \(D_2\) belongs to
\(\mathfrak B_2\), so \(D_2(\gamma-1)\) is an algebraic integer.
The ordinary ray congruence places it in \(\mathfrak m\).
If it is nonzero, its absolute norm is at least \(N\mathfrak m\).
On the other hand
\[
 \big|\sigma\big(D_2(\gamma-1)\big)\big|
      \le 2X^{(T+1)H_0L+\nu/(8d)}.
\]
This uses
\(|N\gamma|=N\mathfrak B_1/N\mathfrak B_2\) and the bound
on \(\mathcal L(\gamma)\). At sufficiently small fixed \(L\)
the norm upper bound is smaller than \(X^\nu\) for large \(X\).
Consequently \(\gamma=1\) and \(\mathfrak B_1=\mathfrak B_2\).
There is no comparison of two distinct field embeddings here.

Each product ideal of polynomial norm has only height-subpower
many ordered factorizations of length \(m_e\) by the divisor bound,
since \(m_e\) is fixed at the height limit. Thus
\[
 \left(\sum_\chi\widehat K_0|F_e|^{2m_e}\right)^{1/(2m_e)}
               \le X^{L\theta_e/2+o_X(1)}
                               |\mathrm{Cl}(K)[e]|^{1/(2m_e)} .
\]
Sum the translated tests over \(r\), then use Fourier inversion
and Hölder with exponent \(2\) for \(S\), exponents \(2m_e\)
for the other moving slots and the constant one for any remaining
exponent. This counts all admissible ideal tuples at
\(R^{1/2}X^{H_0L/2+o_X(1)}|\mathrm{Cl}(K)|^{1/(2T)}\).
Each fixed ideal has at most height-subpower many generators
in the allowed embedding box, by the same unit packing argument.
Minkowski's bound and ideal enumeration give
\(|\mathrm{Cl}(K)|\ll_{d,\xi}|\operatorname{Disc}(K)|^{1/2+\xi}\).
Choose \(T\) large enough that the class factor costs at most
\(X^{\eta L/4}\), then take \(L\) small enough for the norm
comparison. If all slots are trivial, the fixed-ideal generator
bound suffices. Including the tail, mesh and enumeration costs
proves (AT1).

The same choice of moment orders applies to the projective
comparisons in (PR1) and (PR1'). Their bound on the product norms
changes from \(X^{2H_0L}\) to \(X^{(T+1)H_0L}\), and the diagonal
step still holds at sufficiently small fixed \(L\) after choosing
\(T\). Thus the \(g/4\) term in (PR2) can be absorbed in any
positive coefficient allowance when \(g\) is fixed. In (PR7) this
gives the combined defect threshold \(f<1/6\) for every fixed
\(q\ge2\). Indeed the small input lock can still be chosen from
the positive nonatomic mass since \(f<1/3\); the coefficient
dictionary costs \(f\), and the third-letter cap costs \(1/6+o(1)\).
The split opposite-square case still uses the cyclotomic count.
The order \(q\), all field-degree bounds and the chosen moment
multiplier \(T\) are fixed before the small-budget threshold.
[Projective count and weighted powers (PR2, PR6, PR7)](ABC_PROJECTIVE_RAY_CLASS_NORM_MOMENTS.md#weighted-powers-and-the-balanced-endpoint)

## Divisors of binary norms with an absolute input residue

Let \(F(U,V)\) be an irreducible homogeneous integral form of
degree \(d\ge2\) with polynomial coefficient and input heights as
in (PR1), and let \(K\) be its root field. Fix a divisor record
\(M\mid F(u,v)\) with \(M\ge X^\nu\), coprime to
\(2F(1,0)F(0,1)\operatorname{Disc}(F)\). Assume
\(|\operatorname{Disc}(K)|\ll X^{gL}\) and cap the radical of
the form value outside \(PM\) by \(X^{\ell L}\), where \(P\)
is prescribed and of polynomial height. Suppose \(u\bmod M\)
belongs to a set of size \(R\). The number of primitive inputs
is
\[
                          \ll R^{1/2}X^{(\ell/2+\eta)L+\varepsilon}.
                                                               \tag{AT4}
\]
We may use a divisor \(M_*\mid M\) for the absolute residue,
provided \(M_*\ge X^\nu\). Only \(M_*\) enters the ray modulus,
and all primes of \(M\) are still excluded from the remaining cap.

Write \(a_0=F(1,0)\), let \(\theta\) be a root of the monic
polynomial \(a_0^{d-1}F(Z/a_0,1)\), and put
\(\alpha=a_0 u-\theta v\). Then
\(N\alpha=a_0^{d-1}F(u,v)\), and \(\alpha\) determines the
primitive input, since \(\theta\notin\mathbb Q\). At any good
prime \(p\mid M_*\), \(u,v\) are units. The residue of \(a_0u/v\)
modulo the whole \(p^k\parallel M_*\) is one of the simple
unit roots of the monic polynomial by Hensel's lemma. Fix the
root list at height-subpower cost. Let \(\mathfrak p_p\) be its
corresponding unramified degree-one prime of \(K\). For
\[
              \mathfrak m=\prod_{p^k\parallel M_*}
                     \prod_{\substack{\mathfrak p\mid p\\
                                      \mathfrak p\ne\mathfrak p_p}}
                            \mathfrak p^k
\]
we have \(N\mathfrak m=M_*^{d-1}\ge X^\nu\). For any specified
unit \(u\bmod p^k\), the fixed root determines \(v\bmod p^k\)
as well. Hence the residue of \(\alpha\bmod\mathfrak m\) is a
unit and ranges over at most \(R\) classes for the whole list.
Place the primes of \(M\), the optional support \(P\), and the
bad coefficient and polynomial discriminant primes in the
fixed support of (AT1). At all remaining primes dividing
\(F(u,v)\), the ideal of \(\alpha\) uses just the one degree-one
prime specified by \(a_0u/v\bmod p\), so the ideal radical cap
equals the remaining rational radical cap. This proves (AT4),
including the quadratic case. A variable divisor dictionary
can be summed at its actual cardinality.

Without a scalar residue cap the quadratic conclusion with
\(R=1\) fails even for \(F(U,V)=U^2+V^2\) and one fixed growing
divisor. Multiplication by \((2+i)^k\) gives many primitive
representations of norms divisible by \(5^k\) and lets the residue
on the inactive conjugate prime vary. This is the quadratic
obstruction to the projective target-divisor comparison.
[Quadratic counterexample for the projective condition](ABC_PROJECTIVE_RAY_CLASS_NORM_MOMENTS.md#the-target-divisor-statement-fails-in-degree-two)

## Consequences for abc residue blocks and minimal target repair

Suppose \(a=A u^2,\ b=B v^2\) for coprime positive squarefree
kernels \(A,B\), with \(a+b=c\) primitive. Fix \(A,B\) and a
good divisor \(M\mid c\) of height at least \(X^\nu\), coprime to
\(2AB\). The form \(A U^2+B V^2\) is irreducible. Its field is
\(\mathbb Q(\sqrt{-AB})\), with discriminant of absolute value
at most \(4AB\). An actual dictionary for the absolute residues
\(a\bmod M\) of size \(R\) gives one for \(u\bmod M\) of size at
most \(2^{\omega(M)}R\), since all roots are units and Hensel
lifting gives at most two choices at each prime. Thus the
bound (AT4) holds with the same square root dictionary cost,
up to a height subpower. Fixing the kernels also fixes the
weighted form, with no need to prescribe the complete \(c\).
This applies to a genuine deep residue block in (RB1).
[Deep residue blocks (RB1)](ABC_DEEP_RESIDUE_BLOCK_INVERSE.md#the-marked-deep-residue-block-theorem)

For example, in the balanced graph of that result put
\(c_0=1/6,\ n=X^L\), and cap \(AB\le n^f\) for a fixed \(f\).
On each surviving prime \(p\ge n^\beta\) one has
\(M=P_p\), at most \(n^{\beta+\kappa+\eta_0}\) summand residues,
and a target radical outside \(p\) at most \(n^{2c_0+\kappa-\beta}\).
Here \(\eta_0\) denotes the tolerance of (RB1).
For fixed kernels with a row at \(p\), primitivity makes \(p\)
prime to \(AB\), and \(p>2\) at large height. The cancellation
of \(\beta\) between the residue cost and the removed target
radical gives the bound \(n^{c_0+\kappa+\eta_0/2+\eta}X^\varepsilon\)
per prime and kernel pair, allowing all its surviving residue
blocks at once. Summing over at most \(n^{\beta+\gamma}\) primes
and \(n^{f+o_X(1)}\) kernel pairs gives coefficient
\(c_0+f+\beta+\gamma+\kappa+\eta_0/2+\eta\), in addition to the
exceptional term in (RB1). This is a strict estimate on the
marked graph when \(f+\beta<c_0\) by a fixed margin. Its new
input is the absolute target-divisor count in a packet where
no separate divisor on a summand root was supplied.

For minimal target repair let \(h\) be a sufficiently large
fixed repair order, let \(H_h(c)=\prod_{p\mid c}
p^{\max(v_p(c),h)}\le Y\), and set \(L=1/h,\ X=Y\).
In the same weighted-square sector assume \(AB\le Y^{f/h}\).
Fix \(f,\mu,\sigma\ge0\) and \(0\le\beta\le1\).
Suppose each row has a good target divisor \(M\ge Y^\nu\)
from a dictionary of size at most \(Y^{\mu/h}\), with
\(\operatorname{rad}(M)\ge Y^{\beta/h}\). For each fixed
\(M\) and kernel pair, require an absolute residue dictionary
for \(a\bmod M\) of size at most \(Y^{\sigma/h}\). The fixed
target repair height implies \(\operatorname{rad}(c)\le Y^{1/h}\),
so (AT4) gives
\[
              \ll Y^{\{f+\mu+(1-\beta+\sigma)/2+\eta\}/h
                           +\varepsilon} .                  \tag{AT5}
\]
This counts primitive quotients. The minimal multiplier is
determined by \(c\) and contributes no extra choices.
The strict inequality \(f+\mu+(1-\beta+\sigma)/2<1\)
gives a genuine improvement over \(Y^{1/h}\) in this sector.

The norm argument uses absolute residues and an irreducible
bounded-degree carrier. For arbitrary diffuse profiles the
balanced packet theorem supplies the former. The remaining
extraction problem is to represent enough rows by such carriers
at a coefficient dictionary cost that leaves the indicated
margin. The square-kernel encoding alone need not have that cost.

## Reducible binary forms and composite target packets

**AT6.** Let \(F\) be a nonzero separable homogeneous integral
binary form of fixed degree, with polynomial coefficient height
and with nonzero first and last coefficients. Suppose \(F\)
has one fixed irreducible factor \(F_0\) of degree \(d\ge2\)
whose root field \(K\) has \(|\operatorname{Disc}(K)|\ll X^{gL}\).
Let \(M\mid F(u,v)\) be a fixed divisor of polynomial height of
a nonzero value at a primitive input of polynomial height.
It suffices that \(M\) have a divisor \(M_*\ge X^\nu\)
coprime to \(2F(1,0)F(0,1)\operatorname{Disc}(F)\).
If the radical of \(F(u,v)\) outside \(PM\) is at most
\(X^{\ell L}\), and the \(u\)-residues modulo \(M_*\) run through
at most \(R\) classes, the number of inputs satisfies (AT4)
with these \(\ell,R\). Here \(P\) is a prescribed integer of
polynomial height. All fixed parameters have the same order of
limits as in (AT1).

For each \(p^j\parallel M_*\), the reduction of \(F\) is squarefree
and has both leading and trailing coefficients nonzero.
Primitivity makes \(u,v\) units at \(p\). The ratio \(u/v\)
lies in one of a bounded number of residue classes of roots
of \(F(T,1)\bmod p^j\). Record these classes for every \(p\)
at height-subpower cost. With \(a_0=F_0(1,0)\) and \(\theta\)
as in (AT4) for \(F_0\), set \(\alpha=a_0u-\theta v\).
If the recorded root reduces to a root of \(F_0(T,1)\),
omit only its corresponding degree-one prime of \(K\) over
\(p\). If it reduces to a root of a different factor, keep all
primes of \(K\) over \(p\). Use the exponent \(j\) on every
kept prime. In both cases the resulting contribution to
\(N\mathfrak m\) is at least \(p^{j(d-1)}\). The chosen
root and the specified absolute residue of \(u\) determine
\(v\) and the unit residue of \(\alpha\) at all retained
primes. For primes outside \(PM\) and the separability and
coefficient pool, the ideal of \(\alpha\) contributes no
more radical norm than \(F_0(u,v)\), which divides \(F(u,v)\).
Include the extra polynomial-height pool in (AT1).
The resulting count is (AT4), since the irrationality
of \(\theta\) recovers \(u,v\) from \(\alpha\).
In particular the divisor \(M_*\) may lie entirely on
linear factors of \(F(u,v)\). It still supplies an ordinary
ray modulus for the fixed nonlinear factor.

**AT7.** In the composite packet graph of (AP1)–(AP2), take
\(a=A u^q\), \(b=B v^q\) with canonical positive \(q\)-free
coefficients and fixed \(q\ge2\), and suppose
\(\operatorname{rad}(AB)\le n^f\).
After the exceptional loss in (AP2), the marked count in a
cohort of empirical packet price \(s\) is at most
\[
                   O\left(n^{c_0+f+(3s-\beta)/2+O(\eta)}\right)
 \tag{AT7}
\]
with an arbitrarily small additional normalized height allowance,
where \(\kappa\le\eta/100\) as in (AP2).
An upper bound \(\bar s\) for the prices of near-full cohorts
therefore gives a strict saving if
\[
                     f+\tfrac32\bar s-\tfrac12\beta<c_0
\]
with fixed margin. This application permits composite packets
all of whose primes have subpower size, with no restriction
on the allocation of the packet among factors of the
binomial.

Indeed, the same-side form \(F(U,V)=A U^q+B V^q\) has a
nonlinear irreducible factor, since \(F(T,1)\) has a nonreal
root. Remove the primes of \(2q\) from the ray modulus. The
loss is bounded at fixed \(L\) by the compact-order hypothesis
\(v_p(M)\le S/L\), and the ray modulus still has numerical
height at least \(X^{w/2}\) for sufficiently large height.
Primitivity gives \(p\nmid AB\) on the target. The other
good-prime conditions are automatic off \(2qAB\), since the
binomial has simple, nonzero roots there. The field
discriminant of the fixed factor is bounded by
\(O_q(n^{(q-1)f})\), by tame ramification outside \(q\)
and bounded valuations at the fixed primes of \(q\) for
canonical \(q\)-free coefficients. This is the same
discriminant bound as in (PR6). The choice of factor
can be made once for each pair \(A,B\).

The residue list for the actual \(a\bmod M\) in (AP2)
has size at most \(n^{s+\kappa+\eta}\).
At fixed \(A,B,M\), this gives at most
\(q^{\omega(M)}n^{s+\kappa+\eta}\) unit residue classes
for \(u\) on the stripped modulus, by simple-root
lifting of \(u^q\) away from \(q\).
The residual target radical outside all of \(M\)
has cap \(n^{2c_0+\kappa-\beta}\).
Use (AT6), then sum over at most \(n^{s+\kappa}\) packets
and \(n^{f+o_X(1)}\) coefficient pairs. The exponents are
the sum of the packet price, the coefficient price,
half the residue-list price, and half the residual
target radical price. This gives (AT7).

If the coefficient representation depends on the target
repair problem, the same proof fixes a good target divisor
dictionary before counting. In (AT5) one may replace
squares by any fixed \(q\ge2\), count the canonical
\(q\)-free defect radical at price \(f\), and use a divisor
and residue dictionary on the original letters. The
same exponent holds when a fixed positive numerical part of
each target divisor remains after excluding \(2qAB\).
The estimate still concerns least-repaired primitive
quotients.
[Composite packets and their empirical price (AP1–AP2)](ABC_AGGREGATE_DEEP_PACKET_INVERSE.md#composite-packet-inverse)
· [Weighted binomial field discriminants (PR6)](ABC_PROJECTIVE_RAY_CLASS_NORM_MOMENTS.md#weighted-powers-and-the-balanced-endpoint)
