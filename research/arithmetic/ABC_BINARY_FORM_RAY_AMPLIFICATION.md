# Radical counts for binary forms in deep ray classes

A large congruence modulus lets high even moments absorb the
ordinary class group factor. The resulting half coefficient on
the radical of a binary-form value holds for any fixed form
degree. The projective congruence can come from a large divisor
of the value itself when the degree is at least three.

## BR1. Primitive inputs with prescribed projective residues

Put \(n=X^L\). Fix \(d\ge2\), \(B,C,\nu>0\), \(r,b\ge0\)
and \(\eta>0\). Let \(P(U,V)\) be a degree-\(d\) integer
binary form of nonzero discriminant and nonzero leading
coefficient \(a=P(1,0)\), with coefficients bounded by
\(n^C\). Let \(M\le X^B\) have \(M\ge X^\nu\) and be prime
to \(a\operatorname{disc}(P)\). Count coprime positive
inputs \(x,y\le X^B\) in a list of at most \(n^b\)
projective residues modulo \(M\). Assume \(P(x,y)\ne0\)
and \(\gcd(P(x,y),M)=1\), and impose
\[
              \operatorname{rad}_{\mathrm{out}}(P(x,y))\le n^r.
\]
Here the radical can omit a prescribed prime support given
by an integer of height at most \(X^B\). Projective equality
means proportional primitive coordinate pairs by a unit
modulo \(M\). For all sufficiently small fixed \(L>0\)
the count is
\[
               O_{d,B,C,\nu,r,b,\eta,L,\varepsilon}
                     \bigl(n^{(r+b)/2+\eta}X^\varepsilon\bigr).
                                                        \tag{BR1}
\]
The threshold on \(L\) and constants are uniform over the
form, modulus, residue list and omitted support.

The same estimate applies with a large integral ideal
\(\mathfrak m\) in place of the modulus, as follows. Use the
splitting field \(E\) of \(P(U,1)\), and any two distinct
integral roots \(s,t\) of \(a^{d-1}P(U/a,1)\). Require
\(X^\nu\le N\mathfrak m\le X^B\), and require that the quotient \((a x-s y)/(a x-t y)\)
is a unit at every prime of \(\mathfrak m\) with a list of
at most \(n^b\) unit residues there. One can change the fixed polynomial height
bound \(B\) to accommodate the field degree. For the
rational version use \(\mathfrak m=M\mathcal O_E\):
the primitive projective class determines the quotient
residue. It also ensures both factors are units under
the coprimality condition on \(P(x,y)\).

## The quotient of two linear factors

In \(E\) the integral roots \(s_j\) give
\[
       \prod_{j=1}^d(a x-s_jy)=a^{d-1}P(x,y).
\]
Choose two roots \(s,t\) and set
\(\gamma=(a x-s y)/(a x-t y)\). The quotient uniquely
recovers \(x/y\), and hence \(x,y\) by positivity and
primitivity. The field degree \(D=[E:\mathbb Q]\) is at
most \(d!\), and its discriminant is \(n^{O_{d,C}(1)}\).
For example the integral roots have all conjugates of size
\(n^{O_{d,C}(1)}\). An integer linear combination of the
roots with coefficients bounded in terms of \(d\) separates
all \(D\) embeddings. This follows by avoiding the finitely
many equality hyperplanes on a fixed integer grid. The
discriminant of its power basis then has the required
bound. Minkowski's theorem and the fixed-degree ideal
divisor bound give
\[
                        |\mathrm{Cl}(E)|\le n^{C_1}
\]
for an effectively bounded \(C_1=C_1(d,C)\) after increasing
it to absorb fixed constants. The bad rational primes
dividing the discriminant of the monic polynomial or
\(a\) are among those dividing \(a\operatorname{disc}(P)\).
Outside that support the field is unramified and the
roots are pairwise distinct modulo every prime of \(E\).

At a good rational prime \(\ell\mid P(x,y)\) we have
\(\ell\nmid y\). At any prime over \(\ell\) at most one of
the \(d\) factors vanishes. Its full valuation equals
\(v_\ell(P(x,y))\). The primes contributing positively and
negatively to the ideal \((\gamma)\) are therefore disjoint
subsets of the primes over \(\ell\). Some \(\ell\) have no
contribution. For each \(\ell\), fix a pair of disjoint
subsets by averaging over their finitely many possibilities.
Some assignment retains \(X^{-o_X(1)}\) of the input pairs,
since each prescribes at most
\(O_{d,B,C}(\log X/\log\log X)\) choices. Restrict to that
assignment. For a rational prime with a nonempty pair
of subsets, let \(\mathfrak j_\ell\) be the product of the
positive primes divided by the product of the negative primes.
On the integers supported on these rational primes the
map \(J(u)=\prod_\ell\mathfrak j_\ell^{v_\ell(u)}\) is
injective. Every ideal used after the modulus and support
records will be prime to \(\mathfrak m\).

Fix the exponents of \((\gamma)\) on the bad and omitted
supports, and on the rational primes below \(\mathfrak m\).
These supports have polynomial height. Each has
\(O(\log X/\log\log X)\) prime slots, and the sum of the
absolute ideal valuation orders is \(O_{d,B,C}(\log X)\),
so all such signed exponent vectors cost \(X^{o_X(1)}\).
The latter assertion follows by stars and bars and the
bounded degree. Fix the contribution of the remaining
primes at which \(v_\ell(P(x,y))>J_0\) as well. Its radical
is at most \(X^{O_{d,B,C}(1)/(J_0+1)}\). Choose \(J_0\)
large at the fixed budget so this costs \(n^{\eta/8}\)
up to height-subpower factors. The remaining ideal is
\[
          (\gamma)=\mathfrak d\prod_{e=1}^{J_0}J(u_e)^e ,
\]
with fixed \(\mathfrak d\) and \(u_e\) the squarefree product
of the active rational primes of exact order \(e\). The
\(u_e\) are disjoint and their product is at most \(n^r\).
On a fixed finite cap mesh write
\(u_e\le n^{\theta_e}\) with \(H_0=\sum_e\theta_e\le r+\eta/4\).
Zero slots can be omitted. All nonzero conjugates of each
linear factor are \(X^{O_{d,B,C}(1)}\), with reciprocals
bounded by the same type of power because their norms
are nonzero integers. Thus the logarithmic coordinates
of \(\gamma\) have size \(O_{d,B,C}(\log X)\).

## Compact ray group and Archimedean cutoffs

Let \(V\) be the real vector space of vectors on the
real places and complex conjugate pairs of \(E\) with
weighted sum zero, the weights being \(1,2\) respectively.
Write \(r_E(u)\) for its normalized log vector with entries
\(\log|\sigma u|-D^{-1}\log|N_{E/\mathbb Q}u|\).
Use the additive compact group
\[
 G=\bigl(I(\mathfrak m)\times V\bigr)/
   \{((u),-r_E(u)):u\equiv1\bmod\mathfrak m\}.
\]
Here \(I(\mathfrak m)\) is the group of fractional ideals
prime to \(\mathfrak m\), and the congruence is in each
local ring at its primes. Write \(g(\mathfrak a)\) and
\(h(v)\) for the images of \((\mathfrak a,0)\) and
\((1,v)\). The kernel of \(h\) is the full Dirichlet
lattice of ray units. Hence \(h(V)\) is a compact torus
and \(G/h(V)\) is the ordinary finite-modulus ray class
group. No signs at real places are imposed. The ray exact
sequence and the local unit root count give, for \(e\)
fixed at the fixed budget,
\[
          |G[e]| \le e^{\dim V}|\mathrm{Cl}(E)[e]|
                         |(\mathcal O_E/\mathfrak m)^\times[e]|
                 \le n^{C_1}X^{o_X(1)}.
\]
At a local factor of residue characteristic \(\ell\)
prime to \(e\), all \(e\)-torsion comes from residue-field
units and is at most \(e\). For each of the finitely many
\(\ell\mid e\), the logarithm identifies \(1+\ell^2\mathcal
O_{E,\mathfrak p}\) with its additive logarithm lattice.
The index of this subgroup in the unit group is at most
\(\ell^{2D}\), and the kernel of multiplication by \(e\)
in the reduced logarithm lattice is at most
\(\ell^{D v_\ell(e)}\), including for ramified primes and
\(\ell=2\). The total is height-subpower at fixed \(e,D\),
since there are \(O_D(\log X/\log\log X)\) local factors.
Quotients of finite abelian groups do not
increase the number of elements killed by \(e\).
[Ray class groups and their exact sequence (Stevenhagen)](https://websites.math.leidenuniv.nl/algebra/Stevenhagen-CFT.pdf)
· [Arakelov class groups (Schoof)](https://library.slmath.org/books/Book44/files/14schoof.pdf)

Set \(W=\nu\log X/(100D)\). Choose one of the boundedly
many log boxes for \(r_E(\gamma)\) with center \(v_0\)
and radius \(W/2\) in the maximum norm on \(V\).
The autocorrelation \(K_0\) of the indicator of
\(h(\{v:\|v\|_\infty\le W\})\), normalized by the Haar
measure of that image, is supported in the projected
box of radius \(2W\). Its Fourier coefficients are
nonnegative and sum to \(K_0(0)=1\). It is bounded below
by a constant depending only on \(D\) on the projected
box of radius \(W/2\). Indeed the overlap contains the
image of the half-size box, and the large box is covered
by a bounded number of translates of that half-size box
in the fixed-dimensional space. This measure comparison
holds also when the projection winds around the torus.
For \(V=0\) take the indicator of the identity.

For each actual quotient residue \(c\) in the list, choose
an anchor \(\gamma_c\) and put
\(z_c=g((\gamma_c))-h(r_E(\gamma_c))\). For any quotient
\(\gamma\) with this residue,
\(g((\gamma))=z_c+h(r_E(\gamma))\) by the ray relation.
The required ideals therefore meet the translated majorant
\(\sum_c K_0(\,\cdot-z_c-h(v_0))\).
With probability weights \(\widehat K_0\) on the dual group,
\[
   \left\|\sum_c\chi(z_c)\right\|_2^2
          =\sum_{c,c'}K_0(z_c-z_{c'})
          \ll |S| X^{o_X(1)},                            \tag{BR2}
\]
where \(S\) is the quotient residue list. To prove the
last inequality, if a term is nonzero then
\(z_c-z_{c'}=h(v)\) for some \(\|v\|_\infty\le2W\).
The definition of \(G\) then gives a unit \(u\) with
residue \(c/c'\) and log vector \(-v\). The number of units
of log norm \(O(\log X)\) is \(O_D((1+\log X)^D)\) uniformly
in the field. For this uniform upper bound partition a
log box into unit boxes. The quotient of any two units
in one box is an algebraic integer of bounded degree with
all conjugates bounded by an absolute constant, and there
are only boundedly many such algebraic integers. This also
counts the roots of unity.

## Amplified moments and the class number price

Let \(F_e(\chi)=\sum_{u_e}\chi(e\,g(J(u_e)))\). Take any
small fixed \(0<\xi\le1/2\) and set
\(m_e=\lceil H_0/(2\xi\theta_e)\rceil\). Then
\(\sum_e1/(2m_e)\le\xi\) and
\(m_e\theta_e\le (1+1/(2\xi))H_0\).
Under the same probability weights \(\widehat K_0\), expand
the \(2m_e\)-th moment. After fixing the first product,
the second products \(u\) must have \(e\,g(J(u))\) in
one given translate of the projected box of radius \(2W\).
Partition them into at most \(|G[e]|\) classes by taking
the bounded log vector of that box divided by \(e\).
This is valid because \(e v\in h(B_{2W})\) implies
\(v\in G[e]+h(B_{2W/e})\).

Two second products \(u_1,u_2\) in the same part satisfy
\(J(u_1/u_2)=(\zeta)\) with \(\zeta\equiv1\bmod\mathfrak m\)
and log norm at most \(4W/e\). The integer \(Q=u_1u_2\)
clears the denominator of this ideal. Also
\(|N\zeta|\le Q^D\), because every signed ideal order is
bounded by the corresponding rational prime order in \(Q\).
It follows that every conjugate of \(Q(\zeta-1)\) has modulus
at most \(2Q^2\exp(4W/e)\). This algebraic integer is in
\(\mathfrak m\). By taking \(L\) sufficiently small in terms
of \(\xi,r,\eta,\nu,D\), we have
\((2Q^2\exp(4W/e))^D<N\mathfrak m\). A nonzero integral element
of \(\mathfrak m\) has absolute norm at least \(N\mathfrak m\).
Therefore \(\zeta=1\) and \(u_1=u_2\) by injectivity.
A fixed \(u\) has height-subpower many ordered factorizations
at fixed \(m_e\). The resulting moment norm obeys
\[
                    \|F_e\|_{2m_e}
                \le n^{\theta_e/2+C_1/(2m_e)}X^{o_X(1)}.
\]

Use Fourier inversion and Hölder with the probability
weights \(\widehat K_0\). The factor from (BR2) is bounded
in \(L^{1/(1-\xi)}\) by its \(L^2\) norm. The \(F_e\)
contribute the displayed norms, with a constant function
of norm one filling any unused Hölder weight. The result
counts ideal tuples at \(n^{(r+b)/2+\eta/4+C_1\xi}X^\varepsilon\)
after the mesh allowance. In each fixed tuple the quotient
has only a height-subpower number of possible generators
by the unit log count. Each generator recovers \(x,y\).
Choose \(\xi\) so that \(C_1\xi<\eta/4\) before choosing \(L\)
and then \(J_0\). The tail price, finite bins and
height-subpower factors prove (BR1). In particular the
class group factor can have arbitrarily small normalized
price before the height limit, without any hypothesis
on the exact valuations retained in the form value.
The moment amplification is the same mechanism as in the
rational residue count. [Relative expansion (RX1)](ABC_RESIDUE_EXPANSION_AND_JOINT_RAY_COUNTS.md#the-relative-expansion-price-for-a-deep-absolute-modulus)

## BR3. Placing the ray modulus on the form value

Suppose \(d\ge3\) and a fixed rational \(M\ge X^\nu\) divides
\(P(x,y)\), with \(M\) prime to \(a\operatorname{disc}(P)\).
The bound (BR1) with \(b=0\) holds with \(r\) measuring
only the radical outside the primes of \(M\) and any
prescribed support. To see this, \(y\) is a unit modulo
\(M\) and Hensel lifting gives at most \(d^{\omega(M)}\)
residue choices for \(x/y\),
since \(P(U,1)\) has only simple roots at its good primes.
Fix a choice at height-subpower cost.

For each prime ideal over a rational prime of \(M\),
at most one of the \(d\) linear factors vanishes. Among
unordered pairs of factors at least the fraction
\((d-2)/d\) are both units there. Assign each prime ideal
the weight of its full power in \(M\mathcal O_E\), measured
by log norm. Some pair therefore has both factors units
on an ideal \(\mathfrak m\mid M\mathcal O_E\) of norm
at least \(M^{D(d-2)/d}\). These choices depend only on
the fixed residue and modulus. The quotient residue on
this ideal is fixed by the projective rational residue.
Record the powers on \(M\)'s rational support at subpower
cost and apply the ideal version of BR1. In this argument
the complete prime-power depths, rather than the support
of \(M\), supply the Archimedean norm cutoff.

## BR4. Two weighted powers and the original three caps

Fix any integer \(k\ge2\). Write two letters as \(U=A x^k\)
and \(V=B y^k\) with \(k\)-free positive coefficients
\(A,B\) and let \(T\) be the third original letter.
Assume \(A,B\le n^C\) for fixed \(C\), as supplied by
bounded coefficient radical caps at fixed \(k\).
Suppose the actual coefficient-pair catalogue costs
\(n^{f_0}\) up to height-subpower factors and the actual
divisor catalogue costs \(n^\mu\) after fixing the
coefficients. With a fixed positive-depth divisor on one
root or a product of divisors on both roots, BR1 counts all three original
letters at exponent at most
\[
                    f_0+\mu+r_T/2+\eta .                 \tag{BR4}
\]
Here \(r_T\) is the cap on the third-letter radical,
with the same option of free prescribed support.
For \(k\ge3\) the same exponent holds with a positive-depth
divisor on \(T\), prime to \(kAB\), using BR3 and charging
only the third radical outside the purchased support.
Take \(P\) to be the corresponding signed binomial. For the root case choose two roots \(s,t\) of
\(Z^k+\epsilon A B^{k-1}\), where \(T=B y^k+\epsilon A x^k\).
At a prime of \(M_x\mid x\), \(B,y\) are units, so
\((By-sx)/(By-tx)\equiv1\bmod M_x\). At a prime of
\(M_y\mid y\), \(A,x\) are units, and \(B/s\) is integral
because its valuation is \(v(B)/k\). Dividing the two factors
by \(sx\) gives quotient residue \(s/t\bmod M_y\), a unit
also at primes dividing \(kAB\). The two moduli are coprime
by primitivity. The ideal version of BR1 therefore applies
to their product with one residue. In both cases the root ratio determines the
entire original primitive triple within the record. This
works also for reducible binomials and for opposite exact
squares. The elementary shared-power slope reconstruction
already gives a cheaper count when the third divisor
exceeds the square of the root height.
[Primitive slope reconstruction (CPW5)](ABC_COPRIME_SHORT_POWER_WINDOWS.md#an-adaptive-common-power-sector-with-a-third-letter-lock)

In particular, for any fixed \(f<1/6\) and fixed \(k\ge2\),
the original small-budget sector in which two letters
have combined \(k\)-defect radical at most \(c^{fL}\)
has a strict coefficient below \(1/3\).
The \(k\)-defect radical is the product of primes whose
positive valuations on the two letters are not multiples
of \(k\). To prove the assertion, if this sector approached
the endpoint, the one-letter quarter bound would force all
three raw caps to \(1/3+o(1)\). The nonatomic masses on
the chosen pair sum to at least \(1/3\). Splitting their
measures by exact divisibility by \(k\), the retained
nondefective masses therefore have positive nonatomic part.
A compact positive band of arbitrarily small positive radical
mass supplies a full root divisor with fixed positive numerical
depth. Its raw price can be made arbitrarily small relative
to the saving margin. Coefficient pairs cost at most
\(f+o(1)\), because the combined defect support determines
the finitely many exponent residues at each prime.
Thus BR4 gives exponent at most \(1/6+f+o(1)\), a contradiction.
All widths and prices in this argument are chosen in the
limiting profile before the inner height limits. The
choice of \(k\) and of the ordered pair can be made from
any fixed finite menu. The same proof bounds the
minimally repaired primitive quotients by applying their
original height and radical restrictions before the repair.
[Balanced masses and the pair nonatomic bound](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)
· [Minimal target repair and its height](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

More generally every balanced saturating cohort must have
actual joint \(k\)-free coefficient price at least \(1/6\)
on every pair of roles and for every fixed \(k\). In fact
each pair has total nonatomic mass at least \(1/3\).
The full powers of the corresponding roots retain a fixed
fraction of the numerical weight of any compact positive
band as \(L\) tends to zero, even if all chosen primes also
divide the \(k\)-free coefficients. Thus some root divisor
of fixed positive depth has arbitrarily small price in each
pair, and (BR4) applies before the class catalogue is summed.
Divisor enumeration conditional on the coefficient record
costs only height-subpower if the divisor is supported
entirely on the coefficients. This conclusion uses actual
coefficient dictionaries, which can be cheaper than the raw
sum of defect prices.

This strengthens the earlier two-square \(2/15\) defect
threshold and leaves no class group price in the leading
coefficient. It still requires an arithmetic extraction of
a cheap exact coefficient record. Weak profiles alone
allow every valuation to be \(1\) modulo every order in
any fixed menu. For the original binary form \(UV(U+V)\)
with a lock on one factor, the other two factors at balanced
caps have combined radical \(n^{2/3+o(1)}\). The half count
at that input price reaches \(1/3\). For the leading coefficient
condition one can use \(x=U+V,\ y=U+2V\), which gives the form
\(x(2x-y)(y-x)\) with leading coefficient \(-2\). A universal saving
needs a discount on this input or a different use of the
simultaneous valuations.


## BR5. Several coprime factor groups

Factor \(P\) over \(\mathbb Q\) into primitive integer
factors, after taking out its content. Each nonlinear
factor is an eligible group, and disjoint pairs of
rational linear factors are eligible groups. Let \(j\)
be the number of groups, with unused linear factors
allowed. Away from the content and the resultants
of the factors, their values on a primitive input
have disjoint supports. These are fixed polynomial-height
bad supports. Use two of the full form's integral roots corresponding
to the chosen group. Away from the fixed bad support,
the active rational primes for their quotient all divide
the value of that group. This is the only place in the
proof of BR1 where the radical cap is used.
For the count with a fixed projective
residue list and unit form value modulo \(M\), apply BR1
to the group of least outside radical. It costs at
most \(r/j\) in the cap, giving exponent
\[
                            r/(2j)+b/2+\eta.
\]
Each group still determines the input from its quotient.
The loss for the choice of group is bounded in terms of
\(d\).

For a modulus on the form value, a group of degree at
least three is eligible at any positive limiting share
of \(M\)'s log depth supported on its own value by BR3.
Any group is eligible if a fixed positive proportion of
the modulus depth lies on other factors, by BR1 applied
at those primes. For \(j\ge2\), at most one group can
have more than half the \(M\)-depth on its own support.
Thus one always has a bound \(n^{r/(2(j-1))+\eta}X^\varepsilon\)
with \(r\) measured outside the fixed \(M\)-support.
This uses only good primes, as in BR3. If all groups are
eligible with a uniform positive depth, the denominator
can be \(2j\). The projected residue of the input remains
fixed after the Hensel decomposition on \(M\).

## Comparison with binary-form value counts

The primitive input count here uses a bound on the
product of distinct primes and a positive-height
projective or ray condition. Balog, Blomer, Dartyge
and Tenenbaum study the largest prime factor of a binary
form value in their friability counts. Fouvry and
Waldschmidt count represented integers for families of
binomial forms. These theorems give benchmarks for friability and represented
values. BR1 uses the product of the distinct prime factors and
the numerical depth of the imposed ray condition. The analytic predecessor used here is the rational
residue moment and its quadratic ray extension. The
value-divisor argument above uses the number of good
root pairs in each local factor of a general splitting field.
[Friable values of binary forms, Theorems 1–2](https://ems.press/content/serial-article-files/43293)
· [Fouvry–Waldschmidt, Theorems 1.1 and 1.4](https://webusers.imj-prg.fr/~michel.waldschmidt/articles/pdf/BinomialFormsAA-EF-MW.pdf)
· [Projective norm-form ray moments (PR1–PR5)](ABC_PROJECTIVE_RAY_CLASS_NORM_MOMENTS.md)
· [Quadratic ray moments (WR1–WR5)](ABC_WEIGHTED_SQUARE_RAY_CLASS_MOMENTS.md)
