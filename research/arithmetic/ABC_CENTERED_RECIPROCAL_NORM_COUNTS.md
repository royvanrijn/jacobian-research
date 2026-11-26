# Joint radical bounds for linear forms and an Eisenstein norm

For a positive primitive triple \(a+b=c\) with \(a<b\), put
\[
 D=abc,\qquad C=(b-a)(2a+b)(a+2b),\qquad Q=a^2+ab+b^2.
\]
For every \(\eta>0\) and every sufficiently small fixed \(\lambda>0\),
\[
 \#\{(a,b,c):c\le X,\ \operatorname{rad}(DCQ)\le c^\lambda\}
 \ll_{\eta,\lambda,\varepsilon}
 X^{(1/8+\eta)\lambda+\varepsilon}.
 \tag{CN1}
\]
The sixth linear factor supplies a useful additional constraint even
when the norm is counted by a separate cap. Specifically, impose
\[
 \operatorname{rad}(D)\le c^{B_D\lambda},\qquad
 \operatorname{rad}(C)\le c^{B_C\lambda},\qquad
 \operatorname{rad}(Q)\le c^{B_Q\lambda}
\]
for fixed nonnegative \(B_D,B_C,B_Q\). With the same quantifiers, the
coefficient in the count is at most
\[
 \min\!\left\{
 \frac{B_D}{3},\frac{B_C}{3},B_Q,\
 \max\!\left(\frac{B_D+B_C}{8},\frac{B_Q}{2}\right)
 \right\}+\eta .
 \tag{CN2}
\]
In particular, two unit radical caps on \(D,C\) and any
fixed \(B_Q<2/3\) give a strict improvement below \(1/3\).

## The centered reciprocal companion

For \(z=a/c\), set \(F(z)=z-z^{-1}+(1-z)^{-1}\) and
\(H(z)=2F(z)-1\). Direct substitution gives
\[
 H(z)=-C/D,\qquad
 R(z)=\frac{1-2z}{2-z},\qquad H(R(z))=27/H(z).
 \tag{CN3}
\]
Indeed the positive primitive triple of \(R(z)\) is obtained by
dividing \((b-a,2a+b,a+2b)\) by their common gcd, which divides
three. Applying this linear transformation twice multiplies the
two summands by three. In particular a cap on the primitive
companion's radical is equivalent to a cap on \(C\) with a
bounded factor. The Eisenstein identity is
\[
 C^2+27D^2=4Q^3.
\]
For example, if \(\omega^2+\omega+1=0\) and
\(\alpha=a-b\omega\), then \(\alpha^3+\bar\alpha^3=-C\) for
\(\omega=(-1+\sqrt{-3})/2\), and
\(\alpha^3-\bar\alpha^3=-3(\omega-\bar\omega)D\).
Also \(N_{\mathbb Q(\omega)/\mathbb Q}(\alpha)=Q\).

For each pair of the six linear forms, the absolute determinant
of the two coefficient vectors belongs to \(\{1,2,3\}\).
On their respective primitive kernel vectors, the values of
\(Q\) are \(1,1,1,3,3,3\), in the order
\(a,b,a+b,b-a,2a+b,a+2b\). Thus no prime outside \(\{2,3\}\)
can divide two different values among the six forms and \(Q\)
at a primitive input. The sum of their logarithmic radicals
differs from the logarithm of the joint radical by \(O(1)\).

The identity (CN3) is a concrete branch-preserving exception
for the fractional-linear recovery theorem. Its sole fixed
point in the sorted interval is \(2-\sqrt3\), so it pairs
distinct rational arguments. The arithmetic constraints in
(CN1) and (CN2) concern the actual three numerator factors
and the norm.
[Fractional-linear recovery and its branch-pair exceptions](ABC_FRACTIONAL_LINEAR_PROJECTIVE_RECOVERY.md#recovery-and-the-branch-pair-exceptions)

## A norm count with a full input divisor

Let \(L\) be a fixed primitive integral linear form in two
variables. For fixed \(\nu>0\) and \(q,\mu\ge0\), count primitive
\((a,b)\) of height at most \(X\) with \(\operatorname{rad}(Q)\le
X^{q\lambda}\) and some \(M\mid L(a,b)\) satisfying
\(M\ge X^\nu\) and \(\operatorname{rad}(M)\le X^{\mu\lambda}\).
Require \(L(a,b)\ne0\). For every \(\xi>0\) and sufficiently small
fixed \(\lambda\), the count is
\[
 \ll_{L,\nu,q,\mu,\xi,\lambda,\varepsilon}
 X^{(\mu+q/2+\xi)\lambda+\varepsilon}.
 \tag{CN4}
\]
For fixed \(M\) omit the \(\mu\) term. This is the imaginary
quadratic input case of the projective ray moment argument,
with an extension to divisors at fixed bad primes. These
primes matter because a large power of \(2\) or \(3\) can
carry numerical height at zero limiting radical cost.
[Projective ray class moment bound](ABC_PROJECTIVE_RAY_CLASS_NORM_MOMENTS.md#counting-primitive-inputs-of-an-irreducible-binary-form)

Write \(O=\mathbb Z[\omega]\), the norm-Euclidean ring of
Eisenstein integers. Extend \(L\) to a unimodular pair of
integral coordinates \(L,W\). Then
\(\alpha=a-b\omega=\gamma W+\delta L\) with fixed
\(\gamma,\delta\in O\), \(\gamma\ne0\). If \(p\mid M\), primitivity
forces \(W\) to be a \(p\)-adic unit. At primes \(p\nmid N(\gamma)\)
one has \(\alpha\equiv\gamma W\bmod p^k O\) for \(p^k\parallel M\),
and \(\alpha\) is a unit there. At the finitely many other
primes choose integers \(c_p\) with
\(c_p e(\mathfrak p/p)>v_{\mathfrak p}(\gamma)\) for every
\(\mathfrak p\mid p\). Replace \(p^k\) in the modulus by
\(p^{\max(k-c_p,0)}\). At every retained prime above such a \(p\),
the valuation of \(\alpha\) equals that of \(\gamma\). After
removing the known ideal factors above that \(p\), division of
the displayed congruence preserves precision \(p^{k-c_p}\).
The resulting modulus \(M'\) is at least \(M/O_L(1)\).

The integer \(M\) has at most \(X^{\mu\lambda+\varepsilon}\)
possibilities, by numbering its radical and using the
fixed-support count at polynomial height. Fix it. Freeze the
part of \((\alpha)\) on primes above \(3N(\gamma)\), the known
part at primes of \(M'\), and the full part on rational primes
whose valuation in \(Q\) exceeds a cutoff \(J\).
The fixed-prime portions contribute subpower loss. Since
\(Q\ll X^2\), choose \(J\) at the fixed budget so that the
large-valuation portion has radical at most \(X^{\xi\lambda/10}\).
Its norm and ideal factor can be fixed at cost
\(X^{\xi\lambda/10+\varepsilon}\). This uses the fixed-support
estimate followed by the divisor bound for ideals of specified
norm. Denote the frozen ideal by \((\eta_0)\). All factors
left after this freezing are coprime to \(M'\). At retained
prime powers above \(M'\), the element \(\alpha/\eta_0\)
is congruent to \(W(\gamma/\eta_0)\), with both sides units.
Here the congruence is in each completion, so denominators
of \(\gamma/\eta_0\) away from \(M'\) are harmless.

A remaining rational prime dividing \(Q\) splits in \(O\) and
exactly one of the two primes above it divides \(\alpha\).
Indeed an inert prime, or both split primes together, would
divide both integer coordinates of \(\alpha\). The only
ramified prime was frozen. Choose one prime above each split
prime as the allowed orientation. Some orientation retains
an \(X^{-o_X(1)}\) fraction of the rows, by independent
uniform choices and the uniform bound
\(\omega(Q)=O(\log X/\log\log X)\). On the retained rows write
\[
 (\alpha)/(\eta_0)=\prod_{e\le J}\mathfrak a_e^e.
\]
The squarefree ideals \(\mathfrak a_e\) use only allowed primes,
so each is determined by its norm. Bin their norms at upper
caps \(X^{\theta_e\lambda}\) on a fixed mesh at this budget,
with \(\sum_e\theta_e=H_0\le q+\xi/10\). Let
\(\mathcal A_e\) be the dictionary of allowed ideals in the
\(e\)-th cap.

The congruence gives a fixed target for
\(\sum_e e\,g(\mathfrak a_e)\) in the finite abelian group
\[
 G=(O/M'O)^\times\big/
       \big((\mathbb Z/M'\mathbb Z)^\times\,O^\times\big),
\]
where \(g\) sends an ideal prime to \(M'\) to the class of any
generator. For each fixed \(e\), \(G[e]\) has subpower size.
For the ambient unit group, at unramified \(p\ge5\) the
residue-field \(e\)-torsion has size at most \(e^2\), and the
logarithm on principal units bounds the remaining torsion
by \(p^{2v_p(e)}\). At the fixed primes \(2,3\), use the
logarithm on \(1+p^2O_p\) and its additive filtration, at a
bounded index. The logarithm identifies the subgroup modulo
\(p^k\) for \(k\ge2\) with \(p^2O_p/p^kO_p\).
The product of the local bounds is \(X^{o_X(1)}\) since
\(e\) is fixed before the height limit and
\(\omega(M')=O(\log X/\log\log X)\). Passing to the quotient
preserves this upper bound since for a finite abelian group
the size of its \(e\)-torsion equals the size of its quotient
by \(e\)-multiples.

For a positive \(\theta_e\), take
\(m_e=\lceil H_0/\theta_e\rceil\), so
\(m_e\theta_e\le2H_0\) and \(\sum_e m_e^{-1}\le1\).
Zero-cap slots are trivial. On products \(\mathfrak B\) of
\(m_e\) ideals in \(\mathcal A_e\), the map \(g\) is
injective when \(\lambda\) is sufficiently small in terms
of \(\nu,q,\xi\). In fact equal classes for two products
mean that their generators, after one is multiplied by a
unit, are proportional by a rational unit modulo \(M'\).
The determinant of their two integral coordinate vectors
is a multiple of \(M'\). Its absolute value is
\(O(X^{2H_0\lambda})\) since the generator of an integral
Eisenstein ideal of norm \(N\) has coordinates \(O(\sqrt N)\).
For \(2H_0\lambda<\nu\) the determinant is zero. The quotient
of the generators is then rational. A nonzero valuation of
a rational number at a split prime uses both primes upstairs.
The orientation therefore forces equality of the two ideals.

For a uniform character \(\chi\) of \(G\), put
\(S_e(\chi)=\sum_{\mathfrak a\in\mathcal A_e}\chi(e\,g(\mathfrak a))\).
Orthogonality and the injectivity just proved give
\[
 \mathbb E_\chi |S_e(\chi)|^{2m_e}
 \le |G[e]|X^{m_e\theta_e\lambda+o_X(1)}.
\]
There are at most \(X^{m_e\theta_e\lambda}\) first ordered
products. For each, the second products have at most
\(|G[e]|\) possible classes and each ideal product has
subpower many ordered decompositions. Fourier inversion
for the fixed target and Hölder bound the number of ideal
tuples by \(X^{H_0\lambda/2+o_X(1)}\). Each ideal allows
only the finite group of units as generator choices, and
a generator recovers the integer inputs. Including the
frozen factors and the initial modulus count proves (CN4).

## Combining the mixed triples

The proof applies more generally to \(k\ge3\) fixed primitive
pairwise nonproportional integral linear forms
\(L_1,\ldots,L_k\) and the norm \(Q\) on primitive
\((a,b)\) with \(a,b>0\) and \(X/2<a+b\le X\).
Exclude zero values of the forms, which only add finitely many
primitive pairs. Take radicals of the absolute values.
A finite exceptional prime set contains the prime divisors of
the pairwise determinants and the values of \(Q\) on the
primitive kernel vectors. Outside it all these supports
are disjoint. The general joint cap on \(\operatorname{rad}(Q\prod_i L_i)\)
gives limiting coefficient at most
\[
                     \frac{1}{k+2}.
 \tag{CN5}
\]
Here and below a limiting coefficient \(s\) gives the count
\(O_{\eta,\lambda,\varepsilon}(X^{(s+\eta)\lambda+\varepsilon})\)
for every \(\eta>0\) and all sufficiently small fixed \(\lambda\).
Constants may depend on the forms and cap coefficients.

Any three forms \(L_i,L_j,L_l\) satisfy a linear relation
with nonzero fixed integer coefficients given by the pairwise
determinants. Fix the signs of the three values, put the
opposing term on the other side of the additive equation,
and divide by the gcd. This gcd is bounded independently of
\(a,b\), since the gcd of two weighted values divides the
product of the two fixed weights and their determinant at a
primitive input. There are only finitely many sign and gcd
records. For a fixed record the derived primitive positive
triple determines \((a,b)\), and its largest term is comparable
to \(X\) by independence of any pair of forms. At all primes
outside the fixed exceptional set, each valuation on a selected
form is exactly the valuation on its corresponding letter.
The fixed multipliers and bounded gcd affect only the
exceptional set.

We give the argument with explicit ordering of limits. If a
claimed exponent with allowance \(\eta>0\) fails, choose fixed
budgets \(\lambda_j\to0\) at which larger height exponents
occur. A failure of an \(X^{s\lambda_j+\varepsilon}\) estimate
has upper height exponent at least \(s\lambda_j\).
Dyadic decomposition in \(a+b\) contributes a subpower
factor and the mixed heights are in a bounded number of
dyadic shells. At each fixed budget dyadically bin the
complete form values, their valuation roots, and the
fixed-prime and large-degree portions below. These add
\(X^{o_X(1)}\) records at the fixed budget. Refine the
radical of \(Q\) on a fixed logarithmic mesh with width
\(o(\lambda_j)\). Choose a box preserving the count exponent
up to \(o(\lambda_j)\). Take the inner height limit first,
then a subsequence in \(j\). Let \(\beta_i,q,d\) denote the
limiting exponents of \(\operatorname{rad}(L_i)\),
\(\operatorname{rad}(Q)\), and the count, respectively,
all divided by \(\lambda_j\). A bin uses upper endpoints
which differ from the actual normalized masses by a vanishing
amount. Since exceptional prime overlap has logarithmic cost
\(O(1)\), the joint cap passes to
\(\sum_i\beta_i+q\le1\) for (CN5). The separate caps in (CN2)
similarly pass to \(\sum_{D}\beta_i\le B_D\),
\(\sum_{C}\beta_i\le B_C\), \(q\le B_Q\).
An enlargement of a height-relative budget by an arbitrarily
small constant multiple of \(\lambda_j\) absorbs the bounded
changes in height and the bounded exceptional radicals.

The limiting one-third bound on any derived triple gives
\[
      \beta_i+\beta_j+\beta_l\ge3d,\qquad
      \sum_i\beta_i\ge k d .
 \tag{CN6}
\]
The sum inequality follows by averaging over all triples.
Likewise \(q\ge d\). To see this, enumerate its possible
complete norm values by the radical and the fixed-support
count. Each value is the norm of only subpower many integral
Eisenstein ideals, and each ideal has finitely many generators.
[Limiting small-budget coefficient \(1/3\)](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

Suppose some form admits divisors of arbitrarily small
positive normalized radical cost with positive limiting
numerical height. Precisely, for arbitrarily small fixed
\(\mu>0\) there is a subsequence and a fixed \(\nu>0\)
such that each row in its boxes has an available divisor
with \(\operatorname{rad}(M)\le X^{\mu\lambda_j}\) and
\(M\ge X^\nu\). Formula (CN4) at that fixed \(\mu,\nu\)
gives \(d\le \mu+q/2\), taking its allowance to zero.
The thresholds on \(\lambda_j\) then hold along the
subsequence. Sending \(\mu\) to zero gives \(q\ge2d\).
Together with (CN6) this yields
\(\sum_i\beta_i+q\ge(k+2)d\).

## The alternative with no cheap large divisor

Here every form has a uniform positive lower threshold
for the cost of a divisor with positive limiting height,
in the preceding sense. At budget \(\lambda_j\), put all
primes in the fixed exceptional set and all powers with
degree \(e>\lambda_j^{-2}\) into a separate portion of each
form. Its normalized radical cost tends to zero by the
height bound. Hence its numerical height tends to zero
in this alternative. All remaining portions use the exact
squarefree roots \(x_{i,e}\) of the prime powers of degree
\(e\). Let \(\beta_{i,j,e}\) be their radical costs
normalized by \(\lambda_j\). Their total costs are bounded
and their numerical contributions are
\(\lambda_j e\,\beta_{i,j,e}\), with an error vanishing
in the inner height limit.

Choose a fixed \(\xi>0\) such that no form admits a
positive-height divisor of cost at most \(3\xi\) along
any subsequence. Set aside all coordinates having
\(\beta_{i,j,e}\ge\xi\). There are boundedly many for
each form, uniformly in \(j\). If the remaining
coordinates carry positive limiting height on a
subsequence, discard those with \(\lambda_j e<r\),
for a fixed \(r>0\) small enough to leave positive
height. The discarded contribution is at most \(r\)
times the uniformly bounded total cost. Group the
remaining coordinates greedily into a bounded number
of blocks of cost at most \(2\xi\), filling each
nonfinal block to cost at least \(\xi\).
One block carries positive limiting height on a further
subsequence. Its product of full prime powers is
an available divisor for every row at cost at most
\(2\xi+o(1)\), a contradiction.

Thus boundedly many **exact degree coordinates** carry the
whole limiting height of each form. Pass to a subsequence
where the number kept and each ordered coordinate's mass
and numerical height converge, discarding those whose
limiting numerical contribution is zero. Each of the
remaining coordinates has order comparable to
\(1/\lambda_j\). Write \(P_i\) for their total limiting
mass. Every scaled degree measure of a form is atomic.
Indeed the selected coordinates have subsequential limits
of their scaled degrees in a compact positive interval,
and the complement has vanishing numerical weight, so
its weak mass can be supported only at zero.

If one form has limiting numerical height zero, take a
mixed primitive triple containing that form. Its other two
letters \(U,V\) have height comparable to \(X\) by linear
independence and the relation, with \(U<V\) and the short
letter equal to \(V-U\). Accordingly
\[
                         |U/V-1|\le X^{-1+o(1)} .
\]
In the ratio \(U/V\), retain their finitely many positive-weight
roots as separate coordinates with their actual signed
exponents. The ratio of their complementary factors,
including the fixed coefficient and gcd adjustments,
has height \(X^{o(1)}\). An exact value of \(U/V\)
recovers both primitive terms, the third term, and the
original pair for the fixed mixed-triple record. It
therefore fixes all retained roots and the complementary
ratio. Apply the real strip estimate (RUW1) with \(\tau=1/2\),
each narrow height constant fixed, and arbitrary \(s>0\).
For each \(s\) the admissible wider-height threshold is fixed
and positive, so the complementary ratio meets it for all
large \(j\). The narrow exponents are fixed in the inner
height limit, and the estimate's small-budget threshold is
independent of them. It follows that \(d\le (m^2+m)s\),
where \(m\) is the bounded number of retained roots.
Sending \(s\) to zero gives \(d=0\). Thus every form has
positive limiting numerical height when \(d>0\) in this
alternative.
[Real strip with unequal root heights (RUW1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

Every positive-weight coordinate on any form must have
limiting mass at least \(d\). Apply (GE1) to a mixed
triple containing that form, using the entire coordinate
as modulus. The other two scaled degree measures have
zero nonatomic mass. To match the normalization in
(GE1), the mixed triple's total root budget relative to
\(\lambda_j\) tends to \(\beta_i+\beta_j+\beta_l\).
This limit is positive by (CN6) when \(d>0\).
The common exponent and the modulus cost are both divided
by this same factor. The compression may be taken above
all the retained orders, and bounded coefficients and gcd
only affect the fixed exceptional prime set. The modulus
is prime to that set, has positive limiting weight, and
has fullness comparable to the reciprocal of the mixed
budget. These are the hypotheses of (GE1) before its
balance assumption. The bound therefore gives
\(d\le\beta_{i,e}\) in the original normalization.
[Atom grouping with a cheap full modulus (GE1)](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

For any mixed triple freeze the complementary factors of
each form at total cost at most
\(\sum_{i\in I}(\beta_i-P_i)\), where \(|I|=3\).
This follows by numbering the complementary radicals in
the exact degree boxes, then using the subpower
fixed-support count. The portion at the fixed exceptional
primes and above the cutoff costs \(o(1)\) in normalized
units. The resulting additive equation has fixed
coefficients multiplying monomials in boundedly many
pairwise coprime high-degree roots. These roots have
height \(X^{O(\lambda_j)}\), and their degrees tend to
infinity. The bounded-root release (30) counts the
moving roots at vanishing normalized cost, uniformly
in the frozen coefficients. Therefore
\[
 \sum_{i\in I}\beta_i\ge d+\sum_{i\in I}P_i
 \ \ge\ 4d.
 \tag{CN7}
\]
Averaging (CN7) over all triples gives
\(\sum_i\beta_i\ge (4k/3)d\). Since \(q\ge d\), the
total cap is at least \((4k/3+1)d\ge(k+2)d\) for
\(k\ge3\). Combining the two alternatives proves (CN5),
including (CN1) for \(k=6\). In particular the same
method gives \(1/5\) for the original three forms
and \(Q\), and \(1/6\) upon adding any fourth direction.
[Bounded-root release (30)](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#releasing-a-bounded-number-of-high-degree-coordinates)

For (CN2), the two specific triples of linear forms give
\(d\le B_D/3\) and \(d\le B_C/3\) by (CN6).
The norm fibre gives \(d\le B_Q\). In the first alternative
\(d\le B_Q/2\), and in the second (CN7) gives
\(d\le(B_D+B_C)/8\). This proves the displayed
maximum and completes the exponent comparison.

## The norm obstruction for a single original triple

With just \(\operatorname{rad}(D)\le c^\lambda\) and
\(\operatorname{rad}(Q)\le c^{B_Q\lambda}\), the coefficient
in the count is at most
\[
       \min\!\left\{\frac13,B_Q,\
           \max\!\left(\frac14,\frac{B_Q}{2}\right)\right\}
       +\eta .
 \tag{CN8}
\]
Use only \(a,b,a+b\) in the two alternatives. The norm
fibre and the one-third bound give \(d\le B_Q\) and
\(d\le1/3\). The cheap-divisor alternative gives
\(d\le B_Q/2\), and (CN7) gives \(d\le1/4\) in the
other alternative. Thus an unrestricted endpoint family
can have only a strictly smaller contribution from the
sector \(\operatorname{rad}(a^2+ab+b^2)\le c^{(2/3-\rho)\lambda}\)
for any fixed \(\rho>0\). This consequence needs no radical
condition on the companion.

The same argument admits \(r\) fixed quadratic forms
\(Q_j(a,b)=Q(U_j(a,b),V_j(a,b))\) for fixed unimodular
linear changes of variables, provided that the \(Q_j\)
are pairwise nonproportional. Their prime supports are
disjoint from each other and from the linear forms away
from the finite set of divisors of the corresponding
resultants. Apply (CN4) in each changed coordinate pair
with the same original divisor. If a cheap divisor is
available, each norm mass is at least \(2d\), and the
linear masses sum to at least \(kd\). Otherwise each
norm mass is at least \(d\) and (CN7) applies to the
linear masses. A joint cap on the norms and the linear
forms gives the coefficient
\[
             \frac{1}{\min(k+2r,\,4k/3+r)}+\eta .
 \tag{CN9}
\]

## Power classes of the centered reciprocal value

For a prime \(\ell\) and \(h\in\mathbb Q^\times\), write
\[
 {\cal D}_\ell(h)=\prod_{\ell\nmid v_p(h)}p .
\]
For every fixed prime \(\ell\) and \(\varepsilon>0\), the number
of rational \(z\in(0,1/2)\) with
\({\cal D}_\ell(H(z))\le Y\) is
\[
                       O_{\ell,\varepsilon}(Y^{1+\varepsilon}).
 \tag{CN10}
\]
The height of \(z\) is unrestricted. Away from \(2,3\) the
defect primes are exactly the primes on the six linear
forms whose individual valuation is not divisible by
\(\ell\). Thus a cap \(X^{(1/3-\rho)\lambda}\) on their
joint defect gives a strict saving in the original
exceptional set for fixed \(\rho>0\), even when the
companion radical itself is unrestricted. Omitting the
fixed primes \(2,3\) in the defect changes the cap by
at most six.

For the proof let \(E=\mathbb Q(\sqrt{-3},\zeta_\ell)\),
\(A=3\sqrt{-3}\) and \(U(z)=(z-\eta)/(z-\bar\eta)\)
where \(\eta=(1+\sqrt{-3})/2\). The reciprocal identity
gives
\[
 U(z)^3=\frac{H(z)-A}{H(z)+A}.
\]
Each class of \(\mathbb Q^\times/(\mathbb Q^\times)^\ell\)
has an integral representative \(d\ne0\) with all positive
valuation exponents less than \(\ell\). For odd \(\ell\)
one may choose \(d>0\), and for \(\ell=2\) the sign is
fixed by the class. If \(H(z)=d\,t^\ell\), write \(v=U(z)\)
to obtain an \(E\)-point on
\[
             t^\ell=(A/d)(1+v^3)/(1-v^3).
 \tag{CN11}
\]
The parameter \(v\) determines \(z\) since \(U\) is invertible.
There are exactly six distinct split branch points with
orders \(1\) or \(-1\). The smooth projective curve has
genus \(2(\ell-1)\). At primes outside \(6\ell d\), the
sections remain disjoint and the scalar \(A/d\) is a unit,
so the usual tame-cover charts give good reduction.

For completeness let \(\sigma\) act on \(t\) by multiplying
by \(\zeta_\ell\), and let \(\phi=1-\sigma\) on the
Jacobian. The deck norm is zero as the base is a projective
line. Consequently \(\phi^{\ell-1}\) is multiplication by
\(\ell\) times an automorphism and \(\deg\phi=\ell^4\).
The rational branch-point differences have degree-zero
coefficient vectors modulo \(\ell\) in dimension five.
Any principal relation between them is the divisor of
a deck eigenfunction \(t^j f(v)\). Its branch orders
modulo \(\ell\) lie in the one-dimensional span of the
six branch orders of the defining rational function.
Thus the branch differences give all \(\ell^4\) points
of \(J[\phi]\) over \(E\). The Kummer injection for \(\phi\)
then takes values in four copies of the \(\ell\)-power
classes unramified outside \(6\ell d\). This follows by
good reduction and the finite étale \(\phi\) map at the
remaining places. The \(S\)-unit rank and the fixed-field
class group bound the logarithm of the number of these
classes by \(O_\ell(1+\omega(d))\). The rank of \(J(E)\)
has the same bound by the relation between \(\phi\) and
multiplication by \(\ell\). Embedding with a rational
branch point and applying uniform Mordell–Lang gives
\(\exp(O_\ell(1+\omega(d)))\) points for the fixed class.
For each possible \(\operatorname{rad}(|d|)\le Y\),
there are at most \(2(\ell-1)^{\omega(d)}\) representative
choices. Summing with the divisor bound proves (CN10).
For squares, the class is the signed squarefree class
of \(DC\), so this case also follows from split genus-two
descent for six linear factors.
[Centered Cayley identity (GR2)](ABC_PROJECTIVE_RECIPROCAL_GENUS_TWO_RECOVERY.md#affine-pair-recovery-over-a-fixed-quadratic-field)
· [Split cyclic-cover descent](ABC_BINARY_POWER_TWISTS.md#projective-covers-from-conjugate-ratios)
· [Gao–Ge–Kühne, Theorem 1.1′](https://arxiv.org/html/2105.15085v4)

## Scope

Without the ray moment, averaging the mixed three-form bounds
and using the individual norm count gives only
\(\sum_i\beta_i+q\ge(k+1)d\). The joint theorem improves
\(1/7\) to \(1/8\) for the six centered-companion forms.
The unconditional \(abc\) endpoint remains open. The original
radical cap gives no bound here on the companion or the
Eisenstein norm radical. Even with both three-form caps, an
endpoint sequence can survive this argument only with norm
radical mass at least \(2/3\) in units of the original budget.

For the Frey model \(y^2=x(x-a)(x+b)\), one has
\(c_4=16Q\), \(|c_6|=32C\) and \(\Delta=16D^2\).
Thus \(\operatorname{rad}(Q)\), away from two, is the
radical of the numerator of \(j=256Q^3/D^2\).
Counting by conductor alone constrains the denominator
support. The conductor families of Shankar, Shankar and
Wang, Theorems 1.1, 1.4 and 1.5, impose a size restriction
on the discriminant-to-conductor ratio or require it to
be squarefree. The high powers in an exceptional primitive
triple can violate both hypotheses.
[Shankar–Shankar–Wang, conductor families](https://arxiv.org/html/1904.13063)

The fixed-\(S\) binary-form results of Bugeaud, Evertse and
Győry, Theorems 2.4–2.6, estimate the contribution of a
prescribed finite prime set and allow constants depending
on that set. The caps here are on the total support and
the support varies with each row. The unrestricted
count of Bernert, Browning, Lichtman and Teräväinen,
Theorem 1.3, has exponent \(3/5+\varepsilon\) for
\(0<\lambda<1\). Browning–Verzobio's Theorem 1.2 is a
coefficient-uniform fixed-order count for generalized
Fermat equations. The mixed cap proof uses the internal
small-budget endpoint and a separate norm recovery,
with no requirement of exact primitive power orders.
[Bugeaud–Evertse–Győry, binary form \(S\)-parts](https://arxiv.org/pdf/1708.08290)
· [Bernert–Browning–Lichtman–Teräväinen, Theorem 1.3](https://arxiv.org/html/2410.12234v2)
· [Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1)
