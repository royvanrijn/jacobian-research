# Primitive norm sums for every negative prime discriminant -p, p >= 7

Fix a prime \(p\ge7\) with \(p\equiv3\pmod4\), put
\(K=\mathbb Q(\sqrt{-p})\), and write \(G=\mathrm{Cl}(\mathcal O_K)\).
The group \(G\) automatically has odd order by PrimeDiscriminantOdd below.
No class-number restriction is imposed, and \(G\) need not be cyclic.
Write
\[
 Q_p(x,y)=x^2+xy+\frac{p+1}{4}y^2,\qquad
 P_K=\{Q_p(x,y)>0:\gcd(x,y)=1\},
\]
and let \(R_K(N)=\#\{(a,b)\in P_K^2:a+b=N\}\) count ordered value pairs.
For the quadratic character \(\chi=\chi_{-p}\), put
\[
 J_K(N)=
 \prod_{\substack{\ell\mid N\\\ell\nmid2p,\ \chi(\ell)=-1}}
                    \frac{\ell-1}{\ell-2}.                    \tag{OC.1}
\]
The exclusion of 2 is essential when 2 is inert.

**Theorem OddClassCounts.** There are positive constants \(c_K,C_K,N_K\)
such that, for every \(N\ge N_K\) satisfying
\[
 p^2\nmid N,\qquad
 \text{and }2\mid N\text{ if }p\equiv3\pmod8,                 \tag{OC.2}
\]
one has
\[
 c_K\frac{N J_K(N)}{\log N}
 \le R_K(N)\le
 C_K\frac{N J_K(N)}{\log N}.                                \tag{OC.3}
\]
Condition (OC.2) is precisely primitive local admissibility of the sum
at all finite primes. Every positive target failing it has count zero.
The constants and threshold depend on the fixed field and fixed auxiliary
primes; they are neither numerical nor uniform as the field varies.

This is a deduction from the existing quadratic splitting-semigroup
machinery with finite exact-prime profiles. It extends the explicit
[class-number-three construction](PRIMITIVE_NORM_CLASS_NUMBER_THREE.md)
to every negative prime discriminant in the stated range. Individual primitive
principal norms can still have genuine class-group obstructions. The
additive lower bound only needs a sufficient family on which finitely
many fixed primes cancel those obstructions. No nonabelian additive
distribution estimate is required. No historical novelty or external
verification is claimed.

## The class group automatically has odd order

**Lemma PrimeDiscriminantOdd.** If \(p\equiv3\pmod4\) is prime,
the ideal class group of \(\mathbb Q(\sqrt{-p})\) has odd order.
This classical fact is proved here to remove a redundant hypothesis;
the lemma also holds for \(p=3\), although the additive theorem
retains its separate restriction \(p\ge7\).

**Proof.** For a nonzero fractional ideal \(\mathfrak a\),
\[
          \mathfrak a\overline{\mathfrak a}=(N\mathfrak a).
\]
This follows prime by prime: a split conjugate pair multiplies to
the rational prime ideal, an inert prime ideal is already rational,
and a ramified prime ideal squares to the rational prime ideal.
The identity also allows negative ideal exponents. Conjugation thus
acts by inversion on ideal classes.

If \([\mathfrak a]^2=1\), then
\(\mathfrak a=(\alpha)\overline{\mathfrak a}\) for some
\(\alpha\in K^*\). Taking absolute ideal norms gives
\(\lvert N_{K/\mathbb Q}\alpha\rvert=1\); imaginary quadratic norm
positivity gives \(N\alpha=\alpha\bar\alpha=1\).
If \(\alpha\ne-1\), put \(\beta=1+\alpha\); then
\(\beta/\bar\beta=\alpha\). If \(\alpha=-1\), take
\(\beta=\sqrt{-p}\), with the same quotient. Hence
\[
 \mathfrak b=\mathfrak a/(\beta)
    =\overline{\mathfrak a}/(\bar\beta)=\overline{\mathfrak b}.
\]
There is no integrality requirement on \(\beta\): fractional ideals
are the appropriate setting for this step.

Factor the invariant fractional ideal \(\mathfrak b\).
Split conjugate primes have equal exponents and give rational
principal factors; inert prime ideals are themselves rational.
Modulo principal ideals, only ramified prime classes remain, with
their exponents reduced modulo two. This reduction works for negative
exponents as well. The field discriminant is \(-p\), so there is only
one ramified prime ideal, and it is principal:
\[
                  \mathfrak p_p=(\sqrt{-p}).
\]
The generator is integral and has norm \(p\). In particular 2 is
unramified and contributes no additional class. It follows that
\(\mathfrak b\), and therefore \(\mathfrak a\), is principal.

The finite class group has no nonidentity element fixed by inversion.
Pairing every other class with its inverse shows its order is odd.
Finiteness follows from the retained Minkowski theorem, Milne v3.08,
Theorem 4.3. The elementary norm-one construction above imports no
Hilbert--90 or genus-theory hypothesis. \(\square\)

## Primitive ideals and a finite signed class correction

In a maximal quadratic order, a primitive norm has no inert prime
factor, at most exponent one at each ramified prime, and an orientation
of every entire split prime power whose ideal-class sum is zero. Whole
prime powers must be oriented together: assigning their copies separate
orientations can produce an imprimitive generator. These facts and the
local norm descriptions are proved in
[Primitive norm local structure](PRIMITIVE_NORM_LOCAL_STRUCTURE.md).

**Lemma OddClassAbsorber.** Let \(G\) be a finite abelian group of odd
order and choose nonidentity generators \(g_1,\ldots,g_r\), of orders
\(h_i\). They need not be independent. Put \(m_i=h_i-1\).
A block of \(m_i\) independent signs with classes \(\pm g_i\) can
contribute any element of \(\langle g_i\rangle\). All blocks together
can therefore cancel any prescribed element of \(G\).

If, in a maximal imaginary quadratic order with full class group \(G\),
there are distinct split rational primes in each of these blocks of
classes \(\pm g_i\), assign each prime exponent exactly one. Their product
times any cofactor coprime to them, having no inert factor and ramified
valuations at most one, is a primitive principal norm.

**Proof.** With \(j\) negative signs in block \(i\), its class sum is
\[
             (m_i-2j)g_i,\qquad 0\le j\le h_i-1.             \tag{OC.4}
\]
These are all \(h_i\) elements because multiplication by 2 is invertible
modulo odd \(h_i\). The cyclic subgroups of the generators sum to \(G\),
even if the generator list has relations.

Orient the entire split prime powers of the cofactor arbitrarily. The
resulting ideal has some class \(c\). Ramified prime ideal classes have
order dividing two and hence are trivial in \(G\). Choose the reserved
prime orientations to contribute \(-c\); the product ideal is principal.
It is not contained in \(\ell\mathcal O_K\) for any rational prime
\(\ell\), by the absence of inert primes,
the ramified exponent bound and the exact split orientations. An integral
generator is therefore primitive. Distinct reserved rational primes and
cofactor coprimality are necessary for the independent exact-exponent-one
argument. \(\square\)

For the trivial group use an empty generator list and empty prime pools.
This recovers the class-number-one semigroup situation.

The exponent condition cannot be replaced by divisibility alone. In the
class-number-three example, \(351=3^3\cdot13\) has only one active
nonprincipal prime power and is not a primitive principal norm, although
both 3 and 13 divide it. The residual coprimality guard must persist
through every later cofactor restoration.

## Existence of the finitely many class primes

For each \(i\), reserve a pool \(\mathcal F_i\) of \(2m_i\) distinct
odd rational split primes of classes \(\pm g_i\). Choose all pools
disjoint and avoid \(p\). Such finite lists exist.

The precise classical input is Milne, *Class Field Theory*, v4.03,
VIII Theorem 7.2, printed p.258/PDF p.267: prime ideals in every
congruence class have Dirichlet density the reciprocal of the index.
Apply modulus 1 and the principal ideal subgroup. Prime ideals in each
ordinary ideal class have density \(1/|G|\). Restricting to residue degree
one over \(\mathbb Q\) preserves this density, as in VI Proposition 4.5
and Corollary 4.6, printed pp.196--197/PDF pp.205--206. Equivalently,
the degree-at-least-two prime ideals have bounded Dirichlet sum at 1,
since their norms contribute a convergent sum over rational prime
squares. After deleting ramified primes, degree-one primes in a
quadratic field lie over split rational primes. Each rational prime
accounts for at most two such ideals. Every class thus supplies infinitely
many distinct rational primes, and every finite bad set can be avoided.

The exact retained primary edition is
[sources/milne-cft403.pdf](sources/milne-cft403.pdf); its checksum and the
inspected statements are recorded in the source ledger. This input chooses
only a fixed finite pool. The additive estimates below involve the
quadratic Dirichlet character \(\chi\), with no varying class-field
prime-distribution requirement.

Write
\[
                      F=\prod_i\prod_{q\in\mathcal F_i}q.
\]
For any target \(M\), give each pool prime an exact valuation pair
\((0,0),(1,0),(0,1)\), or \((1,1)\) in its two core summands.

**Lemma OddClassProfiles.** Every target has a profile giving each
summand exactly \(m_i\) primes from pool \(\mathcal F_i\), with common
primes dividing the target and exclusive primes not dividing the target.
At every pool prime, that profile has positive local mass on
\(a+t b=M\), for \(t=1\) or 2, uniformly in all allowed target valuations.

**Proof.** Let \(k_i=\#\{q\in\mathcal F_i:q\mid M\}\).
If \(k_i\ge m_i\), select \(m_i\) dividing primes common to both summands.
If \(k_i<m_i\), select all \(k_i\) dividing primes common and add
\(m_i-k_i\) nondividing primes exclusively to each side, using distinct
primes for the two sides. There are
\[
              2m_i-k_i\ge2(m_i-k_i)
\]
nondividing primes available. Set every remaining pool prime absent from
both summands. This gives exactly the required block on each side.

For a pool prime \(q\), hence odd and split, both coefficients of
\(a+t b=M\) are units.
For exact exponent pairs 0 or 1, the weighted local masses for
\(V=\tau\mathbf1_C\), where \(C\) is the odd split semigroup, are
\[
\begin{array}{c|c|c}
(\epsilon,\delta)&\text{condition}&\text{mass}\\ \hline
(0,0)&q\nmid M&(q-2)/q\\
(0,0)&q\mid M&(q-1)/q\\
(1,0),(0,1)&q\nmid M&2(q-1)/q^2\\
(1,1)&v_q(M)=1&4(q-2)/q^2\\
(1,1)&v_q(M)\ge2&4(q-1)/q^2 .
\end{array}                                                \tag{OC.5}
\]
These follow from the complete residues modulo \(q^2\): each exponent-one
factor has weight 2. Unweighted counts have the same positivity. The
assignment never selects a zero-mass case. Since \(q>2\), every listed
mass is positive and bounded below independently of target valuation.
The unit linear derivative gives all higher lifts; in particular common
exponents one do not require the target to have valuation exactly one.

Let \(a_0,b_0\) be the assigned products. Their gcd \(c\) divides \(M\).
In the normalized equation
\[
 (a_0/c)x+t(b_0/c)y=M/c,\qquad (xy,F)=1,                    \tag{OC.6}
\]
exclusive coefficient primes do not divide the target, and common
coefficient primes have been canceled. Other pool primes are unit
conditions on both residual variables. These local equations all retain
positive mass, at arbitrary remaining target depth. \(\square\)

For \(t=1\), or \(t=2\) when 2 splits, every resulting value
\(a_0x\) or \(t b_0y\), with split-supported residuals coprime to
\(F\), is a primitive principal norm by OddClassAbsorber. If \(t=2\), the extra split prime can have any ideal
class; its contribution is part of the cofactor class canceled by the
odd reserved blocks. This does not require 2 to be principal.

## The fixed-period and fixed-profile analytic transfer

Here is the analytic interface needed from the existing theory. Set
\(L=2p\), let \(H=\ker\chi\) on units modulo \(L\), and let \(C\)
contain the integers supported on odd split primes, excluding \(p\).
Use normalized shifted weights
\[
 V=\tau\mathbf1_C,\qquad F_B(m)=V(Bm)/V(B),\qquad
 G_D(l)=V(Dl)/V(D).
\]
For \(t=1\), the uniform weighted line estimate is the quadratic-kernel
specialization of (H2) in
[Asymmetric subgroup splitting](ASYMMETRIC_SUBGROUP_SPLITTING.md).
For \(t=2\), used only when \(\chi(2)=1\), apply the same explicit
transfer as [CN2Filter](PRIMITIVE_NORM_CLASS_NUMBER_TWO.md) and
CN3.5--CN3.5a in the class-number-three note, with \(L=2p\).
The first-variable filter is
\[
 \mathbf1_{2\mid M-um}\,
 \mathbf1_{((M-um)/2,\,2p)=1,\ \chi((M-um)/2)=1}.             \tag{OC.7}
\]
It has fixed period \(4p\). Its support removes the 2-part of
\(um\equiv M\pmod{2vka}\), leaving only the good varying modulus
\(vka\). Expand this filter on the first unit variable.

The old two character factors, coefficient reduction, retained smooth
scalar phases and signed bilinear bounds apply with this fixed period.
Kernel indices need only be units with \(\chi=1\), not numbers with
split prime support. Principal twists reassemble by the pointwise
identity \(F_B\xi=F_B\) for \(\xi\in\{1,\chi\}\). The fixed-period
box theorem and (H17)--(H18) error hierarchy have no special requirement
on \(p=23\); \(4p\) is fixed before the target grows. Thus there are
fixed \(\eta,\theta>0\), with \(\theta\) chosen last, such that
\[
 \sum_{um+tvl=M}F_B(m)G_D(l)
 =
 \frac{M}{utv}\Gamma_{t;u,v}(M)\rho_{F_B}\rho_{G_D}
 \prod_{\ell\nmid2p}\frac{D_\ell}{M_{A,\ell}M_{B,\ell}}
 +O_K(M^{1-\eta}).                                       \tag{OC.8}
\]
The good-prime factors are exactly (H3)--(H4) on \(ux+tvy=M\),
including coefficient/gcd zero factors. The estimate is uniform in
good shifts and coefficients \(B,D,u,v\in C\) at most \(M^\theta\),
all their gcds and target valuations, and the small divisibility
predicates used in unweighting. Fixed coefficients eventually fit
the power budget. At good primes \(t\) is a unit.

The actual finite density is
\[
 \Gamma_{t;u,v}(M)=
 \frac{tL}{|H_{tL}|^2}
 \#\{(a,b)\in H_{tL}^2:ua+tvb=M\pmod{tL}\}.                \tag{OC.9}
\]
Allowed coefficients permute \(H_{tL}\). On the appropriate parity,
for \(p\nmid M\),
\[
 \Gamma_t(M)=\frac{8p}{(p-1)^2}s_p(M),\qquad
 s_p(M)=\frac{p-1-2(M/p)}4>0.                             \tag{OC.10}
\]
Here \(s_p(M)\) counts ordered nonzero square residues whose sum is
\(M\pmod p\). For a unit target,
\(\sum_a\chi(a(M-a))=1\): completing the square reduces this to
\(\sum_x\chi(x^2-d)=-1\), \(d\ne0\), which follows by counting
\((x-y)(x+y)=d\). Expanding the two square indicators, excluding
\(a=0,M\), gives \(4s_p(M)=p-1-2\chi(M)\).
Since \(p\ge7\), it is at least \((p-3)/4\ge1\).
For \(t=2\), \(\chi(2)=1\); the extra lifts at period \(4p\) give the
same normalization as for \(t=1\). For \(p\mid M\), the density is zero.

We also need the positive unweighted conclusion, rather than a weighted
asymptotic alone. Fix one of the finitely many profiles in (OC.5).
After factoring \(a_0,b_0\), its weights factor exactly:
\[
 V(a_0x)=V(a_0)V(x),\qquad
 V(b_0y)=V(b_0)V(y),\qquad (xy,F)=1.                       \tag{OC.11}
\]
Keep \(F\) in the fixed local block throughout section 6 of the subgroup
theorem and the positive restoration argument in
[Asymmetric semigroup scale](ASYMMETRIC_SEMIGROUP_SCALE.md).
Exclude its primes from every variable roughness, directive and
common-cofactor pool. This restriction concerns actual restored factors;
analytic kernel indices may still contain pool primes and their local factors
must be retained. Consequently every restored factor is coprime
to \(F\), and none changes an assigned exponent. Merely imposing the
initial guard would be insufficient.

Deleting this fixed finite set changes only constant Euler factors.
The remaining weights have rough dimension two, restoration exponent
one, and the same cofactor size budget and recoverable patterns as the
existing proof. The shifted estimates (OC.8) cover the same small-power
tests; the finitely many new local masses have a common positive lower
bound by (OC.5). Removing common coefficients as in (H15) preserves
parity and \(\chi(M)\), because they are odd split products. The finite
density is unchanged, and \(J_K(M)\) is unchanged by these factors.
Hence every chosen positive profile contributes
\[
                         \gg_K M J_K(M)/\log M            \tag{OC.12}
\]
unweighted ordered value pairs for every sufficiently large allowed
target, uniformly in all target valuations. This is exactly the
CN3Profile proof with a larger fixed pool, not a new assertion that
weighted asymptotics automatically unweight. Taking a maximum cutoff
and minimum constant over the finite family of target-divisibility
templates makes the conclusion pointwise for every admissible target.

## Local conditions, ramification and completion

At an odd inert prime, a primitive norm is a unit and every unit is a
local norm. Two units can sum to any target over an odd residue field,
and the norm maps and the linear equation lift. At a split prime, a
primitive local norm may have any valuation: use one unit component of
the split algebra. These primes give no additive local exclusions.

At \(p\), complete the square:
\[
                     4Q_p=(2x+y)^2+p y^2.
\]
A primitive unit norm has nonzero square residue modulo \(p\). A
primitive norm divisible by \(p\) has valuation exactly one and a
nonzero square unit quotient. These conditions are sufficient by
lifting unit squares. Since \(-1\) is nonsquare, two unit norms cannot
sum to zero modulo \(p\); if a sum is divisible by \(p\), both summands
must have valuation one. Their unit quotients cannot sum to zero
modulo \(p\) either. Thus \(p^2\mid N\) is impossible. For a unit target,
or a target of valuation one after dividing by \(p\), positivity of
\(s_p\) and unit-square lifting give local representations. This uses
\(p\ge7\); \(p=3\) has an additional ramified residue obstruction and
is not included here.

At 2, \(p\equiv3\pmod8\) gives the unramified quadratic algebra with
primitive norms exactly the odd units, and their sums cover exactly
the even targets. For \(p\equiv7\pmod8\) the algebra splits; one unit
component permits every valuation and every additive target. This
proves exactly (OC.2) as the complete primitive local sum condition.

Write an admissible target as \(N=p^eM\), where \(e=0\) or 1 and
\(p\nmid M\). For even \(M\), use odd cores on \(a+b=M\), with \(t=1\).
For odd \(M\), admissibility implies 2 split; use \(a+2b=M\), with
\(t=2\). Assign the profiles above. The absorbers make both constructed
values primitive principal norms, so (OC.12) gives the lower bound
for \(M\). These are sufficient subsets only: no decomposition
\(P_K=C\cup2C\) is asserted.

If \(e=1\), multiply both norm values by \(p\), using multiplication
of their generators by \(\sqrt{-p}\). For a \(p\)-unit primitive norm
this preserves primitiveness: at primes other than \(p\) the map is
invertible, and at \(p\) the new ideal has exact ramified exponent one.
It injects ordered value pairs into the target \(N\). The factor \(p\)
only changes a field-dependent constant, and \(J_K(N)=J_K(M)\).

For the upper bound ignore all class constraints, ramified valuations
and the prime 2; retain only the absence of odd inert prime factors
in each summand. The direct upper sieve from the existing quadratic
semigroup proof has dimension one and gives
\(\ll_K N J_K(N)/\log N\), uniformly in the target. This applies also
to odd targets when 2 splits, without applying an odd-core theorem
outside its parity hypotheses. Together with the local zero cases and
the lower construction, this proves OddClassCounts.

## Evidence and remaining scope

[The finite absorber checker](verify_primitive_norm_odd_class_absorbers.py)
checks complete abstract class supports, target profiles, class-pair
compensation and corruption controls on its recorded finite groups.
It also checks the square-pair and finite-density formulas at the recorded
finite list of primes, including rejection of the inert-2 coefficient mode.
It does not certify prime lists for every field or any infinite analytic
estimate. The C3 note separately supplies explicit class certificates
for its four rational primes and an exhaustive individual norm audit.

The [source ledger](primitive_norm_odd_class_source_ledger.json) pins the
retained primary source and exact canonical proof dependencies. Reproduce
the finite receipt with:

    python3 research/arithmetic/verify_primitive_norm_odd_class_absorbers.py --output research/arithmetic/primitive_norm_odd_class_absorber_checks.json

The general proof uses the primary class-prime
density theorem and the existing analytic theorems as imported inputs.
Internal agents audited the finite-group argument, local conditions,
fixed-period transfer and restoration guards. These reviews are not
formal verification or external review.

No explicit cutoff, exact list of finite admissible exceptions, leading
constant, shrinking-window estimate, uniformity in \(K\), theorem for
even class groups, general discriminants, nonmaximal orders, mixed fields
or simultaneous norm restrictions is established in this note. Those
extensions require their own local, genus and analytic arguments. In
particular, fixed class absorbers do not address the two arbitrary
subgroup convolution problem represented by \(C_{12}+C_{12}\).
