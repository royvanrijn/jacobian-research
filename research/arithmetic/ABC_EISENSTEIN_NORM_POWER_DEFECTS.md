# Exact Eisenstein norm powers under the original radical caps

Use \(n=X^L\) with height first at fixed \(L\), and
primitive positive rows \(a+b=c\in(X/2,X]\).
Write \(N=a^2+ab+b^2\). For a fixed order \(q\ge2\)
and a total original cap \(\operatorname{rad}(abc)\le n^B\),
the count with \(N=w^q\) has limiting coefficient
\[
 \frac{B}{(2a_0+1)\tau(s)},\qquad
           3q=3^{a_0}s,\quad 3\nmid s.                  \tag{A}
\]
In particular the coefficients at orders \(2,3,4\)
are \(B/6,B/5,B/9\). Uniformly over adaptive orders \(q\ge2\),
the proper-power sector has limiting coefficient \(2B/9\).
These bounds use a total original cap. The further defect
bound below uses three individual original caps.

For a limiting coefficient \(D\), the convention throughout
is a bound \(O_{B,L,\eta,\varepsilon}(n^{D+\eta}X^\varepsilon)\)
for each \(\eta>0\) and sufficiently small fixed \(L\).
In fixed-order claims, the threshold on \(L\) and constants
may depend on the fixed order. All exponent allowances are
chosen before taking \(X\to\infty\).

## EN1. Exact factorization and counting

Let \(\omega^2+\omega+1=0\), put
\(\alpha=a-b\omega\), and use the Euclidean property of
\(\mathbb Z[\omega]\) as in CP4. The ideals of \(\alpha\)
and \(\bar\alpha\) are coprime outside 3: a common prime
ideal would divide both \(a-b\omega\) and \(a-b\bar\omega\),
hence \(b\) and \(a\) at an unramified prime, contrary
to primitivity. If 3 divides \(N\), then \(a\equiv b\bmod3\)
and both are units at 3. The identity
\(N=(a-b)^2+3ab\) gives \(v_3(N)=1\). Thus the
exact-power condition excludes ramification, and
\(\alpha=\epsilon\beta^q\) for an Eisenstein integer
\(\beta=x-y\omega\) and a unit \(\epsilon\).
The conjugate ideals of \(\beta\) are coprime, so
\(\gcd(x,y)=1\). Also \(|\beta|\ll X^{1/q}\) with
an absolute implied constant.

The identity
\[
 \prod_{j=0}^2
  (\omega^j\alpha-\omega^{-j}\bar\alpha)
       = \alpha^3-\bar\alpha^3
       = \pm(\sqrt{-3})^3abc .
\]
Each difference divided by \(\sqrt{-3}=\omega-\omega^2\)
is one of the complete letters with a sign.
Since a unit and its conjugate have the same cube,
\[
              \beta^{3q}-\bar\beta^{3q}
                          =\pm(\sqrt{-3})^3abc.         \tag{B}
\]

First use (B). The \(\Phi_d(\beta,\bar\beta)\), for
\(d\mid3q\), give distinct homogeneous factors, each
fixed at the chosen \(q\). If \(3\nmid d\), the
whole factor is rational after division by a nonzero
constant where necessary (in particular \(d=1\)
is proportional to \(y\)). If \(3\mid d\), split the
roots of \(\Phi_d\) into the two classes
\(j\equiv1,-1\bmod3\). Each product of
\(\beta-\zeta_d^j\bar\beta\) has coefficients in
\(\mathbb Q(\omega)\). Conjugation replaces each factor
by \(-\zeta_d^{-j}\) times itself, so the conjugate of
the product is a constant multiple of the product.
Divide by any nonzero coefficient of that product to
give a rational binary form, and clear denominators.
The subgroup fixing \(\mathbb Q(\omega)\) is transitive on
each set of roots and on all roots when \(3\nmid d\).
Thus these rational forms are irreducible. Their
geometric root sets are disjoint, as are those for different
\(d\), since the change from \((x,y)\) to
\((\beta,\bar\beta)\) is invertible over the algebraic
closure. Thus all resultant primes and normalization
primes can be put in a fixed finite pool. Away from
that pool, primitivity of \((x,y)\) makes the supports
of the factor values pairwise disjoint, and the union
of the supports divides the original \(abc\). The
number of these factors is \(T=(2a_0+1)\tau(s)\).

There are three linear factors for odd \(q\), from
\(d=1,3\), proportional to \(x,y,x+y\). For even \(q\)
the indices \(2,6\) add \(2x+y,x-y,x+2y\). Every
three of the resulting distinct directions yield an
additive primitive triple after applying their fixed
linear relation, signs, and bounded common gcd.
Its height is comparable to \(X^{1/q}\), since any
two of the fixed independent directions control the
height of \(\beta\). It recovers \((x,y)\) at bounded
cost and its radical cap is the sum of the three
form caps up to fixed primes. Choose exact radical
bins at each fixed budget, and apply MCB2 or the
one-third bound in the relative budget \(qL\). Thus
each three linear forms have combined normalized
radical cap at least \(3d\) at row count exponent \(d\).
Average over all such triples to get \(k_{\rm lin}d\)
as the lower bound for the total linear cap.

A fixed exact value of each nonlinear factor determines
\((x,y)\) at \(X^\varepsilon\) cost. After a fixed
change of coordinates and clearing fixed denominators,
the irreducible form is a constant times the norm of
a linear expression in a fixed field. At a fixed value,
its integral ideal divides a fixed multiple of that value
and has subpower many possibilities. Every embedding
of a candidate generator has polynomial height. The
norm is a nonzero integer, giving reciprocal polynomial
lower bounds for its embeddings. For each ideal the
logarithmic unit lattice therefore gives polylogarithmically
many generators in this range. An exact generator
recovers the two coordinates since the degree is at
least two. Consequently each
nonlinear factor costs normalized radical at least
\(d\), by supported-power enumeration of its exact
polynomial-height values. Summing the linear and
nonlinear charges gives \(B\ge Td\) and proves (A).

## EN2. Uniformity over exact power orders

If \(q\) has a fixed divisor two or three, take
that divisor. Otherwise take a prime \(p\ge5\)
dividing it. The unramified \(p\)-th power case has
\(T=6\), including exactly three linear forms. Even
without applying an abc count at relative budget
\(pL\), fixing any two linear values recovers
\((x,y)\), so the sum of the three linear caps is
at least \(3d/2\). Add the three nonlinear charges
\(3d\) to get \(d\le2B/9\). This uses no smallness
assumption on \(pL\). At each fixed budget one
need apply this argument only for \(p\le P(L)\),
with a finite cutoff large enough in terms of
\(B,L\) and the desired allowance. For larger
orders, enumerate \(x,y\ll X^{1/p}\) directly
at exponent \(2/p\), which is at most the desired
fixed positive multiple of \(L\) for this cutoff.
There are \(O(\log X)\) possible orders for the
nontrivial norms. For the count up to \(X\), sum the
shell estimates with any fixed allowance in the exponent.
The constants for the finite
order list are taken at the fixed-budget height
limit.

## EN3. Retaining a small radical of norm defects

For \(q\ge2\), put
\[
 D_q(N)=\prod_{\substack{p\mid N\\q\nmid v_p(N)}}p .
\]
Under the three individual caps \(n^{1/3+o(1)}\) at the
balanced small-budget limit, a bound \(D_q(N)\le n^f\)
with fixed \(f<1/6\) eliminates saturation, uniformly
over all choices of the order \(q\). For each fixed
\(f<1/3\) there is also a fixed \(Q_*(f)\) such that
\(D_q(N)\le n^f,\ q\ge Q_*(f)\) eliminates saturation.
Concretely each statement provides positive constants
\(\delta,\kappa,L_0\), depending on \(f\), such that the
count with three caps \(n^{1/3+\kappa}\) and
\(0<L<L_0\) is \(O_{f,L,\varepsilon}
(n^{1/3-\delta}X^\varepsilon)\).

Here is the arithmetic count used in both assertions.
A prime dividing \(abc\) cannot divide \(N\), by
primitivity and the congruence obtained by setting
the corresponding letter to zero. Split the ideal
of \(\alpha\) into quotient and remainder exponents
modulo \(q\). The UFD gives
\[
                          \alpha=\delta_0\beta^q,
\]
with an exact integral coefficient \(\delta_0\)
including the unit choice. Its support is over
\(D_q(N)\), since each split rational prime selects
just one of the conjugate primes and the ramified
exponent is at most one. For a fixed defect radical
there are height-subpower choices of the orientations
and all exact coefficient powers, uniformly up to
polynomial height. Thus coefficients cost
\(n^fX^\varepsilon\). The root is primitive and
\(|\beta|\ll X^{1/q}\).

Choose a full divisor \(D\) from any one original
letter, with \(X^w\le D\) and radical at most \(n^\mu\).
Its choices cost \(n^\mu X^\varepsilon\).
The quotient \(\alpha/\bar\alpha\) is a specified
cube root of unity modulo \(D\mathbb Z[\omega]\).
For the original positions \(a,b,c\) these are
\(\omega^2,1,\omega\), respectively. Both
\(\delta_0\) and \(\beta\) are units modulo \(D\),
so \((\beta/\bar\beta)^q\) is fixed modulo \(D\).
For bounded \(q\), the number of possible unit ratio
classes is height-subpower. At unramified primes it
is bounded locally by a constant depending on \(q\)
and the field, using cyclic finite-field units
and the principal unit filtration. At the finitely
many primes dividing \(3q\), the number of \(q\)-torsion
classes stays bounded as the modulus depth grows, by
the local logarithm on a subgroup of fixed finite index.
If two roots have the same class, their cross-product
difference \(\beta\bar\beta'-\beta'\bar\beta\) belongs to
\(D\mathbb Z[\omega]\) and has absolute value
\(O(X^{2/q})\). For \(qw>2\) by a fixed margin the
difference vanishes. The roots are rationally proportional
and primitive, hence equal up to sign. This proves the
count with exponent
\[
                             f+\mu+\eta .                \tag{C}
\]
The arbitrary normalized allowance \(\eta>0\) absorbs
adaptive orders above a finite cutoff at each fixed
\(L\): enumerate \(\beta\) in a box of height
\(O(X^{1/q})\) when \(2/q\le \eta L\). Nonunit roots
require \(q=O(\log X)\). A unit root has its whole
original norm support in the coefficient, and can
be counted directly at cost \(n^fX^\varepsilon\),
without enumerating unbounded orders.

There is a second count for each fixed \(q\), without
requiring \(qw>2\). Write the signed original letters
as degree-\(q\) homogeneous integer forms in the two
coordinates of \(\beta\). Their root sets are the three
distinct sets obtained by requiring
\((\beta/\bar\beta)^q\) to equal a fixed cube root times
\(\bar\delta_0/\delta_0\). Their roots are simple and
the three sets disjoint outside primes dividing
\(3q\,N(\delta_0)\). The latter norm is automatically
coprime to the bought divisor. Omit the fixed primes
dividing \(3q\) when choosing \(D\).
Also \(|N(\delta_0)|\le n^{(q-1)f}\), so the form
coefficients satisfy the polynomial-in-\(n\) bound
at fixed \(q\).

Let \(S\) be either other role and let \(P_S\) be
its form. Of the \(q+1\) charts with first column
\((i,1)\), \(0\le i\le q\), and second column \((1,0)\),
at least one has leading \(P_S\) coefficient a unit
at each prime of \(D\) larger than \(q\). Omit the
other fixed small primes at purchase. At least one
chart retains product of assigned full prime powers
of size at least \(D^{1/(q+1)}\). For fixed \(D\) and the
coefficient, the assignment is deterministic using
the first eligible column. On the retained submodulus
the bought form has at most \(q^{\omega(D)}\)
projective root classes, by simple-root Hensel lifting
in the two local coordinate charts. The tested value
\(P_S\) is a unit by original primitivity. Apply BR1,
after fixed sign changes if needed to make the
root coordinates positive. Primitive inputs with a
zero coordinate have only bounded possibilities.
Writing \(r_S\) for the original radical cap on \(S\)
in base \(n\), the normalized bound is
\[
                             f+\mu+r_S/2+\eta .          \tag{D}
\]
The choice of \(D\) is charged once. Its retained
submodulus depends on the same records.
[Binary-form ray estimate BR1](ABC_BINARY_FORM_RAY_AMPLIFICATION.md#br1-primitive-inputs-with-prescribed-projective-residues)

To deduce the uniform claims, argue by a sequence
at balanced saturation with vanishing positive
cap tolerance and vanishing exponent loss. On a
limiting profile GE1 supplies a positive compact
nonatomic window on some original letter. It
supplies full divisors of arbitrarily small fixed
price \(\mu\) and positive fixed depth \(w\) after
refining. Use the window away from the finitely
many bad primes at bounded \(q\). For \(f<1/6\),
any sequence of bounded orders has a constant-order
subsequence to which (D) applies with positive
slack. Any divergent order subsequence satisfies
\(qw>2\) with slack and (C) applies. For \(f<1/3\)
and orders required to be at least \(Q_*\), failure
for every \(Q_*\) would give a divergent order
sequence and the same contradiction. In this
argument one can first pick the window away from
any fixed finite prime support in the bounded
order case: its atom-free mass can be retained
outside the boundedly many large-price primes,
and small-price primes with fixed rational
names have vanishing total raw price. Counting
orders at each fixed \(L\) uses the finite cutoff
in (C) before the height limit.
[GE1 and the positive nonatomic window](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

## EN4. The cyclic pair count

For a fixed original row the primitive cyclic determinant
output recovers a second row at subpower cost. The
unreduced determinant pair is the image of the two
primitive coordinates of row two under an integer matrix
of determinant \(N_1\). Its content \(g\) divides \(N_1\)
by the adjugate and has subpower choices. The bounded
signed-output choices and exact \(g\) then recover the
second row. Exclude the degenerate zero-coordinate
output, which here means identical oriented rows and
costs at most the first-row count. The primitive output norm is \(N_1 N_2/g^2\).

If the first-row catalogue costs \(n^{d_0}\), the
primitive output has height in \([X^\nu,X^B]\) for
fixed positive \(\nu,B\), its radical cap is \(n^s\)
and its norm is an exact fixed \(q\)-th power, the
pair exponent is at most \(d_0+s/T+\eta\), with \(T\)
from (A). To apply EN1 at each fixed input \(L\),
divide the output into height shells and round their
relative budgets upward to a finite mesh of fixed
budgets covering the possible relative budgets from
shell depths in \([\nu,B]\). Take the mesh fine enough
to absorb its error in \(\eta\) before the height limit.
For adaptive exact power output the corresponding
bound is \(d_0+2s/9+\eta\). Without a lower height
floor one can use \(s/(T-k_{\rm lin}/2)\) in the fixed
order bound, since any two independent linear values
recover a root pair. The criterion for beating the
original full pair exponent is that the stated
pair price, including any additional records,
sits strictly below \(2/3\). An endpoint extraction
forcing this norm condition remains open.
[Cyclic determinant transfer CP4](ABC_CYCLIC_PROJECTIVE_PAIR_TRANSFER.md)

## EN5. Joint original divisors, a real strip and a norm ideal

Keep the notation of EN3 and the original three caps. Fix the
order \(q\), or allow it to vary up to the finite-order cutoff
in EN3 at fixed \(L\). The exact coefficient
\(\alpha=\delta_0\beta^q\) costs \(n^{f}X^\varepsilon\) when
\(D_q(N)\le n^f\). The following records can be supplied jointly
with the coefficient and priced by the number of joint choices.

* Prescribe one real interval of length \(O(X^{-\sigma})\) for
  \(a/c\), with \(\sigma\ge0\). An available interval cover is
  charged at its actual price. In particular one interval
  suffices when \(a\le X^{1-\sigma}\), when \(b\) satisfies
  this bound instead, or when \(|a-b|\le X^{1-\sigma}\).
* Prescribe actual numerical divisors \(D_a\mid a\), \(D_b\mid b\),
  \(D_c\mid c\), and put \(D=D_aD_bD_c\). This costs
  \(n^\mu X^\varepsilon\) if the joint radical is at most
  \(n^\mu\), assigning its primes to their roles as part of the
  record at height-subpower cost.
* Prescribe a rational divisor \(E\mid N(\beta)\). Once
  \(D_q(N)\) is fixed, this has cost \(n^\xi X^\varepsilon\)
  if its radical outside that support is at most \(n^\xi\).
  Retaining the entire full power on an already recorded
  prime costs only height-subpower.

If for a fixed positive margin \(\nu\) the records satisfy
\[
       q\bigl(\log_X D+\sigma+\log_X E\bigr)
                           +\log_X N(\delta_0)\ge 2+\nu,
                                                        \tag{EN5}
\]
then the count is at most \(n^{f+\mu+\xi+\eta}X^\varepsilon\)
for any fixed allowance \(\eta>0\). Add the price of a moving
interval cover if used. More generally use the actual joint
record price in place of \(f+\mu+\xi\). The height error bounds
are uniform in record values at the fixed \(L,q\). Orders
above \(2/(\eta' L)\), with \(0<\eta'<\eta\) fixed, can be counted
at \(n^{f+\eta}X^\varepsilon\) by the EN3 direct root enumeration.
This makes the bound valid for adaptive \(q\) with the usual
fixed-budget order of quantifiers.

The three moduli are pairwise coprime by primitivity. The
corresponding residue of \(\alpha/\bar\alpha\) is
\(\omega^2,1,\omega\), respectively, and each is a unit
because \(\gcd(N,abc)=1\). Fix their joint CRT unit class
for \((\beta/\bar\beta)^q\) after dividing by
\(\delta_0/\bar\delta_0\). The possible unit ratio classes
of \(\beta/\bar\beta\bmod D\) cost height-subpower by EN3,
including primes dividing \(3q\).

On the real interval, \(\alpha/\bar\alpha\) varies over an
arc of length \(O(X^{-\sigma})\). In fact it is a fixed smooth
fractional linear function of \(a/c\) with denominator bounded
away from zero on \([0,1]\). Dividing by \(\delta_0/\bar\delta_0\)
and extracting the \(q\)-th root restricts \(\beta/\bar\beta\)
to at most \(q\) short arcs of length \(O(X^{-\sigma}/q)\).
For \(\sigma=0\) one can just use the unrestricted unit circle.

An ideal \(\mathfrak e\) dividing \((\beta)\) with norm \(E\)
is specified by the orientations of the split primes on \(E\),
at height-subpower cost. The prime over 3 has zero exponent
on \(\beta\), and an inert prime would divide both coordinates.
For two roots in the same class and real arc, with the same
\(\mathfrak e\) and coefficient, the determinant
\[
                           \beta\bar\beta'-\beta'\bar\beta
\]
belongs to \(D\mathbb Z[\omega]\) by the unit congruence. It
also belongs to \(\mathfrak e\overline{\mathfrak e}
=E\mathbb Z[\omega]\). Since \(D,E\) are coprime, a nonzero
determinant has absolute value at least \(DE\). Its real
size is bounded by
\[
                 O\!\left(X^{2/q-\sigma}N(\delta_0)^{-1/q}\right)
\]
using \(|\alpha|\ll X\) and the arc bound on the unit ratio.
Condition (EN5) forces it to vanish. The primitive positive
row is recovered at bounded cost since the two roots are
rationally proportional and primitive. This argument uses the
original primitive roots throughout, also on families selected
by a least-repair height cap.
[Related slope count for weighted powers CPW7](ABC_COPRIME_SHORT_POWER_WINDOWS.md#finite-height-slope-counts-for-unequal-orders)

## EN6. The optimized depth and the already recorded norm part

Let \(N_{\rm def}\) be the full part of \(N\) supported on
\(D_q(N)\). Fix its exact values on that already recorded
support and take \(E=(N_{\rm def}/N(\delta_0))^{1/q}\)
in (EN5). This is an integer and is an available divisor
of \(N(\beta)\). It has \(\xi=0\). The depth condition is
\[
                    q(\log_X D+\sigma)+\log_X N_{\rm def}
                                      \ge 2+\nu.        \tag{EN6}
\]
Any further full part of \(N\) away from the recorded support
with valuations divisible by \(q\) can be used through its
\(q\)-th root in \(E\); add its raw price if it has to be bought.

Here is a consequence uniform over balanced profiles. For
every fixed \(f<1/3\), there exist \(\kappa,L_0,\delta_1,\theta>0\)
such that, at all fixed \(0<L<L_0\), the count with three
caps \(n^{1/3+\kappa}\), \(D_q(N)\le n^f\) and
\(N/N_{\rm def}\le X^\theta\) for some \(q\ge2\) satisfies
\(O_{L,f,\varepsilon}(n^{1/3-\delta_1}X^\varepsilon)\).
Otherwise take a height-first saturation sequence with
\(\theta\to0\) and \(\kappa\to0\). By GE1, a limiting balanced
profile gives a positive compact nonatomic band on some
original letter. Select a full block there of arbitrarily small
fixed price \(\mu<1/3-f\) and fixed positive depth \(w\).
Since \(N\asymp X^2\), inequality (EN6) holds with slack, already
for the lowest possible order \(q=2\) and \(\sigma=0\).
Count by EN5 and use a fixed allowance smaller than \(1/3-f-\mu\).
This contradicts saturation and proves the claimed uniform
threshold, with no maximum price condition on individual
original primes.
[Positive nonatomic windows GE1](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

For a second consequence, sort the available original full
prime powers by decreasing \(L v_p(abc)\), after excluding any
prepaid powers. Write \(R\) for their total normalized raw
price and \(W\) for their total numerical weight in base \(X\).
If the maximum price of one power is \(\varkappa\), any
target depth \(0\le v<W\) can be reached at price at most
\(Rv/W+\varkappa\). Indeed the top fraction of the divisible
knapsack has average efficiency at least \(W/R\). Round
the last power of its initial segment up to a whole power.
Prepaid original and norm weights and the known real width
reduce the target \(v\) through (EN5). The actual joint
divisor catalogue can improve this raw price.

In particular suppose the three balanced original caps are
\(1/3+o(1)\) in base \(n\), and the largest individual price
on any original letter tends to zero in the fixed-budget
height-first limits. Then \(W=2+h+o(1)\) and \(R\le1+o(1)\),
where \(h\in[0,1]\) is the limiting numerical short-summand
height. Use \(\sigma=1-h\) to an arbitrarily small fixed
allowance, and discard the free norm weight for this estimate.
The required price is at most
\[
        \frac{(2/q+h-1)_+}{2+h}+o(1)
                                  \ \le\frac{2}{3q}+o(1).
\]
Thus for fixed \(q\ge5\), \(f<1/3-2/(3q)\) with a fixed
margin gives a strict saving on this sector, for every height
profile. At order two the more precise expression improves the
\(1/6\) threshold of EN3 for \(h<2/5\), and at order three for
\(h<4/5\). The strictly short range \(h<1-2/q\) by a fixed
margin needs no modulus and has count at price \(f\) with no
individual prime assumption. At varying heights, bin the
normalized short height into finitely many intervals of
sufficiently small constant width before taking the height
limit. At adaptive orders use the bounded-order arguments for
a finite range at fixed \(L\) and the direct enumeration for
large orders, as above.

A fractional version that buys part of a single prime's
support can give an unrealizable depth. As a formal example,
on each role place one prime of price \(1/6\) at a scaled
degree close to 6, and diffuse mass \(1/6\) at positive
scaled degrees close to zero. This has balanced total masses,
three heights tending to one, and nonatomic total \(1/2\).
A price strictly below \(1/6\) can then have arbitrarily small
numerical depth because it purchases none of the three
macroscopic powers. These are formal valuation data. On
genuine families (EN5) uses whole prime powers, actual
individual prices and the actual norm defect support.

## EN7. Fresh cyclic radicals for pairs with matching defects

Fix \(q\ge2\) before the inner height limit. Work with an
actual family of positive primitive abc rows in \(X/2<c\le X\)
and \(n=X^L\). Suppose its size is at most \(n^{d+o(1)}\).
The original three radical caps, and any additional
eligibility restrictions, remain in force. Put
\(\alpha_i=a_i-b_i\omega\), \(\omega^2+\omega+1=0\).
Use the exact cyclic determinants and their common content
\(g\) of CP1–CP4, so their Eisenstein product is
\(\alpha_1\bar\alpha_2=-\Delta_1+\Delta_0\omega\).
Restrict to distinct ordered pairs with the same exact
defect coefficient \(\delta\) in
\(\alpha_i=\delta\beta_i^q\), using the same generator
and unit convention. It is enough in what follows to
assume directly that the primitive cyclic output has
Eisenstein norm an exact \(q\)-th power, and hence is
a \(q\)-th power up to a unit by EN1.

Set
\[
 k(q)=
 \begin{cases}
    \tau(q)-1,&3\nmid q,\
    \tau(q/3^{v_3(q)}),&3\mid q .
 \end{cases}
\]
Fix one of the three cyclic positions, or allow its choice
to vary over the three. Let \(s\) bound the radical price
of that primitive signed cyclic value outside all six
original letters and any prescribed polynomial-height pool.
Then the pair count exponent \(E\) satisfies
\[
                          E\le d+\frac{s}{2k(q)}+\eta
\]
for every fixed positive normalized allowance \(\eta\).
There is no numerical lower bound on the primitive
cyclic output. This count is for fixed orders or a fixed
finite list at a fixed budget, with inner constants
depending on the list.

First verify that a same-\(\delta\) pair supplies the
required exact cyclic power. Put
\(\gamma=\beta_1\bar\beta_2\), and let \(g_0\) be the
ordinary content of its two coordinates. The Eisenstein
product is \(N(\delta) g_0^q(\gamma/g_0)^q\).
The \(\beta_i\) are units at the prime over 3, because
3 has norm exponent at most one on each primitive original
row and all of that exponent is supported in its defect
record at order \(q\ge2\). Away from 3, primitiveness of
\(\gamma/g_0\) means that it is nonzero in at least one
factor of the semisimple residue ring modulo the rational
prime; raising to a power preserves that property.
At 3 it is still a unit. Thus
\((\gamma/g_0)^q\) is primitive, so the cyclic content is
exactly \(N(\delta)g_0^q\).

For a primitive \(\beta=x-y\omega\), the three coordinates
of \(\beta^q\) in the cyclic additive representation are
proportional to
\[
                   \beta^q-\omega^j\bar\beta^q,
                \qquad j=0,1,2,
\]
after the corresponding bounded choices of units and common
sign and multiplication by a fixed scalar in the field.
The scalar normalization does not affect the radical
outside a fixed finite set. More explicitly work with
\(\omega^i\beta^q-\omega^{-i}\bar\beta^q\) and divide by
\(\sqrt{-3}\) to get an integer signed coordinate.

In the unshifted phase, the factors
\(\Phi_e(\beta,\bar\beta)\), \(e\mid q\), are disjoint in
geometric roots. All \(e\ge3\) give squarefree factors
of degree at least two. If \(q\) is even, combine the
two independent linear factors from \(e=1,2\) into one
squarefree quadratic. If \(q\) is odd, discard its
single polynomial linear factor. The result is
\(\tau(q)-1\) reconstructing factors.

If \(3\nmid q\), a shifted phase contains one of the
two rational halves of \(\Phi_{3e}(\beta,\bar\beta)\)
for each \(e\mid q\). Its degree is \(\varphi(e)\).
There is one linear factor for odd \(q\) and two
for even \(q\); discard the former or combine the
latter into a squarefree quadratic. This again gives
\(\tau(q)-1\) reconstructing factors. If \(q=3^v q'\)
with \(v\ge1\), either shifted phase contains one
rational half of \(\Phi_{3^{v+1}e}(\beta,\bar\beta)\)
for each \(e\mid q'\), all of degree at least three.
This gives \(\tau(q')\) reconstructing factors, and
the unshifted phase has at least that many.
The rational half factor construction is the two
\(j\bmod3\) conjugation argument of EN1. In particular
\(\Phi_3(\beta,\bar\beta)=3x(x+y)\) is allowed as a
squarefree reconstructing quadratic although its
individual linear factors are rational.

Fixing a nonzero exact value of a reconstructing factor
recovers \((x,y)\) at height-subpower cost. Use a fixed
norm representation count for the irreducible factors,
and, for a product of independent linear factors,
factor the exact value into their two values at divisor
cost. The fixed signed primitive cyclic output then
follows. All selected factor supports are pairwise
disjoint outside their fixed resultants and normalization
primes by primitiveness.

To price the old support globally, label each pair by
its exact cyclic content \(g\) and the squarefree radical
\(S\) of its shared original prime support. The identity
and adjugate calculation of EN4 show \(g\mid N_1\)
and \(g\mid N_2\). Also \(S\) divides both original
products \(abc\). Thus every original row is incident
to only divisor-many possible joint \((g,S)\) labels
at polynomial height. If a chosen pair family has
exponent \(E\), discard labels with either projection
smaller than \(n^{E-d-\eta'}\); their pairs number at
most \(n^{E-\eta'+o(1)}\) by summing opposite
incidences. Remaining joint labels have catalogue
price at most \(2d-E+\eta'\).

If an original prime divides a cyclic determinant
then it divides some original role of both rows.
Indeed a prime of any role is one of the three
projective boundary directions of that primitive
row modulo the prime; vanishing determinant aligns
the two directions, and the cyclic rotation
permutes the three boundary directions. This holds
also at 2. Division by content does not introduce
primes. Thus all old determinant primes are already
in the fixed shared support \(S\) after the joint
record is selected. Under the cap \(s\), one of
the \(k(q)\) reconstructing factors has radical
outside the now fixed pools at most \(n^{s/k(q)}\).
Enumerate its exact polynomial-height value at that
price, recover the primitive cyclic output as above,
and use the fixed content \(g\) with CP4 to recover
both rows. Thus \(E\le2d-E+s/k(q)+\eta\), as claimed
after shrinking subsidiary allowances. The argument allows every positive primitive output
of polynomial height. Count distinct oriented input rows.

As one consequence, in a family of exponent \(d\)
with exact defect coefficient dictionary of price
\(f<d\), Cauchy–Schwarz supplies \(n^{2d-f-o(1)}\)
off-diagonal same-defect pairs. Any subfamily
retaining this exponent therefore requires each
primitive cyclic value's fresh radical price
at least \(2k(q)(d-f)-o(1)\). For every fixed
margin the subfamily below that floor is power-smaller
than the guaranteed equal-defect pair count. The original
coefficient dictionary can be priced by its norm
defect radical as in EN3, without also paying it
inside the joint \((g,S)\) pair record.

The same record argument also gives a joint version.
If the radical of the product of all three primitive
outputs outside the six original supports and prescribed
pool is at most \(n^s\), replace \(k(q)\) by
\(T-k_{\rm lin}/2\) in the pair bound. In each fixed
joint label, the factor accounting of EN1 works
outside the now prescribed pools. Each fixed nonlinear
factor value gives height-subpower fibres. Any two
independent linear factor values suffice, so the
linear charge is \(k_{\rm lin}/2\) times the common
count exponent, by averaging over pairs. This avoids
requiring a positive numerical lower bound on the
output height. A proposed pair extraction must still
supply a sufficiently small actual fresh radical.
[Rich-record pricing for product transfers](ABC_FIELD_FREE_CROSS_PRODUCT_TRANSFERS.md#pricing-shared-support-by-its-actual-pair-incidence)
