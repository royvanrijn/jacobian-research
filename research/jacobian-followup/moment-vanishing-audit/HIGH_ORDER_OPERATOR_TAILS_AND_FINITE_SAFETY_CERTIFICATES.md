# High operator tails and finite certificates for linear safety

15 November 2026. Internal written theorem, independent review pending.
This is an all-degree transfer for constant-coefficient operators. Unlike
the [principal-pair theorem](FILTERED_PRINCIPAL_PAIR_GVC_TRANSFER.md), it
does not require positive weights balancing the operator and phase.
It assumes a linear mixed cutoff for the starting pair. The
[calculator and small checker](verify_operator_tail_safety.py) verify
rounding, identities and boundary controls; they do not establish the
universal argument from samples. No novelty or formal verification claimed.

## 1. Exact assumptions and terminology

Let A be any commutative Q-algebra, let P in A[x_1,...,x_n] have degree
at most D>=1, and let B and R be constant-coefficient operators whose
symbols may be arbitrary formal power series in derivatives. Such an
operator acts locally finitely on every polynomial. Assume an integer
C>=1 satisfies

    B^j(MP^j)=0 for every polynomial M and j>C deg M.          (1)

Constants M are included, so (1) includes all positive pure identities.
Call this a linear cutoff C. This is stronger than separately knowing
eventual vanishing for each multiplier without a degree bound.

Suppose every monomial of the tail R has derivative order at least s,
and put a=CD+1. The sufficient hypothesis is

    s>D(CD+1)=Da.                                           (2)

The zero tail satisfies every such lower bound. The phase is held fixed;
no phase-perturbation assertion is made.

## 2. The two-inequality argument

Write Lambda=B+R. For a marker Q of degree delta, commutativity gives

    Lambda^m(QP^m)
      =sum_(l=0)^m binomial(m,l) R^l B^(m-l)
           [Q P^l P^(m-l)].                                (3)

The inserted marker QP^l has degree at most delta+Dl. Hence the bracket
after B^(m-l) vanishes if

    m>C delta+a l.                                         (4)

Every constant-coefficient operator is degree nonincreasing. Also
ord(R^l)>=sl, including over rings with zero divisors. Thus the same term
vanishes by degree if

    sl>delta+Dm.                                           (5)

Any nonzero term must satisfy both opposite inequalities. Eliminating l,

    m<=C delta+a l<=C delta+a(delta+Dm)/s,
    m<= [(sC+a)/(s-Da)] delta.                              (6)

Consequently a sufficient full mixed cutoff is

    C_tail=ceil((sC+CD+1)/(s-D(CD+1))).                      (7)

There is no restriction on the number or coefficients of tail monomials.
The formal-symbol case follows from the same finite binomial expansion:
only finitely many derivatives act on the polynomial in each term.
No division occurs in A; the rational arithmetic in (6)--(7) is only
an inequality for integer degrees.

## 3. Pure equality and a convenient bound

For Q=1, any term with l>=1 vanishes. Indeed either m>a l, when (4)
applies, or m<=a l, when sl>Da l>=Dm and (5) applies. Therefore

    (B+R)^m P^m=B^m P^m for every m>=1.                     (8)

Under(1) both sides are zero. Identity(8) is asserted under the starting
cutoff assumption, not as an identity on the full unrestricted coefficient
space. Mixed outputs need not agree.

A convenient integral specialization is

    s>=2D(CD+1) implies C_tail<=2C+1.                       (9)

Indeed the ratio in(7) is at most2C+1/D. The exact rational expression
may provide a better sufficient constant; no sharpness is claimed.

The strict order condition matters. Merely asking s>deg P would be false.
Take

    B=partial_x^2, P=x y^3, D=4, C=1, R=partial_y^5.

The base cutoff follows from x-degree: B^m(QP^m)=0 for m>deg Q. Yet

    (B+R)P=0, (B+R)^2 P^2=2880 y!=0.                       (10)

The tail order5 exceeds phase degree4 but does not meet(2), which requires
s>20. With R=e partial_y^5 over Q[e]/(e^2), the second output is2880ey,
still nonzero. Nilpotent tail coefficients do not remove this obstruction.

## 4. Literal coefficient-ideal consequence

For a finite coefficient family for B,P, let I_C(B,P) be its full ordinary
cutoff ideal. Freely adjoin coefficients of a tail of order>=s satisfying(2),
with no bound on its upper order, and let T denote the resulting polynomial
coefficient ring. Apply the theorem in A=T/I_C(B,P)T. It proves

    I_Ctail(B+R,P) subset I_C(B,P)T.                         (11)

This is an ordinary ideal containment, retaining nilpotents and embedded
components. It is not equality of mixed ideals. Over the original C-cutoff
scheme, every higher tail coefficient is free and the perturbed pair has
cutoff C_tail. Formula(8) holds in that quotient too.

Thus every certified positive finite core generates an infinite family of
positive operators, including tails mixing all spatial directions. A
coordinate grading or principal decomposition is not required. This does
not turn a pure unsafe core into a safe one, give a converse for arbitrary
low jets, or preserve an optimal cutoff.

## 5. Formal linear safety has a finite truncation certificate

Let lambda be a formal symbol and lambda_<=r its ordinary jet through
order r. A formal pair has some integer linear cutoff if and only if there
are finite integers r,C such that

    (lambda_<=r(partial),P) has cutoff C,
    r+1>D(CD+1).                                           (12)

Sufficiency is(7). For necessity, suppose the original formal pair has
cutoff C0. Apply the same theorem to the original operator and the negative
of its tail. For every r with r+1>=2D(C0D+1), the finite truncation has
cutoff at most C1=2C0+1. Choosing

    r=D((2C0+1)D+1)                                       (13)

satisfies that inequality and r+1>D(C1D+1). Thus(12) holds. The proof
does not assume that an arbitrary formal-symbol safe pair has a linear
cutoff in the first place.

Each condition in(12) involves finitely many phase/symbol coefficients
and the ordinary finite-order cutoff ideal I_C. The
[retained algorithm](HOLONOMIC_SOURCES_AND_CUTOFF_LOCI.md) makes that
ideal theoretically computable. Therefore the functor of formal pairs
having a linear cutoff, on arbitrary Q-algebras, is a union of closed
finite-jet cylinders defined by these ordinary ideals. No Noetherianity of
the infinite symbol coefficient ring is used, and no single closed
representative or equality of prescribed cutoff ideals is asserted.

The forward certificate can have a larger cutoff than the actual formal
pair. This finite-certificate theorem does not compute the least formal
mixed slope or decide whether a specified exact slope is valid.

## 6. A theoretical positive semidecision for computable formal symbols

Let all coefficients lie in a fixed effective number field K. Assume a
procedure supplies the coefficient of each formal symbol monomial, and
let P be a supplied polynomial. The following is a positive semidecision
for the existence of a linear cutoff, with purity included:

1. For r=0,1,..., form the finite operator lambda_<=r.
2. Decide its all-order purity and all-multiplier safety using the retained
   [finite-operator algorithms](CONSTANT_COEFFICIENT_OPERATOR_DECISION.md).
   If it is not pure-safe, continue.
3. For a pure-safe truncation, find its least integer cutoff C_r>=1 by
   testing I_1,I_2,... at its coefficient point. The full ordinary cutoff
   ideal algorithm terminates at each test; finite-order safety supplies
   some linear cutoff, so this search terminates.
4. If r+1>D(C_rD+1), return r,C_r and the transferred bound(7).

Every returned certificate is sound by(12). If the formal pair has a
linear cutoff C0, all sufficiently long truncations have cutoff<=2C0+1
by(9). Thus eventually C_r<=2C0+1 and the test in step4 succeeds. Each
earlier finite-operator decision terminates. Outside this positive class,
the procedure may run forever. Coefficients require exact arithmetic and
zero tests; an arbitrary numerical complex oracle would not suffice.

This algorithm is written, not implemented or run. The small implementation
only calculates the tail bounds and validates supplied numerical degree
data. It does not provide the retained general finite-operator/ideal solvers.
An individual nonzero pure output can separately be detected by enumerating
finite expressions, but no total formal-symbol decision follows.

## 7. The remaining formal-symbol question

**OP-GVC-FORMAL-LINEAR.** Over characteristic-zero fields, for a formal
constant-coefficient symbol and a polynomial phase with all positive pure
outputs zero, does eventual vanishing separately for every polynomial
multiplier imply some common linear degree cutoff?

The answer is positive for finite-order symbols by the retained holonomic
source theorem, and for the formal quadratic/univariate classes by the
[finite-seed results](FINITE_SEEDS_AND_FORMAL_OPERATOR_UNIFORMITY.md).
The current proof only characterizes the linear-cutoff class. It does not
settle this implication for general formal symbols, assume their source
series holonomic, or upgrade the positive semidecision to full eventual
safety. Keep this gap separate from the already closed finite-order theorem.

For context, [de Bondt1206.2836v2](https://arxiv.org/html/1206.2836v2)
Theorem1/Corollary2 provide a finite-operator commutator reduction to
shifted powers; the corollary retains a separate eventuality assertion for
each shift. Their proof does not supply the missing formal-symbol uniform
bound used in(1). The latest listed revision remains6February2013.
The current tail proof is elementary and imports no new theorem from it.

The [Hexagon two-variable preprint](https://hexagonmath.org/2610.00580v2)
was checked for current metadata and provenance only: its v2 page is dated
27October2026 and its history records current28October. Its statement is
for finite polynomial operators and carries an explicit incomplete-review
provenance. No new proof or assurance is adopted from that landing page.
The retained two-variable/three-variable endpoints remain unchanged.

## 8. Scope and next structural gates

The result supplies an all-degree, all-ring operator perturbation transfer
and finite certificates for linear formal safety without a weight grading.
It does not deform the phase, prove positivity of a new unrestricted GVC
dimension, repair a retained failure, cover variable coefficients, or cover
formal nonpolynomial phases/multipliers. Independent review remains pending.

Exact controls retain(10), its nilpotent version, binomial operator
expansions, a changed mixed output under an admitted high tail, exact
rational/rounded bounds and finite-core inequalities. No large marked-tail
campaign or general solver is run. Positive-core classification and the
formal eventuality-to-linearity question are structural next gates; refining
these sufficient numerical constants is secondary.
