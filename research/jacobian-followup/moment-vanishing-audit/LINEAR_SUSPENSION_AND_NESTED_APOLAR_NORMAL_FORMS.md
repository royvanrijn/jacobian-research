# Linear suspension and nested apolar normal forms

16 November 2026. Internal written all-degree proofs; no novelty,
independent replay, formal verification or external review is claimed.
See the [checker](verify_linear_suspension.py),
[receipt](linear_suspension_verification.json) and
[source audit](linear_suspension_source_audit.json).

The full formal eventuality-to-linearity problem can be represented on
phases whose ACTUAL power-annihilator ideals are nested. Add one passive
linear variable; this preserves the full qualitative and linear-gap
questions, gives exact ordinary pure-prefix ideals, and preserves existence
of the newly retained scalar time equation by an explicit gauge.

It is false that the original power-annihilators are always nested. We give
an integer family in EVERY degree d>=3, already in TWO variables and with
a PURE SAFE operator of sharp cutoff1/d. Universal nesting holds for all
nonzero phases exactly in the unary or degree<=2 ranges. Thus the lift is
a genuine replacement, not a proof of the false original nesting premise.
The positive-gap question for the resulting nested family remains OPEN.

## 1. One linear variable preserves the entire formal problem

Let A be ANY commutative Q-algebra, P in A[x], and
T=lambda(partial_x) for a formal symbol with constant coefficient EXACTLY0.
Introduce one independent variable u, let T ignore u, and put

    P_plus=P+u.

The symbol, its order when finite, and its original coefficients do not
change. For degP>=1 the phase degree does not change either. The exact
zero-constant assumption is automatic for a nonzero pure phase over a
field, but is an explicit assumption in the all-ring assertion here.

**GVCLINEARSUSP20261116.** The following hold over every such A.

1. EVERY ordinary pure-prefix coefficient ideal for(T,P_plus) is EXACTLY
   the original one, at the SAME prefix. In particular purity is equivalent.
2. All-polynomial-marker eventuality is equivalent in the two pairs,
   WITHOUT a purity premise.
3. A sufficient common real cutoff C>=0 for(T,P) gives cutoff C+1 for
   (T,P_plus). A cutoff for the lifted pair is a cutoff for the original.
   Thus linear-cutoff existence is equivalent.
4. In the universal zero-constant coefficient ring, for ordinary mixed
   cutoff ideals, with ALL output coefficients retained,

       I_(C+1)(plus) subset I_C(base) subset I_C(plus).       (1)

These are ordinary containments, with all nilpotents retained. No infinite
coefficient-ring ACC or numerical finite pure gate is asserted.

Proof. The identity at every m>=1 is

    T^m(P+u)^m=sum_(j=0)^m binom(m,j)u^(m-j)T^mP^j.        (2)

The term j=0 is0 because lambda(0)=0. For1<=j<=m,
T^mP^j=T^(m-j)(T^jP^j), whose output coefficients belong to the ordinary
pure-prefix ideal through m. Conversely the coefficient of u^0 in(2) is
T^mP^m. This proves the exact prefix-ideal assertion, including over rings.

For eventuality fix an ORIGINAL marker Q(x). There is N with
T^j(QP^j)=0 for every j>=N. These terms in the marked version of(2)
also vanish after the extra T^(m-j). For j<N there are finitely many fixed
polynomials QP^j, each killed by all sufficiently high T powers, since T
lowers ordinary polynomial degree by at least1. Hence the lifted marked
sequence eventually vanishes. A general marker in x,u is a finite sum of
u powers times original markers. Reflection follows by taking its u^0
coefficient, using markers independent of u. This proves assertion2.

For assertion3 let delta=degQ(x) and m>(C+1)delta. Terms with j>Cdelta
in(2) vanish by the original cutoff followed by T^(m-j). For1<=j<=Cdelta,
purity and the formal commutator theorem give

    T^(j+delta)(QP^j)=0,

and m>j+delta, so these terms vanish too. For j=0 local nilpotence gives
T^mQ=0 because m>delta. Multiplication of the marker by u^b does not
affect the action and only increases its total degree. This proves C+1.
Reflection is again the u^0 coefficient. In the universal quotient by
I_C(base) the same proof applies, because that quotient imposes the entire
pure ideal as well as the specified cutoff. This gives the first ordinary
containment in(1); coefficient extraction gives the second.

We do NOT identify least marked cutoffs. The C+1 transfer is sufficient.
The zero-constant assumption cannot be dropped: for T=1,P=0 the original
pair has all positive outputs0, but its suspension has output u^m.

## 2. Actual nested ideals and exact running heights

Now let K be a characteristic-zero field, P!=0, and lambda(0)=0. Define
h_j=max{k>=0:T^kP^j!=0}, including h_0=0, and

    H_m=max_(0<=j<=m) h_j.

No purity is needed for these heights: T is locally nilpotent on each
polynomial. Distinct u powers in(2), with T^k replacing T^m, show that

    h_m(plus)=H_m.                                         (3)

The ideal identity in the ACTIVE symbol subring is also exact:

    Ann_(K[[Z_x]])((P+u)^m)=intersection_(j=0)^m Ann(P^j).   (4)

In the FULL lifted symbol ring, the actual ideals are nested too, since

    partial_u(P+u)^m=m(P+u)^(m-1).                         (5)

Every formal constant-coefficient operator commutes with partial_u, so
Ann((P+u)^m) is contained in Ann((P+u)^(m-1)). The active contraction(4)
is a genuine running intersection, not an assumed property of the original
single-power ideals.

On the original pure locus put b_j=j-h_j for j>=1 and b_0=0. The lifted
deficit is

    b_m(plus)=min_(0<=j<=m)(m-j+b_j).                       (6)

It follows that b_m tends to infinity IFF b_m(plus) does. For the forward
direction, sufficiently large j have b_j above any prescribed bound, and
the finitely many smaller j are handled by m-j tending to infinity.
Reflection follows from b_m(plus)<=b_m.

Moreover the EXACT infimum gap and the asymptotic gap agree:

    inf_(m>=1) b_m(plus)/m=inf_(m>=1) b_m/m,
    liminf b_m(plus)/m=liminf b_m/m.                        (7)

Indeed H_m/m>=h_m/m and H_m/m<=sup_j h_j/j. Taking the supremum proves
the first equality. For the second, a fixed finite initial segment has
height divided by m tending to0, while the remaining j<=m satisfy
h_j<=j(rho+epsilon)<=m(rho+epsilon), where rho=limsup h_j/j.
Together with H_m>=h_m this proves equality of the two limsups, hence(7).

Thus ANY hypothetical pure eventual-but-nonlinear pair gives one with
nested actual power-annihilators, nondecreasing heights, and a phase with
a constant derivative. A linear-gap theorem for this lifted class would
settle the full original implication. No such theorem is supplied here.
Monotone heights alone do not prove the required positive gap.

## 3. The actual reversed source changes by one elementary gauge

Use the ACTUAL proper reversed source from
[GVCREVERSEDHEIGHT20261116](SCALAR_TIME_EQUATIONS_AND_ACTUAL_APOLAR_GAPS.md),
retaining all independent outputs. On purity,

    Psi(t,z,x)=sum_(m,k) T^kP^m t^m z^(m-k)/(m!k!).

Direct binomial expansion gives the literal formal identity

    Psi_plus(t,z,x,u)=exp(t z u) Psi(t,z,x).                 (8)

At a fixed t^m coefficient all k sums are finite, and purity makes every
z exponent nonnegative. Rearrangement in(8) is therefore proper, not a
substitution into an undefined Laurent series. The gauge and its inverse
both have the linear support cone time degree=source degree.

**GVCREVERSESUSP20261116.** Equations(3)-(8) preserve the entire actual
gap alternative. Existence of ONE scalar time equation with coefficients
formal in z and polynomial in time is EQUIVALENT for the two reversed
sources, over K(x) and K(x,u), respectively. Thus the weaker finite-equation
gate may also be pursued on the nested lifted normal form without changing
the question.

In the forward direction replace partial_t in a base operator by
partial_t-zu. In the reverse direction replace it in a lifted operator
by partial_t+zu. Both are invertible scalar Weyl gauges, keep finite time
powers/derivative orders, and give a nonzero equation for Psi over
K(x,u)[[z]][t].

We justify descent of this last equation rather than assuming u=0 is a
valid specialization. Formal z coefficients can have unbounded orders of
u poles. Expand each rational function of u coefficientwise at u=0.
Choose an integer q such that the K(x)-linear functional

    ell:K(x,u)->K(x), r(u)->[u^q]r(u)

does not kill one specified nonzero coefficient of the gauged operator.
Apply ell coefficientwise in z to all its finitely many operator positions.
Since Psi is independent of u, ell commutes with multiplication by every
coefficient of Psi and with the finite coefficient sums in the equation.
This yields a NONZERO equation in K(x)[[z]][t]. No common lower Laurent
u bound over all z degrees is required. This proves reflection.

Under eventuality a supplied equation gives a positive height fraction by
MVTIMEINDICIAL20261116. The same HEIGHT-derived sufficient constant works
before and after suspension by(7). This does not assert equality of least
marked cutoffs. Nor does(8) produce a general finite equation: scalar-time
finiteness is preserved, not assumed or established for arbitrary symbols.

## 4. A linear-phase normal form with a precisely constrained operator

The polynomial shear v=u+P(x) makes the suspended phase the coordinate v.
The active constant derivations become

    delta_i=partial_(x_i)+(partial_(x_i)P)partial_v.         (9)

They commute, and each is locally nilpotent. More strongly, every formal
lambda(delta) acts locally finitely on each polynomial: conjugate back to
f(x,u+P(x)), a polynomial with finite total x degree. The entire joint
derivative action above that degree is0. Hence there is no undefined formal
infinite differential action in this normal form.

Purity and all-marker eventuality are exactly transported by this
polynomial automorphism. If D=max(1,degP), both the shear and its inverse
increase total marker degree by at most D. A sufficient cutoff C therefore
transports in either direction with sufficient coefficient DC. Heights are
unchanged by output conjugation.

This is a LINEAR phase for a specific COMMUTING polynomial derivation
system, conjugate to constant derivatives by ONE shear. It is NOT a
constant-coefficient linear-phase theorem, and it does not promote GVC to
arbitrary variable-coefficient operators. The full original gap is thereby
represented on nested power-apolar families, or equivalently on these
constrained linear-phase operators. This gives a global normal form without
shrinking the original phase/symbol class.

## 5. Exact boundary of universal original-power nesting

**GVCAPOLARNESTCLASS20261116.** For nonzero polynomial phases over
characteristic-zero fields, universal adjacent-power apolar nesting holds
in EVERY degree in ONE variable and in EVERY dimension at degree<=2.
For EVERY n>=2 and EVERY EXACT degree d>=3 it fails, already between
powers2 and3, even for a PURE SAFE operator with sharp real cutoff1/d.
Thus these are the exact universal dimension/degree ranges.

Unary proof. If degP=d, the derivatives of P^m have distinct leading
degrees0,...,dm and span all univariate polynomials of degree<=dm. Hence
Ann(P^m)=(Z^(dm+1)); these ideals are nested. Nonzero constant phases are
included with d=0.

Quadratic proof. Extend the field to an algebraic closure; ideal
containments and polynomial identities descend by faithful scalar extension.
Complete the square in the nondegenerate quadratic directions. If a linear
term remains in a radical direction, a linear derivative takes P to a
nonzero constant, and the argument(5) applies. Otherwise translation and
linear coordinates give P=q+c, q=sum_(i=1)^r x_i^2. Translation commutes
with constant derivatives and preserves each annihilator. For r>=1,

    Delta(q^j)=2j(r+2j-2)q^(j-1).

All leading coefficients are nonzero in characteristic zero. The m
polynomials Delta^k(P^m),1<=k<=m, have distinct q degrees m-1,...,0,
so they span ALL polynomials in q of degree<=m-1. In particular P^(m-1)
is a constant linear combination of them. Commuting any formal annihilator
through that finite combination proves nesting. The remaining rank0 case
is a nonzero constant or a linear phase, already handled.

Counterfamily. In TWO variables let d>=3 and put

    P_d=x^d+y^d,
    R_d=3(X^(2d)+Y^(2d))-(3+binom(3d,d))X^dY^d,
    T_d=R_d(partial),
    H_d=binom(3d,d)-3binom(2d,d)+3.                        (10)

The integer H_d counts d-element subsets of three disjoint d-element
blocks meeting ALL THREE blocks: inclusion-exclusion gives precisely(10).
It is strictly positive for every d>=3. Direct differentiation gives

    T_d(P_d^3)=0,
    T_d(P_d^2)=-2(d!)^2 H_d !=0,
    T_d(x^d P_d)=-(d!)^2 H_d !=0.                          (11)

For completeness, write U=(3d)!/d!, V=(2d)!, C=binom(3d,d).
The two degree-d outputs on the cube are

    (partial_x^(2d)+partial_y^(2d))P_d^3=(U+3V)P_d,
    partial_x^d partial_y^d P_d^3=3V P_d.

Since U=CV, their combination in(10) is0. On the square the respective
constant outputs are2V and2(d!)^2, giving(11). On x^dP_d they are V
and(d!)^2, giving its last equality. This is the full all-d proof, not a
fit to the small controls.

The operator has homogeneous order2d and phase degree d, so degree counting
gives T_d^m P_d^m=0 for EVERY m and

    T_d^m(QP_d^m)=0 whenever m>degQ/d.

The last equality in(11) proves cutoff1/d is SHARP. Its actual heights
are h_1=0,h_2=1,h_3=0: a height decrease inside a pure safe pair.
For d=3, R=3X^6-87X^3Y^3+3Y^6, square output -1944 and cube output0.
Padding unused coordinates preserves the failure in every n>=2.

Consequently general Ann(P^m) nesting, and even monotonicity of h_m for
pure safe pairs, are false from the smallest possible dimension2 and
degree3. These pairs are NOT eventual-but-nonlinear counterexamples.
Their suspended heights are running maxima and their lifted actual ideals
are nested; original nonnesting does not obstruct that exact replacement.

## 6. Literature context and the remaining decisive obligation

The retained Flavi--Jelisiejew--Michalek
[2408.02754v2](https://arxiv.org/html/2408.02754v2), published IMRN2025,
was reread at Section3's power-apolar subquotient diagram and all of
Section4 through the encompassing construction. Encompassing phases have
a constant linear derivative in the proof of their maximal-dimension
proposition. General power algebras are not unconditionally symmetric
tensor powers. We use neither such an identification nor a smoothability
theorem here: one simple linear suspension provides the needed nesting.
It does NOT in general make the phase encompassing or imply the paper's
maximal-dimension hypothesis.

The remaining full target can now be pursued on the ACTUAL nested lifted
families, with the same infimum/asymptotic gaps and an equivalent scalar
time-equation gate. A nested family alone need not have a finite Rees
presentation or bounded-error quasi-additivity. The retained actual safe
nilpotence-index counterfamily still rules out automatic subadditivity.
No finite recurrence,positive-gap theorem,or eventual-but-nonlinear witness
for the general lifted class has been supplied.

The checker uses small exact binomial/height/gauge identities,two members
of the all-d counterfamily,quadratic radial controls,the shear formula,
and guards. It does not perform a growing catalecticant,degree or prime
campaign. Both full exceptional MZ cases stay open;concurrent boundary
reductions and the FOUR-Gaussian bridge are preserved. Closure-first and
deferred external review/outreach persist.
