# The quadratic remainder must balance at every actual pole

PF2RULED178. 17 November 2026. Actual complex polynomial plane Keller maps,
within the complete marked quadratic source model.

Every survivor of the ruling-six obstruction has a genuinely quadratic
remainder

    P=H3^2+d(u)z+e(u)w+f(u), d!=0, e!=0.

At every actual A supporting base value, the quadratic and linear
remainder terms have exactly equal polar weight. If its slope is L,
then

    a=ord d, b=ord e, a=b+L,
    L<=a<=2L-2, 0<=b<=L-2.

Their leading constants satisfy one of two explicit resonance ratios,
using the SAME reduced-companion odd index at every supporting base.
These are necessary conditions on the full quadratic branch; that
branch and HC4/JC2 remain open.

## 1. Actual source and retained companion resonance

Use the original source fields and full normalizations of PF171--176:

    w^2=N+M*z, dP wedge dQ=du wedge dw/M,
    P=H3^2+R2, R2=d*z+e*w+f,
    leading_w H3=h/M,
    H3,d*z,e*w,f are regular source-image functions.       (1)

Here d,e,f are polynomials in u. No lower-ruling function is put in
the target field or full finite normalization. The actual companion
is globally reduced, with

    deg_w Q=3k, k>=3 odd,
    leading_w Q=beta*(h/M)^k, beta in C*.

Write m_n=[w^-1]P^(n/2) for positive odd n. The cyclic-constant
argument of [PF176](RULING_SIX_COMPANION_RESONANCE_AND_LINEAR_REMAINDER_CLOSURE.md)
gives an actual rational-function identity

    sum_(1<=n<=k,n odd) c_n*m_n=0,
    c_n in C, c_k=beta!=0.                               (2)

This lemma in PF176 is proved before its linear-remainder assumption.
It uses the noncube leading equation from PF173 and the whole actual
Jacobian, not a converse to the separate exactness of each moment.
The index k has no upper bound and is fixed globally for the chosen Q.

PF176 excludes d=0. We show that e=0 is impossible as well, and
obtain the stronger local restrictions when both are present.

## 2. Whole scaled functions at an actual A support

Fix ANY actual supporting base value c of an A pole cluster. Put
v=u-c and let its integer slope be L. The retained square face is

    H3(v,v^-L*Y)=v^-L*(F(Y)+v*H1(v,Y)),
    F(Y)=C*Y*(Y-gamma)^2, C,gamma in C*, M(c)!=0.         (3)

All scaled coefficients are in C[[v]]. This also follows directly
from the regular square form: PF174 gives a=ord_c d>=L, while
ord_c e>=0. Thus the remainder has weight at most L, strictly
below the weight 2L of the actual P face, and the chosen square
root initial is (3).

Set a=ord_c d and b=ord_c e, using infinite order for a zero
polynomial. PF175 gives

    a<2L OR b<L.                                        (4)

Hence R2 has some pole in this scale. Define

    q=max(2L-a,L-b)>0.

Since a>=L and b>=0, one has 1<=q<=L. The entire scaled
remainder satisfies

    R2(v,v^-L*Y)=v^-q*(T(Y)+v*R1(v,Y)),                 (5)

where T is one of

    A*Y^2,          if a-2L<b-L,
    B*Y,            if b-L<a-2L,
    Y*(A*Y+B),      if a-2L=b-L.

Every displayed coefficient A or B is nonzero. The constant term
f-d*N/M is strictly higher because it is regular off M.
The powers and lower terms in (3),(5) are integral in v; all the
original branches and normalizations remain retained.

## 3. A complete odd moment has one possible leading coefficient

Let n=2r-1>=1, so r>=1. The full finite binomial expansion from
PF175 gives

    m_n=sum_(j=r)^floor((3n+1)/4)
           binom(n/2,j)*[w^-1](H3^(n-2j)*R2^j).         (6)

Terms of nonnegative H3 exponent are polynomials and have zero
coefficient residue. Coefficient residue scales by

    [w^-1]A(w)=v^-L*[Y^-1]A(v^-L*Y).

Consequently the j term of (6) has order at least

    -L*(n+1)+(2L-q)*j.                                  (7)

The first index j=r has order at least -r*q. Every later j has
strictly higher lower bound, since 2L-q>=L>0.
The leading coefficient at the first index is exactly

    binom(r-1/2,r)*[Y^-1](T(Y)^r/F(Y)).                  (8)

Reduction modulo v commutes with this coefficient operation:
the scaled cubic has a unit leading Y coefficient, so its inverse
powers have integral Laurent coefficients. In particular

    ord_c m_(2r-1)>=-r*q,                               (9)

with equality whenever (8) is nonzero. These statements use the
entire polynomial and all binomial terms, not a truncated inverse.

## 4. Strictly dominant monomials contradict the actual companion

If T=A*Y^2, the coefficient in (8) apart from its binomial factor is

    (2r-1)*A^r*gamma^(2r-2)/C,

which is nonzero for EVERY r>=1. If T=B*Y, it is

    (r-1)*B^r*gamma^(r-2)/C,

which is nonzero for EVERY r>=2. These identities follow by
expanding at Y infinity, or differentiating the numerator over
the sole finite double pole (Y-gamma)^2.

For the top odd character n=k>=3 in (2), equation (9) is then
an equality with r=(k+1)/2. Every smaller odd character has a
strictly smaller possible pole order. Its leading pole therefore
cannot cancel in (2). Both strict dominance cases are impossible.

It follows that both d and e are nonzero and their weights balance
at EVERY actual A supporting base value:

    a-2L=b-L, a=b+L.                                    (10)

The support inequality (4) gives a<2L and b<L. The rational
character d*du/(Mh) of PF173 excludes a=2L-1, since ord_c h=2L
and M(c) is a unit. Together with PF174, this proves

    L<=a<=2L-2, 0<=b<=L-2.                             (11)

No claim that e divides d, that h divides d or e, or that the
remaining source functions lie in the target is made.

## 5. Exact leading ratio in the balanced case

Now q=L-b=2L-a>=2. Define the actual nonzero leading constants

    A=initial(v^-a*d)/M(c), B=initial(v^-b*e),
    T(Y)=Y*(A*Y+B).

For every r>=1 the coefficient in (8) before its binomial factor is

    [Y^-1](T^r/F)
     =C^-1*gamma^(r-2)*(A*gamma+B)^(r-1)
                         *((2r-1)A*gamma+(r-1)B).       (12)

Indeed T^r/F=C^-1*Y^(r-1)*(AY+B)^r/(Y-gamma)^2.
Its coefficient residue equals the derivative of the numerator at
gamma. Formula (12) is valid also for r=1, with zeroth power
interpreted as one.

If the top-character coefficient (12) were nonzero, the same
strict pole comparison in (2) would contradict the companion.
For r=(k+1)/2>=2 it must therefore vanish. Since all other
factors are nonzero, the exact alternatives are

    A*gamma+B=0,
    OR B/(A*gamma)=-2k/(k-1), k>=3 odd.                 (13)

The first alternative says that the leading remainder vanishes
at the nonzero double root gamma of the cubic face. The second
is a discrete ratio determined by the actual reduced companion.
The SAME k works at every actual A base for that Q. Neither
alternative is asserted sufficient for an actual companion.

For example a slope-five support must have (a,b) among
(5,0),(6,1),(7,2),(8,3). This is a local consequence of (11),
not a closure or a complete global coefficient classification.
When (12) vanishes, the leading pole estimate for that moment
can improve. This theorem makes no claim about its next term,
the full balanced branch or cancellation in the entire identity.

## 6. Proof and exact controls

The checker verifies the unrestricted coefficient formula algebraically,
the finite residue sums on pure and balanced initial models, the complete
binomial order margins, the strict top-pole comparison, both possible
ratio cancellations and their non-sufficiency as leading tests.
The written finite-range and valuation arguments prove the result for
all odd companion indices and all coefficient degrees.

Replay:

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_balanced_quadratic_resonance.py

The arbitrary-degree result depends on the actual companion resonance,
PF171 regular square form, PF173 rational character, PF174 half-order
bound and PF175 pole-support theorem. Full ruling-six closure, target
integrality, normal-tail termination and HC4/JC2 remain open.
Independent whole-proof replay, formal verification and external human
review remain unclaimed.
