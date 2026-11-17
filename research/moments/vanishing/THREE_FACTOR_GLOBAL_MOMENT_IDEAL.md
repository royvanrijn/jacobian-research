# Global ideal generation for three Gaussian factors

14 November 2026. Characteristic zero. Internal written invariant-ring
proof and exact symbolic certificates; no external review, formal
verification or priority is claimed. This settles the global ideal question
left open in the [three-factor classification](NODAL_CUBIC_REDUCTION.md),
not the arbitrary ternary cubic GMC problem.

## Statements: eight detects the set, ten generates the scheme

Let R=Q[a,b,c,d,e,f] be the ring of a symmetric three-factor Gaussian Gram
matrix

    Sigma=[[a,d,e],[d,b,f],[e,f,c]].

Let M_n=E((L1*L2*L3)^n) be the universal Wick polynomial, and put

    I8=(M2,M4,M6,M8),
    I10=(M2,M4,M6,M8,M10),
    I_inf=(M_n:n>=1).

Odd M_n are identically zero. The exact global result is

    I_inf=I10,  M10 not in I8,
    M_n in I8 for every n>=12.                    (1)

Thus **ten is the least initial ideal-generation cutoff**, whereas eight
remains the least universal set-theoretic pure-moment cutoff from the
previous theorem. I8 and I10 have the same radical J, but they are different
schemes. In particular agreement of all six generic local algebras did not
justify global equality.

Writing J=(af,be,cd,abc,def), let nu(I) be the least k with J^k contained
in I. We additionally prove the global bounds

    9<=nu(I8)<=19,  9<=nu(I_inf)<=17.               (2)

The nonzero module I10/I8 is zero at every generic minimal component.
Consequently I8 has embedded associated structure. The subsequent [colon geometry theorem](THREE_FACTOR_COLON_GEOMETRY.md)
now gives exact indices14/13 and the complete annihilator primary
decomposition: three diagonal axes and embedded origin. A full primary
decomposition of I10 itself is supplied by the subsequent
[constructive theorem](THREE_FACTOR_PRIMARY_DECOMPOSITION.md).

The polynomial generation identities in (1) survive every substitution of
Gram entries by coefficients of actual Gaussian linear forms, in any ambient
dimension. The sharp nonmembership assertion is proved in the full universal
Gram ring; it is not silently asserted for every constrained coefficient
parameterization. The new proof and bounds here are in characteristic zero;
the older excluded-prime radical certificate does not automatically certify
these stronger scheme statements in positive characteristic.

## Balanced Wick monomials and a five-generator ring

Here uppercase A,...,E denote covariance monomials, not the linear factors:

    A=abc, B=af^2, C=be^2, D=cd^2, E=def.

Give each covariance its factor-count multidegree, for example
wt(a)=(2,0,0), wt(d)=(1,1,0). Each M_(2m) has weight (2m,2m,2m).
The subring R0 containing precisely these balanced weights is

    R0=Q[A,B,C,D,E]/(BCD-AE^2),                   (3)

with all five generators of balanced degree one.

For completeness, this is an all-degree semigroup statement. In a balanced
monomial the three off-diagonal exponents have the same parity. Remove E
if they are odd, and put m'=m-odd. Write the remaining off-diagonal
exponents as 2p,2q,2r on d,e,f respectively. The diagonal exponents are
m'-p-q, m'-p-r, m'-q-r, all nonnegative. Put

    k=max(0,p+q+r-m').

The pair inequalities give k<=min(p,q,r). The exponents of A,B,C,D,E in
an explicit nonnegative factorization are

    (m'-p-q-r+k, r-k, q-k, p-k, odd+2k).          (4)

They sum to m and pull back to the original monomial. The five exponent
vectors have rank four and primitive integer kernel (-1,1,1,1,-2).
Every binomial relation is therefore a multiple of BCD-AE^2, after removing
a common monomial and a geometric-series factor. Grouping coefficients
of equal monomials proves (3) has no additional polynomial relation.
This does not discard the E=0 covariance faces.

The factor permutations preserve A,E and permute B,C,D. Averaging a
representative and applying the fundamental theorem of symmetric polynomials
gives the symmetric balanced subring

    S=R0^(S3)=Q[A,E,s,q],
    s=B+C+D, q=BC+BD+CD.                          (5)

The third elementary symmetric function is BCD=AE^2. The displayed four
variables are algebraically independent: after localizing at E, the
standard algebraic independence of the three elementary symmetric functions
in B,C,D applies with their product AE^2. This localization checks
independence, rather than replacing the full polynomial ring in (3).
Give A,E,s weight one and q weight two. Every M_(2m) is a polynomial F_m
in S, weighted homogeneous of degree m.

## The exact two-variable quotient

The first two polynomials are

    F1=A+2s+8E,
    F2=27A^2+216As+1728AE+72s^2+720q+1152sE+1728E^2.

Modulo F1,F2 one can eliminate

    A=-2s-8E,
    q=(7s^2+88sE+288E^2)/20.                     (6)

Under (6) the next two moments become

    F3=-1620*F, F4=381024*G,
    F=107s^3+2196s^2E+11280sE^2+17152E^3,
    G=337s^4+12736s^3E+109968s^2E^2
                     +359424sE^3+406272E^4.      (7)

The [certificate](three_factor_global_ideal_certificate.json) gives a complete
Groebner basis of (F,G). Its standard monomials are

    1,E,...,E^5; s,sE,sE^2,sE^3; s^2,s^2E.

Their Hilbert function is 1,2,3,3,2,1, of total length12. In particular
all homogeneous degrees at least six vanish. The basis keeps the E=0
face; this is not an inference from a resultant computed only on E=1.

The tenth moment F5 has the exact nonzero remainder

    F5 = -(225834448896000/92791)*E^5 mod (F,G).   (8)

It kills the unique top-degree class. The quotient after adding F5 has
Hilbert function 1,2,3,3,2, length11. For every m>=6, weighted homogeneity
and (6) turn F_m into an ordinary homogeneous polynomial of degree m in
s,E. The degree-six annihilation therefore proves

    F_m in (F1,F2,F3,F4) for every m>=6.

Pulling this polynomial identity back to R proves membership in (1) for
all higher moments. This all-order step is a grading argument, not an
extrapolation from a computed long moment prefix.

## Why the nonmembership is global in six Gram variables

Nonmembership in an invariant subring alone need not survive extension to a
larger coefficient ring. Here it does, for an explicit reason. Projection
onto balanced weight spaces is R0-linear: adding a balanced weight cannot
turn an unbalanced weight into a balanced one. Averaging over the six
factor permutations is S-linear. Their composite

    pi:R -> R0 -> S

is S-linear and restricts to the identity on S. If F5 belonged to I8 in
R, applying pi to a cofactor expression would put it in (F1,...,F4) in S,
contradicting (8). This proves the strict inclusion in the full ring and
is precisely the contraction step absent from a generic normalization.

There is also a separate direct linear certificate. In multidegree
(10,10,10), all possible homogeneous cofactor multiples of M2,M4,M6,M8
form119 columns in the111-dimensional Gram monomial space. An exact
integer linear functional in the certificate annihilates every column but
has nonzero value on M10. Since the ideal and the target are multihomogeneous,
this exhausts every possible cofactor expression, not just a bounded-degree
ansatz. The same construction proves a*M10 and d*M10 are not in I8.
This last observation does not locate the full annihilator; in particular
it does not by itself prove that the support has positive dimension.

## Balanced nilpotence and full Gram radical powers

There is an explicit global free S-module basis

    1,B,B^2,C,BC,B^2C.

Eliminate D=s-B-C. The remaining relations are the two monic equations

    C^2+BC+B^2-sC-sB+q=0,
    B^3-sB^2+qB-AE^2=0.

For lexicographic order C>B their leading monomials C^2 and B^3 are
coprime, so they form a monic Groebner basis over the polynomial
coefficient ring S. This proves freeness globally, without replacing S
by its fraction field. The basis degrees are0,1,2,1,2,3. Multiplying the
preceding Hilbert functions by1+2z+2z^2+z^3 gives

| Balanced quotient | Hilbert function | Length | Last nonzero degree |
|---|---|---:|---:|
| R0/I8 | 1,4,9,14,16,14,9,4,1 | 72 | 8 |
| R0/I10 | 1,4,9,14,16,13,7,2 | 66 | 7 |

The checker also computes both quotients directly in all five generators
with relation (3), independently checking these Hilbert functions. These
lengths are invariant-ring lengths; they are not the generic Gram-component
lengths24/6 or the full top-dimensional scheme degree90.

Let H=(A,B,C,D,E)R. The table gives H^9 contained in I8 and H^8 contained
in I10. To pass to J, write u=af,v=be,w=cd. The identities

    u^2=aB, v^2=bC, w^2=cD, uvw=AE

imply J^(2k+1) is contained in H^k. Indeed factor as many squares of
u,v,w as possible; at most one of each remains. If all three remain,
replace their product by AE. Counting these H factors proves the bound
for every monomial in J^(2k+1). Thus J^19 is contained in I8 and J^17
is contained in I10.

For the lower bound, use the exact generic null-row algebra already
computed in the preceding note. Set b=c=1,f=z over Q(z), eliminate a by
M2, and reduce by the homogeneous fourth and sixth moments. The eighth
and tenth moments reduce to zero, but e^8 survives as a standard monomial.
The global element (be)^8 in J^8 maps to e^8, so it is in neither I8
nor I10. This proves both lower bounds in (2).

## Embedded structure: what is now closed and what remains

The earlier generic local algebras of I8 and I_inf agree on all six
minimal components (row length24/index9, star length6/index4). Their
nonzero difference I10/I8 is therefore a submodule of R/I8 supported in
strictly smaller dimension. A nonzero Noetherian module has an associated
prime, and associated primes of a submodule are associated primes of the
ambient module. None here is a minimal prime of R/I8, so I8 has embedded
associated primes. Generic agreement did not see this structure.

This proves existence, not a primary decomposition or an exact support
calculation. The [subsequent exact colon calculation](THREE_FACTOR_COLON_GEOMETRY.md)
now completes that annihilator and sharpens the indices to14/13. The
[subsequent associated-prime theorem](THREE_FACTOR_ASSOCIATED_PRIMES.md)
now gives all ten primes and the six isolated primary components of I10.
The [constructive primary theorem](THREE_FACTOR_PRIMARY_DECOMPOSITION.md)
now supplies its four embedded components and intrinsic origin torsion77.
The remaining scheme gates concern constrained linear-form pullbacks. Universal
ten-moment generation survives all pullbacks. The subsequent
[coefficient theorem](THREE_FACTOR_COEFFICIENT_PULLBACK.md) now gives
sharp generation gates2/6/10 and complete flat primary transfers for N>=5;
full N=2,3,4 coefficient thickness and further constraints remain separate.
The arbitrary cubic GMC problem and its nodal Gaussian-width/infinity
step are unchanged by this scheme refinement of an already safe class.

## Reproduction and the literature boundary

Run `.venv/bin/python research/moments/vanishing/verify_three_factor_global_ideal.py`.
The [receipt](three_factor_global_ideal_verification.json) binds the preceding
checker by hash, reconstructs the first five symmetric moment polynomials,
checks their pullbacks against full six-variable Wick polynomials, and
checks M10 independently by pair-count enumeration. It records both exact
invariant quotients, three full-Gram dual nonmembership witnesses and the
radical-power bounds. The semigroup formula and weighted homogeneity prove
the quantified infinite tail. These are internal exact algebra results,
not proof-assistant formalization or independent specialist review.

The [source audit](three_factor_ideal_source_audit.json) records a fresh
product-normal literature pass. The October2026 variance-gamma paper treats
CDFs and absolute moments for real distributions, including two correlated
normal factors. The related Stein-characterization paper treats two factors
with nonzero means and moment recurrences. Neither supplies the complex
three-factor ideal equality or justifies replacing holomorphic pure moments
by absolute moments. No formula from these probability papers is used to
prove (1); no producer code was executed.
