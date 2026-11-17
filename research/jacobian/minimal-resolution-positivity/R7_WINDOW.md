# The normal-r7 global attachment has an 8--11 Picard-rank window

14 November 2026. This applies the minimal-resolution positivity theorem
and its general-line genus corollary from [PROOF.md](PROOF.md) to the exact
single-ramified-row hypothesis in PF2R7BND1. It does not close normal r7.

The [subsequent adjoint-fibration result](../r7-adjoint-exclusion/README.md)
sharpens this necessary window to rank eight and excludes both nonimmersed
coefficient classes under the explicit hypothesis and source-model inputs.
The window calculation below is retained evidence; its old attachment
continuation is superseded by that scoped exclusion.

## Hypotheses

Let F:A2 -> A2 be a hypothetical complex Keller map of geometric degree
six. In its affine finite normalization, assume the entire reduced
ramification divisor is one prime E_ram of generic (e,f)=(2,1), whose
image C is one of PF2R7CUSP1's two nonimmersed coefficient classes.

These are exactly the hypotheses of
[PF2R7BND1](../jc2-boundary-attachment/PROOF.md). In particular this is
NOT an assumption that E_ram is the entire omitted boundary. That result
forces t=2 or 3 additional unramified omitted primes, with images different
from C. Here t denotes the exact number of such affine-normalization
primes, not the number of distinct target curves.

Normalize target P2 in the source function field; let Y be its minimal
resolution and g:Y -> P2 the resulting morphism. Put H=g*L_infinity,
N=H^2=6, and gamma=genus of a smooth general target-line pullback. Write
rho=rho(Y), and R_tot for the full projective determinant divisor.

## 1. A general affine-ramification budget

For any plane Keller map with the preceding minimal model, define

    B = sum over affine ramification primes (e_i-1)*f_i*deg(C_i).

Only ramified primes occur in this sum; omitted unramified primes
contribute zero and are not discarded from the finite cover.

For each prime over the generic target line at infinity let e_j be its
normal index and f_j its residue degree. Set s=sum f_j. Finite degree
conservation gives sum e_j*f_j=N, while s>=1. A general projective target
line avoids the finite singular and exceptional-image loci. Its pullback
is smooth, irreducible and of genus gamma, and the restricted finite map
to P1 has affine ramification degree B and infinity ramification degree

    sum (e_j-1)*f_j = N-s.

Riemann--Hurwitz therefore gives the exact identity

    2*gamma-2 = -2*N+B+(N-s),
    s = B-N+2-2*gamma.                                      (1)

The general-line genus theorem gamma>=2 implies

    B >= N+3,
    2 <= gamma <= floor((B-N+1)/2).                          (2)

This is a weighted projective degree of the affine branch data; it is not
the number of omitted primes, the degree of F, or an infinity contact sum.

## 2. r7 forces gamma=2 and sum f_j=2

The two classified parametrizations p,q have degrees (6,10) and are
birational. A generic target line pulls back on their normalization to a
degree-ten polynomial, so deg(C)=10. The single ramified affine row has
e=2,f=1. Consequently B=10, regardless of how the two or three omitted
unramified companion images intersect C.

With N=6, (1) is s=6-2*gamma. Since s>=1 and gamma>=2,

    gamma=2, s=2.                                         (3)

The prime decompositions over generic target infinity are now exactly:

| number h of primes | unordered generic (e,f) rows |
|---|---|
| 1 | (3,2) |
| 2 | (1,1)+(5,1) |
| 2 | (2,1)+(4,1) |
| 2 | (3,1)+(3,1) |

The one-prime case has two geometric points on a general target-line
pullback over its infinity point. It must not be misread as s=1.
These are necessary generic decomposition types, not realized boundary
models, and contracted components of the resolution are not on this list.

## 3. Picard rank and full ramification

Adjunction gives K.H=2*gamma-2-N=-4. Since K^2=10-rho, Hodge index gives

    6*(10-rho) <= (-4)^2,

so rho>=8. Positivity gives

    Q^2=16+24*gamma-4*rho-3*N=46-4*rho>0,

so rho<=11. Hence the exact numerical window is

    8 <= rho(Y) <= 11.                                   (4)

| rho | K^2 | Q^2 | R_tot^2 | p_a(R_tot) |
|---:|---:|---:|---:|---:|
| 8 | 2 | 14 | 32 | 12 |
| 9 | 1 | 10 | 31 | 11 |
| 10 | 0 | 6 | 30 | 10 |
| 11 | -1 | 2 | 29 | 9 |

These arithmetic genera concern the full determinant divisor, including
infinity and exceptional multiplicities. They are not the genera of C
or of its rational parametrization.

## 4. A bound on the contracted part of the boundary

All boundary primes of the smooth Y are of three disjoint types:

1. the strict transforms of the omitted affine primes: E_ram and the
   t=2 or 3 unramified primes from PF2R7BND1;
2. the h=1 or 2 primes over generic target infinity in the table;
3. the c primes contracted by g, possibly over affine or infinity points.

The whole-A2 boundary basis gives their total number rho. Therefore

    c = rho-(1+t+h),  hence 2 <= c <= 7.                   (5)

Each such curve is exceptional for the resolution of the finite normal
surface. No smooth rational (-1)-curve occurs on this model. The bound
does not bound their weights, intersection multiplicities, singular
boundary contacts, or the number of further blowups needed for SNC.

## Consequence for the next attachment calculation

Instead of allowing an arbitrary projective completion, a complete
candidate under the stated r7 hypothesis must supply one of four infinity
decompositions, rho in {8,9,10,11}, and between two and seven contracted
boundary primes, together with the exact 1+t+h+c ledger. It must still
pass every curve intersection test from PROOF.md, the global Euler
conditions, and the actual morphism and affine-plane-source requirements.

This does not exclude either coefficient class, the immersed r7 family,
normal r5, nonnormal slices, or (75,125). Additional branch images are
retained. The rank window is proved only with the explicit degree-six,
single-ramified-affine-row hypotheses. The symbolic and finite integer
checks in verify.py support the written proof; no formal verification,
independent external review or novelty priority is claimed.
