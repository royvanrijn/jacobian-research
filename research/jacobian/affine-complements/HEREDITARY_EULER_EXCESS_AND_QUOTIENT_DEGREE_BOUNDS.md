# Hereditary Euler excess and bounds in the intermediate degree

16 November 2026. PF2RULED129. Written deductions in EVERY original
degree and for EVERY genuine finite intermediate, smooth or normal.
The full plane Keller/Wright and shared marked trace/energy gates remain
OPEN. Published inputs are credited; no novelty or independent review
is asserted. The finite controls do not prove the unbounded statements.

## 1. Keep the entire intermediate etale locus

Let F:A2->A2 be an actual nonautomorphic complex plane Keller map,
N its full finite normalization, and

    N --h--> M --p--> A2, deg(h)=t, deg(p)=b>1.

Use PF123/127: the critical primes E_i of p are pairwise disjoint smooth
A1 curves, with smooth-axis or small cyclic-quotient-axis local pairs.
Their image normalizations are A1. Write e_i>=2 for transverse index,
f_i for degree E_i->C_i^nu, and s for their number. If several E_i have
the same image C, keep ALL of them. Let

    A=sum_i e_i*f_i=sum_j a_j,

where a_j is the moved support of the j-th quotient branch meridian on
b blocks. These are MOVED LETTERS, not permutation indices or e_i-1.

Let V be the ENTIRE p-etale locus. By HC4MYKU0/1 it is smooth affine,
and p has no isolated critical points outside the E_i. Thus V=M\union E_i.
The finite preimage Z=h^-1(V) contains the WHOLE original A2 openly.
PF128 Section1 proves pi1(Z)=1 and that pi1(V) is finite of order m|t.
In particular b1(V;Q)=0. The Stein two-dimensional CW bound gives

    c=chi(V)=1+b2(V;Q)>=1.                            (1)

This does not assume that M has a plane chart, that h|Z is etale, or
that a generic deck quotient is the entire V.

## 2. Exact Euler integration, including singular source points

For a point z on E_i put q_z=deg_z(E_i->C_i^nu) and let ell_z be the
local finite degree of p at z. Define

    eta_z=ell_z-e_i*q_z,   Eta_M=sum_(i,z on E_i) eta_z.

Here C_i^nu is the normalization of the chosen image branch, even when
several analytic branches meet at p(z). The finite flat map p has a
constant total scheme fibre length b: normal surface M is Cohen--Macaulay
and the target is regular (the retained KU1 miracle-flatness argument).
Near a general point of that chosen image branch, the local factor at z
has q_z critical points of multiplicity e_i, so ell_z>=e_i*q_z. Hence
eta_z is nonnegative, and is zero away from finitely many points.

Proper finite Euler integration of the local-degree function gives

    b=chi(M)+sum_i integral_(E_i)(ell_z-1).

No singular point of M is dropped: it is on one E_i and its full local
degree appears in ell_z. A degree-f_i polynomial A1->A1 has total
finite ramification f_i-1. Therefore integral_(E_i) q_z=f_i and

    integral_(E_i)(ell_z-1)=e_i*f_i-1+sum_z eta_z.

It follows EXACTLY that

    chi(M)=b+s-A-Eta_M,
    chi(V)=b-A-Eta_M=c,
    A+Eta_M+(c-1)=b-1.                              (2)

This is the original Euler ledger's intermediate analogue, with the
nonnegative term c-1 in place of an assumed intermediate plane-boundary
count. It uses the genuine finite map and whole-source topology. It
does not transfer the original omitted-prime Picard basis to M.

Consequently EVERY genuine quotient satisfies the stronger inequality

    A+Eta_M<=b-1.                                   (3)

If B=b-A is the quotient support slack, then

    0<=Eta_M<=B-1,    c=B-Eta_M.                    (4)

All arithmetic uses the INTERMEDIATE degree b, even if the original
degree bt is arbitrarily large. Replacing b by a representation degree
of an unrelated abstract quotient would be invalid.

## 3. Equality identifies more minimal-counterexample cases

Whenever A+Eta_M=b-1, equation (2) makes c=1. PF128 now applies:
the published universal-cover plane chart gives an ACTUAL new plane
Keller map of degree b*m, with m=|pi1(V)| dividing t. If m<t this is
a strictly smaller counterexample. For an original map of MINIMAL
generic degree, m=t and

    h|Z universal cyclic etale, Z=U_N,
    FULL generic rational deck D=C_t, b0=t-1.        (5)

No primitivity is needed. There is no residual original ramification
over V. The full surviving cyclic finite-target case is not excluded.

In particular maximal support A=b-1 AUTOMATICALLY forces Eta_M=0 and
c=1. A local critical germ with a SMOOTH chosen image branch is a
connected finite normal cover branched only over a smooth axis. The
axis complement has pi1=Z, so normal completion is the smooth cyclic
map (x,z)->(x^e,z), with q=1 and eta=0. Conversely a SINGULAR chosen
image branch has eta>=q>=1 by PF115, whose hypotheses PF123 supplies.
Thus Eta_M=0 implies M smooth, every f_i=1, and every image normalization
is immersive everywhere. Self-identifications, contacts between smooth
analytic branches, and multiple image components are still allowed.

For a MINIMAL counterexample in this maximal-support case, N is smooth
and EVERY original critical residue is one. Indeed (5) removes residual
critical primes; the local smooth-axis cover above each E_i is smooth,
and PF127 gives n_D*f_D=n_E*f_E=1. The relative cyclic Galois group
makes the inertia k_i over E_i constant, with k_i|t, so

    over E_i: t/k_i critical primes,
              all original indices e_i*k_i and residues1.

Adding the disjoint smooth E_i to V kills their meridians. Hence here

    pi1(M)=C_(t/lcm_i k_i).                         (6)

No equality k_i=t or simple connectedness of M is assumed. In general
G is contained in C_t wr S_b; solvability follows only when the block
action is solvable. This does not extend PF128's S4 solvability claim
to arbitrary b.

## 4. Slack at most two: an exact all-index local classification

When B<=2, (4) gives Eta_M<=1. In this range EVERY intermediate
critical residue f_i is one. For if f_i>1, PF127 puts its sole finite
ramification at a quotient point, with local q=f_i>=2. Its image
branch must be singular by the smooth-axis argument above. PF115 then
charges eta>=q>=2, impossible.

If Eta_M=0, Section3 describes the smooth immersive case. If Eta_M=1,
there is exactly ONE singular image-normalization parameter in the
entire critical image, and exactly ONE critical prime above its image
component. Otherwise another preimage would charge another positive
excess. PF115 makes eta=1 force q=1, and its COMPLETE equality
classification gives an ordinary cusp with exactly two local profiles:

| e | local q | quotient n | local degree ell | source |
|---|---|---|---|---|
| 2 | 1 | 1 | 3 | smooth |
| 3 | 1 | 2 | 4 | A1 quotient |

The separate smooth e3,q2,n1,ell8 equality germ in PF115 has eta=2,
so is already outside this budget; it is not a third eta=1 case.
All other image parameters are immersive, and all other
points of M are smooth. The remaining source singularity, when
present, is exactly one A1 point on the index3 critical curve.

For B=2,Eta_M=1, c=1, so (5) holds in a MINIMAL counterexample. Every
original critical residue is then also one: n_D*f_D=n_E*f_E is either
1 or2; the factorization 2=1*2 is forbidden because n_D=1 forces f_D=1.
In the smooth cusp case the full N is smooth. In the A1 cusp case all
singularities of N are A1 points above that cusp, with relative inertia
k necessarily ODD (gcd(2,3k)=1), original index3k and local degree4k.
There are t/k such points. No formula (6) for singular M is inferred:
at a quotient point, adding the axis kills more than the generic
meridian relation alone.

For B=2,Eta_M=0 instead c=2, and the compression conclusion (5) is NOT
automatic. No finite-index or residue-one profile is declared realized.

## 5. The all-degree torus-cusp bound descends in degree b

Fix a quotient branch image C with moved support a_C>b/2. PF122's
disjoint-support theorem at multibranch values makes its normalization
injective; C is contractible. It is an original nonproper image component,
so Chau's retained no-smooth-A1-component result excludes a smooth C.
The retained Palka/Lin--Zaidenberg classification from PF113/114 makes
C a one-Puiseux-pair torus cusp in polynomial target coordinates.
These are PRIOR explicit imports, not new classifications.

The LOCAL PF114 solid-torus inequality applies to each critical factor
of M by PF123. Over the cusp, write q_z for normalized local degrees.
It gives e_i*q_z<=5*eta_z-2*q_z. Since sum_z q_z=f_i, summing gives

    a_C<=5*Eta_C-2*F_C,    F_C=sum_(E_i over C) f_i>=1.

The FULL quotient branch cannot be this one contractible C, by PF122.
Another quotient branch image therefore costs at least TWO moved
letters. Use (2), c>=1 and all remaining nonnegative excesses:

    b>=a_C+2+Eta_C+1,
    6*a_C<=5*b-15-2*F_C.                           (7)

In particular for EVERY b>=9, EVERY quotient branch image and EVERY
intermediate critical prime satisfy

    a_C<=floor((5*b-17)/6),
    e_i<=floor((5*b-17)/6).                         (8)

For a_C<=b/2 the same bound follows arithmetically when b>=9; otherwise
use (7). This transfers the actual LOCAL knot/filling obstruction via
the newly proved Euler ledger, not by pretending M is a plane.

PF106's extremal quotient-order argument, now with PF127's proved
hereditary power-map location, gives f_i=1 for indices1,2,3,4,6 in
EVERY intermediate. Therefore f_i>1 would require e_i>=5 and e_i*f_i>=10.
For b<=10 equation (3) already forbids this. For 11<=b<=15, a_C>=10>b/2
contradicts (8), whose right side is at most9. We conclude

    EVERY genuine intermediate of degree b<=15
    has f_i=1 at EVERY critical prime.              (9)

The original degree may be arbitrarily large. Equation (9) is about
critical residues in M; it does not say all original residues are one.
For an original prime D above E_i, PF127 still requires the exact
invariant n_D*f_D=n_E*f_E. In particular the useful degree-only bound

    n_D*f_D=n_E*f_E<=(e_E-1)*f_E<=b-2              (10)

holds in EVERY degree b>1, whether or not (9) applies. Original
invisible primes are outside (10). The smooth/nontrivial-index case
n=1 is included. No intermediate omitted-prime coverage is assumed.

## 6. Scope and evidence

The new universal input is the exact calculation (2) and its use with
the whole-source finite pi1 statement. PF114/115/123/127/128 and their
published source manifests are retained rather than replayed. No new
classification, low-sheet inversion theorem or theorem about solvable
monodromy is imported. The Domrina degree-four route was screened:
its source boundary comes from blow-ups of the plane and is unimodular;
a pseudo-plane boundary is not supplied with that marking. Its published
plane theorem is not applied to a pseudo-plane source.

Exact controls check Euler integration at coincident values, the binary
profiles, the degree-only inequalities, and the strict scope boundaries.
No degree/depth census, actual plane Keller realization, full natural
A/S or binary affine exclusion, formal verification or external review
is claimed. The cyclic finite-target and full trace/energy/empty-original-
fibre obligations remain OPEN.
