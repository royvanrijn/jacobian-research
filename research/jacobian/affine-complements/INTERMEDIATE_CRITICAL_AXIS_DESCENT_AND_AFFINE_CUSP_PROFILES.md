# Critical axes descend to intermediates; binary affine quotients have two profiles

PF2RULED123, 16 November 2026. Work with an ACTUAL complex plane Keller
map, its ORIGINAL finite normalization N, the WHOLE U=A2 and ALL omitted
primes. For EVERY nontrivial sheet-block quotient, local monodromy is
intransitive at every finite target point. Every intermediate normalization
has disjoint smooth A1 critical curves; its singular points are cyclic
quotient-axis points of marked type 1/n(1,-e), with n<=e-1.

For a binary affine quotient the intermediate surface is entirely SMOOTH.
The visible original primes have residue degree1 and even indices bounded
by twice the block size. The common branch image is either immersive or
has exactly ONE ORDINARY cusp. This closes higher/multiple cusp and
visible singular-source/residue-cover alternatives for ALL such quotients,
without assuming an intermediate plane. Natural/affine compatibility and
the full finite algebra remain OPEN; the full goal stays ACTIVE.

## 1. No intermediate has locally transitive full sheet monodromy

Use PF122's notation Omega->Gamma, |Gamma|=b>1, block size t, d=tb,
and the finite factorization N->M->A2. Fix a finite target value y.
Let K be original local monodromy. Suppose its image is transitive on
Gamma. Every original local factor O is a K-orbit whose nonempty block
image is therefore ALL Gamma. Thus |O|>=b>1; O comes from an omitted
critical point p, whose local branch image is just ONE analytic branch
beta. Every meridian about another branch acts trivially on O and hence
on Gamma. All original local factors have the same visible branch beta:
a different one would force the beta meridian trivial on Gamma as well.

The intermediate full branch union B_Q has exactly this one branch at y.
Otherwise another visible branch meridian would act nontrivially. Its
irreducible component C has only ONE normalization parameter over y,
since all branches of a globally visible component are visible. All
points of the complete original fibre lie on omitted D_i mapping to C.
For each such D_i, the local residue degrees above this parameter sum
to f_i. As in PF122, the exact original normalized ledger gives

    d=sum_(p over y) ell_p
      =sum_(D_i over C)(e_i*f_i+sum_(p over y)eta_p)
      <=sum_(ALL D_i)e_i*f_i+Eta=d-1,

a contradiction. Normal meridian generation suffices; no chosen global
meridian generating list is assumed. This applies at EVERY finite y,
not only at a weighted origin. It allows f_i>1 and singular original N.
Consequently every local degree of M->A2 is at most b-1 and every
complete M fibre has at least TWO points. These are normalization fibres,
not source fibres of F; source surjectivity is not inferred.

## 2. New literature imports, with the exact hypotheses used

Pedersen, *Images of maps from (C2,0) to (Cn,0)*, arXiv2401.14941v3
(29 October2025), Theorem1.1: a NORMAL surface germ that is the image
of a finite holomorphic map from C2 is a finite quotient surface germ.
We import this implication, not a new proof of the link classification.
The corrected v3 statement and its proof context/Lemmas3.1--3.2 were
read in primary HTML; the earlier v2 proof is not our source version.

Giraldo--Roeder, *Pulling back singularities of codimension one objects*,
Proc. AMS148 (2020),1207--1217, Theorem A (arXiv1911.00628v1): if the
REDUCED inverse image of a hypersurface under a finite map between smooth
equal-dimensional germs is smooth, the target hypersurface is smooth.
The theorem and its explicit set-theoretic convention were read in primary
HTML, with Proposition4.1 and the dimension2 proof context. No nonvanishing
Jacobian hypothesis on the hypersurface is required. We use only dimension2.

Both are explicit theorem IMPORTS. The primary HTML bytes are pinned.
No fresh PDF visual audit, independent proof of either imported theorem,
withdrawn general pullback claim, or independent CFSG proof is asserted.

## 3. The critical-axis geometry descends to EVERY finite intermediate

Fix q in M and p in N above it. The original source germ is smooth or
has a small cyclic uniformization rho_N:C2->(N,p). Its composite finite
map to (M,q) is surjective as a germ. By Pedersen, M is a quotient germ.
Choose its standard small uniformization rho_M:C2->(M,q), with finite
group G free away from0. The map C2 minus0 -> M minus q lifts through
rho_M, since C2 minus0 is simply connected. Hartogs extends the two
coordinates across0. The lift h:C2->C2 is finite: h^-1(0)={0}, as
the original composite is a finite germ. Thus

    pi_N composed rho_N = pi_M composed rho_M composed h.

Let R be the reduced critical divisor of pi_M composed rho_M. If q is
a critical point, R is nonempty. Its reduced inverse image under h is
contained in the ONE original smooth critical axis, by the Jacobian
chain rule. Finiteness makes that inverse image a nonempty curve, hence
exactly the axis. Giraldo--Roeder now makes R smooth, including its
irreducibility. The group G preserves R. Equivariant coordinates
straightening R and linearizing the finite action exist by averaging,
as in PF106. The representation has an invariant line, hence two
one-dimensional summands. Smallness makes each character faithful;
therefore G is CYCLIC. No noncyclic quotient or multiple critical-axis
germ can occur in M.

At a noncritical q, purity makes pi_M etale near q, hence M smooth.
At a critical q let e be the generic intermediate ramification index.
With R={X=0}, the invariant target pullbacks satisfy

    Jac(P,Q)=unit*X^(e-1).

The EXACT local argument of PF106 Section2 applies to these actual
invariant functions: determinant character gives alpha=-e mod n,
gcd(n,e)=1; a nonzero coefficient X^(e-1) uses a monomial pair X^i,
X^jY with i+j=e and n|i. Hence

    marked type 1/n(1,-e),   2<=n<=e-1.              (1)

This is a verified use of that LOCAL lemma after proving its hypotheses,
not application of a whole-plane boundary theorem to M. It proves in
particular that every intermediate index2 critical point is smooth.

The reduced critical locus of M is locally one smooth quotient axis,
so its irreducible curves are smooth and pairwise disjoint. Each is
the finite image of an original omitted A1 normalization, hence itself
has normalization A1 and is an A1. Their image normalizations are A1.
No statement about a whole-plane open in M, its entire omitted boundary,
or global higher-index residue power maps is inferred.

## 4. Binary affine quotient: smooth intermediate and exact local rank

Suppose Q is a primitive binary affine quotient on b=2^a blocks.
PF122 gives ONE full quotient branch image C, transvection inertia,
noncontractible C, and ORIGINAL moved support w_C=d/2, with zero
extra movement inside fixed blocks. In particular a>=2, and d is a
multiple of4. The intermediate indices are all2, so Section3 makes
the entire M smooth. Its critical curves are disjoint A1s.

At each critical germ the Jacobian has a simple zero on a smooth axis.
PF119's LOCAL rank1 calculation therefore applies: at an image branch
of multiplicity m, the intermediate local rank is m+1, its normalized
residue degree is1, and its local group is S_(m+1). The image orders
are m,m+1. Thus each critical A1 maps UNRAMIFIED to C^nu=A1. A finite
unramified polynomial A1->A1 has degree1, so each maps isomorphically.
Transvection inertia then gives exactly b/4 critical primes in M.

Consider N->M over one of these critical curves E. Every inverse-image
curve is originally ramified over C, since its index is at least2.
In the original local germ its reduced preimage is the single critical
axis. Away from E this germ of N->M is unramified: another ramification
curve would make pi_N critical off that axis. The target pair (M,E)
is now a SMOOTH axis pair. Its local complement has cyclic fundamental
group, so a connected finite normal cover branched only over E has
the ordinary cyclic form

    (s,z)->(s^k,z).                                (2)

Consequently ORIGINAL N is smooth along every prime above C, and the
normalized residue map there has local degree1 everywhere. Each such
prime therefore has GLOBAL f_i=1, even though original indices need
NOT be2. Writing e_i=2k_i, degree conservation over each E gives

    sum_(D_i over E) k_i=t,    1<=k_i<=t.           (3)

Thus every visible original index is even and e_i<=2t=2d/b. There
are b/4 packets (3), and their union is exactly the original critical
primes above C. Invisible original branch images and their singularities
are not excluded.

## 5. The two ALL-INDEX profiles and the exact excess budget

At a chosen normalization parameter of C of branch multiplicity m,
the original local factors above each E have rank k_i*(m+1), while
e_i=2k_i. Hence the total original excess over this parameter is

    eta_C(parameter)=sum_i k_i*(m-1)=d/4*(m-1).

Put L=sum_parameters(m-1). Summing, and retaining ALL other primes,

    eta_C=(d/4)*L,
    outside_cost=d/2-1-(d/4)*L>=0.                 (4)

The outside cost is the sum of other original branch supports, ALL
unramified omitted-prime costs, and excess off C. It is not assumed
to be just a number of unramified primes. Equation (4) forces L=0 or1:

| Case | Image normalization | Original excess over C | Outside cost |
|---|---|---:|---:|
| I | Immersive everywhere | 0 | d/2-1 |
| II | Exactly ONE ordinary cusp (2,3) | d/4 | d/4-1 |

All higher cusps and multiple cusps are excluded for EVERY binary
affine block quotient, even when the original monodromy is imprimitive
and original visible indices are larger than2. This transfer is now
PROVED using Section3, rather than assumed as was forbidden in PF122.

At any double normalization value of C, the two branch fibres cost at
least2w_C=d. Thus both branches are smooth with zero excess on all
visible original primes; there are exactly two parameters, arbitrary
contact is allowed, and the original open source fibre is EMPTY.
No other omitted prime can meet that fibre. Since C is noncontractible
such double values exist. No unramified omitted prime can have image C,
as its finite map would supply an extra point over such a double value.
The ordinary cusp cannot share its value with another C parameter.

In Case II the intermediate local orbit profile is 3^(b/4)1^(b/4).
At the cusp the total original rank on visible critical primes is
3d/4; the remaining rank d/4 may include other omitted points. No
exact retained-source cusp count or surjectivity is claimed.

## 6. Genuine local controls keep the remaining distinction visible

For EVERY k>=1 the actual finite map

    g_k(x,y)=(x^k-3y^2, x^k*y-2y^3)

has Jacobian k*x^(2k-1), one smooth critical axis, index e=2k, ordinary
cusp image, local rank3k and excess k. It factors through (s,z)=(x^k,y)
and the ordinary simple cusp germ. Subtracting y times its first fibre
equation from the second gives y^3; the local algebra is free of rank k
over C[y]/(y^3), so has length3k. This is an UNBOUNDED construction.
It is NOT an actual Keller normalization with the whole-plane marked
open: its unramified open is C* times C. It confirms that e_i>2 is
compatible with the new local geometry and must not be declared excluded.

Every intermediate has the local cyclic-axis restriction (1), but a
binary affine intermediate is only proved SMOOTH, not polynomially
isomorphic to A2. Published quartic PLANE inversion therefore still
does not retire a degree4 intermediate. Natural alternating/symmetric
quotients, binary affine global compatibility, invisible original branch
images, global trace/multiplication and full B=C[x,y] remain OPEN.
All original primes, sealed predecessors, shared priorities and quartic
retirement are retained. No historical priority, independent proof of
the imported geometric theorems, formal verification, external review
or full goal completion is claimed.
