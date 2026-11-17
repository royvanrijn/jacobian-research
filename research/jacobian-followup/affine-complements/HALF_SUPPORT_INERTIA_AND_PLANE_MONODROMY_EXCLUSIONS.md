# Half-support inertia and plane monodromy exclusions

PF2RULED118, 16 November 2026. Every hypothetical nonautomorphic complex
PLANE Keller map has a nonidentity divisorial inertia element moving at
most half its generic sheets. This excludes whole faithful permutation
actions, including both exceptional eleven-sheet PSL2(11) actions,
natural projective-linear actions, and nonbinary affine-linear actions.
The binary affine case has a precise surviving boundary profile below.
No full Jacobian-conjecture proof or exclusion in higher source dimension
is asserted. The full goal remains ACTIVE.

## 1. A global support budget

Use the ORIGINAL finite normalization pi:N->A2 of an actual polynomial
Keller map with whole marked open U=A2. Its generic field degree is d>1.
Retain EVERY omitted prime, including unramified primes and primes with
the same target image. Let e_i and f_i be transverse index and normalized
residue degree, and eta_i the total normalized local excess along that
prime. PF2OMIT6 gives

    d-1 = sum_i e_i f_i + sum_i eta_i,   eta_i>=0.      (1)

Let C_1,...,C_k be the DISTINCT irreducible components of the FULL finite
map branch curve. These are images of ramified omitted primes only;
extra unramified nonproper images need not belong to this branch curve.
For a general geometric meridian sigma_j around C_j put

    w_j = |supp(sigma_j)| = sum_(i: image D_i=C_j) e_i f_i.

Indeed its nontrivial cycles have length e_i, repeated f_i times.
Generic residue characteristic is zero. Unramified omitted primes
contribute fixed sheets, not moved sheets. Consequently the exact
monodromy form of (1) is

    sum_(j=1)^k w_j + t0 + Eta = d-1,                 (2)

where t0 is the number of unramified omitted primes (each has f=1 by
PF105/PF2OMIT2) and Eta=sum_i eta_i. In particular sum_j w_j<=d-1.
The full branch curve is nonempty: otherwise the finite normalization
would be etale over simply connected A2 and have d=1.

If k>=2, some w_j <= floor((d-1)/2). If k=1, the full-branch mass
theorem PF113 gives w_1<=floor(d/2). Its proof is worth recalling:
2w_1>d forces the image normalization A1->C_1 to be injective, since
two distinct normalization points at one finite value would consume
at least 2w_1 sheets, including their actual local excess. The image
is then contractible. The Lin--Zaidenberg theorem makes it polynomially
equivalent to x^a=y^b. The entire finite branch complement contracts by
a positive weighted radial flow onto a small local complement. Its
globally transitive cover would therefore have transitive LOCAL
monodromy, contrary to PF113's at-least-two complete-fibre factors.
This uses the FULL branch locus and all omitted primes; it does not
assume that the original open Keller map is a finite cover outside it.

Thus in ALL cases some branch meridian satisfies

    2 <= |supp(sigma_j)| <= floor(d/2),
    |Fix(sigma_j)| >= ceil(d/2).                       (3)

This is a corollary of the retained global ledger and PF112--113,
not a new derivation of those geometric inputs. The prior HC4MYKP2
node support screen is much more local; it did not exclude the global
exceptional actions. No new linear trace-image refinement is used here.

For a faithful transitive action G on d letters, define its minimal
degree mu(G)=min_(g!=1)|supp(g)|. Equation (3) proves the all-degree
obstruction

    mu(G)>floor(d/2)  =>  NO plane Keller realization. (4)

More quantitatively k*mu(G)+t0+Eta<=d-1. If d is even and mu(G)=d/2,
the only possible case is

    k=1,   w_1=d/2,   t0+Eta=d/2-1.                  (5)

Every branch inertia then realizes the minimum support. Equation (4)
does not exclude the equality case, alternating/symmetric actions in
general, semilinear actions not checked by the bounds below, or every
map of any particular surviving field degree.

## 2. Both exceptional eleven-sheet PSL2(11) actions are excluded in the plane

The retained complete KMFI95 mark calculation gives, in BOTH exceptional
actions, the following cyclic fixed-sheet counts:

| Element order | 2 | 3 | 5 | 6 | 11 |
|---|---:|---:|---:|---:|---:|
| Fixed sheets | 3 | 2 | 1 | 0 | 0 |

These exhaust the nonidentity element orders of PSL2(11). Hence the
minimum support is 11-3=8>floor(11/2)=5. Equation (4) excludes BOTH
exceptional degree11 actions for every ACTUAL polynomial PLANE Keller
map, without assumptions on knots, target singularities, omitted-prime
count, residue degrees or normalization smoothness.

This is a plane exclusion. KMFI95's geometric and arithmetic local
mark statements remain valid; its no-absolute-action-exclusion flag
concerns arbitrary affine dimension and is not overwritten. The natural
PSL2(11) action has degree12 and is a DIFFERENT action, covered by Section3.
Jones--Zvonkin, arXiv2104.12015v2, Section3, supplies primary context for
the two A5-stabilizer actions and their dessin passports. The complete
cyclic mark input is the sealed KMFI95 census, reused without rerunning
its subgroup enumeration. No claim of historical priority is made.

## 3. All natural projective-linear actions

Let G be any faithful transitive subgroup of PGL_n(F_q) on the points
of projective (n-1)-space, with n>=2. Put

    d=1+q+...+q^(n-1),    H=1+q+...+q^(n-2).

A projective point is fixed by a matrix precisely when its line lies
in an eigenspace for an eigenvalue in F_q^*. For a nonidentity projective
matrix every eigenspace has dimension at most n-1. If q=2, the only
possible nonzero rational eigenvalue is1, so at most H points are fixed.
Then d=2H+1 and every nonidentity element moves at least H+1>d/2.

For q>=3, distinct eigenspaces are independent. Their dimensions sum
to at most n. The number of eigenlines is largest for dimensions n-1
and1, hence at most H+1. This elementary convexity bound follows by
merging positive dimensions; for a,b>=1,

    q^(a+b-1)+q-q^a-q^b
      = q (q^(a-1)-1)(q^(b-1)-1) >=0.

If there is only one rational eigenspace the stronger H bound holds.
Since d=qH+1, the H+1 bound is strictly below d/2 exactly when
(q-2)H>1. This holds for every q>=3,n>=2 except q=3,n=2.
The same argument applies to hyperplanes by the contragredient action.

Thus equation (4) excludes ALL these faithful transitive projective-linear
point and hyperplane actions for

    q=2,n>=2;    q=3,n>=3;    q>=4,n>=2.               (6)

In particular every natural PSL_n(q) and PGL_n(q) action in these ranges
is excluded, with no prime-degree or congruence hypothesis. For PSL2(q)
on q+1 points this extends beyond the already retained ALL-DIMENSION
KMFI67 obstruction for odd q=3 mod4; that older result is credited.
It also excludes the natural projective examples of degrees7,13,17,31,73
in Jones--Zvonkin, without changing their existence as genuine dessins.
The sole exceptional parameter q=3,n=2 has degree4: its plane Keller
realizations are already closed by published Domrina2000/HC4MYKR0.
The half-support test itself does not exclude S4 in that action.

## 4. Nonbinary affine actions, and the binary case which survives the test

For a nonidentity affine map z->Az+b on F_q^n, its fixed-point set is
empty or an affine subspace of dimension at most n-1. Thus in degree
d=q^n it fixes at most d/q points. If q>=3, every nonidentity moves
at least d-d/q>d/2. Equation (4) excludes EVERY faithful transitive
subgroup of AGL_n(F_q) in its natural affine action for q>=3,n>=1.
Here q is the field SIZE. Fields F_4,F_8,... are included; binary means
the base field F_2, not all fields of characteristic two.

For q=2 every nonidentity moves at least d/2. If a transitive subgroup
has strictly larger minimum support it is excluded. Otherwise (5)
forces an irreducible full branch curve C and half-support inertia.
For n>=2 a half-support affine transformation has a fixed hyperplane;
translation to a fixed point conjugates it to I+v*l, where l is a
nonzero linear functional, v!=0, and invertibility over F_2 forces
l(v)=0. It is an affine transvection, an involution with cycle type

    2^(d/4) 1^(d/2).

In an ACTUAL plane realization PF105's e=2 theorem gives f_i=1.
There are therefore EXACTLY d/4 ramified omitted primes, all of index2
and all with the same image C. PF106's n_quotient<=e-1 rules out a
singular quotient point on them. Unramified omitted primes are etale
at all finite points. Since U=A2 is smooth and these exhaust the
boundary, the original N is SMOOTH everywhere. Each ramified curve
is a smooth A1 with normalized image map an isomorphism, and

    t0+Eta=d/2-1.                                    (7)

The curve C cannot be contractible by PF113. Its normalization is A1,
so it has a finite value with at least two normalization preimages.
If b_C(y) is their number at a finite value y, the d/4 disjoint ramified
curves contribute at least 2*(d/4)*b_C(y) sheets. Hence b_C(y)<=2.
At ANY two-branch value all d sheets are consumed: the finite fibre
has exactly d/2 points, each of local length2 on a ramified curve,
and the open source fibre is EMPTY. There are no other normalization
points over y. Local excess is zero on both branches, so PF115's
positive excess for a singular image branch forces both branches to
be smooth. Their contact order is unrestricted; a node is not assumed.

No omitted unramified prime can have image C: it would add points over
these double values. Every other unramified image must avoid ALL such
values. This does NOT prove source surjectivity; the empty source
fibres at these values are a necessary part of the surviving profile.

AGL_5(2) in degree32 illustrates the boundary of the new exclusion:
its minimum support is16, so it is NOT excluded. Any actual plane
realization would have an irreducible C, eight smooth ramified A1
curves of index2/residue1, smooth N, and t0+Eta=15, with the above
double-value restrictions. The n=2 degree4 case is already excluded
by the published quartic theorem. No exclusion of the binary affine
family in general or of all degree6 maps is claimed.

## 5. Verification and the next global gate

The universal support theorem uses the sealed PF112/113 geometry and
PF2OMIT6 ledger, not finite group samples. The projective and affine
bounds are elementary dimension arguments written above, not a finite
field enumeration or a classification of finite simple groups. The
checker reuses the sealed KMFI95 cyclic marks, verifies small exact
matrix actions and the binary transvection boundary, and retains
countercontrols for S_d, A_d and the equality case. It does not rerun
the 620-subgroup census or start a parameter campaign.

Primary HTML and retained source hashes are pinned in
sources/half-support-inertia/manifest.json. PF117's local trace comparison
remains closed and automatic for finite covers. The NEW conclusion is
an all-degree class exclusion from WHOLE-plane geometry, not global
trace-moment regularity or equality of the finite normalization with
C[x,y]. The binary affine profile, actions with smaller minimum support,
general target branch arrangements, and the global finite-algebra
finish line remain OPEN. Shared trace priorities, published quartic
retirement, all omitted primes and closure before external review
are preserved. No formal verification or external review is claimed.
