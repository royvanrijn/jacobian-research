# A central-cover obstruction: absolute A4 and natural PSL2(q)

Work geometrically over an algebraically closed field k of characteristic0.
For a polynomial Keller self-map F:A^d_X -> A^d_Y write
K=k(Y),L=k(X),M/K its generic Galois closure, G=Gal(M/K), and
H=Gal(M/L). Coefficients C_l below have trivial G action.

**KMFI65 (central prime-cover obstruction).** For every prime l with
l not dividing |H|,
\[
                         H^2(G,C_l)=0.                     \tag{1}
\]
Equivalently, no non-split central extension
1 -> C_l -> Gtilde -> G ->1 can exist for such a prime.
This is a necessary condition on the faithful transitive action
itself, in every affine dimension. It is stronger than the
normally-generating-stabilizer condition in these examples.

**KMFI66 (absolute A4 exclusion).** A4 is not the geometric generic
monodromy group of any polynomial Keller self-map of ordinary affine
space in characteristic0, in **any faithful transitive action**.
In particular the natural degree4 action is excluded absolutely,
not merely through a chart or a curve-induced presentation.

**KMFI67 (absolute natural PSL2(q) exclusion).** For every odd prime
power q congruent3 modulo4, PSL2(q) on P1(F_q) cannot be geometric
generic monodromy of a polynomial Keller self-map in any dimension.
This includes natural PSL2(11) of degree12 and PSL2(7) of degree8.
The exceptional degree11 PSL2(11) actions remain undecided:
their A5 point stabilizers have even order, so(1) at l2 does not apply.

**KMFI68 (complete action spectrum through degree4).**
The possibilities are the identity action in degree1, none in degree2,
natural S3 in degree3, and natural S4 in degree4.
Both nontrivial possibilities have actual polynomial affine-space
witnesses from KMFI11,14,23. For a realization over any characteristic0
coefficient field, the arithmetic group is also S3 or S4 respectively
and there are no additional normal-closure constants.

[Exact group and cohomology helper](central_lift.py) ·
[checker](verify_central_lift.py) · [model](central_lift_model.json) ·
[receipt](central_lift_receipt.json) ·
[pinned Brauer source](sources/central-lift-manifest.json).

## 1. Inputs and the scope of the argument

KMFI1 gives two global facts: every affine divisorial inertia of M/K
is contained in a conjugate of H, and A^d over k has no nontrivial
connected finite etale cover. Purity extends a finite normalization
that is unramified in codimension1 to such an etale cover.

For the Brauer step we use Kummer
H^2(K,mu_l)=Br(K)[l], codimension-one purity for Br on a regular
scheme, and Br(A^d_k)[l]=Br(k)[l]=0.
These precise inputs are Section1.3.4, Corollary3.7.3 and
Theorem5.1.1 of the pinned 31 October2019 author draft of
[Colliot-Thelene--Skorobogatov](https://www.imo.universite-paris-saclay.fr/~jean-louis.colliot-thelene/BGgroup_book.pdf).
The local unramified interpretation uses their Theorem3.4.2.
These are cited geometric inputs, not results certified by the
finite group checker.

The central-extension, local splitting, Kummer twist and source
contradiction are proved below. No assertion concerns only a chosen
construction, and no condition at infinity is silently imposed.

## 2. A central extension over prime-to-l inertia is unramified

Assume for contradiction that a non-split central extension
\[
             1\longrightarrow C_l\longrightarrow Gtilde
                         \longrightarrow G\longrightarrow1          \tag{2}
\]
exists, with l prime and l not dividing |H|.
Choose mu_l=C_l, possible since k is algebraically closed.
Its class a in H^2(G,C_l) pulls back along the generic Galois
representation rho:Gamma_K -> G to
\[
                        b\in H^2(K,mu_l)=Br(K)[l].
\]

At an affine prime divisor D of Y, let D_G and I_G be a decomposition
group and its inertia subgroup in G.
KMFI1 places I_G in a conjugate of H, so |I_G| is prime to l.

A central extension of a finite group I of order prime to l by C_l
has a unique complement to its central kernel.
Here is an elementary cohomological proof. In additive C_l notation,
let c(g,h) be a normalized2-cocycle and define
\[
                b_I(g)=|I|^{-1}\sum_{t\in I}c(g,t).
\]
Summing the cocycle identity gives
c(g,h)=b_I(g)+b_I(h)-b_I(gh), so it is a coboundary.
Two splittings differ by a homomorphism I -> C_l, necessarily zero.
Thus the complement is unique, and is normal in every overgroup
normalizing the full preimage of I.

Apply this to I_G inside the preimage of D_G in Gtilde.
Quotienting by that unique complement gives a central extension
of D_G/I_G by C_l. Pulling it back recovers the original extension
over D_G. The local Galois representation into D_G/I_G is unramified.
Consequently b is inflated from the residue field at D, and has
zero Brauer residue there.

This holds for **every** affine target prime D. Brauer purity and
the vanishing on A^d give b=0 globally.
No claim that the Brauer group of K itself vanishes is being made.

The vanishing is exactly the weak embedding-problem condition.
Indeed a set-theoretic section of(2) gives its2-cocycle c; a continuous
cochain killing rho^*c changes that section into a continuous
homomorphism rhotilde:Gamma_K -> Gtilde lifting rho.
Its image projects onto G. If it missed the prime-order kernel,
its image would be a section of(2), contradicting nonsplitting.
Therefore this lift is surjective.

## 3. Remove all new central affine ramification

The first lift may introduce l-ramification not present in M.
This must be removed before using the source; omitting this step
would leave the obstruction unproved.

At an affine D, the preimage of I_G is the product of its unique
prime-to-l complement with C_l. The central component of rhotilde
on the procyclic tame inertia is a character to mu_l.
It is a power of the Kummer inertia character of a local uniformizer.
There are only finitely many D where it is nontrivial, because
the finite lifted extension has finitely many ramified divisors.

The target coordinate ring k[Y] is factorial. Choose a reduced
polynomial d_D for each such prime and exponents a_D modulo l so
that the Kummer character of
\[
                              f=\prod_D d_D^{a_D}
\]
has the opposite central inertia at every listed D.
Twist rhotilde by the global character of the lth root of f.
At each selected D its central part cancels; at every other affine
prime v_D(f)=0, so no new divisorial inertia is introduced.

The resulting lift still projects onto G and remains surjective
by nonsplitting of(2). Let N/K be its Gtilde-Galois field.
All its affine target inertia groups now have order prime to l.
Ramification at infinity is irrelevant and has not been suppressed.

## 4. The affine source forces a contradiction

Let Htilde be the preimage of H in Gtilde. Since |H| is prime to l,
the same unique-complement argument gives
\[
                           Htilde=C_l\times H.
\]
Hence it has a quotient C_l nontrivial on the central kernel.
The corresponding intermediate field E of N/L is a nontrivial
cyclic extension of degree l.

Every prime divisor of the affine source X maps dominantly to
an affine prime divisor of Y: F is etale and quasi-finite.
The inertia of N/L there is a subgroup of an inertia group of N/K.
By Section3 it has order prime to l, so its image in the C_l
quotient is trivial. Thus E/L is unramified at every affine
source prime.

Normalize X=A^d in E. This is finite and unramified in codimension1.
Purity makes it finite etale everywhere. Simple connectedness of
affine space then gives E=L, contrary to [E:L]=l.
This proves(1), with no finiteness assumption on F.

The argument survives all algebraically closed characteristic0
constant extensions, identity stabilizations and arbitrary field
presentations of the same action. It is an absolute action obstruction.

## 5. The natural projective PSL2 family

For odd q there is a central extension
\[
       1\to\{I,-I\}\to SL2(F_q)\to PSL2(F_q)\to1.              \tag{3}
\]
It does not split. In odd characteristic a determinant-one matrix
whose square is I is diagonalizable with eigenvalues in{1,-1};
determinant1 forces it to be I or -I.
Thus -I is the unique nonidentity involution of SL2(F_q), and it
maps to the identity in PSL2(F_q).
But the projective class of [[0,-1],[1,0]] is a nonidentity involution.
A section would lift it to a noncentral involution, impossible.

The stabilizer of infinity in the natural action has order
q(q-1)/2. This is odd exactly for q congruent3 modulo4.
Equation(1), with l2 and the non-split cover(3), excludes every
action in KMFI67, for all odd prime powers in this congruence class.
No formula for a general Schur multiplier is required.

For q3 the projective group is A4 on four letters; the checker
identifies the entire twelve-element permutation group.
This proves the natural A4 exclusion.

To cover **all** faithful transitive A4 actions, its possible
core-free stabilizers are1,C2,C3. Stabilizer1 cannot normally
generate A4, and the normal closure of C2 is only V4, so KMFI1
already excludes their degree12 and6 actions. The C3 case,
degree4, is excluded by the central cover SL2(F3).
The remaining proper subgroup V4 has nontrivial core and does
not define a faithful A4 action. This proves KMFI66 in full.

## 6. Checks against actual occurrence and hypothesis boundaries

The actual binary eight-letter group SL2(F3) also has an odd
point stabilizer C3. It is **not** excluded: the exact normalized
cochain computation gives H^2(SL2(F3),C2)=0, whereas
H^2(A4,C2) has dimension1. This is a useful control against
confusing a group with its central quotient.

Natural A5 on5 is likewise untouched: its stabilizer A4 has even
order. The exceptional degree11 PSL2(11) stabilizers have order60,
so the l2 hypothesis fails. Their preimages of order120 in SL2(11)
are perfect, contrasting with the quadratic character on the natural
stabilizer preimage of order110. The checker verifies both facts;
it does not claim an exceptional-action exclusion.

Coprimality is essential locally: in C4 -> C2 neither lift of
the nonidentity element has odd order. A quadratic twist cannot
convert these order4 lifts into order2 lifts. The corresponding
order-l^2 examples give the same warning for general l.

Ordinary affine space is essential globally. Smooth affine boundary
complements can have nontrivial etale cyclic source covers, and
their target Brauer group or available divisor equations need not
have the properties used above. The earlier A4 chart construction
therefore does not contradict the theorem.

## 7. Complete small-degree spectrum

For degree4, a transitive subgroup of S4 has order4,8,12 or24.
The respective possibilities are C4 or V4, D8, A4 and S4.
The2-groups are excluded by KMFI2; A4 is excluded above.
S4 has actual witnesses in KMFI14 and23.
In degree3, C3 is excluded by KMFI2 and S3 is realized by KMFI11.
Degree2 is the excluded regular C2 action; degree1 is identity.
This proves the complete table in KMFI68.

Over an arbitrary characteristic0 coefficient field, geometric
S3 or S4 already fills its symmetric ambient group, so the
arithmetic group is the same and has no extra constant quotient.
The existence statements here allow the coefficient fields of the
retained actual witnesses; no new descent to a smaller field is asserted.

The written proof uses the pinned Brauer and purity inputs.
Exact finite matrices, cohomology ranks and subgroup enumerations
check consequences and boundaries. Formal verification, novelty
priority and independent whole-proof review are not claimed.
The exceptional degree11 PSL2(11) occurrence and full action
classification remain open.
