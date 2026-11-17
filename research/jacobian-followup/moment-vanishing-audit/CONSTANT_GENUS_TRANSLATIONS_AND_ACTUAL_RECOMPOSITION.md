# Constant-genus translations and actual recomposition

16 November 2026. Internal written proofs. Independent replay, formal
verification and external review are not claimed. Both COMPLETE exceptional
Mathieu--Zhao cases remain **OPEN**; review stays deferred.

Subsequent result,16 November 2026: the [intrinsic safety theorem](INTRINSIC_SOURCE_SPLITTING_AND_FULL_HAAR_SAFETY.md)
closes the residual K gate under EVENTUAL PURITY. Mixed fixed part makes
the second SOURCE quotient admissible; its finite LOCAL monodromy
gives an actual Borel-unit retraction without a Hodge-compatible actual
kernel or a finite GLOBAL monodromy inference. Source unit splitting
is now necessary AND sufficient. The account below records the prior
reduction and its then-open limitation; its abstract controls were not
actual eventually pure phases. Both FULL exceptional cases stay OPEN.

Every phase now has two exhaustive alternatives. A nonsplit GLOBAL unit
in its canonical boundary SOURCE certifies failure whenever the phase is
eventually pure. If that source unit splits, all nonconstant genus can be
removed while preserving the full actual unit problem. The residual
source has CONSTANT bottom space and FINITE top monodromy. On a finite
base cover it is an abelian translation system, with radical square zero.
Only the constant boundary-genus fixed part can carry residual motion.

The actual Cauchy FUNCTION connection, not just its geometric source,
also preserves and reflects global semisimplicity through EVERY polynomial
recomposition. A semisimple actual primitive system of rank>1 admits NO
eventually pure polynomial composite, in any degree.

See the [checker](verify_haar_constant_genus.py),
[receipt](haar_constant_genus_verification.json) and
[source audit](haar_constant_genus_source_audit.json).

## 1. The global intrinsic source-unit gate

Fix f on X=SL2/Gamma, ANY finite Gamma. Retain the FULL regular base and
the canonical boundary quotient E_f=B_boundary from the
[boundary reduction theorem](CANONICAL_BOUNDARY_REDUCTION_AND_LOWEST_WEIGHT_INVARIANTS.md).
It has original weight4 ambient line L=Q(-2), weight5 pure genus system
P_f=GrW5 E_f, and weight6 pure Tate system T_f. The pure grades are
globally semisimple; T_f has finite global monodromy. Let u denote the
ambient generator in a fibre. Call the intrinsic global gate SPLIT if

    some global flat p:E_f -> C satisfies p(u)=1.         (1)

This is a retraction of the rational unit connection; p need not be a
Hodge morphism. It is equivalent to a globally fixed ambient lift in
the FULL dual source, with no physical-orbit condition yet.

If (1) fails, no quotient of E_f retaining u can have such a retraction:
compose one with the quotient map to contradict failure of (1). The actual
boundary period quotient retains u. Therefore

    intrinsic source-unit NONSPLIT + eventual purity
       implies full ALL-original-marker FAILURE.        (2)

This is a global obstruction and can persist even when each chosen local
extension has a fixed lift. If (1) holds, the actual quotient may still
identify u with another radical direction; that is the residual case
handled next. A full-source fixed lift outside the physical orbit is not
silently substituted for the actual Haar lift.

The gate (1) reflects EVERY nonconstant H before purity. Indeed the
reduced composite source is finite induction of E_f modulo trace-zero
ambient L. A primitive retraction induces one by summing its values on
all branches; it kills those trace-zero ambient relations. Conversely a
composite retraction restricts to the injected distinguished primitive
branch, preserving p(u)=1. Thus

    source unit of E_f splits iff source unit of E_H(f) splits. (3)

No branch is discarded, no primitive purity is assumed, and a nonsplit
source alone does not construct an eventually pure H(f).

## 2. Split-source case: remove ALL nonconstant genus

Assume (1). Its restriction splits the weight4 line in W5 E_f. Hence

    W5 E_f is globally semisimple, isomorphic to L direct_sum P_f.

Its canonical trivial-isotypic part and complement are

    F_f=(W5 E_f)^global,
    C5_f=sum_g image(g-I)|W5 E_f,
    W5 E_f=F_f direct_sum C5_f.                         (4)

Write c_f=dim(P_f^global). Then dim F_f=1+c_f. The module C5_f is
semisimple with NO trivial subquotients. This decomposition concerns
GLOBAL local systems. No new Hodge splitting of mixed weights is assumed.

Put

    V_f=E_f/C5_f,   0 -> F_f -> V_f -> T_f -> 0.         (5)

F_f is a constant local system and T_f has finite monodromy. This
regular algebraic connection quotient need not be a newly constructed
geometric MHS or algebraic 1-motive presentation.

Let M_boundary be the actual quotient from the preceding theorem.
Quotient it by the image of C5_f, obtaining N_f. That image cannot
contain u and has no homomorphism to the trivial connection. Therefore

    Hom_D(N_f,C(s))=Hom_D(M_boundary,C(s)),              (6)

with identical evaluation on u. Together with the earlier compact-middle
removal this retains the EXACT full all-marker criterion for EVERY phase
in the split-source branch. No scalar-function realization of N_f is
claimed. The actual Haar period kernel remains part of its definition.

## 3. A constant translation-span criterion, with radical square zero

Take the finite cover killing the monodromy of T_f. For a loop h on that
cover, (T_h-I) kills F_f and maps V_f into F_f, inducing a well-defined
linear map T_f -> F_f independent of a choice of lifts. Define

    R_f=sum_h image(T_h-I) subset F_f.                  (7)

Every product (T_h-I)(T_k-I) is zero on this cover, so these translations
commute. The quotient V_f/R_f has finite monodromy. Conversely the
translation operators generate a square-zero ideal in the full monodromy
algebra. Since finite-image representations are semisimple in
characteristic zero, (7) is EXACTLY the image of its Jacobson radical J:

    J V_f=R_f,    J^2=0.                               (8)

Because the source unit still splits after removing C5_f, u is NOT in
R_f. In particular dim R_f<=c_f. Let

    K_f=F_f intersect ker(V_f -> N_f).                  (9)

The actual unit survives, so u is not in K_f either. Passing to the
actual quotient sends the radical image to the image of R_f. A proof
without any exactness assumption on invariants is: N_f/image(R_f) is a
quotient of the finite-monodromy V_f/R_f, while the images of the cover
translations still generate a square-zero ideal. Thus, under eventual
purity, the complete criterion is

    ALL original markers safe iff u NOT in R_f+K_f.     (10)

This is a condition between CONSTANT finite-dimensional spaces. The
actual kernel K_f cannot be replaced by a convenient Hodge subspace or
arbitrary fixed lift. It need not respect the original mixed weights.
We have not computed it universally or replaced an identity by a finite
moment prefix.

If c_f=0 then F_f=C*u, and u not in R_f forces R_f=0. Therefore, allowing
arbitrary rank and genuinely varying boundary genus,

    P_f has NO global fixed part:
    eventual-pure full safety iff intrinsic source unit SPLIT. (11)

This generalizes the genus-weight-zero criterion. It does NOT say that
every such phase is safe: the nonsplit alternative (2) remains real.
For arbitrary c_f, every residual failure AFTER intrinsic splitting
requires a nonzero CONSTANT genus fixed part and an actual relation
placing u in R_f+K_f. Varying genus no longer supplies an independent
obstruction in that branch. Large constant spaces alone do not imply
failure; the previous pure SAFE examples remain intact.

## 4. The constant data do not grow with polynomial recomposition

Strictness of finite-cover recomposition gives P_H(f)=pi_*P_f. Thus

    c_H(f)=dim H0(T,P_f)=c_f.                           (12)

Deleting extra regular values does not change primitive monodromy image.
In the split-source case the fixed bottom space F_H(f) is canonically
identified with F_f by branch SUM. More precisely

    V_H(f)=pi_*V_f / (F_f tensor trace-zero permutation space). (13)

To see this, induce (4). Induced C5_f has no global invariants. The
nonconstant permutation copies of the fixed GENUS part also belong to
the composite nontrivial isotypic complement. Removing that complement
leaves (13), which identifies ALL copies of F_f, not just the ambient
line. Every primitive branch still injects after restriction.

Under that fixed-bottom identification,

    R_H(f)=R_f.                                        (14)

For one inclusion, V_f/R_f has finite monodromy; finite induction and
quotients preserve finite image, so the radical image of (13) lies in
R_f. Conversely realize any primitive kernel-of-finite-monodromy loop
on the primitive cover and take a finite positive power killing every
composite permutation/top action. On the injected primitive branch the
resulting translation is a nonzero integer multiple of the original
one. Hence R_f lies in the composite translation span over Q or C.

Thus the SAME F_f and R_f govern all H degrees in the split-source
case. The actual relation subspace K_H(f) is retained separately; no
claim that it is automatically constant under H or explicitly produced
is made. Conditions (3), (11) and (12) give an intrinsic all-degree
dichotomy for every fixed f with no fixed genus part.

## 5. The ACTUAL Cauchy connection reflects global semisimplicity

Return to the full scalar-function connection M_f, of rank r_f. For
EVERY H of degree d, before purity, there are a surjection and injection

    pi_* M_f ->> M_H(f),
    M_f -> pi^* M_H(f).                                 (15)

These are about ACTUAL functions and their actual relations, not only
the abstract complement source. Here is an explicit proof.

Write r=H(tau). For deg a<d, partial fractions give

    C_(Q*a(f))^H(f)(r)
       =sum_(H(tau_j)=r) a(tau_j)/H'(tau_j)*C_Q^f(tau_j). (16)

The d elements a(tau)/H'(tau), for a=1,tau,...,tau^(d-1), form a basis
of C(tau) over C(r). Thus trace of the restriction-of-scalars module
lands in exactly M_H(f) and is onto. All actual function relations
remain relations on every branch. For the injection use the retained
divided-difference reconstruction

    C_Q^f(tau)=sum_(j<d) r_j(tau)*C_(Q*f^j)^H(f)(H(tau)),
    sum_j r_j(tau) Z^j=(H(tau)-H(Z))/(tau-Z).             (17)

It embeds the literal primitive scalar functions into the pullback
connection. Both identities hold in convergent physical germs at
infinity and then along compatible continuations. No new amplitude
subalgebra is substituted.

Finite-index induction/restriction preserve semisimplicity. Applying
this to (15) proves

    M_f globally semisimple iff M_H(f) globally semisimple. (18)

Trace on the rank-d rational subconnection in pi_*M_f has kernel d-1.
Consequently

    r_f <= r_H(f) <= 1+d*(r_f-1).                        (19)

Rank cannot drop. A rank-one M_f is precisely the rational connection;
the retained actual-image/indecomposability theorem then gives a unique
safe centre c with f-c all-marker safe. Hence:

* M_f semisimple of rank>1: NO eventually pure H(f) exists, in ANY H
  degree. Otherwise (18) and eventual purity make M_H(f) safe and
  rank1, contradicting (19).
* M_f rank1: eventual-pure H(f), and full safety, are exactly H(c)=0
  for that unique safe centre; no assertion of normalized full purity.
* M_f nonsemisimple: EVERY eventually pure H(f), IF ANY, is unsafe.

These statements keep all original markers. They give no general
producer of actual monodromy or pure composites and do not reflect
purity from H(f) to f-a. In particular a primitive source for a full
eventual-pure failure must already have NONSEMISIMPLE ACTUAL Cauchy
monodromy, before any primitive purity hypothesis.

## 6. A genuine pure failure occupies the intrinsic nonsplit branch

Use Long's known SU2 phase in determinant coordinates

    f=(1+y)(1+z+yz),   Q=-y,   xw-yz=1.

The primary Theorem4.2 gives purity and h(Q f^n)=(-1)^(n-1)/(n+1)
for every n>=1. The following fibre geometry is a new exact check of
where that already-known failure lies; it is NOT an exceptional descent.

Put a=1+y. For s not0,-1, its fibre satisfies

    z=(s-a)/a^2,   a!=0,
    U=a*x, V=a*w,   U*V=(s+1)*a-s.

Therefore

    F_s ~= A2_(U,V) minus {U*V=-s},
    a=(U*V+s)/(s+1).                                   (20)

The inverse is x=U/a,w=V/a,y=a-1,z=(s-a)/a^2. It obeys both equations
and is regular on the indicated open sets. The deleted hyperbola is
C*, so H2(F_s)=Q(-2) by the Gysin sequence in A2. Thus compact-middle
and genus weight3 vanish, while the Tate graph rank is1. The full
complement rank is2. Its actual functions include

    C_1=1/s,
    C_Q=1/s-log(1+1/s),
    d(C_Q-1/s)/ds=1/(s*(s+1)).                           (21)

The last derivative has nonzero residues at0 and-1, hence C_Q is not
rational. These two functions already have rank2, equal to the full
source rank. Their connection is a nontrivial unipotent extension of
the rational unit, with the intrinsic source-unit gate NONSPLIT.

This rules out any shortcut claiming positivity merely from zero genus
weight or from the virtually abelian form. The actual Long phase is
NOT right2O/2I invariant, and prior unary/regularity barriers still forbid
promoting it to an exceptional phase. Zero-constant composites give the
already-proved all-degree failures, without a deeper moment campaign.

## 7. The remaining whole-quotient target

For EVERY actual exceptional eventually pure phase, prove first that its
intrinsic GLOBAL boundary source unit splits, and then that u stays
outside R_f+K_f in the constant-genus reduction, or realize a genuine
failure. These alternatives are exhaustive. Compact-middle and varying
genus no longer need independent campaigns; finite recomposition cannot
create a new intrinsic split status, fixed-genus rank or translation
range. The actual relation space and exceptional regularity remain.

No algebraic1-motive presentation, universal monodromy/relation producer,
coefficient-ring topological promotion or universal machine verification
is claimed. All phases, original markers, branches and H degrees remain
in scope. Both FULL exceptional MZ cases remain OPEN; review deferred.
