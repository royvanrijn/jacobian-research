# Automorphisms preserve the normalized trace on every affine reductive homogeneous space

Let K be any field of characteristic zero, let G be a connected reductive
linear algebraic K-group, and let H be a closed K-subgroup such that
X=G/H is affine. The quotient is the scheme/fppf quotient; no transitivity
of G(K) on X(K) and no splitness assumption is made. Put A=O(X),
g=Lie(G), h=Lie(H), and n=dim X.

There is a unique normalized G-invariant functional Lambda_X:A->K.
The proof establishes the volume lemma in Sections1--13, then uses
the intrinsic determinant cover in Sections14--21. The complete
volume argument is included below.

**Theorem AffineReductiveHomogeneousAllAutomorphismTrace.**
For every K-algebra automorphism sigma of A,

    Lambda_X(sigma f)=Lambda_X(f),  f in A.                 (V1)

For every extension L/K, the functional on A_L is L tensor Lambda_X,
and (V1) holds for EVERY L-algebra automorphism of A_L. There is no
hypothesis on the units of A, the components of H, or the determinant
character of H.

The normalized trace is also independent of the chosen reductive
homogeneous presentation and is preserved by every isomorphism between
two such varieties (Section24). For ANY finite etale morphism q:Z->X
between two such affine homogeneous varieties, with no equivariance
assumption, Section25 proves the stronger marked identity

    Lambda_Z(h q^*a)=Lambda_X((Tr_q(h)/deg q) a).

Thus phase pullback preserves and reflects purity, safety and unsafety,
with one fixed marker transported by normalized algebra trace. Section26
extends these conclusions to geometrically homogeneous K-forms without
assuming a K-point.

Consequently pure functions, safe functions, unsafe functions, and the
complete set of surviving fixed-marker orders are transported in both
directions by arbitrary automorphisms. These are trace statements; no degree filtration, primitive lift,
power descent or unsafe starting phase is produced. An arbitrary phase
is not made invariant under a finite quotient. If it already descends
through a finite etale map between the stated varieties, Section25
pushes any fixed marker down by normalized algebra trace.

## Volume lemma and its applications

The setup and claims until Section13 have the explicitly stated
volume hypothesis (H1). They establish the complete lemma used in
the unrestricted argument.

Let K be ANY field of characteristic zero. Let G be a connected reductive
linear algebraic K-group, and H a closed K-subgroup. Put X=G/H, the
fppf homogeneous quotient, and suppose X is affine. In characteristic
zero all these group schemes are smooth; X is a smooth geometrically
integral variety with the K-point o=eH. No K-split hypothesis and no
transitivity assertion for the set G(K) on X(K) is imposed.

Write g=Lie(G), h=Lie(H), n=dimX and A=O(X). The determinant isotropy
condition in this note is

    delta_H = det(Ad(H) on g/h) = 1                       (H1)

as an algebraic character of H, including its components. Equivalently,
the one-dimensional H-representation Lambda^n(g/h)^* is trivial.

**Volume lemma.**
Under these hypotheses and(H1), there is a unique G-invariant K-linear
functional Lambda:A->K with Lambda(1)=1. EVERY K-algebra automorphism
sigma of A satisfies

    Lambda(sigma f)=Lambda(f),  f in A.                 (H2)

For EVERY coefficient-field extension L/K the same theorem holds on
A_L=L tensor_K A, with G_L,H_L and Lambda_L=1 tensor Lambda. No assumption
that the units of A, or its geometric scalar extension, are constant
is needed.

In particular the theorem applies to:

- every connected reductive G and connected reductive H;
- every connected reductive G and finite closed H;
- G itself, including semisimple groups and algebraic tori;
- the rank-one trace-one homogeneous spaces of the retained all-pair
  proof and their normalized trace;
- SL_2/Gamma for any finite subgroup scheme Gamma, over its field of
  definition, including binary octahedral and binary icosahedral groups.

The volume lemma gives trace invariance. It does not produce a starting
pure unsafe function, identify a degree filtration, or assert that an
arbitrary automorphism lifts to G. For finite quotients, the phase must
actually descend; once it does through a finite etale map between the
stated varieties, Section25 supplies fixed-marker descent by T_q.

The proof uses classical homogeneous-space, complete-reducibility and
Rosenlicht-unit results at explicit scopes, together with the exact
top-form argument below. It is not a claim of a new general method or
of historical priority.

## 1. The abstract exact-form criterion

Here is a useful preliminary statement that separates the mechanism
from homogeneous-space structure.

Let Y be a smooth geometrically integral affine K-variety of positive
dimension n with a nowhere-zero top form omega. Suppose a K-linear
functional l:O(Y)->K has l(1)=1 and

    d Omega^(n-1)(Y) = (ker l) omega.                   (H3)

Assume moreover that

    u in O(Y)^*, l(u)!=0  ==>  u is constant.           (H4)

Then every K-algebra automorphism sigma preserves l.

Indeed Omega^n(Y)=O(Y)omega, so(H3) identifies the top algebraic de Rham
cohomology with K[omega], and [omega] is NONZERO. Because Y is affine,
global algebraic forms compute this cohomology: in top degree it is
simply Omega^n(Y)/dOmega^(n-1)(Y). Pullback by an automorphism is an
isomorphism on this quotient.

Write sigma^*omega=u omega. The multiplier u is an invertible regular
function since both forms generate the same invertible module. Hence

    [sigma^*omega]=l(u)[omega] !=0.

Condition(H4) forces u=c in K*. Exact forms are preserved by pullback,
and therefore sigma(ker l)=ker l. For f-l(f) in ker l, the fact that
sigma fixes K and1 gives l(sigma f)=l(f).

This argument does NOT assume automorphisms preserve omega pointwise.
The nonzero constant c can depend on sigma. It also does not assume
that every unit is constant. The one-dimensional nonzero top class is
the additional input that will eliminate character-valued multipliers.

Dimension zero is separate and tautological: a connected geometric
homogeneous quotient of dimension zero is a point, A=K and l=id.

## 2. Quotient and field issues

The quotient G/H exists as an algebraic K-scheme, the quotient map
pi:G->X is faithfully flat, and G is an H-torsor over X. These are
the usual scheme homogeneous quotients, not the possibly smaller set
of quotients of K-rational points. Since H is smooth in characteristic
zero, pi is smooth. Thus X is smooth, and the derivative of pi at e
identifies T_oX with g/h.

After any field extension the quotient is G_L/H_L. Over an algebraic
closure G is connected and irreducible, and its surjective image X is
irreducible; smoothness makes X geometrically integral. A regular
G-invariant function is constant on the geometric transitive orbit,
so

    A^G=K,  (A_L)^(G_L)=L.                              (H5)

If H is reductive, affineness follows from Matsushima's criterion
over an algebraic closure and faithfully flat descent. Equivalently,
one can use X=Spec O(G)^H: flat scalar extension commutes with this
invariant subspace, since it is the kernel of the coaction minus the
trivial coaction; the geometric quotient is the affine quotient.
Finite generation and the isomorphism descend. This supplies the
asserted affine quotients for connected reductive or finite H.
The volume lemma permits taking affineness as an explicit
hypothesis, without importing the converse part of Matsushima.
The full theorem uses the converse in Section15 to prove that H^0
is reductive.

All uses of smoothness, transitivity, determinants and derivatives
in the proof are geometric or scheme-theoretic. No density of G(K)
and no chosen embedding K->C is required.

## 3. Existence of a global invariant volume form

The tangent bundle of X is the associated bundle G times^H(g/h).
Its top cotangent bundle is therefore associated to the H-line
Lambda^n(g/h)^*. Under(H1), choose any nonzero vector v in this
one-dimensional K-space. Its H-action is trivial, so its translates
by G descend along the H-torsor to a global top form omega on X.
It is G-invariant and nowhere zero:

    Omega_A^n = A omega,  L_(D_Z)omega=0, Z in g.       (H6)

Conversely a G-invariant nowhere-zero form has an H-fixed nonzero
value at o; thus(H1) is EXACTLY the invariant-volume criterion.
This proves the criterion directly by the associated bundle, rather
than appealing to a topological orientation analogy.

If G and H are both connected reductive, (H1) is automatic.
The adjoint determinant of any connected reductive group is trivial:
its connected central torus acts trivially in its adjoint action,
and the determinant character is trivial on its semisimple derived
group. This can also be checked after algebraic closure from the
opposite-root pairs and the trivial torus weights. Consequently

    det(Ad(H)|g/h)
      = det(Ad(H)|g) / det(Ad(H)|h) = 1.

The numerator is the restriction of the adjoint determinant of G;
the denominator is the adjoint determinant of H.

For finite H, h=0, and the numerator is still1 even on disconnected
H. Thus every finite H satisfies(H1). Constant finite groups need
not have all their matrices defined over the prime field; the
statement is over the coefficient field on which the subgroup scheme
is defined.

## 4. Exact top forms are Lie derivatives of functions

Let D_Z be the fundamental derivation associated to Z in g; fix one
consistent sign convention for the action on functions. The fundamental
vector fields span Der_K(A) AS AN A-MODULE. In fact, the bundle map

    O_X tensor_K g -> T_X

is surjective at every geometric point by transitivity and smoothness
of the orbit map. Its coherent cokernel is zero; affineness gives the
claimed surjection on global modules.

Contraction with omega is an A-linear isomorphism

    Der_K(A) -> Omega_A^(n-1),  D |-> i_D omega.

Choose a K-basis Z_i of g. Every (n-1)-form is a finite sum
alpha=sum_i f_i i_(D_Zi)omega. Cartan's formula and(H6) give

    d(i_(D_Zi)omega)=0,
    d(f_i i_(D_Zi)omega)=D_Zi(f_i)omega.                (H7)

The second formula is the exact top-degree identity
df wedge i_Domega=df(D)omega. Thus

    dOmega_A^(n-1) = (gA)omega,
    gA := span_K{D_Z(f): Z in g, f in A}.              (H8)

Here gA is a VECTOR SUBSPACE, not an ideal. This is essential:
the quotient in(H8) is a coinvariant space, not a quotient ring.

## 5. The normalized trace and the nonzero top class

The action of G on A is rational and locally finite. In characteristic
zero reductivity gives complete reducibility of every finite-dimensional
rational submodule over K itself. Thus A is the sum of its simple
rational submodules.

For a nontrivial simple submodule V, gV is a G-stable subspace:
conjugation transports D_Z to D_(Ad(g)Z). If gV were zero, the
differential representation would be zero. Connectedness and
characteristic zero then force the group representation to be trivial.
For example, its image is a connected algebraic group of zero Lie
dimension, hence is trivial. Therefore gV=V by simplicity.

On a trivial submodule every D_Z acts as zero. Using(H5), the entire
trivial isotypic subspace of A is K1. It follows that

    A=K1 direct_sum gA.                                (H9)

Let Lambda be the projection onto the first summand. It is G-invariant,
Lambda(1)=1, and every G-invariant normalized functional must agree on
both summands. In particular it is UNIQUE, and

    ker Lambda=gA,
    Omega_A^n/dOmega_A^(n-1)=K[omega], [omega]!=0.        (H10)

The complete-reducibility argument can be applied inside one finite
G-submodule at a time; it does not presume finite-dimensionality of A.
It also shows directly that the construction is compatible with
arbitrary field extension: derivatives and the split quotient in(H9)
tensor to the corresponding split quotient over the larger field.

## 6. Nonconstant units have zero trace

Let u be a unit of A. Its pullback to G is a unit of O(G), and it
does not vanish at e. Normalize by c=u(o) in K*. Rosenlicht's
unit theorem, in the form of Conrad Corollary1.2, gives

    pi^*u = c chi,  chi:G->G_m a K-group character.     (H11)

Right H-invariance of the pullback implies chi|H=1. Conversely a
character trivial on H descends to a unit on X. In particular any
unit u is a G-semiinvariant. If u is nonconstant, chi is nontrivial,
and G-invariance of Lambda forces

    Lambda(u)=0.                                      (H12)

For precision over a nonalgebraically-closed field, the character in
(H11) is defined over K: the normalized unit is a K-morphism with
value1 at the identity, exactly the hypothesis of Conrad's result.
One can also differentiate. A nontrivial character of a connected
group in characteristic zero has nonzero differential. For a suitable
Z in g, D_Z(u)=a u with a in K*; then Lambda(u)=a^-1 Lambda(D_Zu)=0.

The character can be nontrivial even when G has a central torus.
Hence this step handles rather than excludes torus factors. In the
semisimple case no nontrivial character exists, so all units are
constant, recovering the simpler sufficient hypothesis of the
retained rank-one trace-one proof.

## 7. All automorphisms and all coefficient-field extensions

Combine(H8)--(H12) with Section1. For any automorphism sigma,
sigma^*omega=u omega and the top cohomology isomorphism forces
Lambda(u)!=0. But(H12) makes every nonconstant unit have zero trace.
Therefore u=c is a NONZERO CONSTANT. Exact forms are preserved and

    sigma(ker Lambda)=ker Lambda,
    Lambda(sigma f)=Lambda(f) for every f.             (H13)

No restriction to equivariant, tame, linear, group or lifted
automorphisms has been made.

Every hypothesis persists after an arbitrary extension L/K:
G_L remains connected reductive, H_L remains smooth, X_L is the
smooth affine quotient, and the determinant character in(H1) remains
trivial. Applying the proof over L gives exactly Lambda_L from(H9).
Therefore(H13) holds for EVERY L-algebra automorphism of A_L,
including automorphisms not descending to K.

One could alternatively extend to an algebraic closure, prove the
identities there and descend their coefficients. The explicit proof
above shows that no split-form or constant-field extension is hidden
in the construction of omega or Lambda.

## 8. Moment and fixed-marker transport

Call f pure if Lambda(f^m)=0 for every m>=1. Call f safe if for each
fixed h in A all sufficiently large m satisfy Lambda(hf^m)=0;
an unsafe phase has one fixed h with infinitely many surviving orders.

For every automorphism sigma and fixed pair f,h,

    Lambda((sigma h)(sigma f)^m)=Lambda(hf^m), m>=0.     (H14)

Consequently purity, safety and unsafety are all preserved in BOTH
directions. The transported marker is one fixed polynomial; it is
not chosen anew for each exponent. Its precise set of surviving
orders is preserved.

These notions for a fixed K-phase also persist under coefficient-field
extension. A marker over L is a finite L-linear combination of markers
over K, so finitely many K-cutoffs suffice for safety over L. Conversely
nonzero scalar moments and a fixed unsafe K-marker remain nonzero over
L. No polynomial-power descent, radial degree or original Image-lift
construction is asserted for arbitrary homogeneous spaces.

For SL_2/Gamma with Gamma finite, A=O(SL_2)^Gamma for the right action,
and its normalized trace is the restriction of the left-invariant
normalized trace of O(SL_2). The theorem applies without requiring
an automorphism of A to lift through the finite covering. To obtain an actual failure on such a quotient, one needs an actual
Gamma-invariant phase and an upstairs fixed unsafe marker for that
phase. The finite etale map SL_2->SL_2/Gamma is in the stated class,
so Section25 pushes that marker down by normalized algebra trace,
with exactly the same surviving orders. A preexisting invariant marker
is therefore unnecessary once the phase descends. Averaging an
arbitrary phase is not claimed to preserve purity.

## 9. Determinant and nonreductive examples

The volume lemma permits disconnected H under its displayed determinant
condition. The full theorem has no such restriction; its additional
argument is in Sections14--21. Reductivity of H^0 alone does not imply
that the full determinant character is trivial.

For example take G=SL_2 and H=N(T), the normalizer of its diagonal
torus. In the quotient tangent space g/h, use the off-diagonal root
vectors E,F. The Weyl representative

    w = [[0,1],[-1,0]]

acts by E->-F,F->-E, with determinant-1. Thus delta_H is nontrivial
and X admits no G-invariant algebraic top form. This prevents applying the volume lemma directly to the base.
Sections15--21 instead use its finite determinant cover to prove
automorphism invariance of the normalized trace on this quotient.

Reductivity of the acting group is also a genuine boundary of the
argument. For Ga acting transitively on A1, there is no normalized
Ga-invariant functional: differentiation sends x to1. Equivalently
dOmega^0(A1)=Omega^1(A1), so its top exact-form quotient is zero.
The complete-reducibility/normalized-top-class mechanism has failed.
This is not an assertion that every nonreductive example fails.

Over a fixed field, nonconstant units do not obstruct the theorem:
they are explicitly eliminated as volume multipliers in Section7.
The torus example G=G_m,H=1 has Lambda equal to Laurent constant
term, omega=dt/t and dimension-one top quotient. Its algebraic
automorphisms t->a t or a/t multiply omega by+1 or-1 and preserve
Lambda, in agreement with the general theorem.

## 10. Primary inputs and exact literature boundary

The following actual sources were opened and read at their indicated
statements. Their roles are explicit; the normalized all-automorphism
trace deduction is written above rather than imported from a search
result.

1. Ivan V. Arzhantsev, Invariant ideals and Matsushima's criterion,
   arXiv:math/0506430v1 (2005), Section1 and Theorem2.1/Proposition2.2
   on physicalpp1--2:
   https://arxiv.org/pdf/math/0506430 .
   The paper proves the affine homogeneous quotient criterion in
   algebraically closed characteristic zero; it also displays the
   constant/nontrivial simple-submodule projection on K[G]. We use
   its affine direction and descend as in Section2. Section15 also
   uses its converse: affine G/H^0 implies reductive H^0.

2. J.S. Milne, Algebraic Groups (2017 author PDF),
   Theorem5.28/Proposition5.25/Corollaries5.26--5.27,
   printedpp105--106; Theorem22.42/Corollary22.43 printedpp477--478;
   Section12l, Remark12.53 and the following Reynolds paragraph,
   printedpp248--249:
   https://www.jmilne.org/math/Books/iAG2017.pdf .
   These supply quotient/torsor existence and complete reducibility
   over general fields; they are not an all-automorphism trace theorem.
   Standard associated-bundle descent, Kähler forms and flat scalar
   extension are used explicitly in Sections2--5.

3. Brian Conrad, Units on product varieties, author note,
   Theorem1.1 and Corollary1.2, physicalp1:
   https://math.stanford.edu/~conrad/papers/unitthm.pdf .
   Corollary1.2 proves that a normalized map from a smooth connected
   group to a torus is a group homomorphism. This is the precise
   Rosenlicht-unit input in(H11); the proof works over general fields.
   The source credits Rosenlicht rather than claiming a new unit theorem.

The retained RANK_ONE_TRACE_AUTOMORPHISM_INVARIANCE_ALL_IMAGE_PAIRS.md
already gives the exact-top-form/codimension-one argument for its
specific rank-one trace-one homogeneous spaces; the earlier Hopf
orbit note gives the two-dimensional special case. Those arguments
are antecedents. The present proof applies the same mechanism using
classical reductive homogeneous-space structure, and Section6 removes
the unnecessary constant-unit hypothesis by a nonzero top-class test.

Kaliman--Kutzschebauch, On algebraic volume density property,
arXiv:1201.4769v2, gives directly related classical antecedents.
Proposition8.4 (physicalpp24--25) characterizes invariant volume
forms by the modular-character identity; Example8.5 (physicalp25)
gives SL_2/N(T) as a quotient without such a form. Proposition8.6
(physicalp26) proves over C, assuming H^n(X,C)!=0, that a G-invariant
volume represents a generator, by compact-group Haar averaging.
The algebraic Cartan/coinvariant proof in Sections4--5 establishes
the needed nonzero top class over arbitrary characteristic-zero
fields without that assumption or analytic comparison.

Their Proposition3.6 and Remark3.7 (physicalp9) discuss constant
volume multipliers for a connected automorphism group and explain
that arbitrary automorphisms of general volume varieties need not
have constant multipliers. Our Section7 establishes constancy for
every automorphism of the volume-admitting homogeneous varieties
by testing the multiplier against the nonzero top class and the
Rosenlicht unit theorem. None of the density-property results is
a premise of the proof.

No historical-priority or external-verification claim is asserted.

## 11. The rank-one trace-one application

For n>=2 let G=GL_n and let H be its block-diagonal subgroup

    H={diag(a,B):a in GL_1, B in GL_(n-1)}.

This H is connected reductive. It is exactly the stabilizer, under
conjugation, of e=diag(1,0,...,0). The corresponding homogeneous space
is the affine variety

    X_n={M in Mat_n: rank M=1 and tr M=1}.               (H15)

Here one may define X_n by the two-by-two minors and trM=1. The minor
identities give M^2=(trM)M=M, so these are precisely the rank-one
idempotents. Conversely an idempotent of rank1 has trace1 in
characteristic zero.

Every such M gives V=imM direct_sum kerM. Choosing its image vector
and a basis of its kernel conjugates M to e over a geometric field.
Its stabilizer has no off-diagonal blocks, by the literal equation
ge=e g, so it is H. Scheme-theoretically the same orbit isomorphism
can be seen from the standard charts: write M=v ell, ell(v)=1, with
the redundancy (v,ell)->(c v,c^-1 ell). On M_ij!=0 use the jth column
for v and the ith row divided by M_ij for ell. The idempotent identity
gives ell(v)=1. These descriptions identify the orbit quotient with
the displayed smooth affine rank-one variety, rather than merely
identifying their geometric point sets.

Thus EVERY automorphism of O(X_n), over every characteristic-zero
extension, preserves its unique normalized GL_n-conjugation trace.
No automorphism lifting, chosen parameter chart, or balanced Image
kernel is needed to obtain this application.

The retained427 note
RANK_ONE_TRACE_AUTOMORPHISM_INVARIANCE_ALL_IMAGE_PAIRS.md studies this
same variety via M=v ell and the balanced quotient R=ell(v)=1.
It gives a separate checked monomial formula for its Lambda_n. The
compatibility statement is a direct CONDITIONAL consequence of
uniqueness: whenever that retained GL_n-invariance/normalization
identity is used, its functional equals the trace of this note.
That older identity is context for comparison, not a premise of
the theorem or the new rank-one application. No Image-kernel or
original proper-power assertion is imported.

## 12. A finite-averaging control

Automorphism trace invariance does NOT authorize averaging phases
before testing purity. On the torus A=K[t,t^-1], put f=t and
tau(t)=t^-1. Both Lambda(f^m) and Lambda(tau(f)^m) are zero for
every positive m. Nevertheless

    f_av=(f+tau(f))/2,
    Lambda(f_av^2)=1/2 !=0.                             (H16)

The invariant average is not pure. Therefore applying the theorem
on G/Gamma with Gamma finite still requires an ACTUAL invariant
phase; averaging a known phase is not a valid substitute. Once a
phase descends through the finite etale quotient, Section25 pushes
any fixed marker down by T_q with identical moments. This linear
marker operation is valid and is distinct from phase averaging.
The example tests phase averaging even for an automorphism of
order2 that individually preserves the trace in both directions.

## 13. Declared primary and bounded-control boundary

The retained primary source files for the three actual imports are

    sources/arzhantsev-matsushima-math0506430v1.pdf
    sources/milne-algebraic-groups-2017.pdf
    sources/conrad-units-on-product-varieties.pdf

Their exact roles are those of Section10. The full Milne book is
authenticated as a source artifact; ONLY the indicated quotient and
representation statements are used and were inspected. Neither a
whole-book audit nor an independent theorem establishing all-Aut
trace invariance is asserted.

The source-guarded checker replays the SL_2 adjoint determinant,
its local invariant-volume divergences and Cartan identity, finite
symmetric-power Lie coinvariants, rank-one matrix identities and
stabilizer controls, Laurent-torus trace/marker transport and the
inversion sign, the Ga exact-top-class obstruction, the N(T)
determinant obstruction, and(H16). Its finite models do not prove
complete reducibility, quotient descent, Rosenlicht's theorem or
the all-field conclusion. Those are supplied by the written proof
and the narrow primary inputs. The extension in Sections14--23 adds
the specific group-covering proposition cited there.

A replay is

    python3 verify_reductive_homogeneous_trace.py --repo-root /path/to/repository --output /tmp/reductive-homogeneous-full-trace.json

## 14. Inputs and the trace before choosing a volume

All homogeneous quotients and finite covers below are finite type.
In characteristic zero the groups and their closed subgroup schemes
are smooth. The scheme X is smooth and geometrically integral and has
the K-point o=eH.

The rational G-action on A is locally finite. Complete reducibility
over K gives, exactly as in the volume-admitting trace proof,

    A=K1 direct_sum gA,
    gA=span_K{D_Z f: Z in g, f in A}.                     (V2)

Indeed A^G=K by geometric transitivity; every nontrivial simple
G-submodule V has gV=V because gV is G-stable and a representation
of a connected characteristic-zero group with zero differential is
trivial. The entire trivial isotypic part is K1. Projection in (V2)
defines the unique normalized G-invariant Lambda_X, whether or not an
invariant volume exists. The construction commutes with arbitrary
extension of the coefficient field: (gA)_L=g_L A_L.

We use the already proved volume-admitting statement at precisely this
scope:

**Volume lemma.** If an affine homogeneous quotient of a connected
reductive group has trivial full determinant-isotropy character, it
has an invariant nowhere-zero top form omega,

    d Omega^(n-1)=(ker Lambda)omega,
    H^n_dR=K[omega],  [omega]!=0,                         (V3)

and every algebra automorphism preserves Lambda.

For clarity, the last assertion uses the following mechanism. Fundamental
vector fields span the tangent module, so Cartan's formula makes exact
top coefficients equal gA. For an automorphism, omega pulls back to
u omega with u a unit, and the induced nonzero map on (V3) forces
Lambda(u)!=0. A unit pulled back to the acting group is a scalar times
a character by the Rosenlicht--Conrad unit theorem. Every nonconstant
such unit has trace zero, hence u is constant. Pullback preserves exact
forms and fixes 1, which proves trace invariance. Thus the volume lemma
includes torus factors and does not assume constant units.

Sections1--8 above supply the full all-field proof of (V2)--(V3),
including the unit argument.

## 15. The determinant cover exists and is finite

We first work over an algebraically closed field k of characteristic
zero. All field descent is deferred to Section 21.

Let

    delta_H=det(Ad(H) on g/h).

The neutral component H^0 is reductive. One way to see the precise
affineness input is that G/H^0->G/H is finite etale: its fiber is the
finite component set H/H^0. Thus G/H^0 is affine, and the connected
case of Matsushima's criterion applies. The determinant of the adjoint
representation of a connected reductive group is trivial, so

    delta_H|H^0=1.

Hence delta_H has finite image. Define H^+=ker(delta_H), and let m be
the order of its image. Then

    p:Y=G/H^+ -> X=G/H                                  (V4)

is a connected finite etale Galois cover of degree m, with group H/H^+.
Over k this group is a finite subgroup of G_m, so it is cyclic.
The cover Y is affine, and delta_(H^+) is trivial since its Lie algebra
is h. The volume lemma therefore gives

    H^n_dR(Y)=k and all Aut_k(Y) preserve Lambda_Y.       (V5)

The case m=1 is allowed. Dimension zero is tautological: X is a point
and Lambda_X is the identity. Below one may assume n>0.

**Optional sharp determinant lemma.** In fact delta_H^2=1, so
m is either1 or2. This refinement is not needed below: finiteness
of the determinant image suffices for every argument.

To prove the refinement, work over the algebraically closed field.
The reductive Lie algebra h decomposes into its center and semisimple
derived ideal. Every h0 in H preserves both. On the semisimple ideal,
Ad(h0) preserves the nondegenerate Killing form, so its determinant
has square1. The connected center of H^0 is a torus, and conjugation
by h0 induces an automorphism of its integral character lattice;
the induced map on its Lie algebra has determinant+1 or-1.
Hence det(Ad(h0)|h)^2=1. Since det(Ad(h0)|g)=1 for the connected
reductive group G, the quotient determinant delta_H(h0) also has
square1. Smoothness and density of geometric points give this as
an identity of algebraic characters, not just pointwise values.
The only inputs are the central/derived decomposition of a connected
reductive group (Milne21.60--21.61, printedpp442--443), the character
lattice of a torus (Milne12.23, printedp240), and nondegeneracy of
the Killing form of a semisimple Lie algebra in characteristic zero
(Milne Section22c, printedp475).

## 16. The exact classical covering input

Only the following GROUP covering result is imported:

**Brion--Szamuely, Proposition 1.1(a).**
Over an algebraically closed field of characteristic p>=0, every
connected finite etale Galois cover of a connected algebraic group,
of degree prime to p, carries the structure of a central isogeny.
In characteristic zero this covers every finite degree. The identity
on the covering group can be chosen above the identity of the base.

The actual primary is arXiv:1109.2802v4, physical page 1,
Proposition 1.1(a), with proof in Section 2, physical page 3.
Its Theorem 1.2(a) assumes CONNECTED stabilizers. We do not apply that
theorem to the present possibly disconnected H. The exact extension to disconnected stabilizers is already
Borovoi--Demarche, Proposition5.2, physicalpp17--18, in the author
PDF dated17December2012. Section17 reproduces its characteristic-zero
specialization: the connected pullback component and its deck stabilizer,
Proposition1.1(a), and rigidity of etale lifts. The homogeneous-cover
lemma is a classical input, not a new assertion of this note.

Only Galois covers are needed in this note. No assertion about an
arbitrary non-Galois cover is required for the argument.

## 17. Galois covers of G/H are homogeneous after a central isogeny

This is the specialization of Borovoi--Demarche Proposition5.2 just
cited. The argument is reproduced to make the exact hypotheses and
its dependence on Brion--Szamuely Proposition1.1(a) explicit.

Let q:Z->X be a connected finite etale Galois cover and choose z0 over o.
Form G times_X Z and choose its connected component T containing (e,z0).
The projection T->G is a connected finite etale Galois cover. Explicitly,
the Galois group Gamma of Z/X acts on the pullback; the stabilizer Gamma'
of T acts simply transitively on every fiber of T/G. If two points of
one such fiber lie in T, the unique element of Gamma relating them
preserves T, since its translate meets that connected component.
This proves the stated Galois assertion, including the possibility
Gamma' is trivial.

By Proposition 1.1(a), identify T with a connected algebraic group G'
and write its projection as a finite central isogeny

    pi:G'->G,

choosing e'=(e,z0) as identity. Since G is affine and the cover is finite,
G' is affine. It is reductive: its connected unipotent radical maps to
a connected normal unipotent subgroup of G and hence is contained in
the finite kernel of pi, so is trivial.

Let F:G'->Z be the second projection and set

    P=pi^(-1)(H),   H'=F^(-1)(z0) subset P.             (V6)

The restriction F:P->q^(-1)(o), whose target is finite etale, shows
that H' is an open-and-closed finite-index union of components of P
containing e'. It is a closed subgroup, as follows.

For any h' in H', the morphisms F composed with right translation by
h' and F are lifts through the finite etale map q of the SAME map
G'->X. They agree at e'. The equalizer of two lifts through a separated
etale morphism is open and closed. Since G' is connected, the lifts
agree everywhere:

    F(g'h')=F(g'),  g' in G'.                            (V7)

This gives closure under products; evaluating at g'=(h')^(-1) gives
closure under inverses. The assertions are scheme-theoretic: P and its
open-and-closed union H' are smooth over the algebraically closed field,
and these equalities on geometric points identify the morphisms.
Equation (V7) therefore gives right-H' invariance and descent to

    f:G'/H' -> Z over X.                                (V8)

The quotient G'/H'->X is finite etale. This can be checked after the
smooth faithfully flat base change G'->X, where it is the product
cover with finite fiber P/H'. Consequently (V8), being a map between
finite etale X-schemes, is finite etale as well. Its image is nonempty
open and closed in connected Z, so it is surjective.

Its fiber at z0 is one point: the inverse image under F of z0 is
exactly H', whose quotient by H' is one coset. All schemes in that
fiber are etale and therefore reduced. A finite etale surjection
of connected varieties has constant degree, so f has degree 1 and
is an isomorphism. We have proved

    Z=G'/H',  H' finite index in pi^(-1)H.               (V9)

Differentiating pi identifies g' with g and h' with h. The isotropy
representations are compatible, hence

    delta_(H')=delta_H composed with pi|H'.              (V10)

This proves precisely the disconnected-stabilizer statement needed
here; it is not attributed to Theorem 1.2 of the primary.

## 18. Top de Rham cohomology detects the determinant character

More generally let Z=G'/H' be an affine homogeneous quotient as in
(V9). Put H'_0=(H')^0 and Z0=G'/H'_0. The map Z0->Z is a finite
etale Galois cover with component group Delta=H'/H'_0. Both spaces
are affine. The connected subgroup H'_0 is reductive (alternatively
it is the connected finite-isogeny preimage of H^0), and so Z0 has
a G'-invariant volume omega0. The volume lemma gives

    H^n_dR(Z0)=k[omega0],   [omega0]!=0.                 (V11)

For h' in H', right translation R_h' on G'/H'_0 commutes with left G'.
It sends omega0 to a constant multiple of omega0. Its multiplier is

    R_h'^* omega0=delta_(H')(h')^(-1) omega0.            (V12)

Indeed evaluate at eH'_0, translate the target point h'H'_0 back by
left multiplication by (h')^(-1), and obtain Ad((h')^(-1)) on
g'/h'. Taking its determinant gives (V12). An inverse convention for
the right deck action replaces the character by its inverse and has
the same triviality criterion.

Finite etale descent of Kahler differentials and characteristic-zero
averaging give

    H^n_dR(Z)=H^n_dR(Z0)^Delta.                          (V13)

For a direct algebraic verification, the finite etale torsor gives
Omega^r(Z0)=O(Z0) tensor_(O(Z)) Omega^r(Z), whose Delta-invariants
are exactly Omega^r(Z). The invariant functor is exact by averaging,
so taking invariants commutes with cohomology of the de Rham complex.
Affineness makes these global complexes compute algebraic de Rham
cohomology. Equations (V11)--(V13) prove the exact criterion

    H^n_dR(G'/H')!=0  IFF  delta_(H')=1.                 (V14)

When nonzero it is one-dimensional. The character is trivial on H'_0,
so testing the finite deck group loses no information.

## 19. The determinant cover is intrinsic among the covers needed here

Let q:Z->X be any connected finite etale Galois cover with nonzero
H^n_dR(Z). By Section 17 it is G'/H' with pi a central isogeny.
Equations (V10) and (V14) imply

    H' subset pi^(-1)(H^+).

The quotient map therefore factors through the determinant cover:

    G'/H' -> G/H^+ -> G/H.                              (V15)

All these maps are over X; in particular the first map is a map of
finite etale X-schemes and is finite etale. If q has degree m, equal
to that of (V4), the first map has degree 1 and is an isomorphism.

Thus Y is, up to isomorphism over X, the UNIQUE connected finite
etale GALOIS cover of degree m whose top algebraic de Rham cohomology
is nonzero. This uniqueness is only invoked in the displayed Galois
class. No universal covering space or classification of all finite
covers has been assumed.

## 20. Every automorphism lifts to the determinant cover

Let phi:X->X be any automorphism of varieties over k. Pull back (V4)
along phi and denote the resulting cover by q:Y_phi->X. It is connected,
finite etale, cyclic Galois of degree m. Its other projection
Y_phi->Y is an isomorphism of abstract varieties, since phi is an
isomorphism. Therefore H^n_dR(Y_phi) is nonzero by (V5).

Section 19 gives an isomorphism Y_phi->Y over X. Inverting that
isomorphism and composing with the other projection produces a lift

    phi_tilde:Y->Y,  p phi_tilde=phi p.                  (V16)

Uniqueness of this lift is NOT asserted; deck transformations account
for the usual ambiguity.

The pullback p^*:O(X)->O(Y) is G-equivariant. Consequently the
restriction Lambda_Y composed with p^* is the unique normalized
G-invariant functional on O(X), namely Lambda_X. By the volume lemma,
phi_tilde preserves Lambda_Y. Hence

    Lambda_X(phi^*f)
      =Lambda_Y(p^*phi^*f)
      =Lambda_Y(phi_tilde^*p^*f)
      =Lambda_Y(p^*f)
      =Lambda_X(f).                                    (V17)

This proves the theorem over every algebraically closed
characteristic-zero field. The proof uses no analytic realization
of G/H and no embedding into C.

## 21. Arbitrary fields, extensions and fixed markers

Return to the original K. The normalized trace is defined over K by
(V2). After extension to an algebraic closure it is exactly the
normalized geometric trace, since (gA)_L=g_L A_L and the constant
summand in (V2) tensors to L. Any K-automorphism extends to that
algebraic closure, where (V17) holds. Faithfulness of scalar extension
therefore descends (V1) to K.

The same argument applied to ANY extension L/K proves invariance
under every automorphism of A_L, including ones that do not descend
to K. We do not claim that a geometric lift (V16) necessarily descends
to K; its descent is unnecessary for the trace identity.

For fixed f,h in A, (V1) gives

    Lambda_X((sigma h)(sigma f)^r)=Lambda_X(hf^r),
    r>=0.                                               (V18)

Applying the inverse automorphism gives equivalence of purity,
safety and unsafety. A fixed h is transported to the single fixed
element sigma h, with exactly the same surviving orders. Extension
safety for a K-phase follows because every element of A_L is a
finite L-linear combination of elements of A; the finitely many
corresponding K-cutoffs have a common maximum.

## 22. Boundaries and examples

The result covers G/H with disconnected H and nontrivial determinant
character, including SL2/N(T). The Weyl component in that example
acts with determinant -1, so the base has zero top algebraic de Rham
cohomology by (V14). Its nontrivial double determinant cover has
nonzero top cohomology. The new argument lifts automorphisms to
that cover; it does not pretend the base has an invariant volume.

Torus factors cause no obstruction. They were already handled by the
unit-character/top-class step in the volume lemma; the central group
covers in Section 17 may include isogenies of tori. No constant-unit
assumption or unique unramified torus cover is used.

The simply connected semisimple case is included, rather than needed
as a fallback. The argument here avoids any separate use of the
etale simple connectedness of such groups.

The theorem makes no assertion for nonreductive acting groups. The
Ga action on A1 has no normalized invariant functional because
d/dx sends x to1, so reductivity cannot simply be dropped from this
proof. Nor does the theorem assert a starting unsafe phase or make an
arbitrary phase descend through a finite quotient. When a phase
already descends through a finite etale map between the stated
varieties, Section25 does give fixed-marker descent by T_q.

## 23. Exact primary source boundary

1. Michel Brion and Tamas Szamuely, Prime-to-p etale covers of algebraic
   groups and homogeneous spaces, Bull. London Math. Soc.45(2013),
   602--612, DOI10.1112/blms/bds110; arXiv:1109.2802v4.
   https://arxiv.org/pdf/1109.2802v4 .
   ACTUALLY USED: Proposition1.1(a), physical page1, and its Section2
   proof on physical page3, only over an algebraically closed
   characteristic-zero field. Theorem1.2(a) requires connected
   stabilizers and is NOT imported for disconnected H. The homogeneous
   cover in Section17 is proved below from the group proposition,
   following the credited argument of Borovoi--Demarche.

2. Sections1--8 supply the exact-top-form, complete-reducibility and
   Rosenlicht-unit argument at (V2)--(V3), using the Milne, Arzhantsev
   and Conrad source scopes stated in Section10.
   Matsushima is additionally used at its converse direction for the
   connected subgroup H^0: G/H^0 affine implies H^0 reductive.
   Arzhantsev Theorem2.1 over an algebraically closed characteristic-zero
   field gives precisely that assertion.

Finite etale descent, quotient existence, and averaging are used at
their explicit algebraic scopes above. No connected-stabilizer covering theorem is applied to disconnected H.
No topological orientation analogy, automorphism classification,
historical-priority claim or external-verification assertion is made.


3. Mikhail Borovoi and Cyril Demarche, Le groupe fondamental d'un
   espace homogene d'un groupe algebrique lineaire, author PDF
   dated17December2012, Proposition5.2 and its complete proof,
   physicalpp17--18:
   https://webusers.imj-prg.fr/~cyril.demarche/articles/GroupeFondamental-171212.pdf .
   This is the exact disconnected-stabilizer Galois-cover theorem.
   Section17 reproduces its specialized proof and is explicitly
   attributed to it; only characteristic zero is used.

4. Shulim Kaliman and Frank Kutzschebauch, On algebraic volume density
   property, arXiv:1201.4769v2:
   https://arxiv.org/pdf/1201.4769v2 .
   Precisely inspected: Proposition3.6/Remark3.7 physicalp9,
   Proposition8.4/Example8.5 physicalpp24--25, and Proposition8.6
   physicalp26. Their roles are the antecedent comparisons in
   Section10. These are contextual sources, not additional logical
   premises of the self-contained algebraic argument.

## 24. Presentation independence and every variety isomorphism

**Corollary IntrinsicNormalizedTrace.** If one affine K-variety has
two connected reductive homogeneous presentations, their normalized
traces coincide. More generally, for any K-isomorphism i:X->Xprime
between two affine reductive homogeneous K-varieties,

    Lambda_X(i^*a)=Lambda_Xprime(a), a in O(Xprime).       (I1)

Every assertion holds over every coefficient-field extension.

First work over an algebraically closed characteristic-zero field k.
For two presentations of V, the full theorem makes the normalized
trace Lambda_1 invariant under EVERY automorphism of V, in particular
under translations by all k-points of the second acting group Q.
This implies rational Q-invariance: for each f, a finite-dimensional
rational Q-submodule contains f, and
g |-> Lambda_1(g f)-Lambda_1(f) is a regular function on Q. It vanishes
at every k-point, so is zero because Q is smooth, reduced and has
dense k-points. Thus Lambda_1 is Q-invariant as a coaction identity.
Uniqueness of the normalized Q-trace yields Lambda_1=Lambda_2.

For i:X->Xprime, transport the action on Xprime through i to a second
homogeneous presentation of X. Its normalized trace is
Lambda_Xprime composed (i^*)^-1. The preceding equality proves (I1).
Applying i^-1 proves the separately valid reverse transport, so all
unmarked and fixed-marker moments correspond in both directions.
Over arbitrary K, extend both traces and the isomorphism to an algebraic
closure, use the result there, and descend by faithfulness of scalar
extension. Repeating with an arbitrary extension field proves the
full extension quantifiers. No assertion for a general morphism is
made by this corollary.

## 25. Finite etale morphisms and exact marked trace

**Theorem FiniteEtaleNormalizedTrace.** Let X and Z be affine
homogeneous K-varieties of connected reductive linear K-groups, with
possibly unrelated scheme quotient presentations, and let q:Z->X be
ANY finite etale K-morphism. No equivariance is assumed. Put A=O(X),
B=O(Z), and d=deg q. These geometrically integral varieties are
connected, and the image of the nonempty finite etale map q is open
and closed in X, hence all of X. Thus d>=1 is a constant integer.
For the finite locally free rank-d A-algebra B define

    T_q(h)=(1/d) Tr_(B/A)(h), h in B,

where Tr_(B/A)(h) is the trace of multiplication by h on B. Then

    Lambda_Z(h q^*a)=Lambda_X(T_q(h) a), h in B,a in A.    (E1)
    Lambda_Z(q^*a)=Lambda_X(a).                           (E2)

The statements hold after EVERY extension L/K, for arbitrary elements
of A_L and B_L and the base-changed morphism, not just elements
descending to K.

### 25.1. Compatibility on connected Galois covers

First work over an algebraically closed characteristic-zero field k.
Let V=G/H be affine, and let p:W->V be a connected finite etale
GALOIS cover. Section17, the credited Borovoi--Demarche Proposition5.2
specialization, identifies it with Gprime/Hprime->G/H for a finite
central isogeny pi:Gprime->G and finite-index
Hprime subset pi^-1(H). As in Section17, Gprime is connected affine
reductive, so W has a reductive homogeneous presentation.

The pullback p^* is Gprime-equivariant when Gprime acts on V through
pi. Hence Lambda_W composed p^* is normalized and Gprime-invariant.
Faithfully flat surjectivity of pi makes this the same as G-invariance,
so uniqueness gives

    Lambda_W(p^*a)=Lambda_V(a).                           (E4)

The trace on W is intrinsic by Section24, so this statement is
independent of a chosen homogeneous presentation of W.

### 25.2. A common connected Galois dominating cover

Choose a connected finite etale Galois cover p:W->X dominating q,
with an X-map r:W->Z and p=q r. The precise classical source is
Stacks Project Section58.7, Tag03SF, applying Lemma58.3.8, Tag0BN2:

    https://stacks.math.columbia.edu/tag/03SF
    https://stacks.math.columbia.edu/tag/0BN2

Here is its explicit finite-cover construction. In the d-fold fiber
power Z^d_X, the locus U of ordered d-tuples of pairwise distinct
points is open and closed, since diagonals of finite etale maps are
open and closed. It is a finite etale X-scheme. Choose a geometric
point of U and the connected component W containing it. Projection
to one factor gives r:W->Z. A nonempty connected component of a
finite etale cover of connected X is surjective onto X.

The symmetric group S_d permutes U. Let Gamma0 be the stabilizer
of W. In each geometric fiber of U, S_d acts simply transitively
on the ordered enumerations of that fiber of Z. If two such tuples
lie in W, the unique permutation relating them preserves W because
its translate meets that connected component. Therefore Gamma0
acts simply transitively on every W/X fiber. Its action is faithful;
an X-automorphism of connected W is determined by the image of one
geometric point, by etale lift rigidity. Thus p:W->X is Galois,
with deck group Gamma=Aut_X(W)=Gamma0. This is the finite etale
specialization of the cited Stacks proof, not an imported
classification of non-Galois homogeneous covers.

The X-map r between finite etale X-schemes is finite etale, as can
be checked etale locally after splitting both covers into finite
disjoint unions of the base. Its image is nonempty open and closed
in connected Z, hence it is surjective. It is also Galois. If w1,w2
in a geometric fiber of r are given, the unique gamma in Gamma
with gamma(w1)=w2 has r gamma and r agreeing at w1. These are lifts
of p through separated etale q, so lift rigidity on connected W
makes them equal everywhere. Consequently

    Delta={gamma in Gamma:r gamma=r}

acts freely and transitively on every fiber of r. This does not
assume that Delta is normal in Gamma.

Apply Section25.1 separately to p:W->X and r:W->Z, using the given
reductive homogeneous presentation of each base. This gives two
possibly different reductive homogeneous presentations of W;
Section24 identifies their normalized traces. Thus

    Lambda_W(p^*a)=Lambda_X(a),
    Lambda_W(r^*b)=Lambda_Z(b).                           (E5)

### 25.3. Averaging computes normalized algebra trace

Put N=|Gamma|=deg p and e=|Delta|=deg r, so N=e d. The full theorem
makes Lambda_W invariant under all deck automorphisms. Fiberwise
Galois counting gives the identity of regular functions

    (1/N) sum_(gamma in Gamma) gamma^*(r^*h)
        =p^*(Tr_(B/A)(h)/d)=p^*T_q(h).                  (E6)

Indeed etale locally q has d separate branches, and multiplication
by h is a diagonal matrix of their d values. Algebra trace sums
these values. Each appears e times among the N points of a p-fiber.
This proves (E6) on all geometric fibers, hence globally on reduced
W. Since p^*a is fixed by every gamma, invariance gives

    Lambda_Z(h q^*a)
      =Lambda_W(r^*h p^*a)
      =Lambda_W((1/N)sum_gamma gamma^*(r^*h p^*a))
      =Lambda_W(p^*(T_q(h) a))
      =Lambda_X(T_q(h) a).

This proves (E1). Tr_(B/A)(1)=d gives (E2). Averaging is applied to
the linear trace identity, not used to claim that phase averaging
preserves purity.

### 25.4. Arbitrary fields and fixed markers

For general K, extend to an algebraic closure. X and Z stay connected
and integral because they are geometric homogeneous varieties of
connected groups, and q stays finite etale of the same degree.
The normalized traces commute with extension. Algebra trace on a
finite locally free algebra also commutes with every scalar extension:
locally choose a finite free basis, and the multiplication matrix and
its trace extend entry by entry; these local identities glue.
Thus the geometric identity is precisely scalar extension of (E1)
over K and descends faithfully. The integer d is invertible because
K has characteristic zero. Applying the argument with any extension
field L as base proves the full extension assertion without requiring
a Galois dominator or lift to descend.

For each f in A and ONE fixed h in B, (E1) gives

    Lambda_Z(h(q^*f)^r)=Lambda_X(T_q(h)f^r), ALL r>=0.    (E3)

Purity is therefore equivalent by h=1. If f is safe, the single fixed
marker T_q(h) has eventually zero moments, so q^*f is safe. Conversely
if q^*f is safe, each downstairs marker a can be pulled back: scalar
multiplication on the rank-d A-module B gives T_q(q^*a)=a, so (E3)
proves safety of f. Unsafety is equivalent as well. An upstairs unsafe
marker h pushes down to T_q(h) with the identical surviving orders,
and a downstairs unsafe marker a pulls up to q^*a with identical
orders. A marker with infinitely many surviving orders cannot push
to zero. Pure unsafety is consequently equivalent when purity is
also imposed. No unsafe starting phase is produced.

## 26. Geometrically homogeneous forms without rational points

**Corollary GeometricallyHomogeneousForms.** The trace existence,
uniqueness, all-automorphism invariance, presentation independence,
isomorphism transport and finite etale marked trace theorem above
also hold for any smooth geometrically integral affine K-variety
V endowed with a K-action of a connected reductive linear K-group
that becomes transitive over an algebraic closure. V(K) need not
be nonempty and no quotient presentation G/H over K is required.
The same extension quantifiers and purity/safety/unsafety conclusions
hold. In the finite etale theorem both varieties must have such
reductive geometrically homogeneous presentations.

Indeed geometric transitivity gives O(V)^G=K. Complete reducibility
and connectedness in characteristic zero, exactly as in Section5,
give

    O(V)=K1 direct_sum Lie(G)O(V)

over K, defining the unique normalized G-invariant trace. Its
construction commutes with every field extension. Over an algebraic
closure choose a point of V; the resulting geometric presentation
is Gbar/H and the quotient theorem applies. Every K-automorphism
extends there, so its trace identity descends faithfully. Repeat
over any extension L/K to obtain invariance under all L-automorphisms.

For two such presentations or an isomorphism, perform Section24
over the algebraic closure and descend the trace equality. For a
finite etale q between two such forms, apply Section25 over that
algebraic closure; normalized traces and algebra trace base-change
as already proved, so (E1)--(E3) descend. This uses no rational point
of either original K-variety and no descent of a chosen geometric
stabilizer or Galois cover.
