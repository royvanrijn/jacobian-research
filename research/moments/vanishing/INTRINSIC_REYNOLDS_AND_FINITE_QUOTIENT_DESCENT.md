# Intrinsic Reynolds periods and finite-quotient descent at every degree

15 November 2026. Internal written deductions; no novelty, independent
specialist review or formal verification is claimed. The
[source audit](intrinsic_haar_source_audit.json) records imported inputs.
The [small exact checker](verify_intrinsic_haar_controls.py) and
[receipt](intrinsic_haar_verification.json) check finite-group derived
subgroups and boundary examples. They do **not** machine-verify the
topological or curve-normalization proofs below.

The purpose is to replace further support-by-support enumeration with
general statements about the entire coordinate algebra. There are four
conclusions, with no degree restriction.

* Normalized Reynolds integration on an affine reductive homogeneous
  space is preserved by every algebraic automorphism, including nonlinear
  ones. It is intrinsic to the underlying variety, independent of its
  reductive homogeneous presentation.
* Every finite étale map between connected affine reductive homogeneous
  spaces intertwines their normalized functionals, even without an
  equivariance assumption. Normalized algebraic trace reduces every
  marked moment upstairs to one downstairs.
* On a normal unirational affine variety with only constant units, an
  invariant regular rational reparameterization R(f) forces f, after a
  constant shift, to transform by a character of the finite group. For a
  pure-null f the shift disappears whenever the character is nontrivial.
* In particular, no noninvariant regular function on SL_2 can be made
  right-2I-invariant by **any** nonconstant rational reparameterization
  that stays regular. For right2O, pure-null reparameterizations reduce
  to invariance already present in f or in f^2.

These close transformation and unary-descent questions. They do not
decide MZ on either exceptional quotient. The unresolved issue is the
existence of a bad pure-null function in the full invariant algebra,
or an all-degree positive safety criterion. A degree cutoff on a sparse
union is not such a criterion.

## 1. Compact cycles, units, and the orientation cover

Let G be a connected complex reductive group, H a closed reductive
subgroup, X=G/H, A=C[X], d=dim_C X. Choose compatible compact forms
L=K intersect H, with K connected. The compact homogeneous manifold
M=K/L has real dimension d, and X is a real vector bundle over M.
In particular M->X is a homotopy equivalence. We use the classical
compact-quotient lemma as stated in
[Doebeli, Section1.9, printed509](https://www.numdam.org/article/BSMF_1994__122_4_505_0.pdf).

The determinant character

    epsilon_H(h)=det(Ad(h) on g/h)

takes values in {+1,-1}. Indeed on the Zariski-dense compact group L
this representation is the complexification of its real orthogonal
action on k/l. Thus epsilon_H^2=1 on H. It is the orientation character
of M. If it is trivial, an invariant nonzero top wedge at eH descends
to a nowhere-vanishing algebraic d-form Omega_X. Its restriction to M
is a nonzero constant complex multiple of invariant real volume.

In the orientable case normalize Omega_X by integral_M Omega_X=1.
Then

    h_X(f)=integral_M f Omega_X,
    [f Omega_X]=h_X(f)[Omega_X] in H^d_dR(X;C).          (1)

All holomorphic top forms are closed, and H^d(X;C)=C by the compact
homotopy equivalence and orientation. Equation (1) follows by pairing
with [M]; it is not a claim that the complex manifold has real dimension
d. Its real dimension is2d, whereas M has real dimension d.

Every unit u in A is a constant times a character of G trivial on H.
Pull back u along G->G/H and apply Rosenlicht's unit theorem, in the
form of [Conrad, Corollary1.2](https://math.stanford.edu/~conrad/papers/unitthm.pdf).
Consequently

    u in A^*, h_X(u)!=0  =>  u is constant.              (2)

A nontrivial character has mean zero by left invariance. This step
does not assume that X has no nonconstant units; tori are included.

If epsilon_H is nontrivial, put H^+=ker(epsilon_H). The connected
finite étale double cover

    X^+=G/H^+ -> X                                      (3)

is the orientation double cover of M after the vector-bundle
retraction. This cover is intrinsic: among connected double covers
of a nonorientable closed connected d-manifold it is the unique one
with nonzero integral top homology. Orientability and this description
are invariant under homotopy equivalence. For example the compact core
of SL_2/N(T) is RP^2; H^2(X;C)=0 there. One must use its double cover
SL_2/T, not apply (1) directly to the nonorientable quotient.
The volume-monodromy warning also appears in
[Arzhantsev et al., Remark4.17](https://arxiv.org/pdf/1011.5375v2).

The normalized functional on (3) restricts to h_X: composition with
pullback is normalized and G-invariant. This elementary equivariant
case will be used to descend the general proof.

## 2. A weighted period identity and arbitrary finite étale maps

First let X,Y be orientable spaces as above, of equal complex
dimension d, and let F:X->Y be any regular map. Write

    F^*Omega_Y=J_F Omega_X,

where J_F is a regular function, possibly zero. Pulling back (1) for Y
and pairing with [M_X] gives, for every f in C[Y],

    h_X(J_F F^*f)=h_X(J_F) h_Y(f).                     (4)

Thus the all-degree weighted change-of-variables formula is immediate
from one-dimensional top cohomology. A nonzero **constant** J_F gives
the unweighted identity h_X(F^*f)=h_Y(f), without a finiteness assumption.

Now assume F is finite étale, of degree N. Its differential is an
isomorphism everywhere, so J_F is a unit. Pullback on cohomology is
injective: the topological covering transfer composed with pullback
is multiplication by N. Hence F^*[Omega_Y] is nonzero, so h_X(J_F)!=0
by (1). Equation (2) forces J_F to be constant. Therefore

    h_X(F^*f)=h_Y(f),  f in C[Y].                       (5)

No group homomorphism or equivariance was assumed.

For nonorientable Y, pull back its orientation cover Y^+ to X and
take the connected components X_i of X times_Y Y^+. The deformation
retraction of Y lifts through F, identifying a compact core of X with
the finite covering F^(-1)(M_Y) of M_Y. Its orientation character is
the pullback of that of M_Y. This character agrees with that of M_X:
orientation characters of closed manifolds are homotopy invariant,
as follows also from the characterization of the orientation double
cover by its nonzero integral top homology. Therefore the pulled-back
double cover is split when X is orientable, and otherwise is exactly
its orientation cover. Each connected X_i is consequently either X
itself or its intrinsic orientation double cover. In either case it
is orientable and is an affine reductive homogeneous space: the
identity case, or G_X/H_X^+ from (3). The topological identification
of these covers is algebraic: in their finite étale fiber product a
degree-one connected covering component is an algebraic isomorphism.
In particular no assertion that every finite cover has already been
given a homogeneous presentation is needed.

Apply (5) to X_i->Y^+, then restrict the functionals through (3).
This proves (5) in full generality. If Y is orientable, a connected
finite cover of its compact core is orientable; hence X is orientable
as well and the first proof already applies.

For an algebraic automorphism F of X, N=1. In the orientable case
J_F=+1 or-1 with our normalized volume, since F acts on H_d(X;Z)=Z
by one of these signs. The sign cancels in normalized integration.
For a possibly nonorientable X the cover argument gives the same
functional identity. An isomorphism between two reductive homogeneous
presentations is a special finite étale map; hence Reynolds integration
is intrinsic to the affine variety.

The étale or constant-Jacobian hypotheses cannot be dropped. On SL_2
the square map S(g)=g^2 has

    J_S(g)=det(1+Ad(g^(-1)))=2 tr(g)^2,
    h(tr(g))=0, h(tr(g^2))=-1.                          (6)

The determinant is computed first on a diagonal matrix, then extended
by conjugation invariance and density. It vanishes at tr(g)=0, so S
is ramified. Equation (4) still holds and h(J_S)=2, but unweighted
preservation fails. The checker replays25 weighted trace controls
using h(tr(g)^(2k))=Catalan(k), and the explicit failed unweighted one.
These controls illustrate the hypothesis; (4)--(5) have written proofs.

## 3. Algebraic trace detects every marker upstairs

For finite étale F:X->Y of degree N, let

    E_F=(1/N) Tr_(C[X]/C[Y]):C[X]->C[Y].

The C[Y]-module C[X] is finite locally free of rank N. Then

    h_X(Q)=h_Y(E_F Q), Q in C[X].                       (7)

For the orientable case, algebraic trace on functions also traces
forms by summing local inverse branches and commutes with exterior
differentiation. If F^*Omega_Y=c Omega_X with c!=0, then

    Tr_forms(Q Omega_X)=(1/c) Tr_alg(Q) Omega_Y.

Apply this to [Q Omega_X]=h_X(Q)[Omega_X]. Pairing with [M_Y]
gives h_Y(Tr_alg Q)=N h_X(Q), since Tr_alg1=N. This proves (7).
For a nonorientable target, use the orientation base change from
Section2. Algebraic trace commutes with that base change and is a
sum over its connected components. Their degrees over Y^+ sum to N,
and their normalized h restricts to h_X, proving (7) there as well.

Let P in C[Y], and use P also for F^*P. The trace projection formula
therefore gives the exact all-order identity

    h_X(Q P^m)=h_Y((E_F Q)P^m), every m>=0.              (8)

In particular:

* P is pure-null downstairs iff its pullback is pure-null upstairs.
* P is eventually marked-null for every marker downstairs iff its
  pullback is eventually marked-null for every marker upstairs.
* MZ on C[X] implies MZ on C[Y]; failure on Y pulls back to X.

The second assertion really includes **all** upstairs markers: each
has its fixed trace E_F Q. In the reverse direction use pullback of a
downstairs marker. No extra upstairs marker can reveal a failure of
a pulled-back P that was safe downstairs.

It does **not** follow that MZ on Y implies MZ on X. A pure-null
function upstairs need not be a pullback. Nor does a known arbitrary
failure on X automatically descend. Finite extension of the coordinate
algebras and a compatible trace do not prove that missing assertion.
For SL_2->SL_2/2O and SL_2->SL_2/2I, this is exactly the unresolved
direction needed to reuse an arbitrary SL_2 counterexample.

For a finite Galois quotient with deck group Gamma, (7) is ordinary
group averaging. The marker identity (8) also holds for a finite group
of automorphisms with fixed points: by Section2, h is Gamma-invariant,
and averaging is a C[X]^Gamma-linear projection. No étale assumption
is needed for this latter algebraic averaging statement.

## 4. Unary reparameterizations force characters

Here is a separate algebraic theorem that limits finite-quotient descent
without restricting the degree of f or of the reparameterization.

**Unary theorem.** Let A be a finitely generated normal C-domain with
unirational spectrum and A^*=C^*. Let a finite group Gamma act on A,
and let f in A be nonconstant. Suppose R in C(T) is nonconstant,
R(f) is regular, and R(f) is Gamma-invariant. Then R is a polynomial
and there is a constant c such that

    gamma(f-c)=psi(gamma)(f-c),                         (9)

where psi is a character of Gamma. If h is a Gamma-invariant normalized
linear functional and h(f)=0, then either f is already invariant, or
c=0 and f itself transforms by the nontrivial character psi.

First write R=A_0/B_0 with coprime polynomials. A Bezout identity for
A_0,B_0, together with A_0(f)=B_0(f)R(f), makes B_0(f) a unit of A.
Only constants are units, and f is transcendental over C. Hence B_0
is constant. This proves the assertion about rational functions.

Put u=R(f) and L_0=C(gamma f:gamma in Gamma). It is a finite field
extension of C(u): every gamma f satisfies R(T)-u=0. Let B be the
integral closure of C[u] in L_0. Normalization of an affine curve in a
finite field extension is finite, so B is a normal affine curve ring.
Every element of B lies in A: it is integral over C[u] subset A and
belongs to Frac(A), and A is normal. Also f belongs to B. Gamma
preserves B.

The curve Spec B is unirational, being dominated by Spec A. Its smooth
projective completion has genus zero: compose a dominant rational
parametrization of Spec A with the map to this curve and restrict to
a general line. This gives a nonconstant map P^1 to the completed
curve; Riemann--Hurwitz forces genus zero. Thus Spec B is P^1 minus
a nonempty finite set. It has only constant units because B subset A.
Two or more removed points would give a nonconstant invertible Mobius
function. Exactly one point is removed, and therefore

    B=C[t].                                            (10)

This argument is why normality, unirationality and the unit hypothesis
are retained explicitly.

Gamma acts on (10) by affine maps. A finite affine group has a common
fixed point, obtained by averaging an orbit, so after shifting t it
acts as gamma(t)=chi(gamma)t. Write f=F(t), deg F=e, and R of degree
D. If D=1, f is invariant. Otherwise let c be the unique constant
such that the coefficient of T^(D-1) in R(T+c) is zero, and set U=F-c.
The top coefficient in R(F(chi t))=R(F(t)) gives

    psi=chi^e, psi^D=1.

If V=U(chi t)-psi U(t) were a nonzero polynomial of degree v<e, the
difference of the leading Dth powers would have degree e(D-1)+v.
All lower terms of the depressed R have degree at most e(D-2), which
is smaller. Cancellation is impossible. Thus V=0, proving (9).
The constant term in (9), after applying h and using h(f)=0, gives
(1-psi)c=0. This finishes the proof.

The shift matters without mean zero. For example f=t^3+3t+5 is not
odd under t->-t, while (f-5)^2+1 is invariant. The checker retains this
boundary example; it is not a Haar counterexample.

The theorem applies to SL_2 and to normal unirational affine varieties
with constant units more generally. In particular, connected semisimple
homogeneous spaces have the required unit property by Section1 and are
unirational as images of the rational connected algebraic group G.

## 5. Exact failure descent in the character case

Let h be Gamma-invariant on a C-algebra A, B=A^Gamma, and let f be
pure-null and transform by a character of order d. Then P=f^d belongs
to B and is pure-null. Moreover

    f is marked-safe in A iff f^d is marked-safe in B.   (11)

One direction expands marked powers of f^d as powers of f. For the
other, for any fixed Q write every exponent as dm+r, 0<=r<d. Then

    h(Q f^(dm+r))=h(E_Gamma(Qf^r) P^m).                (12)

There are finitely many residue classes, and every E_Gamma(Qf^r) is
a fixed marker in B. Safety of P proves safety of f. Equivalently, a
bad fixed tail for f has an infinite nonzero subsequence in one residue
class, and (12) gives a fixed downstairs marker for a bad tail of P.
This is an all-order argument; no finite prefix is substituted for it.

The order d divides the exponent of Gamma_ab. Combining (9)--(12),
an invariant regular rational reparameterization of a pure-null f
cannot bypass the character case. If f is bad, its character power
already gives the relevant descended bad function. No assertion that
every polynomial R(f) separately preserves a bad tail is required.

For right2I on SL_2, the derived group is all of2I, so there are no
nontrivial characters. Thus for **every** nonconstant f and rational
R with R(f) regular,

    R(f) right-2I-invariant => f right-2I-invariant.      (13)

For right2O, its derived group is2T and its abelianization is C_2.
The exact finite-group checker constructs the derived subgroups as
normal closures of generator commutators and checks orders120 and24.
The latter contains the central -1 and maps to the derived A_4 in
the projective S_4 quotient, so it is precisely the binary2T subgroup.
For any pure-null f, an invariant regular R(f) consequently forces
either f invariant, or f anti-invariant through that C_2 character.
In the latter case f is right-2T-invariant and f^2 already descends to
SL_2/2O. An arbitrary rational function of the old noninvariant
normalizer witness cannot evade this condition. This extends the old
[particular-circuit power obstruction](NONSPHERICAL_PAIRS_AND_EXCEPTIONAL_SL2.md)
to a degree-free statement about **all** regular unary constructions.
Multivariable orbit averages and norms are not unary reparameterizations
and are outside this exclusion.

## 6. Intrinsic safe classes and the revised all-degree questions

Let Z_X be the pure-null set and S_X the subset safe for every fixed
marker. Every algebraic automorphism preserves both sets exactly:

    h(Q (phi^*P)^m)=h(((phi^*)^(-1)Q) P^m).             (14)

Thus applying any automorphism, even outside the known shear group,
to a safe function cannot manufacture MZ failure. A degree filtration
may change under phi; eventual vanishing in (14) does not depend on it.

If N_G is the original G-nullcone, the entire class

    C_X=union_(phi in Aut_alg X) phi^* N_G

is contained in S_X. More generally every P contracted by any algebraic
C^* action on X is safe: that action has finitely many positive weights
on P and finitely bounded negative weights on a fixed marker, while
h is invariant by Section2. This includes conjugates of the original
cocharacters. The known Haar-shear examples outside N_G lie in C_X,
so their safety is now part of a general intrinsic theorem.

The [all-pair shear proof](HAAR_SHEARS_AND_GENERAL_NULLCONE_CLASSIFICATION.md)
in fact constructs phi^*p outside N_G with p in N_G for every non-toral
connected reductive pair. Hence C_X strictly contains N_G in every
such case: the original G-nullcone is not preserved by all algebraic
automorphisms. In the toral case Z_X=N_G by the constant-term theorem,
and Z_X is preserved by all automorphisms by (14). Thus there is also
an all-pair structural criterion:

    N_G is Aut_alg(X)-invariant iff G_der is contained in H.

This explains why identifying the intrinsic pure-null set with that
particular group nullcone fails outside the toral case. It does not
identify the larger safe class with Z_X or decide MZ.

An all-degree normal-form theorem

    Z_X subset C_X, or even every pure-null P contracted
    by some algebraic C^* action on X,                  (15)

would therefore prove MZ on X. Neither inclusion is proved here, and
they are sufficient conditions, not claimed equivalent characterizations
of MZ. Long's bad models show that no such inclusion holds universally
for all non-toral homogeneous spaces.

For the two exceptional models the main questions should now be kept
at the following level:

1. Is there an all-degree normal form or direct criterion placing every
   pure-null function in S_X? The larger intrinsic safe class in (15)
   gives a concrete positive route, without assuming the answer must
   be negative.
2. Otherwise, can a multivariable invariant construction give a genuine
   pure-null element with an infinite marked tail? Its actual entire
   marked-period rank must exceed one, as in the
   [period theorem](ENTIRE_HAAR_PERIODS_AND_EXCEPTIONAL_DECISION.md).
3. Can the missing finite-cover failure-descent direction be proved
   under a justified additional hypothesis? Trace compatibility (8)
   alone does not settle it; on2I, unary transformations of a
   noninvariant upstairs witness are excluded by (13).

Further sparse enumeration is not the primary open work item. The
short five-weight experiment is retained as exploratory evidence in
the source audit, without promoting a new support-count theorem or
starting a larger enumeration. The already certified degree12 sparse
results remain available as constraints on candidates. Both exceptional
MZ decisions, and therefore the proposed universal negative MZ
classification, remain OPEN. The all-pair moment-nullcone classification
is a separate internally solved question and is not reopened here.
