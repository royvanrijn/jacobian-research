# Haar-cycle variation and fibre boundary linking

16 November 2026. Internal written proofs, with classical topological and
de Rham inputs identified below. Independent replay, formal verification
and external review are not claimed. Review remains deferred. Both COMPLETE
exceptional Mathieu--Zhao cases remain **OPEN**.

The full original-marker question has a single integral topological
obstruction. For every pure phase on either exceptional quotient, safety
is equivalent to the vanishing of the monodromy variation of the actual
compact Haar cycle around PHASE VALUE ZERO. That variation is the tube of
one uniquely determined integral two-cycle in a nearby fibre. No original
multiplier is discarded in this equivalence.

A linking identity constrains the obstruction. If nearby-fibre monodromy
is semisimple at eigenvalue 1, its relevant component must come from the
boundary at infinity of that fibre. Thus semisimplicity at 1 together with
no invariant boundary image gives a complete sufficient all-marker theorem.
Those hypotheses have not been proved for every exceptional pure phase.

This sharpens the geometric meaning of the
[zero-value logarithm criterion](ALGEBRAIC_RESOLVENT_SPLITTING_AND_ZERO_VALUE_LOGARITHMS.md)
and retains the [full source/core target](SEPARATED_HAAR_CORES_AND_GLOBAL_KUMMER_EXTENSIONS.md).
The [checker](verify_haar_cycle_linking.py),
[receipt](haar_cycle_linking_verification.json) and
[source audit](haar_cycle_linking_source_audit.json) delimit the evidence.

## 1. The integral tube sequence, with the distinguished lift retained

Let X=SL2/Gamma, Gamma=2O or 2I, and M=SU2/Gamma. Polar decomposition
retracts X onto the closed oriented rational homology three-sphere M.
In particular, integrally,

    H_2(X)=0, H_3(X)=Z[gamma], H_4(X)=0.                 (1)

Here H_2(M)=0 follows from integral Poincare duality and finite H_1(M).
Use the descended nowhere-zero volume Omega and normalized compact Haar.
Let P be nonconstant and regular. For a regular parameter s put

    F_s={P=s},       U_s=X minus F_s.                    (2)

Choose a finite exceptional set B such that the fibres and complements,
their homology, and the natural maps form locally constant systems over
S=C minus B. Include all bifurcations at infinity and add 0 to B if needed.
This is obtained by a stratification of a proper partial compactification;
it is not an assumption that P is proper or tame. The nearby F_s is a
smooth affine surface. Its H_2 is a finite free abelian group, since it
has the homotopy type of a finite CW-complex of dimension at most two.

Pick a large s_infty in S outside the compact set P(M). The actual compact
cycle M is then a cycle gamma_infty in U_(s_infty). Transport it along a
specified path to a basepoint in a punctured disc about 0. Call the result
gamma_s. It maps to the generator gamma in H_3(X); this fixes the lift
used below. We do NOT replace it by an arbitrary splitting of homology.

The Thom sequence of the pair (X,U_s), with the complex orientation on
the normal line of F_s, includes

    H_4(X) -> H_2(F_s) --Tube--> H_3(U_s)
                          -> H_3(X) -> H_1(F_s).       (3)

The map to H_3(X) is onto at s_infty because gamma_infty lifts gamma.
It is therefore onto everywhere on S. By (1), (3) becomes the exact
sequence of INTEGRAL local systems

    0 -> H_2(F_s;Z) --Tube--> H_3(U_s;Z) -> Z -> 0.     (4)

All groups in (4) are free; a nonzero tube cannot be hidden by torsion.
For clarity, (3) uses ordinary compact-chain homology: excision to a
tubular disc bundle gives H_k(X,U_s)=H_(k-2)(F_s). This homological Thom
identification follows fibrewise from (D2,S1) and glues by the complex
orientation. In fact the normal line here is trivialized by dP, so the
relative Kunneth calculation for F_s times (D2,S1) proves it directly.
It is not an identification of ordinary and Borel--Moore
homology of the noncompact F_s.

Let T_U and T_F be positive-loop monodromy about 0 on the groups in (4).
The quotient Z is constant. There is consequently one unique integral
cycle d_P in H_2(F_s;Z) such that

    (T_U-I)gamma_s=Tube(d_P).                            (5)

It depends on the specified continuation of the actual Haar lift. A
different arbitrary lift gamma_s+Tube(b) changes d_P by (T_F-I)b.
Thus an abstract split extension or a zero class in the coinvariants is
not, by itself, the vanishing of the ACTUAL cycle (5).

## 2. Exact equivalence with ALL original marked tails

For original Q in C[X], the actual continued Cauchy period is

    C_Q(s)=integral_(gamma_s) Q Omega/(s-P)
                    / integral_gamma Omega.            (6)

It is the germ continued from the physical compact resolvent at infinity.
It belongs to a regular-singular algebraic period system, including at 0,
as proved in the retained compact-resolvent theorem. No other solution of
an annihilating equation is substituted.

**Cycle criterion.** If P is pure, or eventually pure, then

    ALL original fixed markers are eventually zero
      <=> T_U gamma_s=gamma_s <=> d_P=0.                (7)

The equivalence holds for any one compatible continuation. In the safe
case it holds for every continuation; a failure has d_P!=0 for each.

First suppose all marked tails vanish. Each (6) is a finite Laurent
polynomial in 1/s, so has no monodromy. Differentiating with a flat cycle
gives all higher-pole periods:

    partial_s^(k-1) C_Q(s)
      =(-1)^(k-1)(k-1)! integral_(gamma_s)
                           Q Omega/(s-P)^k / norm.      (8)

Their monodromy is also zero. These forms span ALL of H^3_dR(U_s):
U_s is smooth affine, so top de Rham classes are represented by regular
top forms; Omega trivializes its canonical module, and its coordinate
ring is C[X][(s-P)^-1]. Terms with denominator zero are included because
Q Omega=s Q Omega/(s-P)-PQ Omega/(s-P). Therefore (T_U-I)gamma_s pairs
to zero with every de Rham class. The perfect de Rham--Betti pairing gives
(T_U-I)gamma_s=0 over C, hence integrally by freeness in (4).

Conversely, suppose gamma_s is fixed. Every (6) is single-valued at 0.
Regular singularity makes it meromorphic there: the local power/log
normal form has only integral powers and no logarithms in this case.
Take its finite principal part at zero. On the entire actual Cauchy
module, modulo polynomials in s, this principal-part map commutes with
s and partial_s. It retracts the point module

    D_0=span{s^-1,s^-2,...}.                             (9)

For a pure P, C_1=1/s. For an eventually pure P its nonzero polynomial
exponential unit period also generates (9). The intrinsic Haar-source
theorem gives the actual exponential image V a one-dimensional quotient
V/tV and makes V indecomposable. The principal-part retraction therefore
forces V=C[t]. This is exactly safety for ALL original markers. This
is the retained zero-log splitting argument, now triggered by the
distinguished geometric cycle, not by a hypothesis on a larger source.

In particular a nonconstant eventually pure phase has 0 as a bifurcation
value of the complement family. If the family and its original rational
forms extended locally trivially across 0, a transported compact cycle
would give holomorphic periods there, whereas C_1 is a nonzero strictly
proper finite Laurent polynomial. Monodromy may still fix gamma_s; a
bifurcation is not a failure certificate.

## 3. Only the actual eigenvalue-1 variation can cause failure

Over Q, let E_1 be the generalized eigenvalue-1 subspace of T_F and pi_1
its rational spectral projector. Put d_(P,1)=pi_1 d_P. The corresponding
projector on H_3(U_s;Q) maps compatibly through (4). Then, for an
eventually pure phase,

    safety <=> d_P=0 <=> d_(P,1)=0.                     (10)

Indeed d_(P,1)=0 means the generalized eigenvalue-1 part of gamma_s is
fixed by T_U. It gives no integer-exponent logarithms in any original
period (6). All remaining eigenvalues give noninteger exponent classes.
The zero-log principal-part retraction again proves safety. The reverse
implications follow from (7). Thus a pure phase cannot have only a
noninteger-class motion of the Haar lift without the corresponding
eigenvalue-1 logarithmic obstruction.

A valid direct middle-fibre sufficient condition follows: if 1 is absent
from the spectrum of T_F on H_2(F_s), the phase is safe. Here T_F is the
ACTUAL geometric monodromy at phase value 0, with the full bifurcation
set retained. This is stronger than semisimplicity at eigenvalue 1.
Semisimplicity of T_F alone does not control the extension in (4).
Generic-rank data, an ordinary direct-image generic connection, and
monodromy at exponential parameter t=0 cannot be substituted for T_F.

## 4. The canonical linking identity in these threefolds

Let B_F be the symmetric integral intersection pairing on H_2(F_s).
There is an integral linking functional attached to each c in H_3(U_s):

    ell_c(beta)=intersection_X(C_beta,c),
                      boundary C_beta=beta.            (11)

Such a compact three-chain C_beta exists by H_2(X)=0. Its ambiguity is
a compact three-cycle Z in X. The intersection pairing on H_3(X) has
odd middle dimension and rank one, so is zero. Thus Z has zero intersection
with c, and (11) is well-defined. Changing beta through a chain in F_s
also changes nothing, because c avoids F_s. These are integral statements;
no linking denominator is needed. Normalize the common orientation sign
so that the tube identity is

    ell_(Tube(a))(beta)=B_F(a,beta).                     (12)

One can check (12) on intersecting surface discs and their normal circle;
it is the usual local intersection definition of the tube linking number.
Changing all these sign conventions changes neither a vanishing statement
nor any sufficient condition below.

Naturality under transport gives
ell_(T_U c)(T_F beta)=ell_c(beta). Insert c=gamma_s, use (5) and (12):

    B_F(d_P,T_F beta)
       =ell_gamma(beta)-ell_gamma(T_F beta).            (13)

This is an exact constraint on the ACTUAL obstruction, not an arbitrary
extension of two vector spaces. In a column basis, with intersection
matrix B, monodromy T, cycle column d and linking column l, it is

    d^T B T = l^T(I-T),       T^T B T=B.                (14)

In particular d_P is orthogonal to every T_F-fixed vector.

## 5. Semisimple fibre monodromy localizes failure at infinity

Assume T_F is semisimple on E_1. Then E_1=ker(T_F-I). Monodromy
invariance of B_F makes E_1 orthogonal to all other generalized
eigenspaces. Equation (13) therefore gives

    d_(P,1) in Rad(B_F) intersect ker(T_F-I).             (15)

This locates the relevant cycle in a specific boundary group. Choose a
large smooth compact core K in the affine surface F_s, with boundary L.
The inclusion K->F_s is a homotopy equivalence. Poincare--Lefschetz duality
identifies the intersection map with H_2(K;Q)->H_2(K,L;Q). Its kernel,
by the exact sequence of the pair, is

    Rad(B_F)=image(H_2(L;Q)->H_2(F_s;Q)).                (16)

Thus (15) is an invariant cycle coming from the boundary at infinity of
the NEARBY FIBRE. This is a precise identification of the intersection
radical in this argument. It is NOT a general identification of the
earlier original-unit V-filtration defect with boundary vanishing cycles.

**Complete sufficient theorem.** An eventually pure phase on either
exceptional quotient is safe for every original fixed marker if

    T_F is semisimple on E_1, and
    Rad(B_F) intersect ker(T_F-I)=0.                     (17)

With semisimplicity still assumed, it suffices that the intersection pairing
restricted to E_1 be nondegenerate. The case E_1=0 is included. Proof: (15) kills d_(P,1),
and (10) proves all-marker safety.

Consequently every full exceptional failure needs either nontrivial
unipotent eigenvalue-1 monodromy on H_2(F_s), or, when that part is
semisimple, a NONZERO ACTUAL d_(P,1) in the invariant boundary image (16).
The existence of such a boundary space alone proves neither purity nor
failure. Conditions (17) are sufficient and have not been established
for every remaining pure phase. Do not redefine the full target by them.

An abstract split extension can still move a specified lift: if B is
nondegenerate and ell_gamma=B_F(b,-), (13) gives d_P=(T_F-I)b. Subtracting
Tube(b) makes a fixed lift, but it need not be the Haar lift. If T_F has
a nontrivial unipotent eigenvalue-1 block, this distinction is substantive.
The checker retains such a linear-algebra control without claiming its
realization by an exceptional pure phase.

## 6. Exact finite-cover transport and the remaining finish line

For pi:SL2->X and the lifted invariant phase, finite-cover transfer
commutes with (4) and with monodromy. Transfer takes [SU2/Gamma] to
[SU2]. Therefore the ACTUAL obstruction upstairs is

    d_(P,SL2)=transfer_F(d_(P,X)),
    pi_* d_(P,SL2)=|Gamma| d_(P,X).                     (18)

Since H_2(F_s) is free, either cycle vanishes exactly when the other does.
The upstairs cycle is already Gamma-invariant. Nontrivial finite-group
characters or averaging do not automatically kill it. This agrees with
the separated-source descent theorem and retains disconnected lifted
fibres if they occur.

The COMPLETE exceptional MZ problem is now equivalently: for EVERY actual
pure invariant phase, prove d_(P,1)=0 for the actual Haar lift, or construct
an actual phase with d_(P,1)!=0. All earlier strict-leading, mixed-weight,
image-interior and automorphic-orbit restrictions remain. The compulsory
non-tame original-unit defect applies to safe phases as well as failures;
the cycle (5) distinguishes them.

These homological arguments also apply to SL2 itself, which already has
known pure failures. Rational-homology-sphere topology alone therefore
cannot force the required vanishing. The remaining proof must use the
actual right-2O/2I invariant regularity or an equally restrictive geometric
property; no bare finite-character argument supplies it.

The new result is the complete topological equivalence and linking/boundary
constraint, not a full quotient decision or an implemented Betti-monodromy
producer. The checker contains selected exact matrix and period controls;
it does not verify the universal topology or compute the obstruction for
all phases. No new degree, support, density or prefix campaign is proposed.
