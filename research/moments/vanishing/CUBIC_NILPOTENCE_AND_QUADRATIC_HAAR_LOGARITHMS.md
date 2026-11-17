# Cubic nilpotence and quadratic Haar logarithms

16 November 2026. Internal written proofs using the audited classical
limit-MHS inputs below. Independent replay, formal verification and external
review are not claimed. Both COMPLETE exceptional Mathieu--Zhao cases remain
**OPEN**; external review stays deferred.

For every regular phase on every finite quotient SL2/Gamma, every local
Jordan block of the FULL complement monodromy has size at most THREE.
Consequently every original Cauchy period has logarithmic degree at most
TWO. These assertions include nonsemisimple fibres, nonisolated critical
loci, bifurcations caused by infinity, and the puncture at infinity itself.
They do not require purity. All size-three blocks come from the original
compact-middle image W2 H2 of the fibre; boundary and ambient extensions
create no additional size-three blocks.

See the [checker](verify_haar_quadratic_logs.py),
[receipt](haar_quadratic_logs_verification.json) and
[source audit](haar_quadratic_logs_source_audit.json).

## 1. The full complement and its Hodge strip

Fix X=SL2/Gamma for ANY finite Gamma and a nonconstant regular f. Over the
regular locus of the FULL algebraic fibre/complement system put

    B_s=H3(X minus F_s;Q),   K_s=H2(F_s;Q)(-1).

Here H3 and H2 denote COHOMOLOGY. The existing tube/Gysin theorem gives
an exact geometric mixed-Hodge sequence

    0 -> L=Q(-2) -> B -> K -> 0.                       (1)

For the smooth affine surface F_s, H2(F_s) has original weights 2,3,4.
Its original Hodge indices p,q lie between 0 and2, with p+q>=2. Thus

    K: 1<=p,q<=3,  p+q>=4;
    L: (p,q)=(2,2).

Mixed-Hodge strictness in (1) implies, on the ENTIRE complement system,

    F^1 B_C=B_C,    F^4 B_C=0.                        (2)

The original W4 piece has Hodge indices1..3; GrW5 has indices2..3;
GrW6 is Tate of type(3,3). These statements are about original W.

At ANY puncture of the base P1, take a finite base change killing the
semisimple monodromy, and let N=log(T_unip). By geometric admissibility,
the Hodge filtration extends as subbundles of the canonical flat extension.
Consequently (2) holds for the limiting filtration F_infty with the SAME
ranks. The geometric limit (B,M,F_infty) is a mixed Hodge structure and
N is a morphism of type(-1,-1). Original W gives sub-MHS filtrations
of this limit. Exactness and strictness are the geometric facts audited
in the [unipotent extension note](UNIPOTENT_HAAR_EXTENSIONS_AND_GENERAL_RECOMPOSITION.md),
not a claim about arbitrary relatively filtered vector-space sequences.

## 2. A Hodge-strip nilpotence lemma

Let (E,M,F) be any mixed Hodge structure with F^a E=E and F^(b+1) E=0,
and let N:E -> E(-1) be a nilpotent mixed-Hodge morphism. Then

    N^(b-a+1)=0.                                      (3)

Indeed the functorial Deligne bigrading has I^(p,q)=0 unless
a<=p<=b, and N maps I^(p,q) into I^(p-1,q-1). After b-a+1 applications
the target has p<a and is zero. This argument uses the type of N on
the LIMIT mixed Hodge structure. Merely writing Griffiths transversality
on nearby Hodge filtrations would not prove (3).

Applying (3) to (2) gives

    N_B^3=0.                                         (4)

Likewise N_K^3=0, N_(GrW5 B)^2=0, and N_(GrW6 B)=0. Since the limit
is rational, conjugation also bounds q by1..3, so

    M_1 B=0,    M_6 B=B.                              (5)

Finite base change preserves every Jordan length: for a block with
eigenvalue lambda, the derivative of z^e at lambda is nonzero. Therefore
(4) bounds EVERY ORIGINAL eigenvalue block by size3, not only the
original eigenvalue1 block. It does not say (T-I)^3=0 on non-1 blocks.
For the original1 part it DOES give

    (T-I)^3=0,   N=(T-I)-(T-I)^2/2,   N^2=(T-I)^2.     (6)

This bound holds at every degree and puncture, including infinity.

## 3. Exactly where the size-three blocks live

We can sharpen the bound without assuming fibre semisimplicity. Use W
for ORIGINAL weights and M for the RELATIVE limit weights. The limiting
GrW5 B has Hodge strip2..3, hence M_3=0 on that piece. The limiting
GrW6 B is pure Tate of weight6. Strictness in the W-truncated geometric
limit sequences therefore gives

    M_2 B=M_2(W_4 B),    M_3 B subset W_4 B.            (7)

The pure original-weight4 piece W_4 B has its monodromy filtration
centred at4. Its defining isomorphism, together with (5), gives

    N^2:GrM6(W_4 B) -> GrM2(W_4 B) is an isomorphism,
    image N_B^2 = M_2 B = image N_(W_4 B)^2.            (8)

The inclusion image N_B^2 subset M2 follows from (5) and the weight
drop2 of N; the reverse inclusion follows from that centred isomorphism.
This proves equality in (8), not just an inequality of ranks.

The ambient line L is pure LIMIT weight4. It has zero intersection with
M2 B, so (1) induces an isomorphism

    image N_B^2  --~--> image N_(W_4 K)^2.              (9)

For surjectivity, lift from W4 K to W4 B using original-W strictness and
apply N^2. For injectivity use L intersect M2 B=0. Since
W4 K=(W2 H2(F))(-1), this proves

    rank N_B^2 = rank N_(W2 H2(F))^2.                  (10)

Both sides count size-three unipotent Jordan blocks. The weight2 piece
W2 H2(F) is the IMAGE of H2 of a smooth proper compactification, not
automatically its entire cohomology. Thus (10) concerns that precise
compact-middle image. Neither boundary-genus monodromy nor the ambient
extension can increase the size-three count. In particular if the right
side is zero, ALL complement logarithms have degree at most ONE, even
if the full fibre monodromy is nonsemisimple.

The image in (8) is pure limit Tate type(1,1); GrM6 B is pure limit
Tate type(3,3), and N^2 maps that top grade ONTO the bottom grade with
the appropriate twist. It need NOT be injective on all of GrM6 B:
genus/boundary contributions can lie in its kernel. The isomorphism
is on GrM6(W4 B), as specified in (8). This LIMIT Tate statement does NOT
make their nearby period functions algebraic or put every actual Haar
motion in the ORIGINAL top Tate boundary space.

## 4. Quadratic logarithms for every original marker

Let gamma be the physically continued Haar cycle and gamma_1 its
ORIGINAL eigenvalue1 projection. On the dual homological system let
N_h be the logarithm with the positive cycle-continuation convention.
It too has cube zero. In a canonical logarithmic frame, every original
simple-pole Cauchy period, restricted to its integer-exponent sector,
has the UNIQUE form

    C_Q,1(z)=a_Q(z)+ell*b_Q(z)+ell^2*c_Q(z),
    ell=log(z)/(2*pi*i),                               (11)

where a_Q,b_Q,c_Q are single-valued convergent MEROMORPHIC germs with
finite principal part. In a compatible canonical row frame h_Q(z),

    a_Q=h_Q gamma_1/Z,
    b_Q=h_Q N_h gamma_1/Z,
    c_Q=h_Q N_h^2 gamma_1/(2Z).                         (12)

Z is the fixed nonzero Haar normalization. For other eigenvalues one
has the corresponding fractional-power sector and still at most
quadratic logarithms; after finite base change every sector has (11).
The coefficients in (11) need NOT be algebraic. This is stronger in
Jordan scope, but different from the earlier semisimple-fibre algebraic
log-coefficient theorem.

For integral m the exact continuation law is

    T_h^m gamma_1=gamma_1+m*N_h gamma_1
                                  +(m^2/2)*N_h^2 gamma_1. (13)

Thus no third or higher log/finite-difference obstruction is missing
from the full local system. No expanding Jordan-size campaign is needed.

## 5. One finite ORIGINAL marker basis controls all markers and punctures

For each FIXED f there exist rho=1+b2(F_s) original regular functions
Q_1,...,Q_rho on X such that their simple-pole classes

    [Q_i Omega/(s-f)]                                  (14)

are a basis of the full generic H3_dR(U_s) over C(s). This existence
assertion is independent of semisimplicity, purity and the puncture.

Here is the single-pole spanning argument. On the generic smooth fibre,
the residue eta=Res(Omega/(s-f)) is a nowhere vanishing regular2-form;
adjunction trivializes its canonical bundle. Every regular2-form is
Q*eta for some Q in C(s)[F_s]. Since F_s is affine and smooth of dimension2,
such forms represent all H2_dR(F_s). The coordinate restriction from
C(s)[X] to the fibre is surjective, so these classes come from the
residues of original simple-pole3-forms. Moreover

    Omega=s*Omega/(s-f)-f*Omega/(s-f)

supplies the ambient line. The exact residue sequence (1) therefore
makes simple-pole original classes span the full complement cohomology.
Clearing denominators in s and extracting a basis selects rho fixed
original Q_i. These are not parameter-dependent amplitudes: each lies
in C[X]. The basis may have poles or determinant zeros at a puncture;
over its field of meromorphic germs it remains a basis, which is what
the following zero test needs. There is no asserted degree bound or
implemented general basis/canonical-frame producer.

It follows from (12) and nondegenerate de Rham pairing that

    all c_(Q_i)=0 iff N_h^2 gamma_1=0;
    all b_(Q_i)=0 iff N_h gamma_1=0.                    (15)

The first condition alone does NOT imply the second. The second already
forces all quadratic coefficients to vanish. Equivalently test absence
of BOTH log powers on these rho markers. Under EVENTUAL PURITY at0,
the existing actual-cycle theorem now gives the finite, all-degree test

    ALL original markers safe
      iff b_(Q_i)=c_(Q_i)=0 for every i.               (16)

Here zero means identity of a meromorphic germ, not a finite numerical
jet or checked moment prefix. The basis depends on f and its rank may
grow with degree. A statement only about the unit marker is insufficient.

## 6. What this does and does not close

The [general intrinsic Tate theorem](UNIPOTENT_HAAR_EXTENSIONS_AND_GENERAL_RECOMPOSITION.md)
still separates xi from the residual actual beta in F1/ker(A-I).
The present result bounds ALL their local Jordan possibilities. In
particular the ambient line cannot lie in image N_B^2, since that
image lies in M2 whereas L has weight4. The intrinsic obstruction can
therefore extend an ambient block to size2, but never to size3. Every
size3 block is already accounted for by (10). This does NOT make the
physical quadratic-log coefficient vanish: even a split extension may
have a specified lift with nonzero N_h^2 gamma_1.

For both full exceptional quotients the remaining task is still to
prove xi0 and actual beta0 for EVERY eventually pure phase, or realize
a genuine failure with original regular phase and fixed original marker.
Equivalently prove that the finite basis (14) has no logarithms under
purity. Neither the Hodge strip nor the LIMIT Tate image alone supplies
that missing implication. All primitive branches, all polynomial
recomposition degrees and every nonsemisimple case remain in scope.

Classical inputs: Deligne Hodge II, Theorem2.3.5 and Corollaries3.2.15,
3.2.17; Brosnan--El Zein, Theorem2.1, Remark2.2, Corollary2.22 and
Definition3.1/Lemma3.2; and the exact complement/source theorem already
audited in this repository. The new written deductions are (4), (8)--(10),
the original-marker spanning proof and their combination in (11)--(16).
Selected exact matrix controls verify normalizations and logical
boundaries, not the universal geometric existence and functoriality.
