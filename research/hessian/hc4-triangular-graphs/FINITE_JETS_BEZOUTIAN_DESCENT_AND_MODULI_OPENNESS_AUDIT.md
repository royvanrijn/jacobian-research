# Finite jets, Bézoutian descent and the missing openness argument

16 November 2026. SOURCE/ROUTE AUDIT ONLY. No new HC4 theorem is
registered, no published theorem is newly adopted into the HC4 ledger,
and no mathematical checker is executed. The elementary reductions below
identify the global assertions still required by three possible routes.
Full HC4/JC2 remains OPEN.

## 1. Inputs and precise primary-source reading

Retain KU0's finite-projective trace tensor, MA0's full residual-defect
base-change calculation, IV1's primitive-helper/normality boundary, and
MK0's distinction between formal energy identities and polynomial closure.
These registered results are credited; their checkers are NOT replayed.

Stephen McKean, [Bézoutians and injectivity of polynomial maps](https://shmckean.github.io/research/bezoutian-jacobian.pdf),
Journal of Pure and Applied Algebra227(2023),107298,
[DOI10.1016/j.jpaa.2022.107298](https://doi.org/10.1016/j.jpaa.2022.107298):
Definitions1.1--1.2, Sections2--3 and their complete proofs, Proposition4.3,
Lemma4.5 and their proofs, Corollaries1.5--1.6, and the relevant negative
example/Question4.12 on page8 were read. The general residue-duality
foundations cited there and the unrelated Drużkowski classification were
NOT independently reproduced. The required étale diagonal calculation is
derived directly in Section2 below.

Borisov--Gabber--Vasiu,
[On endomorphisms of affine spaces and the Jacobian problem, arXiv:2609.05746v1](https://arxiv.org/html/2609.05746v1),
Section26: Proposition26.1, Lemma26.3, Proposition26.7,
Corollary26.9, Lemma26.10 and their complete proofs, together with the
fixed-degree smoothness statement/proof26.17 were read. Section25's
line-bundle construction and its characteristic-zero limit in Remark25.5
were also screened: the positive construction uses characteristic p;
Remark25.5 explicitly leaves its characteristic-zero étale-plane question
open. IV1 already retained the paper's different weak-type-one result.
The whole175-page paper, the cited inverse-degree theorem, EGA and the
general moduli foundations were NOT independently reproduced or newly
adopted. No stronger low-sheet-degree assertion from its references is
silently imported.

For standard jet functor context, Proposition3.3 and its complete proof on
printed pages43--44 in the primary notes
[Ishii, An introduction to arc spaces](https://www.kurims.kyoto-u.ac.jp/~kyodo/kokyuroku/contents/pdf/1374-3.pdf)
were read, including its distinction from proper birational arc lifting
in Proposition3.2. These are different hypotheses.
The elementary Keller coordinate proof needed here is given in Section3;
no motivic integration or arc-space theorem is imported.

## 2. The Bézoutian is the diagonal projector, not automatic finite descent

Let K=(P,Q) be an ACTUAL complex plane Keller map, normalized to det DK=1.
Use the retained notation B=C[P,Q],R=C[x,y],F=Frac B,L=Frac R, full finite
normalization S and trace-dual Lambda. Put

    D_R=R tensor_B R
       =C[x,y,u,v]/(P(x,y)-P(u,v),Q(x,y)-Q(u,v)).

Use the divided-difference matrix M defined by

    M_i1=(K_i(x,y)-K_i(u,y))/(x-u),
    M_i2=(K_i(u,y)-K_i(u,v))/(y-v),
    beta=det M.

These are polynomials, and

    M*[x-u,y-v]^t=K(x,y)-K(u,v).                         (1)

On the diagonal, beta=det DK=1. At a distinct collision pair (1) has a
nonzero kernel vector, so beta=0. The scheme Spec D_R is étale over the
smooth plane by either projection, hence smooth and REDUCED. Its diagonal
is open (étaleness) and closed (separatedness). Therefore beta's class is
EXACTLY the diagonal idempotent e_R: it is1 on the diagonal factor and0
on the entire off-diagonal factor. In particular

    beta^2=beta in D_R.                                  (2)

No absence of off-diagonal points follows. Their complement is precisely
the support of the other projector1-beta. This recovers the usual
Bézoutian criterion in the present étale setting without importing a
general residue-duality theorem. For any target point b with r retained
preimages, the reduced fibre algebra is C^r. The corresponding Bézoutian
is the r-by-r identity matrix under the decomposition into ordered pairs
of retained points. If r>0 it is a scalar in that tensor algebra exactly
when r=1; an empty fibre has zero algebra. Thus constancy on EVERY fibre
is an injectivity condition, not a consequence of the first-jet Jacobian.
This criterion and its source attribution are McKean's prior work.

KU0 already constructs the generic separability tensor

    e_L=sum_i s_i tensor q_i in S tensor_B Lambda,
    Tr(s_i*q_j)=delta_ij                                 (3)

in a supplied free frame (or with projective dual frames). Its image in
D_R is e_R. These tensors agree generically and D_R is B-flat/torsion-free,
so the agreement is global; alternatively compare their clopen diagonal
values. The inclusion S tensor_B S -> D_R is injective: S and R are B-flat,
so tensor the two inclusions S->R successively. In a free frame, membership
of (3) in S tensor_B S says exactly that every q_i lies in S. Projective
dual frames give the same assertion locally. Hence

    e_L in S tensor_B S
       iff Lambda=S
       iff Tr(Lambda*Lambda) subset B.                  (4)

The last equivalence is the retained KU0 duality argument, since
S subset Lambda and Lambda^vee=S. Thus the Bézoutian version of the
original missing implication is FINITE-NORMALIZATION DESCENT of an already
regular polynomial diagonal projector. Polynomial regularity on the actual
double source D_R is automatic; regularity in S tensor_B S is not proved.
MA0 has already shown why splitting off the diagonal does not kill the
residual different. No new diagonal-descent theorem is supplied here.

Raw constancy of beta in C[x,y,u,v] is still a stronger sufficient test;
it must not be confused with constancy after reduction or with (4).
The primary source explicitly exhibits invertible maps with nonconstant
raw Bézoutian. No raw-constancy degree campaign is proposed.

## 3. Every source jet order retains the SAME full trace defect

For m>=1 write z(t)=z_0+z_1*t+...+z_m*t^m modulo t^(m+1), with each z_j
a pair of independent coordinates. Let K_m be the induced jet map, with
output coefficients y_j of K(z(t)). Taylor expansion has triangular form

    y_0=K(z_0),
    y_j=A(z_0)*z_j+H_j(z_0,...,z_(j-1)), j>=1,
    A=DK.                                                (5)

Because det A=1, A^-1 has POLYNOMIAL entries. Consequently

    Theta_m:(z_0,...,z_m) -> (z_0,y_1,...,y_m)             (6)

is a polynomial automorphism: recover z_j successively from (5), using
A^-1 and the previously recovered coefficients. In these source
coordinates

    K_m composed Theta_m^-1=K times id_(A^(2m)).          (7)

This is not an assumed inverse of K: z_0 is retained in (6). Formula(7)
is the explicit polynomial version of the standard étale jet square.

Let T be the2m independent positive-order TARGET jet variables. Then

    B_m=B[T], R_m=R[T], L_m=L(T),
    S_m=S[T], Lambda_m=Lambda[T].                         (8)

Indeed S[T] is normal, finite over B[T] and has fraction field L(T), so
it is the full normalization. Finite locally free trace duality commutes
with polynomial base change, giving the second equality. If s_i is a
free basis for S, its Gram matrix and inverse for S_m are the SAME matrices
viewed over B[T]. More intrinsically, Lambda_m/S_m=(Lambda/S)[T]. Polynomial
extension is faithfully flat, so

    Tr(Lambda_m*Lambda_m) subset B_m
        iff Tr(Lambda*Lambda) subset B.                  (9)

This covers ALL finite jet orders, not merely a bounded calculation.
Formal étale uniqueness likewise says that every target formal arc lifts
uniquely after a source starting point has been chosen. It does not
produce a missing starting point, eliminate an additional retained sheet,
or give valuative properness. An escaping branch in the finite normalization
need not be an arc based in the actual source at all. No inference from
formal arcs to automatic (4), motivic integration, or HC4 inversion is made.

## 4. Parameter jets are DIFFERENT; all of them already have inverses

Normalize an arbitrary plane Keller K by affine changes so that K(0)=0
and DK(0)=I. Define its scaling family over C[t] by

    K_t(z)=t^-1*K(t*z), K_0(z)=z.                        (10)

The division by t gives a polynomial family; det D_z K_t=1 identically.
For t!=0 it is a linear conjugate of K. This is the deformation used in
BGV Proposition26.7(1), in normalized coordinates. It connects EVERY
bounded-degree Keller map to identity, and does not itself prove inversion.

For each r>=1, the map K_t modulo t^r is already a polynomial automorphism
over C[t]/(t^r), irrespective of the unresolved generic case. More
generally, a polynomial endomorphism reducing to identity modulo a nilpotent
coefficient ideal has a polynomial inverse: correct the identity inverse
one power of the ideal at a time; a correction modulo I^(j+1) is obtained
by subtracting the error modulo I^(j+1), since the derivative is identity
modulo I. Nilpotence makes this procedure terminate. It constructs both
compositions, or proves one composition and uses the same uniqueness for
the other. Thus ALL finite PARAMETER jets are invertible automatically.
This statement is distinct from (7), which concerns jets of source arcs.

The precise missing uniformity is visible without a computation. Write
the unique formal inverse of K at zero as

    V(z)=z+sum_(k>=2) V_k(z), V_k homogeneous of degree k.

Then the inverse of (10) as a formal parameter family is

    V_t(z)=t^-1*V(t*z)
          =z+sum_(k>=2) t^(k-1)*V_k(z).                 (11)

Modulo t^r this is a polynomial in z. The degrees in z of these polynomial
inverses are uniformly bounded as r varies IFF all sufficiently high V_k
vanish, IFF K has a polynomial inverse. In that case (11) is also polynomial
in t. Therefore a proposed passage from every nilpotent parameter jet to
a polynomial inverse must provide a UNIFORM SOURCE-DEGREE bound. Individual
existence of inverses in every truncation has already passed for ANY
hypothetical counterexample. No such uniform bound is proved here.

A standard positive-characteristic control makes this distinction concrete
without claiming a characteristic-zero counterexample. Over characteristic
p, for ANY prime p,

    E_t(x,y)=(x-t*x^p,y), det DE_t=1.

Every t!=0 specialization is a finite étale plane map of degree p, whereas
t=0 is identity. Its inverse modulo t^r is the finite truncation of

    x+sum_(j>=1) t^((p^j-1)/(p-1))*x^(p^j), y.           (12)

The identity follows by applying Frobenius and cancelling consecutive
terms; inverse degrees grow without bound with r. This is the standard
Artin--Schreier phenomenon, not a new example or a characteristic-zero
HC4 counterexample. Its role is to show exactly where uniformity is lost.

## 5. Connected moduli do not supply the missing open set

BGV Corollary26.9 proves connectedness of the fixed-degree determinant-one
moduli scheme. Lemma26.10 proves that geometric degree<=m is a CLOSED
condition; degree1 is the automorphism locus. This is a closedness result,
not the openness needed to extend invertibility away from identity.
The complete relevant proofs were read; no contrary semicontinuity is
assumed. Fixed-degree smoothness is also not available in characteristic
zero for n,q>=2, by the same paper's Theorem26.17.

There is an especially sharp scope check. Suppose that, for every degree
bound D, there were a Zariski OPEN neighbourhood of identity in the plane
Keller moduli of degree<=D consisting of automorphisms. Pull back that
neighbourhood along (10). It is an open subset of A1 containing0 and hence
some t!=0. For that t, K_t is an automorphism, and linear conjugacy then
makes K an automorphism. Thus even this apparently local openness assertion
would ALREADY PROVE JC2. The converse is immediate on geometric points if
JC2 holds. It cannot be used as a harmless consequence of connectedness,
of every parameter-jet inverse, or of the implicit function theorem.

No fixed-degree moduli calculation or finite-order deformation campaign
is proposed. A genuine new theorem giving uniform inverse degree or this
open neighbourhood would be substantial progress, but is not supplied by
the cited connectedness theorem.

## 6. What changes in the next research action

The three routes do not add automatic regularity to the finite trace
algebra. Bézoutian source regularity is automatic but finite-normalization
descent is exactly the missing trace condition. All finite source jets
preserve that condition by faithful polynomial extension. Every finite
parameter jet of the scaling family is invertible already; the missing
uniform degree bound or moduli openness would settle the original problem.

Revisit these routes only with an ACTUAL polynomial-plane theorem proving
that finite descent/uniformity, or with an actual noninverse Keller pair.
Do not replace the full original objective by raw Bézoutian constancy,
jet lifting, formal parameter inverses, or moduli connectedness. The
geometric pseudo-plane compression PF128 remains a separate retained
restriction; neither its solvable surviving monodromy nor an abstract
intermediate surface is assumed to give inversion or a plane model.

FULL THEOREM TARGET: prove or break full HC4/JC2.
CURRENT LAST UNSOLVED IMPLICATION:

    ACTUAL polynomial plane source
       ==> Tr(Lambda_L*Lambda_L) subset B.                OPEN

Full marked-energy integrality is also OPEN. No finite controls are used
as an all-degree proof, and no registered mathematical checker is replayed.
No primary cache is created. No independent replay, formal verification
or external review is claimed. Review and outreach remain DEFERRED.
