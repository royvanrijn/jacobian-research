# The complete cyclic-frame quotient and the global twelve-moment obstruction

16 November 2026. The cyclic-open exhaustion criterion and the Hilbert
residue obstruction were supplied by the user. This note verifies their
constants, constructs the entire quotient as an exact scheme, and retains
the full first-twelve-moment arithmetic input. Global exhaustion is NOT
claimed: the dimension and length of A12 remain uncomputed.

**Subsequent gate resolved negatively.** The
[fourteen-certificate and square-system follow-up](CYCLIC_NILPOTENCE_GATE_AND_PREFIX12_OUTSIDERS.md)
proves s0^10 is not in I12 and additional characteristic-zero first12
points on D!=0. Thus the original length66 target below is false;
its [full-length refinement](ANTI_HALFTURN_SCHEME_LENGTH_AND_A14_GATE.md)
now gives dim_C A12>=13818. Global finiteness and exact length remain open.
The quotient presentation and conditional later-prefix exhaustion
criterion in this note are unchanged.

## Results

**SIC2B33CYCLIC20261116.** For the full sixteen-dimensional harmonic
coefficient representation, the SO(q) invariant ring localized at
D=det(v,Bv,B^2v) has the complete presentation in Section 1. No Gram
entry, spectral discriminant, or cubic coefficient is separately inverted.
It includes nilpotent B and isotropic v. Fixing D=1/8 and removing the
constant harmonic gives a smooth integral eleven-dimensional model in
fifteen variables and four equations. Its first three normalized moment
prefix algebras A1,A2,A3 have Krull dimensions 11,10,9, respectively.

For every N>=12 with finite-dimensional A_N,

    dim_C A_N=3d_N
      => every pure phase with D!=0 is a Rodrigues copy.             (1)

Here d_N is the actual local transverse length for the first N moments,
and d_12=22. The original N=12 target was66, now disproved by the follow-up.
The entire H3=0 control
is excluded from the normalized six-moment locus by an exact unit ideal.
The rank/Gram-rank cases in Section 6 exhaust D=0 as a stratification;
their cubic couplings remain unrestricted and are not classified here.

**SIC2B33HILBERT20261116.** The first twelve moments necessarily have
additional semistable solutions outside the Rodrigues component. This
is a finite-prefix statement, not a new all-order pure component. The
leading invariant Hilbert coefficient is exactly 317/230400. Were the
prefix quotient supported only on the Rodrigues cusp, its residue
would have to be both 659043 and 11, an impossibility. This argument
does not determine whether the extra prefix solutions lie on D=0.

The [checker](verify_cyclic_frame_quotient.py) and
[receipt](cyclic_frame_quotient_verification.json) retain exact algebra
and both residue evaluations. External, independent specialist and
formal verification are not claimed. Classical invariant-ring inputs
are credited in Section 5 and the [source audit](cyclic_frame_sources.json).

## 1. A complete SO(q) quotient on D nonzero

Use physical coordinates (z,t,w),

    q=2zw+t^2,   J=[[0,0,1],[0,1,0],[1,0,0]], det(J)=-1,
    F=c0+<v,x>+<Bx,x>+T(x,x,x),

with B self-adjoint and trace-free and T a symmetric trace-free cubic.
Set

    p=tr(B^2)/2,  r=det(B),  s_j=<v,B^j v> (j=0,1,2),
    E=(v,Bv,B^2v),  D=det(E).

Cayley--Hamilton says B^3=pB+rI. Thus, in this ACTUAL cyclic frame,

    C=E^-1 B E=[[0,0,r],[1,0,p],[0,1,0]],
    G=E^T J E=
       [[s0,s1,s2],
        [s1,s2,p*s1+r*s0],
        [s2,p*s1+r*s0,p*s2+r*s1]].                       (2)

One has C^T G=GC and, crucially,

    det(G)=-D^2.                                        (3)

The minus sign comes from the chosen physical metric. It is not an
extra assumption or a replacement of the oriented invariant D by D^2.

For an especially small moment presentation use the dual frame coordinates
xi=E^T J x. Write the cubic as

    T(x,x,x)=sum_(i,j,k=0)^2 U_ijk xi_i xi_j xi_k,
    U_ijk symmetric.

There are ten coefficients U_ijk with i<=j<=k. The only cubic relations are

    sum_(i,j=0)^2 G_ij U_ijk=0,       k=0,1,2.            (4)

They are the full harmonic trace equations, retaining every cubic coupling.
These U are regular invariant functions on D!=0: changing from physical
coordinates uses E^-1 and only powers of D in denominators.

The exact localized invariant-ring presentation is

    S[D^-1] = C[c0,p,r,s0,s1,s2,D,D^-1,U_000,...,U_222]
                    /(det(G)+D^2, the three equations(4)).          (5)

This is a ring presentation, not merely a birational parameterization.
To verify it, adjoin a frame E with E^T J E=G and det(E)=D. Since G is
nonsingular, such frames exist locally in the etale topology and form an
SO(q) torsor. Reconstruct v=E e0, B=E C E^-1 and the complete cubic
using xi=E^T J x. Relations (2),(4) ensure self-adjointness, trace zero
and harmonicity. Conversely the original cyclic frame gives these data.
Two frames with the same G,D differ by a unique SO(q) transformation.
This identifies the original open scheme with that torsor and proves
(5) by descent. It preserves nonreduced coefficient data as well as points.

In particular, no diagonalization of B or nonzero eigenvalue hypothesis
is used. Setting p=r=0 retains the entire cyclic nilpotent case.

## 2. The normalized global algebra and its exact moment input

The first moment is c0. Set c0=0, D=1/8, and write

    R=C[p,r,s0,s1,s2,U_000,...,U_222]
            /(det(G)+1/64, the three equations(4)).                  (6)

There are fifteen variables and four relations. Since D is fixed nonzero,
G^-1=-64 adj(G). The complete phase in dual coordinates is therefore

    f(xi)=xi0-64 xi^T C adj(G) xi+sum U_ijk xi_i xi_j xi_k.           (7)

For any monomial xi^alpha of even degree d,

    h(xi^alpha)=Wick_G(alpha)/(d+1)!!,

and odd degrees have zero mean. Wick_G is computed by removing the first
leg and pairing it with each possible remaining leg, using covariance G.
Consequently the [producer](cyclic_frame_quotient.py) constructs every
mu_n from (7) by exact rational additions and multiplications. This gives

    A_N = R/(mu_2,...,mu_N).                             (8)

The retained [first-twelve circuit](cyclic_frame_moments_12.circuit.json.gz)
has 733595 arithmetic nodes and compresses to 2745190 bytes. It is the
full normalized system, with no coefficient-support restriction. The
checker replays EVERY gate, checks all twelve Rodrigues moments and
compares generic original-phase controls with direct sphere integration.
The [build checkpoint](cyclic_frame_moment_build_checkpoint.json) records
the prefix sizes and a declared 1200000-node cap. This was a bounded
input construction, not an expanded global Groebner or length calculation.

Current certified global information is

| Algebra | Krull dimension | Vector-space length |
|---|---:|---|
| A1 | 11 | infinite |
| A2 | 10 | infinite |
| A3 | 9 | infinite |
| A12 | UNCOMPUTED | UNCOMPUTED; exhaustion target66 |

Here is the dimension proof, rather than an inference from sampled ranks.
The base det(G)=-1/64 is smooth: the weighted Euler field for weights
(p,r,s0,s1,s2)=(2,3,2,3,4) sends det(G) to 12 det(G), which is nonzero
on that level. It is integral. Over the generic splitting field of the
separable cubic X^3-pX-r, the three moments s0,s1,s2 change invertibly
to the squared components of v. The determinant equation becomes a
nonzero constant times their product equal to -1/64, a geometrically
integral torus equation. No vertical factor can occur: the original
polynomial is primitive over C[p,r], with a nonzero constant term.
This proves the base is a smooth integral fourfold.

The cubic trace kernel is a rank-seven vector bundle. The trace map
splits by the identity Delta(q L)=10L for a linear L, using G^-1,
which is available everywhere in (6). Thus R is smooth integral of
dimension 4+7=11. Moreover

    mu_2=s0/3+4p/15+(2/35)||T||^2.                       (9)

The harmonic cubic norm is nondegenerate of rank seven in every fibre.
Every fibre of mu_2=0 over the base is therefore a geometrically integral
affine quadric of dimension six. Locally in the etale topology its
equation can be made monic in one cubic coordinate, proving flatness.
Hence R/(mu_2) is integral of dimension ten. Finally mu_3 is nonzero
there: at p=r=s0=0,s1=1,s2=1/4,U=0 one has mu_2=0,mu_3=2/5.
It is not a unit because Rodrigues has mu_3=0. Its zero set has dimension
nine. These arguments establish the first three rows globally.

As an exact complete lower-degree regression, on the ENTIRE U=0 locus
the five first nonconstant moments through six have the lex basis

    s0+4p/5, s1+4r/7, s2+16p^2/75,
    p^3+375r^2/98, pr, r^3.

The polynomial det(G) reduces to zero in that ideal. Adding
det(G)+1/64 therefore gives the unit ideal. This excludes H3=0 from
the normalized prefix-six locus. It does not classify unrestricted
H3 or imply the target length66.

## 3. The three Rodrigues points and the length criterion

For F0=w+zt-z^3/2, the original v and B give D=1/8. In the presentation
(6), the three normalized quotient points have

    alpha^3=1,
    p=r=s0=s1=0,  s2=alpha/4,
    U_222=-32 alpha,  every other U=0.                   (10)

Indeed D(lambda F0)=lambda^6/8 and s2(lambda F0)=lambda^4/4.
The normalization forces lambda^6=1. Equal values of the degree-four
and degree-six invariants force the ratio of two scalings to be +/-1;
the rotation (z,t,w)->(-z,t,-w) realizes that sign. Thus (10) consists
of exactly three distinct SO(q) quotient points.

The action on D!=0 is free: a rotation fixing v,B fixes v,Bv,B^2v,
and so is the identity. There is no finite-stabilizer correction to
the local quotient length. Also the derivative of D along scalar phase
scaling is 6D!=0. Fixing D is therefore transverse to scaling and
leaves the normal length d_N unchanged. Phase scalings with sixth power
one identify the three local algebras.

For N>=12 the [local theorem](RODRIGUES_FINITE_JETS_AND_FULL_LOCAL_RIGIDITY.md)
makes these algebras finite, with d_12=22. A finite A_N decomposes as
the direct sum of its Artin local algebras. The three known points
already contribute 3d_N. Equality in (1) leaves no other prefix point.
Every phase with D!=0 can be rescaled into D=1/8, proving the stated
global-open pure classification if that equality is certified.

Neither finiteness nor the equality has been established for A12.
An extra prefix point may fail a later moment. Such a point is not
an additional pure component.

## 4. Updating the local target does not bound global stopping

The actual d_N is nonincreasing for N>=12 and stabilizes by N=90.
Once its filtered multiplication matrices are constructed, later
moments act in a 22-dimensional ambient algebra. Close their span under
multiplication by a,b to obtain the new ideal and d_N by exact finite
linear algebra.

There is an implementation distinction here. The local certificate
currently supplies the weighted associated graded algebra and enough
jets for its length22; it does NOT yet supply the actual filtered
multiplication matrices. One cannot substitute the weighted initials
P,Q for the full local residuals when reducing later moments. The
local nilpotence bound makes the required reconstruction finite
(coefficient jets through order nine suffice), but it has not been
performed because the N=12 target already uses the certified d_12.

The bound90 concerns only the local contribution at the three Rodrigues
points. It supplies no bound on the global moments needed to remove
other candidates in A_N.

## 5. Why the first twelve moments have other global semistable points

The invariant rings for SO(3) and the binary-form SL2 cover agree on

    V=Sym^0+Sym^2+Sym^4+Sym^6.

The cyclic open set proves dim S=13. For the Hilbert coefficient,
remove the trivial Sym^0 first, as required by
[Cayres Pinto--Herbig--Herden--Seaton, Theorem1.1(3)](https://arxiv.org/pdf/1710.02606v3).
The positive weights are (2,2,2,4,4,6), sigma=2, and the Schur partitions
are rho=(3,3,3,2,1,0), delta=(5,4,3,2,1,0). Restoring Sym^0 multiplies
the Hilbert series by (1-t)^-1 and leaves its leading coefficient unchanged.
Exact Jacobi--Trudi evaluation and an independent rational residue
evaluation both give

    gamma0=2 s_rho(2,2,2,4,4,6)/s_delta(2,2,2,4,4,6)
      =(1/(4pi)) integral_R x^2/[(1+x^2)^3(1+4x^2)^2(1+9x^2)] dx
      =317/230400.                                      (11)

Suppose the reduced prefix-twelve quotient had only the Rodrigues
scaling curve as support. Its dimension would be one. The invariant
ring is Cohen--Macaulay by Hochster--Roberts; the original author's
[primary lecture notes, pp124--128](https://sites.lsa.umich.edu/hochster/wp-content/uploads/sites/1337/2026/04/615W20.pdf)
state the theorem for polynomial invariant rings and give its proof route.
Twelve homogeneous generators of height twelve then form a regular
sequence. Their quotient Hilbert residue would be

    12! gamma0=659043.                                  (12)

The Rodrigues quotient curve has normalization parameter u=lambda^2
of coefficient degree two. More explicitly its reduced ring is the
cusp ring C[lambda^4,lambda^6]=C[u^2,u^3]. The degree-four invariant s2
and degree-six invariant D are nonzero there. Every odd-degree invariant
vanishes because F0 and -F0 are rotation-equivalent; every quadratic
invariant vanishes by the individual harmonic norms. Thus these two
restrictions generate the whole curve ring. Its Hilbert series is

    (1-t^12)/[(1-t^4)(1-t^6)],

with residue 1/2. The generic transverse multiplicity is 22 by the
local theorem and scalar equivariance. Its entire contribution is
therefore 22/2=11. Associativity of the one-dimensional Hilbert
multiplicity contradicts (12).

Thus another non-origin quotient point exists outside this curve. It
lifts to a coefficient point: if extension of its maximal ideal to
C[V] contained one, applying the Reynolds projection would put one
in that maximal ideal. This is impossible. A non-origin quotient point
has a nonzero positive-degree invariant, hence its lift is semistable.
The additional points satisfy the first twelve moments; no claim of
all-order purity follows. Their location on D=0 or D!=0 is undetermined,
so (12) does not refute the open-set target66.

## 6. The complete remaining D=0 stratification

Cayley--Hamilton makes W=span(v,Bv,B^2v) invariant under B. Self-adjointness
also preserves W-perp. The ambient form is nondegenerate in dimension
three, so D=0 has exactly these rank/Gram-rank cases:

| dim W | rank(q restricted to W) | Retained structure |
|---:|---:|---|
| 0 | 0 | v=0; B and the cubic remain unrestricted. |
| 1 | 1 | A nonisotropic B-eigenline and its invariant orthogonal plane. |
| 1 | 0 | An isotropic B-eigenline and its invariant perpendicular plane. |
| 2 | 2 | A nondegenerate invariant plane and its orthogonal line. |
| 2 | 1 | A degenerate invariant plane with invariant radical line W-perp. |

A totally isotropic plane cannot occur in a nondegenerate three-space.
The checker gives exact representatives of all five cases. These are
strata to investigate, not closed pure classifications. In particular
the cubic is not required to split across any of the displayed invariant
subspaces; such a restriction would discard legitimate couplings.

## Reproduction and next gate

Run from the repository root:

```sh
.venv/bin/python research/moments/vanishing/verify_cyclic_frame_quotient.py
```

All retained exact checks pass. The next global computation is the
Krull dimension and, if finite, vector-space length of the FULL algebra
(8) at N=12, compared with66. The arithmetic circuit is a concrete
replayable input for that computation, not its result. Higher local jets,
additional support exclusions, or finite-field counts cannot replace
the required exact global dimension/length certificate.
