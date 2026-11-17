# The highest next-Jacobian block retains all sixteen cusp candidates

14 November 2026. PF2RULED23 continues the complete cusp-lift calculation
of [PF2RULED22](FULL_CUSP_FIRST_JET_LIFTS.md). On every one of its sixteen
marked noncontact pole(6,9) candidates, five of the eight highest next
Jacobian conditions vanish identically on the entire affine lift space.
The other three admit exact solutions. The high block removes no
candidate. Its rank-three open locus has 43 boundary parameters before
the eleven remaining low-coefficient membership conditions are imposed.
Those conditions, the rank-drop strata, higher normal coefficients and
global finiteness remain unresolved.

## Exact scope and equations

Work in the specified Wright ring and actual omitted-prime chart of
[PF2WRIGHT1](WRIGHT_MARKED_MODEL.md). Write

    F=F0(v)+u*F1(v)+u^2*F2(v)+...,
    h=v-s, lambda=3h/4+delta, L(F)=h*F1-lambda*F0'.

For each marked common nonic H, all first-cusp-jet lifts have

    P=T^2+R5, Q=T^3+R8, L(P)=L(Q)=0,
    R5 in V5, R8 in V8,

where T is any V3 lift of H. Their complete affine dimensions are
31 and 85. Put p=P0'/h, q=Q0'/h, W=p*q'-p'*q. The next necessary
Keller condition is

    J1=2*h*(P2*q-p*Q2)+lambda^2*W,
    kappa=delta^2*W(s), G=(J1-kappa)/h=0.                 (1)

The predecessor gives deg G<=52. Its eight coefficients in degrees
45,...,52 are independent of the 34 second-coefficient controls.

The new exact identities, on every one of the sixteen complete affine
lift spaces, are

    [v^j]J1=0 for j=49,...,53,
    equivalently [v^j]G=0 for j=48,...,52.               (2)

Equivalence uses division by the monic linear polynomial h and the
predecessor degree bound. The remaining high block is exactly

    [v^46]J1=[v^47]J1=[v^48]J1=0,
    equivalently [v^45]G=[v^46]G=[v^47]G=0.             (3)

These are generally quadratic equations in the full pair parameters.
Their solution does not imply the remaining coefficients of (1).

## Complete section coverage and exact identities

The checker uses the full normalized Cox basis of V_m. Its element
indexed by 0<=j<=m, 0<=i<=m+j has actual chart expression

    (-1)^(m-j)*v^(m-j)*(1-v^2-u*v)^j*(u+v)^(m+j-i),

and boundary symbol r^(2(m-j)+i). The sign is required to match this
symbol to the specified marked frame; omitting it changes H. The basis
has (m+1)(3m+2)/2 elements, with no lower-pole coefficients omitted.
The geometry establishing these complete spaces is the existing
[PF2COORD1 input](SOURCE_COORDINATES.md).

For either candidate polynomial P8, first prove P8 irreducible over Q.
Work exactly in K=Q[y]/(P8), with the predecessor's x(y), H, delta and s.
This simultaneously covers its eight complex embeddings. No selected
numerical root or parameter sample is used. Choose one Cox lift T,
solve the complete linear operators L on V5 and V8, and compute full
kernel bases U_i (31 columns) and V_j (85 columns). This represents
all pairs as

    P=P*+sum a_i*U_i, Q=Q*+sum b_j*V_j.

The u coefficient of the bracket is the bilinear expression

    B(A,B)=2*(A2*B0'-A0'*B2)+A1*B1'-A1'*B1.           (4)

Thus verifying (2) on this entire affine family requires its reference,
all P-linear terms, all Q-linear terms and every mixed term. The checker
verifies the reference and both full linear blocks in degrees49..53.
For the mixed block, the Va,Vb degree bound gives

    deg B(U_i,V_j)<=4*(5+8)-3=49.

It verifies coefficient49 for every one of the 31*85 pairs. All higher
mixed coefficients vanish by the displayed degree bound. This proves
(2) on the full 116-dimensional space, rather than only on witnesses.
The finite exact matrix identities are computational proof inputs;
no independent whole-proof implementation is claimed.

## A rank-three solution over every marked algebra

Fix the reference Q*. Let C be the 3-by-31 matrix

    C_(j,i)=[v^(46+j)]B(U_i,Q*), j=0,1,2.

In each candidate algebra the checker produces three columns whose
3-by-3 determinant is nonzero. It retains the complete matrix, its
exact determinant and inverse in K. These are units at all eight
embeddings, not just at one representative. Consequently (3) is a
surjective affine linear system in P for this fixed Q*. Its solution
space has dimension28. The certificate retains exact global-section
coefficients of one solution P and Q at every embedding, expressed as
corrections to T^2 and T^3 in the complete V5 and V8 Cox bases.

For variable Q, the same minor remains invertible on a nonempty Zariski
open set. Solve the corresponding three P controls there. This gives
an open high-block solution chart with 85+28=113 parameters. The full
first-jet space has 46 boundary parameters and 34+36 normal parameters.
The high block ignores all 70 normal parameters by PF2RULED22's degree
bounds, so this chart has 46-3=43 boundary parameters. These dimensions
describe this rank-three chart; no claim about other components or the
rank-drop strata is made. Those strata must be retained in a complete
exclusion.

## The clean singleton cusp can still be imposed

Both exact witnesses have gcd(p,q)=1 and W(s)!=0. More generally the
three modifications used in PF2RULED22 remain available:

    P <- P+alpha*T2,
    Q <- Q+beta*T2+gamma*T3,
    T_k=(D-s)^k-k*A1*lambda(D)*(D-s)^(k-2), k=2,3.

They preserve the leading symbols and the entire first cusp equation.
Their brackets with the corrected pair, and their mutual brackets, have
u coefficient of v degree<=45. Thus they also preserve all eight high
conditions. The prior two resultant arguments and the nonzero
second/third derivative wedge polynomial apply unchanged. They supply
a nonempty open set with derivative gcd exactly h, a singleton ordinary
cusp and birational boundary parametrization, still after solving (3).
This is a boundary statement; it establishes neither surface finiteness
nor the source Keller identity.

## The remaining next-normal gate

After (3), G has degree<=44. For fixed clean boundary data its remaining
normal freedom is precisely

    G <- G+2*(A*q-p*B), deg A<=10, deg B<=22.

The predecessor proves that this map is injective, of rank34, when
deg p=22, deg q=34 and gcd(p,q)=1. Hence the remaining polynomial must
belong to the image of the same 45-by-34 matrix, giving eleven local
consistency conditions. On the rank-three high-block chart these are
conditions on 43 boundary parameters. Their dependence on these
parameters is rational after choosing nonzero minors; they are not
asserted independent, quadratic, soluble or inconsistent.

The certificate witnesses have deg G=44. Their normal-control matrix has
rank34 and its augmentation by -G has rank35, exactly over each candidate
field. Thus these particular clean boundary witnesses cannot pass the
lower test even after all second-coefficient corrections. This does not
exclude other boundary data in the complete high-block families.
The high block is therefore
solved without solving the whole next-normal equation. No candidate is
excluded or certified as a finite map. The 36 directions invisible
through u^2 must return at higher normal order. Keep the high-block
rank-drop strata, the other height-nine pairs, all contact cases,
higher heights, all-height pole-four coordinates and actual compatible
embedding in the original programme. The first open reduced pairs
remain (5,9),(6,9),(7,9),(8,9); the full goal remains active.

## Replay and assurance

Run from the repository root:

```sh
sage -python research/jacobian/affine-complements/verify_cusp_high_jacobian.py
.venv/bin/python research/jacobian/affine-complements/audit_pole_integration.py
```

The [checker](verify_cusp_high_jacobian.py),
[exact certificates](cusp_high_jacobian_certificates.json) and
[receipt](cusp_high_jacobian_receipt.json) retain full-kernel coverage,
all mixed top coefficients, exact rank minors and actual section
witnesses over both degree-eight fields. Sage and Python versions are
recorded in the receipt. Section geometry, the clean-open argument and
the local matrix-chart interpretation use written proofs. No formal
verification, independent full reproof, external review, new external
theorem or novelty priority is claimed. All finite surface and actual
normalization obligations remain open.

