# Central image fibre: Euler, Milnor, deleted points and retained residue covers

14 November 2026. PF2OMIT19 accounts for the special level h=H_C=0
of PF2OMIT18's nonlinear pencil. The results retain the four oriented
degree-six gap21 hypotheses of PF2OMIT13. They provide exact formulas
in the actual, still unknown image incidences, a compulsory infinity
Milnor jump, and retained-residue constraints in a specified clean
stratum. They do not construct a global normalization or prove properness.

## 1. Exact Euler characteristic of the whole source central fibre

Let X0=(h composed with F)^-1(0), a reduced affine plane curve. The
Keller map is etale, so every source germ above C is analytically the
corresponding germ of C. In particular all its singularities lie above
Sing(C), and its pointwise source cardinality is the existing ledger.

Write c(y),b(y),a(y) for normalization branch counts of C,B,H. Put

    kappa=sum_(y in Sing(C))(c(y)-1), chi(C)=1-kappa,
    epsilon_R=1 in row N,0 otherwise,
    epsilon_B=1 if row B's cusp value belongs to C,0 otherwise.

Let d=deg(h restricted to B^nu), as in PF2OMIT18, and
q_B=#(B^nu points over C). Thus q_B counts distinct zeros of the
nonconstant polynomial h restricted to B^nu, and1<=q_B<=d.
For a distinct H let u=deg(h restricted to H^nu), and q_H be its
distinct zero count. Set the generic unramified point count u and its
central Euler charge z_H as follows:

| Grouping | u | z_H |
|---|---:|---:|
| A, H distinct | degree of h on H^nu | q_H |
| A, H=C | 0 | 1 |
| A, H=B | d | q_B |
| B,N,C | 0 | 0 |

The vertical prime H=C contributes1 here, since its normalization
is A1. It is not replaced by zero just because a generic nonzero
level misses it.

Euler integration of the full pointwise source count over C gives

    chi(X0)=6chi(C)-2-epsilon_R-e_S*q_B-epsilon_B-z_H
           =4-6kappa-epsilon_R-e_S*q_B-epsilon_B-z_H. (C1)

Indeed the integral of c(y) over C is chi(C^nu)=1, whereas the
integral of b(y) is q_B. The excess-one cusp contributes exactly its
indicator. The same normalization argument gives the z_H column.
No tangency multiplicity is substituted for q_B or q_H.

By PF2OMIT18 the generic source Euler characteristic is
-138-e_S*d-u. Therefore the exact central Euler jump is

    Delta0=142-6kappa-epsilon_R+e_S*(d-q_B)
                         +(u-z_H)-epsilon_B.        (C2)

The h-level intersection degree d is a multiplicity sum, while q_B
is a distinct normalization-parameter count. Their difference is a
nonnegative projection ramification defect. If epsilon_B1, the cusp
parameter is at least a double zero of h on B, so d-q_B>=1.

## 2. Source Milnor numbers, normalization and the jump at infinity

Let r(y) be the retained source cardinality from PF2OMIT6/13. Etaleness
gives the exact finite central Milnor number and normalization formula

    mu0=sum_(y in Sing(C)) r(y)*mu_y(C),
    chi(normalization(X0))=chi(X0)
                           +sum_y r(y)*(c(y)-1).     (C3)

These distinguish original source points from the separated branches
of its normalized curve. No connectedness of the special fibre or
of its normalization is inferred from the integral generic fibre.
For example x*(x*y-1) has an integral generic fibre, but its zero fibre
is the disjoint smooth union of x0 and xy1. This polynomial is only a
countercontrol, not a residual Keller realization.

Let f=h composed with F have degree D_f. For any reduced degree D_f
projective plane curve, normalization and the arithmetic-genus formula
give

    chi(projective curve)=D_f*(3-D_f)+sum mu_p.

All f-levels have the same physical points on the infinity line,
determined by the leading homogeneous form. Subtracting the general
level, which is smooth affine, proves

    Lambda_infinity(0):=sum_(p at infinity)
                  (mu_p(projective X0)-mu_p(projective generic fibre))
              =Delta0-mu0.                          (C4)

This identity allows reducible reduced special fibres. Reducedness
holds here because an etale pullback of each reduced h-level is reduced.

The complete first-image catalogue makes(C4) quantitative. For an
immersed C, let N,T,A count nodes, ordinary triples and tacnodes.
Then N+3T+2A=12, kappa=N+2T+A. Put t_B=d-q_B, and, when H is
distinct, t_H=u-q_H. All these defects are nonnegative. In row A
with H distinct and in row B, put

    W=sum_(y in Sing(C))(2b(y)+a(y))*mu_y(C) in A distinct,
    W=sum_(y in Sing(C))2b(y)*mu_y(C) in row B.

The cusp of B cannot belong to Sing(C), so it does not modify W.
The exact results are:

| Row/grouping | Lambda_infinity(0) |
|---|---|
| A, H distinct | 46+12T+4A+W+2t_B+t_H |
| A, H=C | 69+6A+2t_B; triples already forbidden |
| A, H=B or row C | 46+12T+4A+3t_B |
| B | 46+12T+4A+W+2t_B-epsilon_B |

Here mutual singular-locus avoidance in the shared and cyclic rows
gives mu0=2N+6A; the shared H=C row instead has no source points
over any singular C value, so X0 is smooth and mu0=0. In row B,
t_B>=epsilon_B keeps the stated lower bound valid.

For row N, class I has8 nodes, one A5 value and one clean cusp;
class II has11 nodes and one clean cusp. Their kappa values are9,11,
and their finite central Milnor baselines, when B avoids Sing(C),
are32,28. Put W_N=2sum_(Sing(C)) b(y)*mu_y(C)>=0. Then

    class I:  Lambda_infinity(0)=55+2t_B+W_N,
    class II: Lambda_infinity(0)=47+2t_B+W_N.         (C5)

Thus **every residual row requires Lambda_infinity(0)>=46**. Shared
H=C requires at least69, and N requires at least55 or47 in the two
classes. This is an exact obligation on the original source polynomial's
projective central fibre, not a claim that its unknown infinity
singularities have been constructed. The generic nonlinear-pencil
infinity count of PF2OMIT18 is not this Milnor-number difference.

## 3. The finite retained cover over C's normalization

Let Z be the normalization of the reduced base change N times_A2 C^nu.
Before reduction the generic R component has multiplicity2. After
normalization it is a degree-one component over C^nu. Remove that
whole omitted component, and also the degree-one D component if H=C.
The remaining normal finite curve E_C over C^nu=A1 has degree

    m_C=4 if H!=C, and m_C=3 if H=C.                 (C6)

These are exactly the generic retained residue components along C.
Their actual residue degrees sum to m_C. The normalization of X0 is
the open of E_C obtained by deleting all its points mapping into the
original omitted boundary. Reduction of the base change has generic
degree5, not6; the rank2 of the removed original R must not be replaced
by its reduced degree1 in the original sheet accounting.

Over a parameter of C^nu, a cyclic omitted block of original index e
and contact multiplicity m has equation w^e=t^m after absorbing a unit.
Put g=gcd(e,m). Its normalization has

    g points, index e/g each,
    different e-g, deleted-point count g,
    (different)+(deleted points)=e.                  (C7)

All g points are deleted because they map to the same original omitted
point. If e divides m, all are unramified and still deleted. The
original N fibre had one point of ambient length e. It has not acquired
g original omitted divisors or g original N points.

For each ordered pair of distinct C branches at a singular value,
the other R branch supplies(C7) with e2 and their contact m:1 for
nodes and ordinary triple pairs,2 for a tacnode,3 for A5. If H=C,
the corresponding other D branch supplies one further deleted
unramified point. The own D component was removed in(C6).
Each immersive B branch over C supplies(C7) with e=e_S. Each
unramified H branch supplies one deleted unramified point, whatever
its tangency order. Distinct original omitted points are kept separate.

Two noncyclic cusp blocks require their full local cubic instead:

* At the own cusp in row N, use (u,w)->(u,w^3+u*w) and
  C^nu:(u,v)=(-3t^2,2t^3). The inverse is
  (w-2t)*(w+t)^2. Removing the whole R component leaves the degree-one
  residual branch w=2t. Its point t0 is unramified but deleted because
  it maps to the original omitted cusp point. This is not an extra
  omitted divisor in row N.
* If row B's cusp lies on C, that C branch is smooth. A branch
  transverse to the cusp tangent gives w^3=t: one deleted index3
  point, different2. A tangent smooth branch gives w*(w^2+t)=0:
  two deleted points of indices1,2, different1. The original ambient
  length is3 in both cases. These are exhaustive: intersection with
  the ordinary cusp is2 in the first case and3 in the second. Using
  the generic e_S2 alone would give the wrong block.

At a finite base parameter all ramification of E_C occurs at deleted
points; the retained map is etale. But the converse fails, as(C7) and
the cusp examples show. The rank at that parameter is the sum of the
retained r(y) single sheets and all the residual omitted block ranks.
It is m_C, including when the original source point has several C
branches which are separated in this curve normalization.

## 4. Retained residue degrees on the clean nodal stratum

Now add the precise hypotheses: C has twelve ordinary nodes; B avoids
them; every C/B intersection is transverse between smooth branches.
The shared-image and cyclic rows force the required mutual avoidance.
No global degree or nodal hypothesis on B away from C is imposed.

Over A1=C^nu there are24 node-preimage parameters, each giving one
transposition on E_C, and d crossing parameters giving one e_S-cycle.
There are no other finite ramification values in the following rows.
Their local groups generate the finite-cover monodromy because the
punctured affine line is generated by its branch meridians. The
small-cycle graph argument of PF2OMIT8, applied to each orbit, yields:

| Row/grouping | Possible retained residue-degree partitions along C |
|---|---|
| A, H=C | 3; or2+1 |
| A, H=B | 4;3+1;2+2; or2+1+1 |
| C | 4; or3+1 |

In the first row a connected cover has S3 monodromy, otherwise S2
on its double component. In the second row the group is the product
of the symmetric groups on the displayed components. In row C,
the3-cycle requires an orbit of size at least3; a transposition cannot
lie in a separate two-letter orbit with only four letters available.
Thus the groups are S4 or S3 on a triple component with a fixed section.
These are necessary alternatives, not constructions of every option.

There are useful concrete consequences in the split cases.

In row C, a degree-one component is A1 over C^nu and is unramified
at every finite parameter. Every possible deleted point in this clean
stratum has index2 or3 by(C7), so none lies on that component. Its
whole A1 maps into the original source A2. Consequently a split3+1
case requires a **polynomial lift of the entire known parametrization**:

    P(x(t),y(t))=p_C(t), Q(x(t),y(t))=q_C(t),
    x(t),y(t) in C[t].                               (C8)

The lift's image need not be an embedded line; this does not invoke
an injectivity-on-a-line theorem or prove a Keller contradiction.

In row A with H=C, a split2+1 case behaves differently. At each node
parameter the unique fixed sheet is the other unramified D point and
is deleted; at every crossing it is retained. The section is therefore
A1 minus exactly24 known node-preimage parameters in the source.
The connected double component has all24+d branch points as its
deleted finite points. If A_C(t) is the monic node-preimage polynomial
and R_CB(t) the squarefree intersection polynomial, its function field is

    C(t)(sqrt(A_C(t)*R_CB(t))),
    compact genus=11+ceil(d/2).                       (C9)

The two factors are coprime under these clean hypotheses. Source
coordinates on the section are rational functions with poles confined
to those24 parameters, not necessarily polynomials. Its Euler value
is-23; the double component's source Euler value is-46-2d. Their
sum-69-2d agrees with(C1), and X0 is smooth in this grouping.

Normalizing E_C does not restore gluing at original source nodes.
Outside the smooth H=C grouping, source branches over a C-node can
belong to one or several residue components. Original retained fibre
lengths still come from the ambient normalization ledger, not merely
the number of points of this auxiliary normalized curve cover.

## 5. Remaining inputs

The [replay](verify_central_image.py) checks the catalogue Milnor and
Euler charges, all finite row rank completions, cyclic/cusp blocks,
the infinity-jump formulas and the residue permutation alternatives.
The finite-base-change, Euler integration, plane genus/Milnor and
monodromy interpretation remain written arguments. No actual incidence
q_B,q_H,W, residue partition, polynomial lift, source infinity singularity
or global cover is chosen by these controls.

The actual B/H equations and intersections, the full normal algebra
and source-coordinate regularity open are still required to fill the
exceptional-value ledger or exclude all remaining cases. This result
does not prove unrestricted singleton-zero properness or whole r7.
