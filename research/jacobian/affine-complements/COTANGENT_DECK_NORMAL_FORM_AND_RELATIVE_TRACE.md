# Cotangent deck normal form and the remaining global trace gate

PF2RULED135. This is a FIELD-LEVEL reduction for every noncyclic subgroup
of the actual generic deck group, in all original plane Keller degrees.
It uses PF134's actual equivariant complete ruling and the actual target
primitive P*dQ. No new classification is imported. It does not remove
mixed bad fibres or prove the full finite-target/trace-energy implication.

Every such deck extension is, in rational symplectic coordinates, the
cotangent lift of its spherical P1 quotient. Its invariant field and
relative trace have explicit forms. A genuine affine etale control below
retains this structure and an exact target coframe but lacks the whole
plane chart, so generic normal form is not a replacement for that hypothesis.

## 1. Normalize the actual invariant primitive along the generic fibre

Let F=C(P,Q) subset L=C(x,y) be the actual plane Keller extension, with
Omega=dP wedge dQ, full etale completion U_N, and ANY noncyclic subgroup
G<=Aut_F(L), of order g. Keep the whole source A2 and all omitted primes.
PF134 gives an actual G-equivariant complete A1-ruling

    rho:U_N->P1_t

with faithful spherical G action on its base. The generic fibre has ring
C(t)[z0]. Since Omega is regular and nowhere zero there,

    Omega=h(t)*dt wedge dz0, h(t) in C(t)^*.

Indeed a polynomial coefficient in z0 of a nowhere-zero form on A1 over
C(t) is a UNIT, hence constant in z0. Replace z0 by z=h(t)*z0. Then

    L=C(t,z),  Omega=dt wedge dz.

These are rational field coordinates, not global polynomial coordinates
on the original plane or on U_N. Every group element acts on the generic
affine-line fibre by

    gamma(t,z)=(gamma(t), z/gamma'(t)+b_gamma(t)),
    b_gamma(t) in C(t),                               (1)

because it acts affinely in z and preserves Omega. Base derivatives are
nonzero rational functions.

The ACTUAL target functions are regular along that generic A1, so
P,Q belong to C(t)[z]. Let H0 in C(t)[z] be a polynomial antiderivative
in z of P*Q_z. Then

    P*dQ-dH0=B(t,z)*dt, B_z=-1,
    B=-z+c(t).

Average H0 over G. The result H is invariant and still polynomial in z,
because every action in (1) is affine in z. Since P*dQ is invariant,

    P*dQ-dH=-p*dt,  p=z+a(t), a(t) in C(t).             (2)

The coefficient of z remains exactly -1 in the average. The rational
form p*dt is invariant. Consequently the SAME rational shift removes
all the affine cocycles at once:

    gamma(t,p)=(gamma(t), p/gamma'(t)),
    Omega=dt wedge dp.                               (3)

This is an elementary finite averaging proof, not an appeal to an
unverified classification of birational symplectic groups. For comparison,
without the target primitive one may average -z*dt directly; the needed
shift is g^-1*sum_gamma gamma'(t)*b_gamma(t). Construction(2) additionally
retains a target-dependent invariant generating function H.
The original polynomial generating primitive is prior HC4MYLC0 work;
here fibre integration and finite averaging adapt it to the ruled quotient.

## 2. Explicit invariant Darboux field and a three-point deck extension

Let r=R(t) be the quotient coordinate of P1/G. The spherical quotient
has degree g and the three branch signatures retained in PF134. Put

    q=p/R'(t).

The identity R(gamma(t))=R(t) gives

    R'(gamma(t))*gamma'(t)=R'(t),

so q is invariant under (3). Moreover

    L=C(t,q),  E=L^G=C(r,q),
    Omega=dr wedge dq,
    P,Q,H in C(r)[q],
    P_r*Q_q-P_q*Q_r=1,
    P*dQ+q*dr=dH.                                    (4)

To check field equality, C(r,q) is invariant and has degree g in L,
since q is transcendental over C(t). For the polynomial statement,
C(t)[p]=C(t)[q] and q is invariant; every coefficient of an invariant
polynomial in q therefore belongs to C(t)^G=C(r).

Thus ALL the degree-g relative deck extension is in the base variable
t, given by the rational equation R(t)=r; q is unchanged. This is a
three-branch-point FIELD model, after a Möbius change on the quotient
base. It is not a normalization chart across the special fibres.
Zeros and poles of R' and of the rational shift in (2) are retained.

## 3. The relative trace separates from the fibre coordinate

For h(t) in C(t) and every integer k>=0,

    Tr_(L/E)(h(t)*q^k)
        =q^k*Tr_(C(t)/C(r))(h(t))
        =q^k*sum_(gamma in G) h(gamma(t)).              (5)

The basis 1,t,...,t^(g-1) over E consequently has relative trace Gram
entries Tr(t^(i+j)) in C(r), independent of q. This is exact field trace,
and it applies coefficientwise to all C(t)[q]. It does NOT imply that
the actual finite normalization or trace-dual lattice equals the span
of these basis vectors over C[r,q]. Nor are r,q asserted regular on U_N,
T=U_N/G, or the full finite intermediate in E. The total original trace
still factors through Tr_(E/F), whose actual integral lattice remains
part of the unresolved finite-target problem.

## 4. A free affine cotangent control with a real etale target

Take the Klein four group on P1_t generated by t->-t and t->1/t, with
cotangent action p->p/gamma'(t). Put

    R(t)=t^2+t^-2, r=R(t), q=p/R'(t).

Remove from P1 all points over r=infinity,0,2,-2. Call the resulting
affine curve B. It lies in the t affine chart; its omitted points include
0,infinity, the critical points +/-1,+/-i, and the roots of t^4+1.
Let W=B times A1_p. Its cotangent action is regular and FREE: all
nonidentity base fixed points have been removed. Its affine quotient is

    T=Spec C[r, 1/(r*(r^2-4)), q],

and W->T is a finite etale Klein-four cover. The derivative R' is a unit
on B, so these are regular coordinates on the indicated open. Both
surfaces have exact nonvanishing volume, pulled back from dr wedge dq.

There is a REGULAR ETALE target map

    e:T->A2, P=r^2, Q=q/(2*r),
    dP wedge dQ=dr wedge dq,
    [C(r,q):C(P,Q)]=2.                               (6)

The generating-function identity is

    P*dQ+q*dr=d(r*q/2).

The finite trace is explicit too: t satisfies t^4-r*t^2+1=0, and

    Tr(t)=Tr(t^3)=0, Tr(t^2)=2r, Tr(t^4)=2r^2-4.

This does not contradict PF133's noncyclic-deck fixed-field index>=5:
W has nonconstant units and infinite pi1, and contains NO whole A2 chart
(a unit on such a dense chart would be constant). The affine quotient
T also has infinite pi1. The control is NOT a canonical completion of an
actual polynomial plane source. It shows that the cotangent field model,
trace formula, generating function and regular etale target on a chosen
affine open still do not supply the missing whole-source geometry.

## 5. Source comparison and the next actual gate

The retained Danielewski audit already prevents dropping the target-map
condition from a whole-plane exact-volume chart. Control(6) independently
prevents dropping the whole-plane condition after constructing a target
on a cotangent open. Neither control excludes the genuine mixed-fibre
profile or supplies a plane Keller counterexample.

Assi's arXiv1807.03991v1 theorem on rational POLYNOMIAL components requires
a polynomial mate in the ORIGINAL source ring. Its primary abstract and
theorem statement were screened as a comparison only. The rational field
coordinates r,q in (4) do not meet that hypothesis. Similarly, the
logarithmic form dx wedge dy/(xy) in Blanc's symplectic birational group
is not the exact polynomial form used here. No new theorem from either
comparison is adopted.

The next useful implication must control actual orders across ALL the
mixed fibres and the original whole-plane chart: regular multiplication
in the finite normalization, the trace-dual pairing, or full marked-energy
integrality. More rational Darboux coordinates, spherical quotient
formulas or field traces alone do not establish it. Natural A/S, full-AGL,
cyclic and remaining noncyclic deck cases, and full HC4/JC2 remain OPEN.
The checker supplies exact algebraic controls for (1)--(6) and pins the
retained proofs. Written arguments establish the all-degree reduction.
No predecessor mathematical checker replay, independent review, formal
verification, historical priority or full-goal completion is asserted.
