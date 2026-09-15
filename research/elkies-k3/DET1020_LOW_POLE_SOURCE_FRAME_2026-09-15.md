# Determinant1020: a source frame with a disjoint section

The admitted determinant1020 NS contains a fibration frame with root lattice
A2+2A3+A8 and torsion-free MW rank1. A generator has height85/36 and is
disjoint from the zero section. This supplies a concrete alternative to
constructing the source through high-discriminant Humbert equations.
It does not yet supply a K3 equation or a physical nef basis.

## Exact lattice witness

Let R=A2+A3+A3+A8 in the standard Cartan bases. Adjoin a vector p of
norm4, pairing1 with the first simple root in the second A3 and the A8,
and pairing0 with all other simple roots. The positive rank17 Gram F1
is completely specified by this rule; its determinant is1020.

The projected norm of p is

    4 - 3/4 - 8/9 = 85/36.

Any vector outside R has squared norm at least85/36>2, by its orthogonal
decomposition with nonzero integer p coefficient. Hence R is the entire
root system. It is primitive since the first sixteen basis vectors span
it. The checker also enumerates exactly102 roots.

The discriminant group is cyclic of order1020. A pinned reduced dual
generator has norm10613/1020. Multiplication by67 matches it to the
T generator of norm-403/1020 modulo2Z. Thus U+F1(-1) has the signature
and finite form of the [admitted rational NS](DET1020_RATIONAL_MARKING_SOURCE_2026-09-15.md).
The same indefinite uniqueness and nef argument used in the
[MW17 existence proof](DET1020_ARITHMETIC_MW17_EXISTENCE_2026-09-15.md)
gives an actual elliptic fibration on that K3 over Q. All components and
sections descend through the full rational marking.

The primitive root quotient is Z, so MW torsion is zero and this is a
saturated generator. Its component classes are0,0,1,1 on A2,A3,A3,A8.
The height formula then gives P.O=0. The reducible fibers with roots
A8 and A3 are I9 and I4. The A2 fiber could be I3 or IV; the lattice
argument alone does not choose between them.

```sh
sage -python research/elkies-k3/scripts/certify_det1020_rank1_source_frame.py --check
```

The [certificate](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rank1-certificate.json)
is an exact witness. Discovery scanned only torsion-free rank1 frames
with root rank16, at most four ADE factors, height in(2,8], and P.O<=2.
All34517 component combinations in that window were tested;122 gave
even integral norms and32 marked witnesses passed the finite-form test.
There is no claim of a complete fibration classification or optimal cost.

## Exact polynomial chart, still unsolved

Put the I9 fiber at t=0, the I4 fiber met nontrivially by P at infinity,
and the other I4 at t=1. On the split I9 chart one can normalize A(0)=1
and write

    y² = x³ + A(t)x² + B(t)x + C(t)²,     P=(0,C(t)),
    A=1+a1*t+a2*t²+a3*t³+a4*t⁴,
    C=t*(c0+c1*t+c2*t²+c3*t³+c4*t⁴),     deg B<=7.

The degree bounds encode the two nonidentity component meetings. Take
the formal solution w(0)=1 of w³-Aw-2C=0 and set

    B = truncation through degree7 of -(w²-A)(3w²+A)/4.

The [symbolic producer and packet](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/README.md)
verify that the cubic discriminant D is divisible by t^9 identically.
This follows from the nodal cubic parametrization
C=w(w²-A)/2, B=-(w²-A)(3w²+A)/4. The residual D/t^9 has degree13.
The leading two coefficients of D must vanish for I4 at infinity;
D/t^9 must have a fourth-order zero at1 and a third-order zero at another
rational lambda for the proposed I3 chart. These are explicit retained
polynomial equations, not solved coefficient conditions.

For split infinity write a4=d² with d nonzero. With b_i the B coefficients,
the two infinity equations can be written

    b7 + 2*d*c4 = 0,
    d²*b6 + 2*d³*c3 + a3*d*c4 + c4² = 0.

The first is linear in c3,c4. Eliminating c3 on the open set where its
coefficient is nonzero leaves a quadratic in c4. The zero-coefficient
branch must be retained separately. No rational parametrization of this
quadratic cover is asserted.

Next solve a nondegenerate rational point, certify exact fiber orders and
splitting, the section intersections, saturation and Picard19. Such a
point need not be the previously specified intrinsic CM-double period;
a direct rational marking and Picard certificate would admit it on its
own. A finite-field point, formal solution, or nonsplit model is not that
endpoint. The possible IV chart is also not excluded by this work.

## Local nodal equations and bounded elimination

The expanded discriminant is avoidable at the two finite fibers. At a
candidate double root x=r over t=v, expand with u=t-v. Put

    H = 3r + A(v),
    g = 2*A_1*r + B_1,       g2 = 2*A_2*r + B_2,
    F_k = A_k*r² + B_k*r + (C²)_k.

Subscripts mean Taylor coefficients at v. Require f(r,v)=f_x(r,v)=F_1=0
and H nonzero. The next two critical-value conditions are

    E2 = 4H*F_2 - g² = 0,
    E3red = 2H²*F_3 - H*g*g2 + 2H*A_1*F_2 - g*F_2 = 0.

The first gives discriminant order at least3; both give order at least4.
The [exact local identity certificate](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/nodal-jet-certificate.json)
checks the expansion at x=r-g*u/(2H)+v2*u²+v3*u³. The original cubic
condition E3 satisfies E3-4H*E3red=(g-2H*A_1)*E2, so removing the factor
is valid only together with E2=0 and H nonzero. Exact orders, tangent
splitting and the global section conditions remain separate checks.

The subsequent [execution receipt](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/nodal-execution.json)
records a complete F11 slice a1=0 with1464100 base tuples and no full-chart
point. It does not exclude other slices or rational models. A60-second
unsaturated elimination and a60-second elimination on a1=1 with explicit
inverse-product open conditions both timed out without a Groebner basis.
Their full sparse inputs are preserved. No dimension or nonexistence
claim follows, and no longer elimination campaign is prescribed by these
runs. Further structural reduction is needed before increasing bounds.
