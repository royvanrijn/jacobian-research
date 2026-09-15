# Determinant1020: an explicit fully rationally marked K3

**Subsequent physical construction:** the [nef MW17 pencil](DET1020_PHYSICAL_MW17_FIBRATION_2026-09-15.md)
now has an exact source-basis transport and seventeen section classes.
The unresolved transport statements below describe this earlier source
certificate; equation and coordinate compilation remain open.

An exact K3 over Q now has geometric Picard rank19 and its entire
Néron–Severi group defined by rational divisors, with determinant-1020.
The [equation certificate](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/explicit-source-certificate.json)
gives every coefficient of y²=x³+f(t)x+g(t), the section P=(A/3,C),
and square roots of all four nodal tangents. Coefficient arrays are in
ascending powers of t. This source fibration has MW rank1. The
[rootless frame theorem](DET1020_ARITHMETIC_MW17_EXISTENCE_2026-09-15.md)
therefore gives an arithmetic MW17 fibration on this explicit surface,
but its equation and seventeen section coordinates remain to compile.

## Rational construction and restored splitting

The [central-I10 chart](DET1020_CENTER_COMPONENT_SOURCE_CHART_2026-09-15.md)
was reduced by temporarily forgetting the orientations of the other
fibers. A smooth nonsplit F11 point yielded a formal branch; rational
relations reconstructed from it were subsequently verified exactly over
Q(r). Formal fitting was discovery only. The retained rational family has

    e=(r²-106)/84,
    lambda=(r-16/5)(r+8)/(r²-8r/3+136/3),
    mu=(r+128/9)(r+8)²/(r³-50r²/9-9740r/27-43520/27).

Its full coefficient arrays are in
[rational-coarse-family.json](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rational-coarse-family.json).
The exact section identity and discriminant factorization are checked by
`verify_rational_coarse_family.py` in that packet. The four fibers are I10,
I3, I3, I4 at0,1,lambda,mu, and the remaining discriminant is a squarefree
quartic. The section meets components5,1,0,0. The I3 at1 is already split.
The three remaining square conditions, after removing rational squares,
are

    (r-8/7)(r+8),  -(3r-32)(13r+160),  (13r+160)(r+8).

Set r=8(7+x²)/(49-x²). They are satisfied by the genus1 cover

    u²+2x²=306,      w²+2x²=50.

The map a=x², b=xuw/2 lands on E: b²=a(a-25)(a-153).
Starting with P0=(9,144) and T=(25,0), points Q=2nP0+T have the
required square classes; the retained producer checks them exactly for
n=1,2. The certified source uses n=2 and

    r=203702/35733,
    lambda=1115965697/2039997255,
    mu=-186544085931511/182748398174889.

No identification with a chosen Shimura CM origin is required for the
Picard certificate below. The first candidate r=22/13 is retained: its
three tested reductions do not prove Picard19 and are compatible with
CM discriminant-163. No characteristic-zero CM assertion is made for it.

## Full integral rational marking

The checker verifies exact discriminant orders10,3,3,4, nonzero square
tangents, and a squarefree residual quartic avoiding those locations.
Infinity is smooth. Thus the minimal elliptic surface is K3, all sixteen
nonidentity components of its four reducible fibers are rational, and
all other fibers are I1. The section has degrees(A,C)=(6,9), with leading
coefficients(3,1), so its only intersection with O is a simple one at
infinity. At0 the pointed cubic has t^5 dividing B=A²/3+f and C;
its I10 component is5. At1 it meets one nonidentity I3 component, numbered1.
At lambda and mu it avoids the node and meets the identity component.
The height is therefore

    h(P)=4+2(P.O)-5·5/10-1·2/3=17/6.

The classes O, fiber, the sixteen components and P generate a rank19
lattice U+F(-1), where F is the exact central-source Gram previously
certified. Its determinant is1020 and its discriminant group is cyclic.
An even proper overlattice could only have index2, since4 is the only
nontrivial square dividing1020. Its unique order2 discriminant element
has odd integral norm, so it is not isotropic. Once Picard19 is known,
this proves full saturation, MW torsion zero and that every geometric
NS class is an actual rational divisor class. It is stronger than local
splitting without a saturated marking.

## Exact Picard bound

Two implementations, packed-field C++ and pair-coordinate Python,
agree on these counts of the resolved surface:

| p | #X(Fp) | #X(Fp²) | quadratic trace a | reduction NS squareclass |
|---|---:|---:|---:|---:|
|41|2436|2860376|-66|-37|
|59|4548|12193016|-114|-58|

Good reduction is checked from the distinct split multiplicative fibers,
nonzero tangents and squarefree residual discriminant modulo each prime.
Counting the projective Weierstrass fibers includes infinity; resolving
the four split singularities adds (9+2+2+3)q=16q over Fq. Thus

    #X(Fq)=(q+1)² + sum_(t in P1(Fq),x in Fq) chi(x³+f(t)x+g(t)) +16q.

The nineteen known rational NS classes contribute19p to Frobenius trace.
The remaining degree3 factor, as a polynomial in T, is
(1-epsilon*p*T)(1-a*T+p²*T²). The first two counts determine its sign
and a; epsilon=1 at both primes. Neither a/p belongs to{-2,-1,0,1,2},
so neither remaining eigenvalue is p times a root of unity. The Tate
theorem gives geometric reduction Picard rank20. Over Fp² all twenty
algebraic eigenvalues equal p²; Artin–Tate gives discriminant squareclass
a²-4p². These are-37 and-58 respectively. If geometric Picard rank in
characteristic zero were20, specialization at both primes would preserve
its discriminant squareclass (the possible indices are squared).
The mismatch forces rank at most19, hence exactly19.

This uses the Tate and Artin–Tate theorems for K3 surfaces over finite
fields and the usual specialization argument. The two counting programs
are distinct implementations of the same formula; no independent review
of the whole theorem or formal verification is claimed.

## Replay and next gate

    sage -python research/elkies-k3/scripts/certify_det1020_explicit_source.py --check

This rechecks the equation, splitting, height-frame arithmetic, saturation
obstruction, both finite-field counts and their Picard consequences.
The packet preserves the discovery branch, failed finite-field windows,
rational reconstruction inputs and exact generic verification.

The remaining task is to transport the proven rootless U into this
physical source NS basis, reduce it to the nef chamber, and compile its
pencil and seventeen saturated sections. The explicit foundry objective
remains OPEN; this is not yet an MW17 equation for specialization studies.
