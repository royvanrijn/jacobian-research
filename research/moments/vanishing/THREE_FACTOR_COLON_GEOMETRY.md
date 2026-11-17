# The embedded tenth-moment class and exact global nilpotence

14 November 2026. Characteristic zero. Internal exact symbolic/CAS proof;
no external review, independent specialist replay, formal verification or
priority is claimed. The complete primary decomposition below is for an
annihilator module, not for the full Gaussian moment ideal.

## Results

Use the universal Gram ring R=Q[a,b,c,d,e,f], the moments and ideals from
[the global generation theorem](THREE_FACTOR_GLOBAL_MOMENT_IDEAL.md):

    I8=(M2,M4,M6,M8), I_inf=I10=I8+(M10),
    J=sqrt(I8)=sqrt(I10)=(af,be,cd,abc,def).

The least global powers with J^n contained in these ideals are exactly

    nu(I8)=14, nu(I_inf)=13.                       (1)

The remaining difference is the embedded module

    K=I10/I8 = (R/Q)(-15), Q=I8:M10.               (2)

Here -15 records the ordinary Gram degree of M10. We compute Q completely,
including its primary decomposition. Its radical is

    sqrt(Q)=(d,e,f,ab,ac,bc),                      (3)

so K is supported on the three diagonal covariance axes. Its associated
primes are the three axis primes

    p_a=(b,c,d,e,f), p_b=(a,c,d,e,f),
    p_c=(a,b,d,e,f),

and the origin maximal ideal m=(a,b,c,d,e,f). The generic length on each
axis is62. The intrinsic finite-length origin torsion has length79.

There is an exact four-component primary decomposition

    Q=Q_a intersect Q_b intersect Q_c intersect Q_0,
    Q_a=Q:a^infinity, and likewise for b,c,
    Q_0=Q+((a+b+c)^4).                            (4)

The axis components are p_a-,p_b-,p_c-primary. The displayed Q_0 is one
chosen m-primary component of length823. This823 is not the intrinsic
origin torsion length79: an embedded primary component is not unique.
All55 generators of the colon standard basis and38 generators of Q_a
are in the [exact certificate](three_factor_colon_geometry_certificate.json).
Q_b and Q_c follow by permuting factor labels. Thus (4) is an explicit
component description, not only a list of candidate radicals.

This closes the colon and global-index questions left after the preceding
invariant calculation. A full primary decomposition of I10 itself,
constrained coefficient-pullback schemes and positive-characteristic
extensions remain separate problems. The [subsequent theorem](THREE_FACTOR_ASSOCIATED_PRIMES.md)
now classifies all ten associated primes of I10 and its six isolated
primary components. The [constructive primary theorem](THREE_FACTOR_PRIMARY_DECOMPOSITION.md)
now supplies the full I10 decomposition and intrinsic origin torsion77.
The arbitrary cubic GMC and nodal
Gaussian-width/infinity questions are unchanged.

## Reproducible characteristic-zero calculation

The retained [v3 Singular input](three_factor_colon_geometry_v3.sing)
constructs the universal Wick moments exactly over Q, computes standard
bases of I8 and I10, and forms a candidate Q by principal ideal quotient.
No coefficient specialization, modular reconstruction or floating arithmetic
is used. The calculation is limited to50 seconds and768MiB; the successful
run finished within that limit. The [terminal run receipt](three_factor_colon_geometry_v3_run.json)
binds the input and full output by SHA-256.

To turn a candidate quotient into a complete result we check both inclusion
and the full Hilbert-series identity below. The checker also reconstructs
the input moments, validates all printed zero reductions and verifies the
polynomial Hilbert identities exactly. A zero exit code is insufficient:
Singular can report input errors while returning zero. The accepted run
must have no error markers, a terminal END marker and all expected exact
conditions.

## Colon completeness by inclusion and Hilbert series

Each of the55 candidate generators q_i satisfies

    reduce(q_i*M10, standard_basis(I8))=0.

This proves Q_candidate is contained in I8:M10, but would not alone prove
that all cofactor possibilities are covered. For completeness, compute
Hilbert series directly from the three standard bases, over the common
ambient denominator (1-t)^6. The exact numerator identity is

    H_(R/I8)-H_(R/I10)=t^15*H_(R/Q_candidate).     (5)

For the actual colon, the same identity follows from the short exact
sequence

    0 -> (R/(I8:M10))(-15) -> R/I8 -> R/I10 ->0.

The forward inclusion gives a graded surjection from R/Q_candidate onto
R/(I8:M10). Equation (5) says their Hilbert functions agree in every
degree. The kernel must be zero, proving the full colon equality. This
verifies completeness without the raw-generator syzygy computation, whose
terminal timeout is retained separately.

The colon has Hilbert series h_Q(t)/(1-t), where

    h_Q=[1,5,15,32,52,63,47,4,-24,-8,-1].           (6)

The coefficient list means ascending powers of t. Negative numerator
coefficients are permitted; the actual Hilbert function is nonnegative.
The first Hilbert-function values are1,6,21,53. In a six-variable ring this
proves that Q contains no nonzero polynomial of degree at most two and
that its cubic part has dimension three, spanned by

    abc, def, af^2+be^2+cd^2.

Thus the low-degree annihilator is completely determined, not inferred
from testing a few multipliers.

## Support, axis primaries and the embedded origin

For (3), all candidate generators reduce to zero modulo
P=(d,e,f,ab,ac,bc). Conversely Q contains d^6,e^6,f^6 and

    5a^2b^2-12ab*d^2,
    5a^2c^2-12ac*e^2,
    5b^2c^2-12bc*f^2.

These force d,e,f,ab,ac,bc into sqrt(Q). The ideal P is the radical
intersection p_a intersect p_b intersect p_c, proving (3).

Let ell=a+b+c. It vanishes only at the origin on the three support axes.
The exact saturation is

    Q_sat=Q:ell^infinity=Q:ell^4,
    Q_sat=Q_a intersect Q_b intersect Q_c.        (7)

Both containments in the intersection equality are checked by standard
basis reduction. The exponent check ell^4*Q_sat contained in Q is also
explicitly verified. The axis component Q_a has radical p_a: it reduces
to zero modulo p_a and contains sixth powers of all five variables in p_a.
Saturation proves a is a nonzero divisor on R/Q_a. Since this graded
one-dimensional quotient has support the a-axis, it has no embedded origin
prime; hence Q_a is p_a-primary. The other axes follow by symmetry.
This also follows by viewing its finite, torsion-free graded module over
Q[a]. Its fiber at a=1 has dimension62, which is its generic local length.

The saturated quotient has Hilbert numerator

    h_sat=[1,5,15,32,52,57,23,1],                 (8)

over denominator (1-t). Its degree186 is3*62. The difference module

    T=Q_sat/Q=H^0_m(R/Q)

has polynomial Hilbert series

    H_T=6t^5+30t^6+33t^7+9t^8+t^9.               (9)

Thus length(T)=79. It is nonzero, proving the embedded origin prime is
present as well as the three minimal axis primes.

The chosen origin component in (4) follows without a guessed primary
ideal. Set Q_0=Q+(ell^4). If g is in Q_sat intersect Q_0, write
g=i+ell^4 h with i in Q. Since Q_sat is ell-saturated, h is in Q_sat;
then ell^4 h is in Q by (7). This proves Q=Q_sat intersect Q_0. Its
radical is m because ell vanishes on the support only at the origin, so
Q_0 is m-primary. The exact finite quotient length is823, consistently
4*186+79. Different embedded-primary choices need not have this length;
the torsion length79 in (9) is canonical.

## Exact global radical-power indices

The computation builds each power of the monomial ideal J by exact ideal
multiplication and standard-basis minimization, then reduces all of its
minimal monomial generators against the full I8/I10 standard bases.
There is no omission of a selected covariance support.

| Power of J | Minimal generators tested | Nonzero reductions for I8 | Nonzero reductions for I10 |
|---:|---:|---:|---:|
| 12 | 819 | 33 | 6 |
| 13 | 1015 | 6 | 0 |
| 14 | 1240 | 0 | 0 |

The upper containments therefore hold globally. For explicit lower
witnesses put u=af and v=be. The recorded normal forms are

    u^12*v=a^12*b*e*f^12 not in I8,
    u^11*v=a^11*b*e*f^11 not in I10.

They are individual members of J^13 and J^12 respectively, and each
remains unchanged under the corresponding full standard-basis reduction.
Together these prove the exact indices (1). The earlier bounds9..19 and
9..17 remain valid but are now sharpened.

## Reproduction and retained failed routes

Run `.venv/bin/python research/moments/vanishing/verify_three_factor_colon_geometry.py`
for a fresh bounded Singular replay and exact receipt validation. Use
`--audit-existing` for a solver-free integrity/Hilbert audit of the retained
run. The [verification receipt](three_factor_colon_geometry_verification.json)
records Singular4.3.2, SymPy1.14.0 and all current hashes. The parser accepts
only explicit polynomial terms and reconstructs the proof certificate.

The first principal-quotient computation finished exactly. A separate raw
syzygy route hit its50-second limit and remains computational UNKNOWN;
no current claim relies on it. Geometry input v1 had an ideal declaration
error, and v2 used sat where sat_with_exp was required to obtain an exponent.
Their inputs, outputs and failed classifications remain preserved in the
[raw-syzygy timeout](three_factor_scheme_structure_run.json),
[invalid declaration run](three_factor_colon_geometry_run.json) and
[invalid exponent-interface run](three_factor_colon_geometry_v2_run.json). The
current no-error v3 result and Hilbert completeness proof supersede these
routes mathematically; they do not change their historical run outcomes.

A solver-free audit verifies the recorded identities and provenance; it
is not an independent rerun of Singular's standard-basis algorithms or a
formal proof-assistant certificate. The accepted exact characteristic-zero
calculation and the written commutative-algebra arguments establish the
scope claimed here. No positive-characteristic or arbitrary-polynomial
Gaussian conclusion is imported from this result.
