# Seven primary components of the generic quaternary null section

Claim `G4SECTIONPRIMARY20261114`, over C in characteristic zero. The
[preceding geometry theorem](SMALL_DIMENSION_COEFFICIENT_GEOMETRY.md)
gives an exact transverse model at a generic rank-two zero-Gram matrix
in dimension four. Its row relation has three nonzero entries. After
independent row scalings the model is

    H=I10 modulo ell,  ell=a+b+c+2d+2e+2f,
    A=C[a,b,d,e,f],    c=-a-b-2d-2e-2f.

The variables here are normalized transverse Gram coordinates. The
coefficient-space chart has additional smooth parameters; this is not a
statement about torsion at the ordinary coefficient origin.

The full section ideal has the following irredundant primary decomposition:

    H = (H:(bc)^infinity) intersect (H:(ac)^infinity)
        intersect (H:(ab)^infinity) intersect (H:(de)^infinity)
        intersect (H:(df)^infinity) intersect (H:(ef)^infinity)
        intersect [H+(a^35,b^35,d^35,e^35,f^35)].             (1)

In particular, the generic section's primary construction is now closed.
Its regularity is at most34, its intrinsic origin torsion is supported in
degrees at most34, and its torsion length lies in **[1635,2576]**.
The regularity bound and exponent35 are sufficient bounds, not sharp claims.
Neither the exact torsion length nor the remaining global N=4 primary
components are determined here.

## 1. A useful bound survives the failed depth approach

Use the universal ring R=C[a,b,c,d,e,f], its full moment ideal I=I10,
and U, the intersection of the six isolated universal primary components.
The previous theorem proves that ell avoids every associated prime of R/U.
Thus ell is a nonzerodivisor on R/U. It also gives an exact sequence

    0 -> F -> A/H -> (R/U)/ell*(R/U) -> 0,

where F=(U/I)/ell*(U/I) is finite of length1635. Its Hilbert function is
the polynomial (N_H-N_U)/(1-t)^5, with N_H the proved section numerator
and N_U the retained numerator of R/U. The exact polynomial has degree29.

The direct rational minimal-resolution calculation of R/U timed out.
Its monomial initial ideal L_U, however, has a complete exact resolution
over Q in less than a second. The Betti table has regularity34. Its
alternating graded sum equals the retained N_U, checked independently.
Its projective dimension is6, so its depth is0; it supplies no useful
depth>=2 conclusion about R/U.

Graded Betti numbers can only increase in a flat Groebner degeneration.
This follows by taking each degree of the Koszul complex: the terms are
free over the deformation parameter, and matrix rank can only drop at
the special fiber. Hence

    reg(R/U) <= reg(R/L_U) = 34.

For background, the initial-ideal Betti comparison is recalled in
[Murai, Section5, printed page18](https://s-murai.w.waseda.jp/papers/BPP.pdf).
No Borel-plus-powers hypothesis or lex-plus-powers conclusion from that
paper is applied to the moment ideal. The finite-matrix argument above
is the actual transfer used here.

Tensor a minimal graded resolution of R/U with R/(ell). Since ell is a
nonzerodivisor, the result remains exact and minimal over the five-variable
ring. Therefore (R/U)/ell*(R/U) also has regularity at most34. From the
displayed sequence and the finite module F of top degree29, local cohomology
gives

    reg(A/H) <= 34,
    end H^0_m(A/H) <= 34,       m=(a,b,d,e,f).               (2)

No depth>=2 assumption about R/U has been made. Additional origin torsion
from the isolated-component section is allowed throughout this argument.

## 2. Why a single modular Hilbert computation gives an upper bound

Write the five section generators with their retained integer coefficients,
and consider the resulting quotient over the DVR Z_(32003). Its generic
fiber is the rational section just described. The exact standard basis
of its F_32003 fiber has **precisely the same entire Hilbert numerator**
as the independently proved characteristic-zero section numerator.

For each degree d, its finite DVR module is a free module plus a finite
sum of cyclic p-power torsion modules. Equality of the generic and special
dimensions excludes every such torsion summand. Consequently each graded
piece is finite free. This is a direct proof of degreewise flatness;
surjectivity, equal dimension alone, and a Hilbert polynomial agreement
would not suffice.

Here is the rank argument needed for origin torsion. For fixed d and k,
take the map from degree d to the direct sum of degree d+k pieces given
by multiplication by **all** monomials of degree k. Its kernel is

    (0:_M m^k)_d.

All its terms are finite free over the DVR. Generic rank is at least
special rank, so the characteristic-zero kernel dimension is at most
the finite-field kernel dimension. Taking k sufficiently large in a
fixed degree proves the same inequality for H^0_m. Apply the identical
argument to a Groebner flat degeneration of the finite-field section to
its monomial initial ideal L. Thus degree by degree

    dim H^0_m(A/H)_d
       <= dim H^0_m((A/H)_F32003)_d
       <= dim H^0_m(F32003[a,b,d,e,f]/L)_d.                (3)

This is an upper-bound argument. It does not transfer a modular primary
decomposition, primary-component colength, or exact nilpotence index.

The monomial saturation L:m^infinity is computed exactly. The checker
also reconstructs it combinatorially as the intersection of the five
ideals L:x_i^infinity, then enumerates the full finite monomial difference
(L:m^infinity)/L. There are2582 monomials in degrees8..37. The degree35,36,37
pieces have dimensions3,2,1. Bound (2) excludes these six possible classes
in the actual section. Combining (2) and (3) therefore gives

    length H^0_m(A/H) <= 2582-3-2-1 = 2576.

The already proved injected module F supplies the lower bound1635.
The interval is not a claim that the extra941 possible classes occur.

## 3. Primary construction by graded separation

The preceding geometry theorem proves exactly seven associated primes
for H: six minimal restricted Gram primes, and m. In the order of the
first six components in (1), the minimal primes are

    (a,d,e), (b,d,f), (c,e,f),
    (b,c,f), (a,c,e), (a,b,d),

where c has the indicated linear expression in the five-variable ring.
Each isolator bc,ac,ab,de,df,ef avoids its selected prime and belongs to
each of the other five minimal primes and to m. Exact linear-ideal
reductions check these thirty-six incidence statements. Saturation by
the isolator therefore selects the corresponding isolated primary ideal.

Let S be their intersection. Since the associated-prime list is complete,
S=H:m^infinity, and S/H=H^0_m(A/H). By (2) this module has no homogeneous
part in degree35 or higher. Let

    Q0=H+(a^35,b^35,d^35,e^35,f^35).

Its radical is m, so Q0 is primary. Every class of Q0/H has degree at
least35, while every class of S/H has degree at most34. Hence their
intersection is zero, giving S intersect Q0=H. This proves (1) without
an expanded seven-way CAS intersection. The six minimal components are
indispensable, and Q0 is indispensable because S/H has positive length.

The normalized chart is a Laurent/polynomial smooth extension of this
section. Extending (1) and localizing gives primary components on that
chart; the origin prime extends to the generic null-component prime.
The orthogonal reflection exchanging the two maximal-isotropic-plane
families gives the construction on the other generic N=4 null component.
Its intrinsic generic torsion length has the same interval [1635,2576].
The degree35 powers refer to the normalized Gram coordinates, not to
global linear coordinates on coefficient space.

## Evidence and exact failure boundaries

The [checker](verify_quaternary_section_primary.py),
[certificate](quaternary_section_primary_certificate.json) and
[receipt](quaternary_section_primary_verification.json) bind the rational
and modular Betti tables, their Hilbert-numerator identities, complete
monomial saturation/difference, isolators and all retained runs.

```sh
.venv/bin/python research/moments/vanishing/verify_quaternary_section_primary.py
```

The F_32003 section run printed its complete standard basis and matching
Hilbert numerator before entering saturation. Its full process then
timed out after30 seconds. Only that completed basis phase is reused;
the run is still marked incomplete and its saturation remains UNKNOWN.
Two good-prime resolution attempts, before and after pruning the input
to33 minimal generators, each timed out after20 seconds. Principal
linear saturation of the retained section basis also timed out after20
seconds. None supplies a resolution, depth bound or saturation result.
The successful new calculations are the rational monomial resolution,
the finite-field section-initial resolution and the monomial saturation.
Refactored runners reconstruct byte-identical replay inputs.

The [source audit](quaternary_section_source_audit.json) records the selected
background reading, inaccessible downloads and assumptions not imported.
These calculations and written arguments are internal; no independent
replay, external review or formal verification is claimed.

The generic section's constructive primary description is closed. Exact
torsion and sharp truncation powers remain open. The other N=4 zero-Gram
strata, the full N=3 embedded scheme, global N=3,4 primary and nilpotence
classifications, affine factors, arbitrary cubic GMC and lower-power GVC
remain separate obligations.
