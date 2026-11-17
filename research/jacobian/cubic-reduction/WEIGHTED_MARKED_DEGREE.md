# Sharp degree seven in a marked weighted slice

Consider three-variable polynomial maps over a characteristic-zero field
with identity linear part, source/output weights(1,-1,-2), first-output
coefficients -3 on a^2b and -1 on a^3c, and a collision between
(0,1,c0) and(-1,1,4), with c0 free. Within this **marked weighted class**,
the minimum ordinary degree is exactly seven.

The established13-term map supplies degree seven, with c0=-8. This is
not a global three-variable degree minimum: other weights, coefficient
normalizations and collision configurations are not covered. The
minimum ordinary degree across varying dimensions remains the separate
classical degree-three result.

## Complete degree-six coefficient model

Every compatible nonlinear monomial of degree at most six is included:
there are17 atoms, comprising four in the first output, six in the
second and seven in the third. The two designated first-output
coefficients are fixed. The first collision equation forces

\[
F_1=a-3a^2b-a^3c+\theta(4a^3b^2+a^4bc).
\]

The other collision equation and identity linear part leave13 free
coefficient parameters after exact linear elimination. No sparsity
bound or sampled coefficient range is imposed.

Treating theta as a coefficient-field parameter gives a generic unit
ideal. Crucially, the retained lifted certificate clears to

\[
\theta^7=\sum_{i=1}^{25} M_i I_i,
\]

where Ii are the polynomial coefficient equations and Mi are rational
polynomials. [The complete certificate](weighted-newfirst6/parameter-field/certificate.json)
contains only51 multiplier monomials, with maximum degree six.
This exact identity forces theta=0 at every field-valued solution,
including algebraic parameter values. It does not infer an exclusion
at a pole of a generic calculation.

At theta=0, every map in the model belongs to the
[complete weighted tail slice](WEIGHTED_TAIL_CLASSIFICATION.md), with
beta=1. That classification requires the known coefficient8 on
F3's a^3b^3c term, whose degree is seven. Such a term is absent from
the degree-six model. This excludes the sole exceptional parameter
and proves the lower bound. Together with the known map, degree seven
is sharp in precisely the stated class.

## Independent replay and boundary of the result

[The independent checker](verify_weighted_newfirst6.py) uses Fraction
sparse polynomial arithmetic. It enumerates all17 allowed atoms,
reconstructs the full three-variable determinant and both coefficient
collision equations, checks completeness of the linear parameterization,
and multiplies the cleared theta^7 identity. It then checks the exact
embedding of the zero slice into the certified degree-seven tail class
and the missing nonzero degree-seven atom.
[Validation](weighted-newfirst6/validation.json).

The original full-ring Groebner and real-algebraic searches timed out.
Six rational theta slices gave primary unit ideals, but sampling is
not used in the proof. The generic parameter certificate has no
exceptional roots beyond zero. A direct unit-certificate lift at zero
also timed out; the proof above instead uses the completed tail
classification. No missing or timed-out calculation is treated as UNSAT.

The subsequent [KCUB13 classification](WEIGHTED_FULL_DEGREE7_CLASSIFICATION.md)
closes the full23-atom degree-seven marked class, including all new
first-output terms: its unique representative has13 terms. The original
degree-seven real query retains UNKNOWN as an experimental status; the
independent proof supplies the mathematical result. Unrestricted raw12-term
existence remains open.
