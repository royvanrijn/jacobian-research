# Rejected formal root reuse in the exploratory specialization code

The first version of specialize.py created distinct number fields with the
same formal root symbol a:

    QQ.algebraic_field(AlgebraicNumber((q(w), a)))

That representation was unsuitable for reuse across different q.
SymPy's algebraic-field caching identified different supplied minimal
polynomials through the shared root expression. For example, requesting
w-2 and then w-14 could leave both fields carrying the first modulus.
The 48-term rational specialization exposed the mismatch.

The old specialization-manifest.json files are tombstones marked invalid.
The original manifests, case files and source are preserved as the
invalidated v1 experiment. None of those specialization results was promoted
to the mathematical ledger.

The corrected implementation uses distinct Dummy generators, canonicalizes
the residue generator, checks the exact monic coefficient-field modulus,
checks the polynomial-ring coefficient modulus, and evaluates the defining
polynomial at the residue generator before using the field.

Both families were rerun in fresh v2 case directories. Singular independently
reconstructed the explicit coefficient field from each minimal polynomial
and verified all 1,293 specialized minima. This independent replay does not
use the SymPy generator objects or their caches.

[The regression](check_field_generators.py) reproduces the rejected reuse
and checks the corrected fields. Its [receipt](field-regression.json) pins
the implementation and environment.

The earlier polynomial-only sparse censuses, the Moh coordinate reduction,
and the 49-term one-parameter frame theorem used no such field construction
and are unaffected.

A separate serialization correction parenthesized polynomial powers before
passing rational polynomial expressions to Singular. In its parser,
w^2/3 can be read as a noninteger exponent; (w^2)/3 is unambiguous. The passing
replay scripts contain the corrected syntax.

