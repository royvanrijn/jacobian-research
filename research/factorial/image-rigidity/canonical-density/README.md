# Can the constructed obstructions become factorial counterexamples?

## Results

**The two retained examples are excluded as canonical densities and their
moment-polynomial pairs are coprime. This does not exclude all possible SFC
counterexamples.** The results are internally reviewed mathematical proofs
with exact computational certificates, not external peer review or formal
verification.

1. [Canonical density exclusion](canonical_density_exclusion.md) proves that
   at represented degree ten neither constructed current density can occur
   at any integer profile, terminal order or depth. Their exact degrees and
   endpoint orders leave only one candidate each; the actual first remainder
   has the opposite leading sign. The two full pairs and their reflections
   are consequently excluded as well.
2. [Moment-polynomial diagnosis](moment_pair_obstructions.md) proves that
   both original degree-eleven/degree-ten pairs are coprime over every
   characteristic-zero field extension. Full modular Bezout vectors are
   retained. A separate positive Beta-mixture family shows exactly why a
   vanishing nominal multiplier and a genuine common-root collision are
   different events. At theta=72/85 the enlarged independent-law class does
   have a common quadratic factor; no linked factorial polynomial is supplied.
3. [Inverse canonical reachability](inverse_canonical_reachability.md) gives
   a finite necessary-and-sufficient membership test for a supplied pair of
   normalized polynomial densities. A bounded list of integer candidates is
   followed backward by uniquely determined coefficient cancellation, ending
   at the two exact Beta seeds. The test permits negative multipliers and
   signed intermediate densities. It does not prove that every actual chain
   avoids a zero mass.

## The precise unresolved issue

A factorial counterexample in the binary monomial-linear family requires
the actual polynomials J_n(z;r,u) and J_(n+1)(z;r,u), with the same
nonnegative integer profile, to have a common root. Equivalently, one fixed
two-term polynomial must have two consecutive zero factorial moments.
The counterexamples to generic positivity do not supply this.

Canonical-density reachability is stronger than finite-moment reachability.
The inverse note constructs two different strictly positive densities with
the same actual represented quadratic, only one of which can be canonical.
The separate polynomial coprimality certificates are therefore essential.

## Reproduce the exact checks

Use Python 3.11 or newer. Only the inverse-reachability checker requires
SymPy; it was replayed with SymPy 1.14.0. The two other mathematical checkers
and the audit harness use the standard library. Run from this directory:

```sh
python verify_density_exclusion.py
python verify_moment_pairs.py
python verify_inverse_reachability.py
python verify_audit_controls.py
```

The first three commands reconstruct the exact mathematical inputs and
write the three verification receipts. Each has an optional `--output PATH`.
The last command compares optimized runs with those receipts and tests
seven deliberate corruptions in disposable copies. It preserves all proof,
checker and certificate inputs. Proof and certificate hashes are checked
before relying on them; algebraic corruption tests refresh only the copied
fixture pin so independent mathematical checks must detect the error.

The coefficient certificate is [moment_pair_certificates.json](moment_pair_certificates.json).
The package contains no imports from an earlier checkout and no dependency
on the previous curvature/transport certificate. Those earlier properties
identify the two input examples but are not premises of these exclusions.

## Review and provenance

The explicit laws were reconstructed from the preceding research discussion
after workspace maintenance removed the former local packet. The new proof
and certificates are self-contained. The current public repository was
rechecked at commit `41f14391a65112c60fd8f7cafa5b8b45d956b8e1` and left
unchanged. No missing earlier receipt is represented as replayed here.

Root and separate research agents reconstructed the degree and endpoint
arguments, backward coefficient formulas, moment normalization and resultant
identities. Two depth-one wording issues, the nondegenerate-law hypothesis,
and an explicit initial-zero-mass rejection guard were corrected before
freezing the supplied sources. All reviews and replays are internal.
