# First height-six Q80 trace: a bounded rational-contact attempt

**Exploratory construction attempt, not a solution or a global exclusion.**
The rational tangent constructor for an arithmetic-genus-two bisection needs
an actual rational point on its trace-halving curve. For the first height-six
trace P3 in the literal direct11952 source (zero-based), O and the34 signed
basis sections supply no eligible contact. Higher generic words and halves
outside the generic section image remain untested.

The [checker](scripts/verify_q80_rational_contact_p3.py) uses only the actual
generic equation and17 sections. It computes all17 cancelled abscissa
comparisons x(P3)-x(2Pj), thereby covering both signs. After removing coordinate
poles, their degrees are18,24,26 or34. Each is factored over Q, the product
is checked exactly, and each factor has a recorded prime in5..997 where
its leading coefficient is a unit and it has no residue root. This proves
absence of rational roots without relying on a reported irreducibility flag.
There are no abscissa equalities in the weighted infinity chart either.
The only common coordinate pole is that of P3 itself, where P3=O.
The [result](../artifacts/generated-results/elkies-k3-q80-rational-contact-p3-v1/result.json)
retains source/checker hashes, factor degrees, all residue witnesses, infinity
checks and the two signed distance heights for each comparison.

For context, an inherited half has an additional geometric restriction.
In the [ruled quotient](SINGULAR_BISECTION_TRACE_GEOMETRY_2026-09-15.md), the
candidate section C for arithmetic genus2 has C²=1. For an inherited R,
write H=h(T-2R); its quotient image C_R has C_R²=H/2-4, hence

```
C.C_R = H/4 - 3/2.
```

Both curves are tangent to the smooth branch curve at an inherited-half
contact, so distinct C and C_R have local intersection at least2. Thus
H>=14 is necessary for an integral candidate constructed at such a contact.
Equality of C and C_R has split inverse image and is not a candidate.
The pole half O has H=6, and is excluded by this criterion. This also agrees
with its location on the negative section of F1 in the plane tangent model.
This is a written intersection argument; the script checks the actual bank,
not a formal proof of the quotient theory.

The prospective input was fixed before computing: trace index3, O and the
signed basis,25 CPU seconds and1GiB. The calculation and replay each finish
below one second. No exceptional specialization, point-fitting target or
larger enumeration is used. Replay:

```
.venv/bin/python research/elkies-k3/scripts/verify_q80_rational_contact_p3.py
```

The [correlated-gain objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md) remains
open. This attempt supplies neither a new singular bisection nor a matching
cover. A new arithmetic input, rather than reusing this basis bank, is needed
for the proposed tangent construction.
