# Confluent transverse slice for the tagged GVC(3) lift

## Status

This is an exact bounded calculation over \(\mathbf F_{101}\), not a
characteristic-zero theorem and not a proof of \(\operatorname{GVC}(3)\).
It applies the local-degeneration pattern used in Arathoon--Ball--Kvalheim,
*The Maxwell Conjecture is False* (arXiv:2607.27197): isolate a tuned leading
face, retain its first transverse term as a parameter, and test the resulting
finite model before drawing any conclusion about the full problem.
The broader Keller-map version of that protocol is recorded in
[confluent degeneration and collision persistence](CONFLUENT_KELLER_COLLISION_PROTOCOL.md).

The calculation contributes one new transversal regression to the tagged-lift
program.  On the normalized \(\varepsilon=1\) slice below, the first five
pure moments have a one-dimensional affine solution scheme, the first six
have a zero-dimensional nonempty scheme, and adding the seventh moment gives
the unit ideal.  Thus no point on this finite-field transverse slice satisfies
the pure premise through moment seven.

## The confluent coordinate

Use the degree-tagged family from
[the tagged-lift analysis](THREE_VARIABLE_GVC_TAGGED_LIFT.md):

\[
 \Lambda=\partial_t\partial_z+B(\partial_t,\partial_y),
 \qquad
 P=z(t-y)+C(t,y),
\]

where \(B\) and \(C\) are binary cubics.  Normalize the coefficient of
\(t^3\) in \(C\) to one.  The factor-compatible profile is the tuned face
\(C=(t-y)Q\).  Its natural transverse coordinate is

\[
 \varepsilon=C(1,1).
\]

Indeed, every normalized cubic has the unique expression

\[
 C=(t-y)(t^2+qty+ry^2)+\varepsilon y^3.
 \tag{1}
\]

If \(C=t^3+c_1t^2y+c_2ty^2+c_3y^3\), this is the invertible coordinate
change

\[
 c_1=q-1,\qquad c_2=r-q,\qquad c_3=\varepsilon-r.
 \tag{2}
\]

Hence \(\varepsilon=0\) is precisely the factor-compatible divisor; (1)
does not impose an ansatz restriction.  This separates the degenerate model
from its first transverse correction, in the same manner as the rescaled
limit in the electrostatic construction.

## Exact finite-field result

The diagonal identity

\[
 \Lambda^m(P^m)=
 \sum_{k=0}^m\binom{m}{k}^2(m-k)!
 \partial_t^{m-k}B^k\bigl((t-y)^{m-k}C^k\bigr)
 \tag{3}
\]

forms the pure moments without an expansion in the three ambient variables.
Over \(\mathbf F_{101}\), set \(\varepsilon=1\), retain the four
coefficients of \(B\) and \(q,r\), and compute successive Gröbner bases.

| Pure moments imposed | basis size | Krull dimension | normal form of \(1\) |
|---:|---:|---:|---:|
| \(m\leq5\) | 398 | 1 | 1 |
| \(m\leq6\) | 193 | 0 | 1 |
| \(m\leq7\) | 1 | -1 | 0 |

Thus the sixth moment leaves isolated candidates, rather than proving a
transverse branch; the seventh kills them all.  This is useful as a
regression because a search that stops at dimension zero would incorrectly
report possible all-moment continuation.

The calculation is deliberately not extrapolated in either direction:

- it does not show that every nonzero \(\varepsilon\) is equivalent to
  \(\varepsilon=1\);
- it does not compute the full seven-parameter normalized cubic ideal; and
- its empty fiber over one finite field makes no characteristic-zero claim.

Generic perturbation, the final step of the Maxwell paper, is not usable
here: the GVC premise is an infinite list of exact equalities, which a generic
coefficient perturbation destroys rather than preserves.

## Reproduction

```bash
.venv/bin/python scripts/research_three_variable_gvc_confluent_tag.py
```

This uses Singular and writes
[`three_variable_gvc_confluent_tag.json`](../artifacts/generated-results/three_variable_gvc_confluent_tag.json).
The optional `--full-chart` requests the significantly harder elimination in
which \(\varepsilon\) is also free; it was not used for the result above.
