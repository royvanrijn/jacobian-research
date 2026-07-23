# Degree-thirty Hessian synchronization: exact frontier and attacks

This note isolates the remaining synchronization work after the global
all-six degree-thirty intersection was closed.  It distinguishes three
questions which should not be conflated:

1. equality of the six canonical linear lifts on the all-cut scheme;
2. equality on every two-cut scheme;
3. classification of every reduced and nonreduced component by the
   coefficient-decorated Ritt 2-complex.

The first is proved.  The second has five remaining pairs.  The third remains
open even though the six complete-decomposition charts of the Dickson braid
are classified scheme-theoretically.

## 1. Closed global statements

For an outer cut \(d\mid30\), write \(H_d\) for its Hessian residual ideal
and \(\lambda_d\) for its canonical missing linear coefficient.  Exact
factor-chart reductions prove

\[
 \lambda_e-\lambda_d\in H_d+H_e
\]

on the five edges

\[
 (6,2),\ (6,3),\ (3,15),\ (5,15),\ (5,10).
\]

Their Groebner-basis sizes are respectively

\[
 11,\ 6,\ 95,\ 6,\ 11.
\]

These edges form the spanning tree

\[
 2-6-3-15-5-10.
\]

Every edge ideal is contained in
\(H_2+H_3+H_5+H_6+H_{10}+H_{15}\).  Transitivity therefore proves that all
six lifts agree scheme-theoretically on the global all-cut intersection.
This argument sees every component, including any component missing all six
complete-decomposition charts.

The nested pair \(\{2,10\}\) is also closed independently.  On the
`10 o 3` source chart, the common-refinement model is `2 o 5 o 3`.  The top
five coefficients of the degree-ten outer factor recover the `2 o 5`
parameters triangularly.  Replacing its four lower coefficients by
deviations from that graph gives a `4 normal | 7 base` ring.  The transformed
Hessian ideal has a Groebner basis of size four and the lift defect reduces
to zero.

Four further factor-chart reductions certify the non-tree pairs

\[
 \{5,6\},\ \{6,10\},\ \{6,15\},\ \{10,15\}
\]

with basis sizes `502, 189, 12, 96`.  Thus ten of the fifteen degree-thirty
pairs have exact characteristic-zero certificates.  The fast default
regression uses only the spanning tree; the four larger reductions have the
separate replay target

```bash
make audit-degree30-hessian-synchronization-pairs
```

## 2. Exact remaining pair frontier

After the nested certificate, the unclosed pairs are

| outer cuts | common right degree | primitive Ritt core | first normal-form attack |
|---|---:|---:|---|
| \(\{2,3\}\) | 5 | \(2\) versus \(3\), degree 6 | transport the complete degree-six collision through a generic quintic |
| \(\{2,5\}\) | 3 | \(2\) versus \(5\), degree 10 | separate power and Dickson core ideals, then transport through a generic cubic |
| \(\{3,5\}\) | 2 | \(3\) versus \(5\), degree 15 | separate power and Dickson core ideals, then transport through a generic quadratic |
| \(\{2,15\}\) | 1 | \(2\) versus \(15\), degree 30 | direct coprime power/Dickson normal forms |
| \(\{3,10\}\) | 1 | \(3\) versus \(10\), degree 30 | direct coprime power/Dickson normal forms |

The common right degree is

\[
 r=\gcd(30/d,30/e).
\]

Stripping the generic degree-\(r\) right factor produces the primitive core
listed in the table.  This gives a finite componentwise problem rather than
five unrelated ambient eliminations.

## 3. Rejected raw attacks

Two failures are now informative performance data, not mathematical
counterexamples.

- The transported Dickson normal chart for \(\{2,5\}\), with five normal and
  four base coordinates, did not finish an exact full-basis computation in
  300 seconds.
- For \(\{3,5\}\), raw factor-chart computations over the good prime `32003`
  timed out in 45 seconds for `dp`, `Dp`, and the block orders
  `5|4`, `4|5`, `3|6`, and `6|3`.

Changing monomial order alone is therefore not the next attack.  A full
Groebner basis is also stronger than the desired membership statement.

## 4. Ranked gap-closing attacks

### Attack A: componentwise graph-normal certificates

For each row of the frontier table:

1. write the normalized power and Dickson core parametrizations;
2. compose on the right with the generic degree-\(r\) polynomial;
3. recover the model parameters from leading source-chart coefficients;
4. replace all remaining coefficients by graph-normal variables;
5. reduce \(\lambda_e-\lambda_d\) modulo the transformed Hessian ideal.

This is the proven pattern behind the final degree-24 pair and the new
degree-30 nested pair.  It removes irrelevant base directions before
elimination.

### Attack B: prove component completeness before scheme membership

A componentwise vanishing calculation is insufficient unless the normal
forms exhaust the two-cut ideal.  On each source chart:

1. saturate by the exact-degree and factor-reconstruction denominators;
2. compute minimal primes modulo several good primes;
3. match their dimensions and degrees with the power/Dickson graph ideals;
4. lift the graph ideals to characteristic zero and prove containment both
   ways;
5. compare primary structure on their monomial overlap.

This separates reduced component classification from nilpotent
synchronization and prevents a dense-family test from being mistaken for a
global theorem.

### Attack C: targeted syzygy certificates

The desired output is one identity

\[
 \lambda_e-\lambda_d=\sum_i q_iR_i,
\]

not a complete basis of the residual ideal.  After Attack A selects normal
coordinates, compute a modular bounded-degree syzygy for the multipliers
\(q_i\), reconstruct its rational coefficients, and verify the identity by
direct expansion over \(\mathbb Q\).  The final artifact should store the
multipliers, making verification cheaper than discovery.

### Attack D: a differential proof on primitive cores

If two synchronized candidates differ only by their missing linear terms,
their derivatives satisfy

\[
 A'(B)B'-C'(D)D'=\text{constant}.
\]

On each primitive power or Dickson normal form, compare ramification
divisors and show that the constant is zero.  Then linearize this identity
normal to the component to control nilpotent thickenings.  Success here would
replace the three transported calculations by one core theorem and would be
the first step toward all degrees.

### Attack E: conductor and derived gluing

Once pair synchronization is known, compute on every commuting square and
braid hexagon:

- the normalization;
- the conductor/annihilator along the monomial divisor;
- the conormal module and cotangent homology;
- the transverse augmentation and its Tor algebra.

The degree-thirty braid already shows that different paths can have the same
normalization but different nilpotent structures.  These data, not only the
reduced relation graph, are the coefficient decoration required on 2-cells.

## 5. Completion criterion

Degree thirty is closed only after all of the following hold:

- the five frontier pair defects have exact characteristic-zero membership
  certificates;
- the power/Dickson graph ideals exhaust every reduced pair component;
- primary structures and conductors are compared on component overlaps;
- the six-chart braid data glue independently of chart;
- any off-chart all-cut components are excluded or classified.

The first global synchronization obstruction has been removed.  The current
frontier is now component completeness and scheme-theoretic coherence, not
set-theoretic braid connectivity.
