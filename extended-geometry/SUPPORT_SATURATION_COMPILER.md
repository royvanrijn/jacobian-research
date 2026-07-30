# Support-saturation compiler

The reusable implementation
[`jcsearch/support_saturation.py`](../jcsearch/support_saturation.py) compiles
one exact module presentation

\[
N\subset F=R^r
\]

together with a boundary ideal \(I\), an optional class \(c\in F/N\), and an
optional normal ideal \(\mathfrak m\).  Its primary output is a
machine-readable presentation of

\[
H_I^0(F/N)\simeq (N:I^\infty)/N.
\]

This is an algebra compiler, not a new theorem asserting that every input is
saturated.  A nonzero output is retained as a local-cohomology module rather
than treated as a failure of the program.

## Exact output

`SupportSaturationCompiler.compile` records:

- a standard basis for \(N:I^\infty\);
- generators, relations, annihilator, and associated primes of
  \((N:I^\infty)/N\);
- associated primes of \(F/N\), when module primary decomposition is
  requested;
- exact colon tests \(N:f=N\) for candidate \(f\in I\);
- the annihilator and least detected boundary exponent of a distinguished
  class; and
- finite normal jets
  \[
  N_n=N+\mathfrak m^nF
  \]
  with the images and surjectivity of the transitions from order \(n+1\) to
  order \(n\).

For large jets, `distinguished_class_restriction` is the narrower exact
mode.  It checks the boundary exponent of the declared class at every
requested order and certifies that the same representative restricts along
the jet tower.  The related `distinguished_class_colon` mode can instead
search for a lift through one declared boundary-element colon.  Neither mode
reports the full finite-jet local-cohomology module.

All polynomial and module operations are exact Singular computations.
Certificates contain a SHA-256 hash of the canonical input, the full input
presentation, the Singular version, and the precise strategy used.
The backend consumes the complete module returned by `sat(N,I)`.  The fast
regression includes a two-generator saturation specifically to prevent
accidental indexing of that returned module as though it were a
`(saturation, exponent)` pair.

There are two associated-prime modes.  `decompose` runs module primary
decomposition (`modDec`, or `primdecSY` in rank one) and lists the primes.
`regularity` does not claim such a list: it certifies only that no associated
prime contains \(I\), using one exact equality \(N:f=N\).
Likewise, the `regularity` saturation strategy returns a zero
local-cohomology certificate only after finding such an \(f\); it is not a
heuristic replacement for saturation.

Two deliberately partial entry points preserve useful exact information when
the untruncated saturation is too expensive.  A distinguished support witness
certifies a nonzero class killed by \(I^e\); this proves that \(H_I^0\) is
nonzero, that some associated prime contains \(I\), and that no element of
\(I\) is regular, without claiming the complete module or prime list.
Finite-jet-only compilation computes either the full requested jet
saturations or only a distinguished-class restriction while recording that
the untruncated base was not computed.

For rank-one codimension-two perfect ideals, `perfect_height` supplies the
standard faster certificate: a length-two Hilbert--Burch resolution and a
one-step dimension drop after adjoining \(f\) prove that \(f\) is regular.
The emitted certificate records this strategy separately from a direct
colon or saturation computation.

## Repository cases

Run all current adapters with:

```bash
.venv/bin/python scripts/compile_support_saturation_cases.py
```

or select `--case cubic`, `--case degree42`, or `--case plane`.

The adapters deliberately preserve three different mathematical scopes.

### Homogeneous cubic cotangent atlas

The cubic adapter constructs the twelve-generator ambient free module for
each of the ten homogeneous ternary-cubic orbit representatives already
used by
[`verify_cubic_symbol_double_saturation.py`](../scripts/verify_cubic_symbol_double_saturation.py).
For every row it computes the saturation, associated primes, and a regular
boundary element.  This is the complete homogeneous-symbol atlas.  It does
not prove the open statement for arbitrary nonhomogeneous higher lifts of a
universal cubic normalization.

### Degree-42 Ritt synchronization

The degree-42 adapter uses the full unit-pivot-reduced core on the base fiber

\[
(e_1,e_2,t)=(1,2,3)
\]

Over characteristic zero it certifies the nonzero untruncated class

\[
h=v w_0^2(2u-5v),\qquad (w_0,w_2)^2h=0.
\]

Thus \(H^0_{(w_0,w_2)}\) is nonzero, a boundary-containing associated prime
exists, and no boundary element is regular.  The full untruncated saturation
module and associated-prime list remain uncomputed.

At the established good prime \(32003\), the same adapter computes the full
support saturations at normal orders six and seven for
\(\mathfrak m=(u,v)\).  Both local-cohomology modules are nonzero, with
detected annihilation exponents one and two.  Their transition
\(H^0_7\to H^0_6\) is surjective.  The stored order-six representative `C6`
is not itself boundary torsion at order seven, but it has a different
order-seven torsion lift.  This corrects the weaker single-colon observation
\(w_0C6=0\): annihilation by one element of the two-generated boundary does
not alone imply membership in \(H^0_{(w_0,w_2)}\).

The finite-jet computation does not lift itself to characteristic zero and
does not imply generic or completed all-order synchronization.  Passing
`--degree42-characteristic 0` requests the much longer rational finite-jet
calculation.

### Plane-JC boundary layer

The plane adapter uses the normalized Poisson-square coefficient ideal and
the cyclic module

\[
R/(I_0:d_3),
\]

which is the source of the first multiplication-kernel layer in the
certified \(d_3,d_2\) filtration.  It computes its support local cohomology
and associated primes.

This cyclic layer is a currently defined boundary-residue proxy.  It is not
the Case-1 conductor/residue matching cokernel proposed in
[`OPEN_PROBLEMS_FOR_MAP_EXTENSIONS.md`](../OPEN_PROBLEMS_FOR_MAP_EXTENSIONS.md).
That module and its distinguished residue class cannot be compiled until
the omitted lower Newton bands are recovered or a truncation-independence
lemma is proved.

## Fast regression

The quick checker exercises a module with closed-point torsion, a
regular-element certificate, associated-prime decomposition, a
distinguished class, and three finite jets:

```bash
.venv/bin/python scripts/verify_support_saturation_compiler.py
```

Generated case artifacts live in
[`artifacts/generated-results/`](../artifacts/generated-results/).  They are
refreshed only by the explicit compilation command above.
