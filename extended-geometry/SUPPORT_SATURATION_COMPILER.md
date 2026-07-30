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

## Shared input schema

New calculations should use `SupportSaturationProblem`, whose stable
JSON form is
[`support_saturation_input.schema.json`](../schemas/support_saturation_input.schema.json).
It declares, in one place:

- the coefficient ring and monomial ordering;
- the submodule presentation \(N\subset R^r\);
- the support ideal \(I\);
- the completion ideal \(\mathfrak m\);
- parameter/base variables and normal variables;
- requested jet orders and an optional distinguished class; and
- whether the exact backend result is being used as an exact certificate or
  as modular evidence for a different characteristic.

The variable roles are checked for unknown names and overlap.  The completion
ideal is deliberately separate from the normal-variable list: a filtration
may use polynomial combinations of the normal directions.  Empty role lists
are allowed for presentations in which one of the two roles is absent.

A small exact example is
[`support_saturation_example.json`](../schemas/support_saturation_example.json).
Compile any conforming input with:

```bash
.venv/bin/python scripts/compile_support_saturation.py \
  schemas/support_saturation_example.json \
  --output /tmp/support_saturation_certificate.json
```

Existing Python callers of `SupportSaturationCompiler.compile` remain
supported.  Repository adapters should prefer `compile_problem` when the
untruncated module is available; the narrower witness and finite-jet-only
entry points remain necessary when it is not.

## Exact output

`SupportSaturationCompiler.compile` records:

- a standard basis for \(N:I^\infty\);
- generators, relations, annihilator, annihilator radical, and associated
  primes or explicitly partial associated-prime candidates of
  \((N:I^\infty)/N\);
- associated primes of \(F/N\), when module primary decomposition is
  requested;
- exact colon tests \(N:f=N\) for candidate \(f\in I\);
- the annihilator, its radical, and the least detected support exponent of a
  distinguished class;
- finite normal jets
  \[
  N_n=N+\mathfrak m^nF
  \]
  with the images and surjectivity of the transitions from order \(n+1\) to
  order \(n\); and
- the least common annihilating exponent on the requested finite tower,
  accompanied by an explicit statement that a finite prefix is not an
  all-order uniform bound.

For large jets, `distinguished_class_restriction` is the narrower exact
mode.  It checks the boundary exponent of the declared class at every
requested order and certifies that the same representative restricts along
the jet tower.  The related `distinguished_class_colon` mode can instead
search for a lift through one declared boundary-element colon.  Neither mode
reports the full finite-jet local-cohomology module.

All polynomial and module operations are exact Singular computations.
Certificates contain a SHA-256 hash of the canonical input, the full input
presentation, the Singular version, and the precise strategy used.  The
`certificate_state` block separates this backend fact from theorem scope:
a calculation over \(\mathbf F_p\) is exact over that field, while its
default assurance for a characteristic-zero target is `modular` and records
that no characteristic-zero lift is claimed.
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

or select `--case cubic`, `--case cubic-frontier`, `--case degree42`, or
`--case plane`.

The adapters deliberately preserve three different mathematical scopes.
They share the compiler contract used by the cubic-normalization,
degree-42, plane-boundary, and filtered-quantization programmes; adding a
new programme should require an adapter that constructs the presentation,
not another saturation implementation.

### Homogeneous cubic cotangent atlas

The cubic adapter constructs the twelve-generator ambient free module for
each of the ten homogeneous ternary-cubic orbit representatives already
used by
[`verify_cubic_symbol_double_saturation.py`](../scripts/verify_cubic_symbol_double_saturation.py).
For every row it computes the saturation, associated primes, and a regular
boundary element.  This is the complete homogeneous-symbol atlas.  It does
not prove the open statement for arbitrary nonhomogeneous higher lifts of a
universal cubic normalization.

### Cubic formal-gauge annihilator frontier

The `cubic-frontier` adapter imports the exact gauge-cokernel atlas from
[`cubic_formal_gauge_cokernel_atlas.json`](../artifacts/generated-results/cubic_formal_gauge_cokernel_atlas.json).
It treats the proved smooth-symbol theorem as a routing decision: the full
24-parameter quartic cotangent-saturation problem is closed, so the next
smooth work is global algebraization and boundary/Keller-open
compatibility.

For the singular squarefree symbols it emits the finite queue

\[
\begin{array}{c|c|c}
\text{symbol}&\operatorname{Ann}(Q_h)&\dim (Q_h)_4\\ \hline
\text{nodal}&(x)&2\\
\text{cuspidal}&(x^2)&4\\
\text{line + transverse conic}&(yz)&4\\
\text{line + tangent conic}&(y^3)&6\\
\text{triangle}&(xyz)&6\\
\text{concurrent lines}&(x^3)&8.
\end{array}
\]

These annihilator types, rather than coordinate planes or sparse support in
the original 24-element basis, are the search keys for the next
deformation-dependent cotangent adapters.  For double-line, triple-line,
and zero symbols the adapter records the prior generically-étale/Keller
compatibility gate.  This routing artifact is not itself a singular-symbol
saturation theorem.

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

The order-six statement now has a separate characteristic-zero certificate.
Let

\[
 R=\mathbf Q[u,v,w_0,w_1,w_2],\qquad
 \mathfrak m=(u,v),\qquad J_6=I+\mathfrak m^6,
\]

where \(I\) is the full unit-pivot-reduced core on the displayed base fiber.
For the explicit stored class \(c_6\),

\[
 c_6\notin J_6,\qquad w_0c_6\in J_6,\qquad w_2c_6\in J_6.       \tag{1}
\]

Consequently \([c_6]\) is a nonzero class in
\(H^0_{(w_0,w_2)}(R/J_6)\), killed already by the boundary ideal itself.
This is an exact embedded-support obstruction on the sixth normal jet.

The certificate avoids global saturation.  After reduction modulo
\(\mathfrak m^6\), the canonical input has 16 primitive integer generators
and 1,002 nonzero terms, with problem hash
`b5c115b3f2efb61c54017f251274cf7f537d650e5333a039949b4d5833fc4218`.
The two membership statements in (1) use one \(646\)-by-\(2400\) sparse
Macaulay matrix of rank \(547\).  A two-sided block Krylov--Wiedemann solve
over eight 31-bit primes reconstructs 157 nonzero rational multiplier terms
for \(w_0c_6\) and 88 for \(w_2c_6\).

Nonmembership is positive proof data rather than a failed reduction.  A
second Macaulay system has 1,249 equations on 210 possible functional
coordinates and rank 134.  Modular reconstruction gives a nine-term
finite-support functional \(\lambda\) satisfying

\[
 \lambda(c_6)=1,\qquad \lambda(qf_i)=0
\]

for every reduced-core generator \(f_i\) and every monomial \(q\).  Extending
\(\lambda\) by zero outside its finite support also kills
\(\mathfrak m^6\), so it proves \(c_6\notin J_6\) without a degree-bound
assumption.  The minimal checker uses only Python's standard library and
replays the rational polynomial identities and all multiplier tests that
can meet the functional support.

Run the independent replay with:

```bash
python3 scripts/verify_degree42_c6_macaulay.py
```

Intentional reconstruction is:

```bash
.venv/bin/python scripts/compile_degree42_c6_macaulay.py
```

The older modular finite-jet artifact still supplies the full order-six and
order-seven saturation modules and their transition.  The new rational
certificate closes only the fixed order-six class.  It does not compute the
full characteristic-zero saturation, lift the order-seven module, vary the
base specialization, or imply generic or completed all-order
synchronization.  Passing `--degree42-characteristic 0` to the older adapter
still requests the much longer full rational finite-jet calculation.

### Cellular consequence

The class \(c_6\) changes the next step in the Hessian--Ritt programme.
Vanishing on generic local strata cannot be promoted to support saturation:
the sixth jet already has a boundary-supported class concentrated at the
deepest monomial degeneration.  Its location makes a cellular explanation
more plausible than another undirected global colon calculation.

The next comparison with the
[cellular cotangent model](HESSIAN_RITT_COTANGENT_DESCENT_COMPARISON.md)
is therefore:

1. construct the image, if any, of \([c_6]\) in the cut-\(6\) filtered
   cellular/Postnikov coefficient module;
2. determine whether that image is the cut-\(6\) non-splitting class
   detected by `HRCELL2`--`HRCELL4`;
3. restrict it to the rotated cut-\(14\) and cut-\(21\) sectors using their
   completed splittings; and
4. compute the resulting braid restriction class and decide whether
   coherence kills it or makes it persist.

None of these four identifications is currently proved.  In particular,
the completed splitting of the rotated first-conormal extensions does not
by itself kill a class coming from a different finite-jet local-cohomology
module.  The modular order-seven saturation transition supplies evidence
that a lift can exist, but the characteristic-zero certificate above stops
at order six.  This is the precise new interface between support saturation
and cellular cotangent non-splitting.

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
