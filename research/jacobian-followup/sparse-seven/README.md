# Seven term Keller collision exclusion

Every three-variable Keller collision in the frame
`F(0)=F(e1)=0`, `JF(0)=I`, with raw coordinate degrees at most seven,
requires **at least eight nonlinear monomial occurrences**, counted with
coordinate multiplicity.

The complete seven-term search leaves 1,876 determinant-balanced support
orbits. Every exact-support coefficient ideal is the unit ideal over Q in
both Singular and SymPy. An independent implementation reconstructs the
support cover and independently compiles every coefficient ideal.
The six-term regression contains another 450 orbits.

- [Proof and scope](PROOF.md)
- [Complete search manifest](results/manifest.json)
- [Replay validation](validation.json)
- [Independent coverage](coverage-validation.json)
- [Original enumerator comparison and checkpoint regression](regression.json)

From the repository root:

```sh
.venv/bin/python research/jacobian-followup/sparse-seven/verify.py --replay-algebra
```

This closes seven-term attainment in the stated normalized degree range.
The [eight-term successor](../sparse-eight/README.md) also excludes eight
terms and raises the current normalized lower bound to nine. Unrestricted
complexity minima remain open.

