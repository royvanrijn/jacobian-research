# Exact verification for the primitive balanced-degree spectrum

The two scripts use only the Python 3.10 or newer standard library. Each
script contains its exact inputs, namely the manuscript's lattice bounds
and face tables or the four-variable formula for the witness family.

From this directory run

~~~sh
python3 -I verify_faces.py
python3 -I verify_witnesses.py
~~~

The isolated mode removes Python-path and user-site settings. Both scripts
reject optimized mode, which disables their assertions.

| Replay | Arithmetic | Domain and expected result |
| --- | --- | --- |
| Lattice census and cyclic covers | Integers and exact fractions | No three-node crossing lines in degrees 1 or 2, two in degree 3, and fourteen in degree 5. The replay checks every allowed occupation of at least three nodes meeting both signs, then every root multiplicity partition of the associated polynomial. It checks 3 cubic and 73 quintic support/partition cases. A support/partition pair may be counted even when intermediate coefficient conditions forbid that partition. |
| Original witnesses | Sparse integer polynomials in \((\xi_1,\xi_2,z_1,z_2)\) | For six pairs \((a,b)\) of witness parameters, check 22 pure and 22 marked contractions, homogeneity, evaluation on \(R=0\), and a simple-root specialization. The list of pairs and order bounds is in the script. |
| Laurent witnesses | Sparse rational Laurent polynomials | In the balanced septic example, check both traces in the Laurent chart for every exponent from 1 through 8. |

The face checker deduces admissible lattice points from the degree and
enumerates all lines through opposite-sign pairs. The branch check computes
normalized eigenspace dimensions from the branch multiplicities and compares
the total with the genus from Riemann–Hurwitz. It checks the conductor
factors, the rank-two signatures, and the overlarge-multiplicity cases
handled in the written collision argument. The witness checker computes
contractions by literal monomial differentiation and compares them with the
marked factorial formula.

The negative controls remove an admissible line from the face table and add
\(R^7\) to the septic witness. The resulting census is rejected, and the
corrupted witness has first contraction \(8!\). A normal replay on CPython
3.12.3 completed with all stated counts and both controls. This was an
internal execution. The manuscript proves the unrestricted statements,
including the endpoint identity for all orders, the reduction to the finite
support census, and the arithmetic and geometric obstructions. The
[dependency inventory](DEPENDENCIES.md) identifies the proof inputs beyond
the exact computations.
