# q=90 MW-rank-7 node

Selected node from the generic discriminant-948 NS lattice:

- q = 90
- factorization = 9*10
- root rank = 10
- signed roots = 34
- root determinant = 192
- MW rank = 7
- inferred ADE type: A3^2 + A2 + A1^2

The exact 17x17 positive frame Gram is stored in:

    elkies-k3/data/fibrations/q90_mw7_frame.txt

## 1. Analyze the node exactly

    sage elkies-k3/scripts/analyze_q90_mw7_node.sage

This verifies the root lattice, ADE decomposition, Shioda-Tate rank, and
regulator assuming trivial torsion.

## 2. Search second-generation neighbors

Start from the q=90 frame itself, not from the original rootless rank-17 frame:

    sage elkies-k3/scripts/search_q90_second_neighbors.sage \
      --qmin 2 --qmax 16 --enum-cap 4000 --report 20

The target is root rank >= 12, i.e. MW <= 5.

If a good node appears, save that frame and iterate once more rather than going
back to huge-q one-step searches from the original frame.
