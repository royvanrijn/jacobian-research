# MW3 node

Current best chain:

    MW17
      -> q=90, ab=9*10 : MW7
      -> q=15, ab=3*5 : MW4
      -> q=4,  ab=2*2 : MW3

Root invariants:

    rank = 14
    signed roots = 120
    determinant = 132

Expected ADE:

    A10 + A2 + A1^2

Run:

    sage elkies-k3/scripts/analyze_mw3_node.sage

Then one final tiny-neighbor pass:

    sage elkies-k3/scripts/search_mw3_fourth_neighbors.sage       --qs 2,3,4,6,8,9 --enum-cap 6000 --report 20

Targets:

    root rank 15 -> MW2
    root rank 16 -> MW1

If this does not beat MW3, switch to explicit reconstruction of this node.
