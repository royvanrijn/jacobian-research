# q=15 MW4 node

Selected q=15 child of the q=90 MW7 node:

    q=15
    ab=3*5
    root rank=13
    signed roots=68
    root determinant=216
    MW rank=4

Expected ADE:

    A5^2 + A2 + A1

Run:

    sage elkies-k3/scripts/analyze_q15_mw4_node.sage

Then search tiny third-generation neighbors:

    sage elkies-k3/scripts/search_q15_third_neighbors.sage       --qmin 2 --qmax 12 --enum-cap 5000 --report 20

Target:

    root rank 14 -> MW 3
    root rank 15 -> MW 2
    root rank 16 -> MW 1
