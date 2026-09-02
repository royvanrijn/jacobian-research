# MW4 node B

Second MW-rank-4 child of the q=90 MW7 node.

Chain:

    MW17
      -> q=90, ab=9*10: MW7
      -> q=9,  ab=3*3 : MW4

Root invariants:

    rank = 13
    signed roots = 60
    determinant = 576

Expected ADE:

    D5 + A2^2 + A1^4

Run:

    sage elkies-k3/scripts/analyze_q90_mw4b_node.sage

Then:

    sage elkies-k3/scripts/search_mw4b_third_neighbors.sage       --qmin 2 --qmax 12 --enum-cap 5000 --report 20

Compare directly with MW4 node A (q=4 child), looking for MW <= 3.
