# MW4 node

Chain found computationally:

    rootless MW17
      -> q=90, (a,b)=(9,10): MW7
      -> q=4,  (a,b)=(2,2): MW4

The exact vectors for these arrows, and the subsequent `q=4` arrow to the
`E6/MW3` frame, are now pinned and verified in
[`E6_NEIGHBOR_CHAIN.md`](E6_NEIGHBOR_CHAIN.md).  The executable certificate is
[`scripts/verify_e6_neighbor_chain.sage`](scripts/verify_e6_neighbor_chain.sage).

Selected MW4 frame has root invariants:

    rank 13
    60 signed roots
    determinant 360

which identifies ADE = A5 + A4 + A2 + A1^2.

Run:

    sage elkies-k3/scripts/analyze_q90_mw4_node.sage

Then search only tiny neighbors:

    sage elkies-k3/scripts/search_mw4_third_neighbors.sage       --qmin 2 --qmax 12 --enum-cap 5000 --report 20

The target is root rank 14-16, i.e. MW rank 3-1.
