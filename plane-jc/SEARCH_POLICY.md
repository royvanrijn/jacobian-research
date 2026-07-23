# Plane Keller search policy

This policy applies only to searches for two-dimensional Keller maps.

1. Do not spend coefficient-solving time on arbitrary counterexample searches
   with \(\max(\deg P,\deg Q)<125\).  The validated published reduction plus
   the locally replayed \((72,108)\) certificate excludes that range.
2. Low-degree work remains useful for automorphisms, regression examples,
   implementation tests, and validation of individual Newton/valuation
   lemmas.  Such runs must not be advertised as counterexample searches.
3. Every proposed plane counterexample must first be reduced to a standard
   \((m,n)\)-pair and pass the possible-corner/admissible-complete-chain filter.
4. Record orientation, first corner, all successor/final corners, gcd, degree
   ratio, and every vertex nonvanishing condition before constructing a
   coefficient ideal.
5. Before Gröbner elimination, supply exhaustive positive local branch
   scales (including nested resonance charts), compile them with the
   [log-boundary compiler](LOG_BOUNDARY_COMPILER.md), and run the chart-aware
   [boundary-lattice prefilter](BOUNDARY_LATTICE_PREFILTER.md).  Use every
   boundary prime, not only dicritical components; distinguish primitive
   divisor classes from pullback multiplicities; and declare whether the
   chart is `A^2` or a Laurent chart.  Inspect the compiled semigroup
   conductors, differents, and residue degrees.  On a complete `A^2`
   resolution, also run the
   [intrinsic adjunction/Noether gate](INTRINSIC_A2_BOUNDARY_GATE.md), supply
   the global target pole vector, and require effective ordinary and log
   ramification plus an intrinsic dicritical.  A nonproper candidate must
   have canonical free depth at least three.  Corners alone are not a
   compiler input.
6. Treat the source boundary tree as only the first layer of the obstruction.
   A closure claim must also resolve the target nonproper curve, give the
   harmonic map of source and target dual graphs with normal/residue degrees,
   and record the three-section linear series on every noncontracted source
   component.  Bare partitions of a field-degree remainder are not boundary
   packages.
7. If a Poisson-square truncation has the three-layer support
   `P=X^3A+X^2B`, `Q=X^2C+XD` with degree bounds `(3,4;2,3)`, replace its
   reduced coefficient system by the three classified components from
   [POISSON_SQUARE_RIGIDITY.md](POISSON_SQUARE_RIGIDITY.md): the forced
   tangent-pencil closure, `C=0`, and `A=0`.  Retain nonreduced structure
   explicitly when the downstream calculation is scheme-theoretic.  The
   exact associated-prime filtration has three minimal components, the
   three intersection branches, and two deeper core/intersection branches.
   Run proposed lower equations through
   [`cas/poisson_square_filtered_modules.py`](cas/poisson_square_filtered_modules.py)
   before constructing a new global coefficient ideal.  It reports whether
   each dense associated chart is preserved, cut, or eliminated and carries
   the certified `d3,d2` transverse Hilbert vectors needed by later
   scheme-theoretic checks.
8. Before a residue coefficient elimination, run the general polynomial
   right-component remainder sieve for every cover degree dividing the gcd
   of its coordinate degrees.  Do not replace a general right component by
   a parity or fixed-critical-point ansatz.  The sieve requires the complete
   residue coefficient vector: if an archived Newton calculation omits lower
   bands contributing to the residue, derive them or prove a truncation lemma
   before running it; never fill them with guessed zeros.
9. Any localization in a coefficient solve must ship with a complementary-
   strata audit.  A basis containing `1` is not an adequate artifact without
   the input generators, field, order, saturation factors, and an explicit
   identity or independently checkable resultant chain.
10. Rank future work by the tables in
   [FRONTIER_CLOSING_ATTACKS.md](FRONTIER_CLOSING_ATTACKS.md) and
   [NEXT_DEGREE_FRONTIER.md](NEXT_DEGREE_FRONTIER.md).
   The pair \((75,125)\) is the first numerical maximum, but multiple chains at
   \((84,126)\), \((90,135)\), and \((96,144)\) may offer more reusable
   structural tests.

The frontier is a lower bound, not an attainability prediction.  Search
documentation must keep JC(2) separate from the repository's JC(3)
construction program.
