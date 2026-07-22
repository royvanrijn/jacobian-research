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
5. Before Gröbner elimination, translate the candidate into a complete
   resolved boundary/proximity certificate and run the chart-aware
   [boundary-lattice prefilter](BOUNDARY_LATTICE_PREFILTER.md).  Use every
   boundary prime, not only dicritical components; distinguish primitive
   divisor classes from pullback multiplicities; and declare whether the
   chart is `A^2` or a Laurent chart.  Then apply semigroup and ramification
   tests.  A failed conceptual filter is cheaper and usually applies to a
   larger family.
6. Any localization in a coefficient solve must ship with a complementary-
   strata audit.  A basis containing `1` is not an adequate artifact without
   the input generators, field, order, saturation factors, and an explicit
   identity or independently checkable resultant chain.
7. Rank future work by the table in [NEXT_DEGREE_FRONTIER.md](NEXT_DEGREE_FRONTIER.md).
   The pair \((75,125)\) is the first numerical maximum, but multiple chains at
   \((84,126)\), \((90,135)\), and \((96,144)\) may offer more reusable
   structural tests.

The frontier is a lower bound, not an attainability prediction.  Search
documentation must keep JC(2) separate from the repository's JC(3)
construction program.
