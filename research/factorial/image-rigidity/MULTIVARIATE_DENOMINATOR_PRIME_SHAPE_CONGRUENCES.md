# Projective shape congruences at its denominator primes

The [joint native shape test](MULTIVARIATE_JOINT_NATIVE_SHAPE_RECONSTRUCTION.md) uses a common quadratic's rational shape and an exact divisor bound for its denominator. The denominator primes themselves can also supply congruences on the same scaled integer. Their contribution comes from residues of the two candidate roots in different Newton slope faces. The test is for the original factorial parents, with arbitrary row counts, gaps, and heights. It applies with other certificates that force the common factor to be quadratic.

## Congruences from split roots of distinct valuations

Let a possible rational common quadratic be \(H=A z^2+Bz+C\), with roots \(\alpha,\beta\). Put
\[
                \xi=\frac{B^2}{AC}
                    =\frac{\alpha}{\beta}+2+\frac{\beta}{\alpha}.
\]
The even ordinary moment has no real roots, so \(0\le\xi<4\). Take an exact denominator budget \(D_0\), meaning that the reduced denominator of \(\xi\) divides \(D_0\), and put \(T=D_0\xi\in\mathbb Z\). Fix a prime \(p\) and set \(d=v_p(D_0)\), \(D_{\rm unit}=D_0/p^d\).

At \(p\), the possible root valuations of \(H\) are drawn from the shared Newton slopes, with multiplicity. A pair of distinct valuations \(s<t\) forces the quadratic to split over \(\mathbb Q_p\), since conjugate roots of an irreducible local polynomial have equal valuations. Thus \(s,t\) are integers. Set \(w=t-s\), and write \(\alpha=p^s u,\ \beta=p^t v\) with \(u,v\in\mathbb Z_p^\times\). Unequal valuations in the trace show that \(v_p(\xi)=-w\); in particular \(w\le d\). The exact identity
\[
                    p^w\xi=u/v+2p^w+p^{2w}v/u                 \tag{1}
\]
gives the necessary congruence
\[
    T\equiv p^{d-w}D_{\rm unit}
                   \bigl(u/v+2p^w+p^{2w}v/u\bigr)
                   \pmod {p^{d-w+L}}                          \tag{2}
\]
once \(u,v\bmod p^L\) are specified by the permitted local branches. Any chosen positive precision \(L\) is allowed. The first digit has exact valuation \(d-w\), so the congruence accounts for the same spread as the denominator budget.

For the equal-valuation allocation, \(v_p(\xi)\ge0\) and \(T\equiv0\bmod p^d\). If both root valuations are nonintegral half-integers, the trace has integer valuation at least the ceiling of the root valuation, so \(v_p(\xi)\ge1\) and \(T\equiv0\bmod p^{d+1}\). These lower bounds include a vanishing trace.

## Computing the possible units from the original parents

For an integer root valuation \(s\), substitute \(z=p^s x\) in each parent \(P_j\) and divide its coefficients by \(p^{\lambda_j}\), where \(\lambda_j=\min_k(v_p([z^k]P_j)+ks)\). The resulting polynomial is integral over \(\mathbb Z_p\). Any common root of valuation \(s\) yields a unit root \(u\) of both scaled polynomials. Hence the possible residues of \(u\bmod p^L\) lie among their joint unit zeros modulo \(p^L\). This set can be computed from exact truncated coefficients and does not require Hensel simplicity. For \(L=1\) it is the intersection of the unit roots of the two parents' stepped residual polynomials. Any simple common unit disc or other exact local branch restriction can further reduce the set.

For a permitted allocation \(s<t\), apply (2) to all pairs of the resulting unit residues from those two scaled branches. Take the union over the permitted allocations, adding the equal-valuation classes when allowed. This is a necessary residue set for the single integer \(T\). It is valid when the two unit sets have been computed independently. Any further correlation between the branches can only reduce the set. The available [truncated Newton face comparisons](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md) and the [gcd-free Smith caps](MULTIVARIATE_GCD_FREE_TRUNCATED_SMITH_SIEVE.md) give independent restrictions on these allocations.

At a good native modulus \(p^e\) at which \(\xi\equiv r_p\bmod p^e\) is already known, retain \(D_0\) as is and record \(T\equiv D_0 r_p\bmod p^{d+e}\). This is equivalent to stripping p from the denominator budget at good primes, and lets every prime condition refer to the same integer without changing scale. Apply all the local constraints to \(T\in D_0 I\), where \(I\) is a certified real interval for the candidate shape, or use \(I=[0,4)\).

In particular, suppose at each \(p\mid D_0\) the permissible allocations give a single residue modulo \(p^d\). It is enough that a unique split allocation gives fixed units to precision \(L\ge w\), or that all allocations have equal valuations. If \(D_0\) is squarefree, a forced pair of distinct simple residual roots needs only first residue precision. With an additional coprime native modulus \(M\) giving one shape class, the modulus on \(T\) is at least \(D_0M\). Thus \(D_0|I|\) is compared to \(D_0M\), and uniqueness follows from \(|I|<M\), independently of the size of \(D_0\). An empty intersection excludes the common factor. With several allowed residues at a prime one can retain the separate compatible classes instead of treating the prime as uninformative. The scope condition is the availability of sufficiently narrow local residue sets at the actual denominator primes.

## Exact checks

The [checker](verify_multivariate_denominator_prime_shape_congruences.py) applies the unit branch construction to the ordinary orders-five-and-six profile of the preceding shape test. The two distinct shared 2-adic valuations are \(-1,1\), and each permitted scaled unit set modulo two is \(\{1\}\). With \(D_0=4\), all unequal-slope candidates have \(T\) odd and equal-slope candidates have \(T\equiv0\bmod4\). The good modulus 7 alone permits \(T=2,9\) in \(0\le T<16\), and intersecting the 2-adic permitted classes leaves \(T=9\). The fixture verifies the new use of a denominator prime, while the separate good modulus 31 already provides another exclusion route as in the preceding note. The checker also verifies (2) at depths one through three for exact rational conjugate-shape data in several split quadratic local configurations.
