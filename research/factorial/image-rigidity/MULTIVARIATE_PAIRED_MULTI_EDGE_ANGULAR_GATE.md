# Angular separation from paired multi-edge Newton polygons

Let \(P_n,P_N\) be the ordinary factorial moment polynomials of
two distinct exponent columns, with \(N=n+1\) and \(n\ge2\).
Fix a prime \(p\). An edge is binomial if only its endpoints
attain its supporting line. Its width is the difference of
their indices. Use the negatives of the slopes as the root
valuations, taking \(v_p(p)=1\).

**Theorem.** Suppose that every edge of both Newton polygons is
binomial. For each polygon form the set of nonzero differences
of the root valuations on its edges. If these two sets are
disjoint and every lower edge width is coprime to every upper
edge width, then the two root argument sets are disjoint.
The conclusion is unchanged by independent positive rescalings.
This criterion allows both polygons to have several edges.

The residue test below can also distinguish quotient types
when the two difference sets meet. A native all-degree family
has two different slopes at each order, with just two shared
rows supplying the valuations at a single prime.

## The width and difference test

Embed both complex splitting fields in a common algebraic
closure of \(\mathbb Q_p\).
A quotient of two roots has valuation equal to the difference
of their root valuations. If the same quotient occurs in both
degrees, the disjoint difference sets force it to be a unit,
so its two roots in each degree lie on a common edge.

At an edge of width \(w\), scale the argument by an element
of the appropriate valuation and reduce the coefficients on
the supporting line. The nonzero roots of the resulting
binomial have quotients whose orders divide \(w\).
The residues of roots from the original edge belong to this
set by the Newton polygon rule. If \(p\nmid w\), distinct
roots from the edge have distinct residues, since the
nonzero residual roots are simple.

Thus a common quotient of distinct roots in both degrees
would reduce to one, because the two widths involved are
coprime. At least one of them is prime to \(p\), where
distinct roots of the edge cannot have residue quotient one.
This is a contradiction. Equal nonreal arguments would give
equal quotients after dividing each root by its complex
conjugate. The even factorial moment is positive at all
real arguments, so the two degrees cannot share a real
argument. This proves the theorem.

## A residue comparison when slope differences coincide

Here is a direct refinement retaining a common fractional
scale. Choose a positive integer \(L_1\) clearing the
denominators of all root valuations and use one
\(\pi\) with \(\pi^{L_1}=p\).
For an edge \(e\) of root valuation \(s_e\) scale the
variable by \(\pi^{L_1s_e}\) and reduce as above.
Write its nonzero residual roots as roots of
\(x^{w_e}=a_e\), where \(a_e\ne0\) in the algebraically
closed residue field. Write \(w_e=p^{b_e}w_e^0\),
\(p\nmid w_e^0\), and denote by \(a_e^0\) the unique
\(p^{b_e}\)-th root of \(a_e\) in the residue field.

For an ordered pair of edges \((e,f)\) in a given degree,
the possible quotients have valuation \(s_e-s_f\).
After dividing a quotient by \(\pi^{L_1(s_e-s_f)}\),
its residue is contained in the coset
\[
       x^\ell=c,\qquad
       \ell=\operatorname{lcm}(w_e^0,w_f^0),\qquad
       c=(a_e^0)^{\ell/w_e^0}(a_f^0)^{-\ell/w_f^0}.            \tag{1}
\]
This follows by taking a ratio of the residual roots.
If \(e=f\) and \(p\nmid w_e\), delete the value one
when considering quotients of distinct roots.

Compare only pairs of edge types from opposite degrees
with the same valuation difference. Cosets
\(x^\ell=c\) and \(x^{\ell'}=c'\) meet before the
optional deletions exactly when
\[
                  c^{\ell'/g}=(c')^{\ell/g},
                            \qquad g=\gcd(\ell,\ell').         \tag{2}
\]
Necessity is immediate. For sufficiency choose integers
\(u,v\) with \(u\ell+v\ell'=g\). The prescribed value
\(x^g=c^u(c')^v\) is compatible with both equations
by (2), and any solution gives a common point.
Since \(p\nmid g\) there are exactly \(g\) common points
when the intersection is nonempty. Thus the optional
deletions also give an exact finite-field test for these
containing cosets. If every compared pair is disjoint
after the deletions, all two-degree quotients of distinct
roots are disjoint and the angular conclusion follows.
The test need not assert that every residual coincidence
lifts to a common quotient.

## A native all-degree two-edge family

Take \(J\ge4\) at least the
[computable absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Fix any retained ordinary rows. Add two increasing shared
rows of gaps \(d_1\ge J-1\) and \(d_2\ge1\). Choose one
prime \(p>J\) larger than every retained factorial argument
in degrees \(m<J\), and larger than
\((J-1)\max(d_1,d_2)\).
For the smaller exponents \(\ell_1,\ell_2\) of the new
rows require
\[
                 \ell_1+1\equiv p\pmod{p^2},\qquad
                 \ell_2+d_2\equiv p^2\pmod {p^3}.             \tag{3}
\]
These requirements can also be imposed on rows that
already occur in an original support. The rows counted
as retained and fixed may have either orientation or be
private. Any further rows whose sampled factorial factors
are units at this prime after their smaller-end scalar has
been removed are harmless.

At each \(2\le m<J\), strip the scalar \((m\ell)!\) from
each row of smaller exponent \(\ell\).
For the first marked row the only factor divisible by
\(p\) in the full ramp \(m\ell_1+h\), \(1\le h\le md_1\),
occurs at \(h=m\). Indeed its residue is \(h-m\), and
the full ramp length is less than \(p\).
The marker has valuation one since \(m<p\), and \(m\le d_1\)
puts it inside the first block. It is therefore sampled
at every positive coefficient index.

For the second row the residues are \(h-md_2\), so the
unique hit is the last factor \(m(\ell_2+d_2)\).
Its valuation is two, and it is sampled only at the
top index \(k=m\). The full ramp length less than \(p\)
in particular excludes any earlier factor with the same
residue. The binomial coefficients are units since \(p>m\).
The normalized valuation array at degree \(m\)
is consequently
\[
                              (0,1,\ldots,1,3).
\]
The two Newton edges have widths \(m-1,1\), and slopes
\(1/(m-1),2\), respectively. The first edge is binomial
because its intermediate ordinates are exactly one
and its slope is \(1/(m-1)\). For \(m=2\) the first
width is one, and the two slopes are \(1,2\).
The nonzero differences in degree \(m\) are
\(\pm(2-1/(m-1))\), which differ from those in degree
\(m+1\). The cross-degree widths are pairwise coprime,
because the only nontrivial pair is \(m-1,m\).
The theorem proves angular separation for every pair
\(n,N=n+1\) with \(3\le N<J\). The first pair and
the absolute cutoff supply all other degrees.

The bases in (3) can increase arbitrarily along their
progressions, so the added curvature can be arbitrarily
small. A further unmarked shared row can prescribe any
desired signed total gap of the addition. Choose the
orientation and gap to make up the difference from
\(d_1+d_2\), take its smaller exponent a positive
multiple of \(p\), and take the prime larger than its
sampled gap lengths too. Its factorial factors have
residues \(1,\ldots,md\), so it contributes no valuations.
For example, with \(d_1=J-1\), \(d_2=1\), the rows
\[
       (p-1,p-1+d_1),\qquad
       (p^2-1,p^2),\qquad
       (p^3+d_1+1,p^3)
\]
give an explicit fully shared family with mixed orientations
and zero signed total gap. Take \(p>(J-1)(d_1+1)\).
Both polygons in the certificate genuinely have two
edges at every relevant pair.
