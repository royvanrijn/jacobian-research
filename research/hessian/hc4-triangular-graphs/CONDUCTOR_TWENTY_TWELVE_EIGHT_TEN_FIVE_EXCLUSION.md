# The three-stage exclusion for sequence (12,8,10,5)

**Theorem ConductorTwentyTwelveEightTenFiveExclusion.** Let
\((P,Q)\) be a complex polynomial plane Keller map in a source
frame with \(P\) monic of positive longitudinal degree. Its generic
directional coordinate order cannot have reduced characteristic
sequence \((12,8,10,5)\).

## Two approximate-root equations

Over \(K=\mathbb C(y)\), take monic reduced generators \(W,V\)
of longitudinal degrees eight and twelve. The one-place approximate
roots \(S,T\) have longitudinal degrees ten and five and are monic
of degrees three and six in \(W\), respectively. This uses the
same coefficient bounds and normalization as the
[three-stage sequence calculation](CONDUCTOR_FOURTEEN_THREE_STAGE_EMBEDDING.md#the-two-approximate-root-identities).
The primary input is Assi–García-Sánchez, Propositions 2 and 4.
[Retained source](sources/assi-garcia-sanchez-1407.0490v1.pdf)

In the relation for \(S\), the three exponents zero, one and two
in \(W\) give distinct degree residues modulo twelve. Thus terms
of weight above twenty-four cannot cancel, and the highest terms
are \(W^3-V^2\). The coefficient of \(V\) is affine in \(W\).
After completing the square in \(V\) and translating \(W\) to
depress its cube, the first relation is
\[
                   V^2=W^3+kW+C+\eta S,\qquad \eta\ne0,    \tag{1}
\]
with suitable coefficients over \(K\). These target changes
preserve source monicity and the original degrees.

Before making the target changes, divide the monic degree-six
polynomial \(T\) in \(W\) by powers of the monic degree-three
\(S\), with coefficients over \(K[V]\). This writes
\[
                 T=S^2+\sum_{i<3,\ j<2} b_{ij}(V)W^iS^j.
\]
The six weights \(8i+10j\) have distinct residues modulo twelve,
so no terms of weight above twenty can cancel. The only terms
in the sum of weight at most twenty are scalar multiples of
\(1,V,W,VW,W^2,S,WS\). The unique second term of weight twenty
is \(VW\), with coefficient \(-1\) by monicity of \(S^2\).
The target changes in (1) translate \(W\) and change \(V\)
by an affine polynomial in \(W\). They preserve the displayed
span, including the coefficient of \(VW\). Consequently \(T\)
has the form \(S^2+\ell(W)S-VW+\text{lower terms}\), with
\(\ell\) affine in \(W\). Complete the square in \(S\) by replacing
it by \(S+\ell(W)/2\). This changes the right side of (1) only
by an affine polynomial in \(W\), absorbed in \(k,C\). The
degree-ten monic leading term of \(S\) is preserved, and the
correction \(\ell(W)^2/4\) has weight at most sixteen.
Renaming \(S\), the two identities therefore coexist with
\[
                 S^2=VW+pW^2+qV+rW+s+\nu T,\qquad \nu\ne0 .
                                                                  \tag{2}
\]

## Laurent parameter comparison

If the target change to \((W,V)\) has determinant \(\lambda\),
take \(x=\rho z\), \(\rho^{20}=\lambda\), and scale
\(W,V,S\) by \(\rho^{-8},\rho^{-12},\rho^{-10}\), respectively.
The same identities hold over the finite parameter extension
with rescaled coefficients. The target determinant is now one.
The
[parameter identity](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#1-common-source-and-target-identities)
gives
\[
                  W_zV_y-W_yV_z=J\rho+\partial_z H(W,V),   \tag{3}
\]
where \(J\ne0\) is the original constant Jacobian, and \(H\)
is a polynomial in the two target variables.

In the monic Laurent branch at \(z=\infty\), equation (1) gives
\[
 V=W^{3/2}+\tfrac12\eta S W^{-3/2}+\tfrac12 kW^{-1/2}
                                                +O(z^{-12}).
\]
Indeed \(CW^{-3/2}\) has degree \(-12\), and the higher
binomial terms have degree at most \(-16\) since \(\deg S=10\).
Putting this expansion in (2) gives \(S^2=W^{5/2}+pW^2+O(z^{12})\).
The terms from \(VW-W^{5/2}\) have degree at most six, and among
the other lower terms of (2) the largest possible degree is
twelve. Hence
\[
                   S=W^{5/4}+\tfrac p2 W^{3/4}+O(z^2).
\]
Substitution into \(V\) now gives
\[
 V=W^{3/2}+\tfrac12\eta W^{-1/4}+\tfrac12 kW^{-1/2}
                   +\tfrac14\eta p W^{-3/4}+O(z^{-10}).   \tag{4}
\]

The derivatives \(W_y,W_z\) have degrees at most seven.
Differentiation in \(y\) preserves the tail degree bound in
(4), so the tail makes no contribution to the polynomial
Jacobian. The chain rule cancels the terms from differentiating
\(W\) inside each power in (4). Thus the polynomial Jacobian
is the nonnegative Laurent part of
\[
 \tfrac12\eta' W_zW^{-1/4}
      +\tfrac12 k' W_zW^{-1/2}
      +\tfrac14(\eta p)' W_zW^{-3/4}.
\]
The three terms have respective leading degrees five, three
and one, with coefficients \(4\eta'\), \(4k'\) and \(2(\eta p)'\).

The pole semigroup is \(\langle8,10,12,5\rangle\). Equation
(3) first bounds the degree of the evaluated \(H\) by six,
in which range the largest pole degree is five. Thus
\(\partial_z H\) has degree at most four, forcing \(\eta'=0\).
The bound now drops to degree four for \(H\), so \(H\) is
constant in \(z\). Comparing the degrees three and one forces
\(k'=(\eta p)'=0\). No constant contribution to the Jacobian
remains. Equation (3) would give \(J\rho=0\), a contradiction.
