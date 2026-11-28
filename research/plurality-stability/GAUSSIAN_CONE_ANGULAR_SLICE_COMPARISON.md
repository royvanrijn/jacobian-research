# Angular slice comparison for Gaussian conical partitions

The first-order loss for a sequence converging to a maximizer can
also be proved without a slice mass condition by removing the
spectator radius.
[First-order comparison for unrestricted cones](GAUSSIAN_CONE_DECONING_ASYMPTOTICS.md)

Split the ambient coordinates of a cone into a fixed block of
dimension \(p\) and a spectator block of dimension \(n\). If the
angular slices of the spectator block all have the prescribed
masses, dependence on the spectator direction can improve stability
over the best radial profile by at most an exponentially small
amount as \(n\to\infty\). If the cone is antipodally even in the
spectator coordinates, the comparison is exact for each \(n\).

The comparison also bounds the advantage from arbitrary spectator
directions in terms of their slice mass variation. At a nonsimplex
planar maximizing profile, that variation has to be of order at
least \(1/n\) to overcome the first-order radial loss.

## Conditional angular comparison

Use the Gaussian pairs \(X,Y\) in \(\mathbb R^p\) and \(Z,Z'\)
in \(\mathbb R^n\), independent of one another and both at correlation
\(0<r<1\). Write \(L=|Z|/\sqrt n\), \(L'=|Z'|/\sqrt n\),
\(\omega=Z/|Z|\), \(\omega'=Z'/|Z'|\), and \(t=\omega\cdot\omega'\).
Represent an arbitrary measurable conical partition by its tuple of
profiles \(h_\omega(X/L)\). Let \(\mu_n\) be the law of \(X/L\),
\(\sigma\) uniform on \(S^{n-1}\), and \(m(\omega)=\int h_\omega\,d\mu_n\).
Since \(L\) and \(\omega\) are independent, the ambient masses are
\(m=\int m(\omega)\,d\sigma(\omega)\). Write \(B_n(h)\) for the
stability of the radial lift \(h(X/L)\), allowing fractional
simplex-valued tuples in this quadratic form. Write
\(\bar h=\int h_\omega\,d\sigma(\omega)\).

For any \(r<\theta<1\) there exist \(c,C>0\) such that, for all
sufficiently large \(n\),
\[
 S_r(F)\le \theta\int B_n(h_\omega)\,d\sigma(\omega)
              +(1-\theta)B_n(\bar h)+Ce^{-cn}.                \tag{AS1}
\]
If \(F(x,z)=F(x,-z)\) then any \(r^2<\theta<1\) works.
The constants can be uniform for \(r\) in a compact positive
correlation interval with a uniform positive margin in the
condition on \(\theta\).

To prove (AS1), condition on \(t\). The angular pair is a uniform
fixed-angle pair, independent of the radii under this conditioning.
The radial density is proportional to
\[
 (\ell\ell')^{n-1}
 \exp\left[-\frac{n}{2(1-r^2)}
          \bigl(\ell^2+(\ell')^2-2rt\ell\ell'\bigr)\right]\,d\ell\,d\ell'.
\]
Let \(K_t\) be the resulting bilinear form on two bounded profiles in
the \(X/L,Y/L'\) coordinates. For \(t\ge0\) it is positive
semidefinite. Indeed expand the cross exponential into powers of
\(rt\ell\ell'\). Each term is the quadratic form of \(T_r\) in the
\(X\) coordinates applied to the profile integrated against the
radial weight \(\ell^{n-1+k}e^{-n\ell^2/(2(1-r^2))}\) for a
nonnegative integer \(k\). The coefficients are nonnegative.
For bounded profiles the series is absolutely integrable by
Gaussian damping.

The fixed-angle spherical averaging operator \(R_t\) on \(S^{n-1}\)
acts on degree-\(\ell\) spherical harmonics by
\[
       \mathbb E\bigl(t+i\sqrt{1-t^2}\,U\bigr)^\ell,
\]
where \(U\) is a coordinate of a uniform vector in \(S^{n-2}\).
This follows from the rotation-invariant harmonic decomposition by
testing the harmonic \((x_1+i x_2)^\ell\) at \(e_1\).
The absolute eigenvalue in degree one is \(|t|\). In degree at
least two, the complex modulus under the expectation is at most one,
so each absolute eigenvalue is at most
\(t^2+(1-t^2)/(n-1)\). For an even profile only even degrees occur.
[Hwang, Neeman, Parekh, Thompson and Wright, Section 4 for the spherical harmonic decomposition](https://hwangy.github.io/files/vector_borell.pdf)

Apply this bound to the angular function \(\omega\mapsto h_\omega\)
in the Hilbert space of the \(K_t\) form. It can equivalently be
applied to each of the \(T_{\sqrt r}\)-smoothed radial integrals
in the series defining \(K_t\). If the mean-zero norm of \(R_t\)
is at most \(\theta\), the conditional agreement is at most
\[
 \theta\int K_t(h_\omega,h_\omega)\,d\sigma(\omega)
               +(1-\theta)K_t(\bar h,\bar h).
\]
The event \(|t-r|>\eta\) has exponentially small probability for
every fixed \(\eta>0\), uniformly on compact positive correlation
intervals. Indeed the two empirical squared norms and the empirical
cross product of the coordinate pairs have exponential concentration
about \(1,1,r\). This follows by their moment generating functions
near zero and the Chernoff inequality. For the cross product one
can use the difference of two independent scaled squares of
standard normals. On \(|t-r|\le\eta\), choose \(\eta\) small enough
and \(n\) large enough that \(t>0\) and the stated norm bound
is at most \(\theta\). All agreements in the displayed inequality
take values in \([0,1]\) even on the exceptional event. Averaging
gives (AS1).

There is also an exact unweighted comparison for even spectator
profiles:
\[
                      S_r(F)\le\int B_n(h_\omega)\,d\sigma(\omega).
                                                               \tag{AS2}
\]
For \(n\ge2\), condition on \(|t|\). The angular averages at \(t\)
and \(-t\) are identified by evenness. After the two radial
densities are combined using their conditional weights, the cross
exponential becomes a hyperbolic cosine with nonnegative
coefficients in even radial powers. Thus the combined \(K_{|t|}\)
is positive semidefinite. Apply the inequality that its mixed
pairing is at most the mean of the two diagonal pairings, then
average. In fact equality requires that \(h_\omega\) be independent
of \(\omega\) almost everywhere. For \(0<|t|<1\), the combined form
is strictly positive definite on bounded profiles. A zero form
would give zero \(X\)-integrals in the \(T_r\) form for every even
radial moment. Strict positive definiteness of \(T_r\) makes these
moment functions zero almost everywhere in \(X\). For almost every
\(X\), the radial function then vanishes by completeness of
polynomials in \(\ell^2\) for the Gaussian radial weight, equivalently
by uniqueness of the Laplace transform after changing to \(\ell^2\).
Fubini and positive Gaussian densities imply that the profile
vanishes. The conditional mixed-pairing equality thus identifies
the profiles of almost every angular pair under the joint angular
law, whose density is positive. For \(n=1\), evenness already
implies independence from \(\omega\).

## Consequences for the optimum and a nonsimplex tie

Write \(A_{p,n}(r,m)\) for the radial profile optimum from (AX3).
The average profile \(\bar h\) is feasible in its relaxed problem,
which has the same optimum by strict positive definiteness and
rounding. Each \(h_\omega\) can be mass-corrected in \(\mu_n\) at
stability cost at most \(\|m(\omega)-m\|_1\). Apply these bounds
in (AS1) to get
\[
 S_r(F)\le A_{p,n}(r,m)+\theta\int\|m(\omega)-m\|_1\,d\sigma(\omega)
                                   +Ce^{-cn}.                 \tag{AS3}
\]
In particular when \(m(\omega)=m\) almost everywhere, the maximum
over all conical partitions with that condition differs from
\(A_{p,n}\) by at most \(Ce^{-cn}\). The entire first-order expansion
(AX4) holds for this class. For even spectator profiles satisfying
the same condition the optimum is exactly \(A_{p,n}\), and all
optimizers are radial by (AS2).
[Radial profile optimum and expansion](GAUSSIAN_CONE_AXIAL_INSTABILITY.md#exact-leading-loss-for-the-entire-radial-profile-class)

For the balanced three-cell case take \(p=2\), let
\(V=V_2(r,m)\), and let \(A\) be a planar global maximizer, possibly
one that ties the simplex. Write its exact coefficient from (CH2)
as
\[
 C(A)=\tfrac12\bigl(\|b_A\|_{a_A}^2
                           -r^2\langle b_A,\mathcal K_A b_A\rangle\bigr).
\]
Suppose a sequence of cones \(F_n\) converges in Gaussian label
distance to the cylinder of \(A\) over this planar block. For every
\(\theta>r\) with \(\theta<1\),
\[
 n\bigl(V-S_r(F_n)\bigr)
 \ \ge\ C(A)-\theta n\int\|m_n(\omega)-m\|_1\,d\sigma(\omega)+o(1).
                                                               \tag{AS4}
\]
In particular if \(A\) is nonsimplex and \(V-S_r(F_n)=o(1/n)\),
then the limit inferior of the scaled mean absolute slice mass
variation is at least \(C(A)/r>0\). An analogous bound with
\(r^2\) in place of \(r\) holds for even spectator profiles.

For the proof, the local version of the upper bound in (AX4)
says that a balanced ratio profile \(u_n\) converging strongly to
\(A\) has \(B_n(u_n)\le V-C(A)/n+o(n^{-1})\). Indeed the price
repair in that proof only uses the given limit \(A\) and the
uniform law expansion (AX2). This applies also to fractional
simplex-valued profiles: strong convergence leaves reservoirs
with source coordinate at least \(1/2\), and one can transfer
fractional masses within the same bounded patches. The transfer
still changes the quadratic terms by \(O(n^{-2})\) and uses the
same limiting prices. Thus the local upper bound is uniform for
profiles within any shrinking Gaussian \(L^1\) neighborhood of
\(A\), by applying it to a sequence violating the bound.

The mean profile distance of \(h_{n,\omega}\) from \(A\) tends to
zero in \(\mu_n\). In ambient coordinates this follows from the
assumed distance to \(A(X)\) and
\(\mathbb E|1_A(X/L)-1_A(X)|\to0\). The latter is continuity of
Gaussian dilation in measure. The mean also tends to zero in
\(\gamma_2\) by (AX2). Hence \(\bar h_n\) converges strongly to
\(A\) and is balanced for \(\mu_n\). For a set of angles with
measure tending to one, the individual profiles lie in a shrinking
neighborhood of \(A\) and have mass deviations tending uniformly
to zero. Their exact mass repairs still lie in a shrinking
neighborhood, with stability cost at most \(\|m_n(\omega)-m\|_1\).
Apply the local upper bound to these repaired profiles, and apply
the global upper bound \(V\) to repaired profiles at the remaining
angles. This yields an upper bound of
\(V-C(A)/n+o(n^{-1})+\int\|m_n(\omega)-m\|_1\,d\sigma\)
for their unrepaired mean stability. The same bound without the
last term holds for \(\bar h_n\). Formula (AS1) proves (AS4).
[Coefficient of the translated lift](THREE_CELL_CONICAL_HOMOGENIZATION.md#stability-loss-of-the-balanced-lift-of-a-maximizer)

Using dimension-free near-maximizer compactness, every sequence
of balanced three-cell cones with stability tending to \(V\)
admits a subsequence converging to a planar maximizer after a
rotation. The planar block in (AS4) can then be the ambient
centroid plane. Slice mass variation at the scale \(1/n\)
is not controlled by this qualitative compactness.
[Dimension-free near-maximizer compactness](GAUSSIAN_NEAR_MAXIMIZER_COMPACTNESS.md)
