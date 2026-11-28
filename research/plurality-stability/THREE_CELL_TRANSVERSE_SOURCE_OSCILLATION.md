# Oscillation of the spectator source along a reference line

Fix a centered regular three-sector simplex \(T\) in the
plane and \(0<r<1\). Put \(\sigma=\sqrt{1-r^2}\) and let
\(p_\sigma\) be the one-dimensional centered Gaussian
density of variance \(\sigma^2\). For the line of the pair
\(i,j\), use coordinates \((t,v)\) along the forward
reference ray and in its normal direction. The spectator
is \(k\). Its weighted transverse mass and source are
\[
 q_k(t)=\int 1_{A_k}(t,v)p_\sigma(v)\,dv,\qquad f_k=1-3q_k .
\]
Assume the unadjusted pair scores are equal on all three
reference lines, as they are for a maximizing compact
simplex discrepancy. Then, on its own pair line,
\[
 h_k(s)=U_i(s,0)-U_k(s,0)
                    =\tfrac12(p_\sigma*f_k)(rs).             \tag{OS1}
\]
Indeed the conditional Gaussian factors in the two
coordinates, and \(U_i=U_j=(1-U_k)/2\) there.

There is a zero count with multiplicities for (OS1). In
particular, the one-sign-change condition on every source
would exclude a nonsimplex right-contact maximizer in the
compact-defect branch. But ordinary balanced compact
partitions satisfying all three slice identities can have
arbitrarily many source sign changes on **every** line,
even when their scores choose the reference simplex
everywhere off its interfaces. The winner equation for
the defect itself remains an essential further condition.

## The count on a full reference line

[Schoenberg's full paper, Theorem 1 and its direct proof](https://acta.bibl.u-szeged.hu/13626/1/math_012_pars_b_097-106.pdf)
includes the Gaussian as a Pólya frequency kernel and
proves the sign-variation bound for bounded measurable
inputs in Section 2. The variant here counts even the
zeros where the sign stays the same. Count the sign
changes of the source on its ordered nonzero Lebesgue
points, so sets of measure zero have no effect.

Let \(f\) be any nonzero bounded source with finitely many
essential sign changes, say \(m\). The real zeros of its
Gaussian convolution agree with those of the transform
\[
 L(z)=\int_{\mathbb R}e^{zt-t^2/(2\sigma^2)}f(t)\,dt,
           \qquad z=rs/\sigma^2
\]
including multiplicities, since the removed Gaussian
factor is everywhere positive. This transform is entire
by Gaussian domination on compact subsets of \(z\), also
after multiplication of the source by a polynomial.
If \(m=0\), the transform has a strict sign on the real
axis. Otherwise choose a cut \(d\) at a sign transition.
For finite essential variation one can separate consecutive
sign runs by such cuts. Indeed an interval between two
points in a longest alternating sequence can only contain
points of its first sign followed by points of its second
sign, since a reversal between them would lengthen the
sequence. End intervals have only the respective end signs.
Multiplying by \(t-d\) removes one sign change, and
\[
       (e^{-dz}L(z))'
            =e^{-dz}\int e^{zt-t^2/(2\sigma^2)}(t-d)f(t)\,dt .
\]
If a finite list of zeros of \(L\) has total multiplicity
\(M\), the derivative on the left has at least \(M-1\)
zeros counted with multiplicity, by Rolle's theorem
between successive zeros and the order drop at each zero.
Induction bounds \(M\) by \(m\). Thus the sum of the
multiplicities of the real zeros of (OS1) is at most the
number of essential changes of its source. The statement
also permits an infinite number of source changes, in
which case it gives no finite upper bound.

For an actual compact-defect maximizer the two tails of
\(h_k\) have opposite signs, negative backward and positive
forward. Its zero at the origin is simple. An antipodal
germ gives one further sign change on each half-ray and
therefore forces at least three source sign changes for
each \(k\). With the reference germ, at a nonsimplex
right contact at least one forward ray has a strict
spectator interval. That ray has two further sign changes,
again forcing at least three source changes for that \(k\).
An additional dormant even zero costs at least two more
in this source count.
[Central germ and the right-contact arm test](THREE_CELL_COMPACT_DEFECT_CENTRAL_TRIPOD.md#genuine-triple-points-on-the-reference-lines)

## Ordinary compact examples with many changes

Let \(m\) be any positive integer, choose \(a_1\ge4\sigma\)
and \(a_{j+1}>6a_j\), and let \(S\) be the union of the
shells \(a_j<|x|<3a_j\), for \(1\le j\le m\). Start with
the fractional tuple \(b\) equal to the reference simplex
outside \(S\) and to \((1/3,1/3,1/3)\) in \(S\). Its
standard Gaussian masses are \(1/3\) by circular balance.
Multiplying each pair contrast of \(T\) by the radial
indicator of \(S^c\) preserves its oddness under reflection
in its own reference line. Hence all three tilted slice
identities hold, and the pair scores agree on those lines.

In the spectator coordinates above write \(u=-t\).
At \(u\ge0\), the reference spectator occupies
\(|v|\le\sqrt3 u\). Thus \(3q_k(-u)-1\) for the fractional
tuple equals twice the transverse Gaussian mass of
unmodified points inside this wedge minus the transverse
mass of unmodified points outside the wedge. It is
strictly negative at \(u=0\). At \(u=a_j/2\) the core
\(|v|<\sqrt3 a_j/2\) lies in the wedge and is entirely
unmodified, since its radii lie between \(a_j/2\) and
\(a_j\) and all earlier shells end below \(a_j/2\).
Its transverse mass is at least \(11/12\) by Chebyshev's
inequality, and the total negative contribution has
magnitude at most one. The result is strictly positive.
At \(u=a_j\), all points in the wedge are in the \(j\)-th
shell except null boundaries. There is positive
unmodified mass beyond the wedge, for example beyond
the last shell. The sign is strictly negative. Finally
the sign is positive sufficiently far backward, where
the shell modification does not enter the projection.
These evaluations are continuous in \(u\), by dominated
convergence away from the finitely many transverse
boundary points of the wedge and the shells. In either
ordering they give at least \(2m+1\) essential sign
changes of \(f_k\), for all three spectators.

The score winner of the fractional tuple is exactly
\(T\). To see this, fix \(i,j\) and unfold the pair
contrast across the line \(L_{ij}\). On the \(i\)-side
the unsmoothed contrast is nonnegative, and positive on
the portion of \(T_i\) outside the shells. The difference
of the two reflected Gaussian densities is strictly
positive for observer and source on that side. Thus
the smoothed pair contrast has the strict sign of the
normal coordinate away from the line. Its normal
derivative on the line is positive by differentiating
this same difference. This sign and derivative result
holds for every pair.

We now use [the weighted-projection purification theorem](THREE_CELL_COMPACT_PROJECTION_PURIFICATION.md#which-constraints-can-be-preserved)
on a disk containing all the shells. Specify the three
directional pushforwards at exponent \(1/(2\sigma^2)\),
label by label. Their densities, after dividing by
the longitudinal Gaussian factor, preserve each
transverse profile \(q_k\) almost everywhere. Add finite
standard Gaussian mass tests for the labels, and for
each label its mass on the unmodified part of its
reference sector in the disk. The latter tests fix the
reference assignment off the shells, since the
fractional mass there filled the entire unmodified
portion of that sector.

One can also keep the exact reference winner for the
purified scores. The fractional pair contrasts on the
appropriate reference sectors have a common strict
margin of the form \(c\min\{\operatorname{dist}(x,L_{ij}),1\}\),
where \(c>0\) may depend on the fixed shells and \(r\).
On compact sets this follows by positivity, the
strict normal derivative on the line and continuity.
At large \(|x|\), use the same margin for the reference
simplex and the vanishing compact-support score and
gradient perturbation. Within unit distance of the line
integrate that gradient from the line, where the
perturbation is zero. The entire segment stays within
unit distance of its large observation point.
[The saturated reference score margin](POSITIVE_CORRELATION_CONTINUATION.md#relative-exterior-control-below-every-fixed-correlation-bound)

Include in the purification a finite net of tests given
by the conditional kernel and its first observer
derivatives evaluated at finitely many points of a
large observation disk. The test functions are bounded
on the source disk. For any prescribed accuracy,
the scores and their gradients after purification
are uniformly close to the fractional ones on the
observation disk, by kernel equicontinuity in the
observer and the finite net. Outside it the error is
uniformly small by compact source support and Gaussian
tail decay, also for first derivatives. The contrast
error still vanishes on its own pair line by the
preserved projections. Taking the accuracy small
compared with \(c\), the error is at most
\((c/2)\min\{\operatorname{dist}(x,L_{ij}),1\}\),
so all reference winner signs hold. The resulting
ordinary assignment differs from \(T\), since the
preserved profiles have more sign changes than the
single change of the reference profile. It therefore
fails the winner equation for its own labels.
