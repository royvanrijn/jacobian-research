# Proper levels of the junction heights in a three-cell graph

Let \(A\) be an actual balanced planar three-cell Gaussian maximizer
at \(0<r<1\), with exclusive vertical fibers and a finite continuous
height \(h\). Thus the graph is the upper boundary of each of the
two base labels. Use the [exclusive-fiber representation and signed
height laws](THREE_CELL_FIBER_EXCLUSIVITY_RIGIDITY.md#sorting-fibers-and-identifying-their-base-labels).
Write \(B\) for the set of switches between the base labels and
\(\epsilon(s)\in\{-1,1\}\) for their alternating sign. The set
\(B\) is locally finite and nonempty.

Fix a height \(t_0\) attained at a base switch. Suppose there is
a closed neighborhood of \(t_0\) with compact preimage under \(h\),
and both outward height slopes are nonzero at every switch
of height \(t_0\). Then:

- The left and right outward slopes agree at every switch of height
  \(t_0\), and both are strictly positive.
- There are finitely many base switches. The switches of height
  \(t_0\) include the first and last. Between successive switches
  of this height, at most one switch of a different height occurs.

In particular a uniquely attained base-switch height satisfying the
hypotheses forces affine separation. Two distinct base-switch heights
cannot both satisfy the hypotheses. If every base-switch height
satisfies them, all switches have one common height and each is a
strict local minimum of \(h\).

The final section gives a stronger finite-switch conclusion by a
different use of the height laws. Every finite switch set has just one
switch, without assumptions on the junction heights or slopes.

## The signed cumulative distributions

Write \(\sigma=\sqrt{1-r^2}\) and \(p_\sigma\) for the centered
Gaussian density of standard deviation \(\sigma\). For every
observation \(b\in B\), vertical deconvolution gives (FX3),
so
\[
           \int_{\mathbb R}
               \epsilon(s)1_{\{h(s)\le t\}}
                  p_\sigma(s-rb)\,ds=0
                  \qquad(t\in\mathbb R).                    \tag{PSH1}
\]
On each open base interval the height is analytic by the
implicit function theorem applied to the contrast against
label 3, whose vertical derivative is positive. The analytic
branches extend to their finite endpoints. None of these
branches is constant. Indeed constancy on a base interval of
label \(i\), at height \(a\), would make its stationary pair
contrast vanish on the whole horizontal line \(y=a\) by
analyticity. The strict vertical derivative would then give
\(h=a\) on every interval of label \(i\). Equation (PSH1) for
the complementary height set forces the other label's heights
to be \(a\) as well, since its sign is constant there and the
Gaussian density is positive. This would make cell 3 a
half-plane, contrary to the [individual-cell dimension
restriction](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#each-cell-uses-both-active-directions).

Consequently the preimages of \(t_0\) are finite. Compactness
of the height band keeps them in a fixed compact interval,
which meets finitely many base intervals. In each such interval
the analytic nonconstant branch and its extensions at the
endpoints have only finitely many preimages of that height.
Choose disjoint small neighborhoods of these points. The
integral over the rest of the real line is locally constant
as \(t\) crosses \(t_0\). Points outside the compact band
cannot contribute a change in a sufficiently small height
interval, and the remaining compact set with the neighborhoods
removed stays a positive height distance from \(t_0\).

Write \(v=t-t_0\). At an interior preimage, use an analytic
local coordinate \(u\) increasing with \(s\) in which the
height is \(t_0\pm u^n\). With \(\chi\) a primitive of the
local signed density in this coordinate, the cumulative
contribution at an even-order preimage is, up to sign and
an additive constant, a one-sided difference
\(\chi(|v|^{1/n})-\chi(-|v|^{1/n})\). It contains no integer
powers of \(|v|\), because the numerators of the powers are
odd and \(n\) is even. At an odd-order preimage the
contribution has the form
\(\chi(\operatorname{sgn}(v)|v|^{1/n})\), with possible
orientation signs in the argument and the value. Its
integer-power terms form the same signed-\(v\) germ on
both sides, since \(\operatorname{sgn}(v)^k|v|^q=v^q\)
for \(k=nq\). Thus every interior preimage makes zero
contribution to the right-minus-left integer-power germ.
This is the [height-density argument for a single
junction](THREE_CELL_AFFINE_GRAPH_REFLECTION_TESTS.md#the-density-germ-at-the-junction-level),
now applied to each of the signed densities in (PSH1).

## The endpoint germs

Write \(b_1<\cdots<b_m\) for the switches of height \(t_0\)
and let \(\epsilon_k\) be the base sign immediately left
of \(b_k\). The sign immediately right of it is
\(-\epsilon_k\). Let \(u_{k,L}(v)\) and \(u_{k,R}(v)\) be the
analytic inverse displacements of the graph branches extending
outward from \(b_k\). Thus the left inverse is
\(b_k-u_{k,L}(v)\), the right inverse is \(b_k+u_{k,R}(v)\),
and both \(u\)'s vanish at zero. They exist because both
outward slopes are nonzero. Put
\(a_{k,L}=u_{k,L}'(0)\) and \(a_{k,R}=u_{k,R}'(0)\).

For either branch the right-minus-left integer-power germ is
its primitive from displacement zero to its inverse
displacement. If the outward slope is positive, this primitive
occurs only on the right in \(v\). If the slope is negative,
its negative occurs only on the left. Applying this at
each switch and subtracting the two cumulative germs in (PSH1)
gives the exact analytic identity
\[
 \sum_{k=1}^m\epsilon_k
 \left[
  \int_0^{u_{k,L}(v)}p_\sigma(b_k-u-rb)\,du
   -\int_0^{u_{k,R}(v)}p_\sigma(b_k+u-rb)\,du
 \right]=0
 \qquad (b\in B) .                                          \tag{PSH2}
\]
All germ operations are on finite sums of convergent one-sided
Puiseux expansions. In particular the identity holds for all
sufficiently small \(v\) at each fixed observation \(b\).
The radius in \(v\) can be chosen independently of \(b\), since
the inverse displacements are defined by the height branches
and the Gaussian densities in (PSH2) are entire.

The linear coefficient at \(v=0\) in (PSH2) is
\[
 \sum_{k=1}^m \epsilon_k
       (a_{k,L}-a_{k,R})p_\sigma(b_k-rb)=0
                         \qquad(b\in B).                  \tag{PSH3}
\]
Use the \(m\) observation rows \(b=b_1,\ldots,b_m\).
After positive row and column factors, their matrix is
\(\exp(r b_jb_k/\sigma^2)\), which has strictly positive
ordered minors by the [Gaussian total positivity
argument](THREE_CELL_FIBER_EXCLUSIVITY_RIGIDITY.md#counting-switches-with-a-height-band).
It is invertible, hence \(a_{k,L}=a_{k,R}=a_k\) for each \(k\).
The upward slope jump (FX5) is the sum of the outward
slopes, now equal to \(2/a_k\). It is strictly positive,
so \(a_k>0\). In particular both inverse displacements
are positive for every sufficiently small \(v>0\).

## Counting with the endpoint germs

For a fixed sufficiently small \(v>0\), take the function
\[
 g_v(s)=\sum_{k=1}^m\epsilon_k
  \left[
     1_{(b_k-u_{k,L}(v),\,b_k)}(s)
      -1_{(b_k,\,b_k+u_{k,R}(v))}(s)
  \right] .
\]
Its intervals are nonempty and disjoint after choosing
\(v\) small, and (PSH2) says that \(T_r^{[1]}g_v\) vanishes
on \(B\). Each pair of intervals contributes one sign change.
Between successive pairs there is one additional sign change
exactly when \(\epsilon_k=\epsilon_{k+1}\). If \(q\) counts
these equal-sign gaps, there are \(m+q\) sign changes.
The Gaussian zero count gives
\[
                            |B|\le m+q.                   \tag{PSH4}
\]
In particular this argument also excludes an infinite switch set,
since it gives a finite bound using only the selected height.

The count can equivalently be made using Gaussian derivatives
at the switch points. Write
\(u_{k,L/R}(v)=a_kv+c_{k,L/R}v^2+O(v^3)\).
The quadratic equations in (PSH2) then say that all observations
\(b\in B\) are zeros of
\[
  \sum_k\bigl[
     \epsilon_k(c_{k,L}-c_{k,R})p_\sigma(b_k-rb)
      -\epsilon_k a_k^2p_\sigma'(b_k-rb)
  \bigr].
\]
This is the transform of a combination of atoms and their
first derivatives. If there were \(m+q+1\) observation
zeros, interpolate a dual function \(L(s)\) from those
Gaussian kernel rows, with simple zeros at the \(m\)
selected switches and one extra zero in every equal-sign
gap. The construction by the determinant with a variable
column in the Gaussian zero-count argument gives the
alternating signs. The zeros are simple because, after
removal of the positive Gaussian factor, a nonzero
combination of \(m+q+1\) distinct real exponentials has
at most \(m+q\) real zeros counted with multiplicity.
This follows by multiplying by the reciprocal of the
lowest exponential, differentiating, and using induction
and Rolle's theorem. A global sign can therefore make
\(L'(b_k)\) agree in sign with \(-\epsilon_k\) for every \(k\).
The atom contributions vanish, and the derivative terms
give a strictly positive sum
\(\sum_k(-\epsilon_k a_k^2)L'(b_k)\), whereas the choice of
observation rows makes it zero. This is another proof of
(PSH4).

The signs immediately left of all switches in \(B\)
alternate. Thus each equal-sign gap among the selected
switches contains at least one omitted switch, and
\(q\le |B|-m\). Equality is required by (PSH4).
Every omitted switch must therefore be the only omitted
switch in one of those gaps. There can be none before the
first selected switch or after the last. This proves the
stated rigidity of the switch heights.

## Critical height accumulation excludes multiple finite switches

Suppose \(B\) is finite with \(N=|B|\ge2\). The
[common height range on every base interval](THREE_CELL_FIBER_EXCLUSIVITY_RIGIDITY.md#counting-switches-with-a-height-band)
forces \(h\) to be bounded, since at least one base interval
is bounded. Both outer graphs are isolated bounded ends,
so each tends to the offset \(H_*\) given by its pair
pressure. The two pressures are equal by (FX2), and both
graphs have cell 3 above and their base label below.
Thus the offsets agree.
[Bounded ends and the pressure offset](THREE_CELL_ISOLATED_END_RIGIDITY.md#bounded-ends-converge-to-the-pressure-offset)

Neither tail of \(h-H_*\) can have an eventual fixed sign.
The [one-sided end straightening theorem](THREE_CELL_ISOLATED_END_RIGIDITY.md#a-one-sided-approach-forces-exact-straightness)
would make that tail exactly horizontal. Analyticity of
the contrast between its base label and label 3 on the
horizontal line then makes all graph heights on intervals
of that base label equal to \(H_*\). As in the proof of
(PSH1), equality of the signed height laws forces the
other label's heights to be \(H_*\) as well. This
contradicts the noncylindrical cell restriction. Each tail
therefore has infinitely many values of both signs relative
to \(H_*\).

Choose a fixed compact interval containing all switches
in its interior. There are only finitely many exceptional
heights of this interval. They are its endpoint heights,
all junction heights, and the critical values of the
analytic graph branches within it. Finiteness follows
from nonconstancy and the analytic extensions across the
endpoints of every branch in the interval.

On the right tail, take arbitrarily far excursions with
a value above \(H_*\) between two values below \(H_*\).
The maximum on the interval between those below-offset
points is attained at an interior critical point and gives
a critical value \(t>H_*\). The convergence to \(H_*\)
gives a sequence of such values converging to \(H_*\)
from above, from which one may extract distinct values.
Fix one of them outside the finite set of exceptional
heights. Its entire preimage is compact by the limits
at both horizontal infinities and \(t\ne H_*\).
The preimage is finite by local analyticity. All its
critical points are in the two outer base intervals,
and none are switches.

Let \(m\ge2\) be the maximum order of a critical point
in this preimage. If a critical point \(s\) has order \(m\)
and \(h(x)-t=\alpha_s(x-s)^m+O(|x-s|^{m+1})\), its
local mass in the symmetric height band
\([t-\varepsilon,t+\varepsilon]\) under the weight
\(p_\sigma(x-rb)\,dx\), divided by \(\varepsilon^{1/m}\),
tends to
\[
                         2|\alpha_s|^{-1/m}p_\sigma(s-rb).
                                                               \tag{PSH5}
\]
This follows by the change of variables
\(x=s+\varepsilon^{1/m}u\) on the local band.
All lower-order preimages, including the regular ones,
make a vanishing contribution after this division.
The exact equality of the two height laws from (PSH1)
therefore implies
\[
  \sum_{\substack{s:\ h(s)=t,\ {\rm order}=m\\\ell(s)=1}}
          2|\alpha_s|^{-1/m}p_\sigma(s-rb)
    =
  \sum_{\substack{s:\ h(s)=t,\ {\rm order}=m\\\ell(s)=2}}
          2|\alpha_s|^{-1/m}p_\sigma(s-rb)
                 \qquad (b\in B).                        \tag{PSH6}
\]
Both sums must be nonempty, since some point has
order \(m\) and every summand is strictly positive.
This is already impossible when the two outer intervals
have one common base label.

If their labels are opposite, label the right one \(i\)
and the left one \(j\). Choose two observation switches
\(b<b'\). The ratio \(p_\sigma(s-rb')/p_\sigma(s-rb)\)
is a strictly increasing function of \(s\). Every
summand on side \(i\) in (PSH6) accrues a strictly larger
ratio than every summand on side \(j\), since the
contributing critical sets lie in the disjoint ordered
outer intervals. The equal positive sums at \(b\)
therefore become different at \(b'\), contrary to (PSH6).
This excludes \(N\ge2\).

Thus exclusive fibers with a finite continuous height
and a finite number of base switches give a single affine
separating line. In particular a height band with bounded
preimage of positive measure suffices, by the
[finite switch reduction](THREE_CELL_FIBER_EXCLUSIVITY_RIGIDITY.md#counting-switches-with-a-height-band).
