# Two ways to exclude the symmetric two-radius triangle witness

Consider the three-label agreement relaxation on the Gaussian plane.
Write \(\vartheta\in[0,\pi]\) for the smaller angle between two
nonzero observations, \(a=2\pi/3\), and
\(F_T(\vartheta)=(1-\vartheta/a)_+\) for the rotation-averaged
simplex agreement. A disk \(D\) of Gaussian measure \(1/2\) divides
the radii into two types. In the symmetric two-type class the
same-type angular function is identical for both types. A separate
angular function describes pairs of different types.

Suppose \(F_\varepsilon\) is a one-sided perturbation of \(F_T\)
in this class, satisfying \(0\le F_\varepsilon\le1\), symmetric
triangle inequalities, and the exact Gaussian row mean \(1/3\).
For definiteness assume the angular functions of each perturbed
kernel are continuous, their difference quotients are uniformly
bounded, and the one-sided derivatives \(x(\vartheta)\) for same
types and \(y(\vartheta)\) for different types exist for
\(\vartheta\ne a\). Suppose

\[
 x(0)=0,\qquad x(\vartheta)=m\vartheta+o(\vartheta),\qquad
 y(0)=-b,\qquad y(\vartheta)=-b+o(1)
 \quad(\vartheta\downarrow0).
\]

The upper bound at zero gives \(b\ge0\). The condition \(x(0)=0\)
keeps the leading angular agreement of observations of the same
type. Adding positive semidefiniteness forces

\[
 ma\le \left(1-\frac{\sqrt3}{4a}\right)b .                  \tag{TS1}
\]

Alternatively impose the four-point three-label pigeonhole
inequality \(\sum_{i<j}F_\varepsilon(x_i,x_j)\ge1\) without
positive semidefiniteness. It forces

\[
                         ma\le \frac45 b .                 \tag{TS2}
\]

Either bound rules out an improvement in the leading near-one
coefficient for this two-type class.

## Active triangles and available correction mass

For \(u,v\ge0\) with \(u+v<a\), place the three points at
angles \(0,u,u+v\). The simplex saturates the triangle that
compares the two shorter angles with the endpoint angle. With
all three points of one type, and then only one endpoint
of another type, its first variations give

\[
 x(u)+x(v)\le x(u+v),\qquad
 x(u)+y(v)\le y(u+v).
\]

The constraints for these configurations follow from the assumed
pointwise inequalities. They also follow from almost-everywhere
inequalities when the angular kernels are continuous, by
restricting to the indicated radial types and then using
continuity in the angular variables. Subdivide \(\vartheta\)
into \(N\) equal increments and let \(N\to\infty\). This gives
\(x(\vartheta)\ge m\vartheta\) and
\(y(\vartheta)\ge m\vartheta-b\) for \(\vartheta<a\). Outside
\([0,a]\) both derivatives are nonnegative by the lower box
constraint. Write \(u_x,u_y\) for these nonnegative corrections
to the corresponding baselines on \([0,a)\), and use \(u_x=x\),
\(u_y=y\) outside \([0,a]\). The two radial types each have
mass \(1/2\), so the derivative of the exact row mean gives

\[
 \int_0^\pi(x+y)\,d\vartheta=0,\qquad
 \int_0^\pi(u_x+u_y)\,d\vartheta=a(b-ma).                   \tag{TS3}
\]

## The radial contrast in angular frequency two

Let \(\chi\) be \(1\) on the disk and \(-1\) off it. The
Gaussian quadratic form of \(F_T\) on \(\chi(r)\cos(2\theta)\)
vanishes, because \(F_T\) is independent of the radii and
\(\mathbb E\chi=0\). Positive semidefiniteness of the perturbed
kernels therefore requires

\[
                   \int_0^\pi(x-y)\cos(2\vartheta)\,d\vartheta\ge0 .
\]

Indeed averaging \(\cos(2\theta)\cos(2(\theta+\vartheta))\)
over \(\theta\) gives \(\cos(2\vartheta)/2\). The radial
products \(\chi(r)\chi(s)\) are \(1\) for the same type and
\(-1\) for different types, each with probability \(1/2\).
The baseline contribution to the integral is
\(b\int_0^a\cos(2\vartheta)\,d\vartheta=-b\sqrt3/4\).
The correction is at most the total correction mass in (TS3),
which gives (TS1).

## Averaging saturated four-point inequalities

Take four cyclically ordered angles with successive gaps

\[
 (3a/4+t,\ 3a/4+s,\ 3a/4-t,\ 3a/4-s),\qquad
                  t,s\ \text{independent uniform on }[-a/4,a/4].
\]

Give the four points alternating radial types. Each neighboring
gap is in \([a/2,a]\), and the two smaller diagonal angles
are \(3a/2-|t+s|\) and \(3a/2-|s-t|\), both in \([a,\pi]\).
The simplex has zero diagonal agreements and its four
neighboring agreements sum to \(4-(3a)/a=1\), so the
four-point inequality is active. The baseline contribution
to its derivative is \(3ma-4b\). Each neighboring gap has
density \(2/a\), while each diagonal has density bounded by
\(4/a\), since the signed sum or difference has a triangular
density of maximum \(2/a\) before taking its absolute value.
The expected correction is thus at most

\[
 \frac8a\left(\int_{a/2}^a u_y(\vartheta)\,d\vartheta+
                    \int_a^\pi u_x(\vartheta)\,d\vartheta\right)
 \le 8(b-ma).
\]

The derivative of the active inequality must be nonnegative,
so \(4b-3ma\le8(b-ma)\). This proves (TS2).

## Consequence for the near-one objective

Let \(X,Y\) have correlation \(\rho I_2\), and let \(E\) be
the event that they have different radial types. The pair
objective has first variation

\[
 \mathbb E[x(\vartheta)1_{E^c}+y(\vartheta)1_E]
       = m\arccos\rho-b\,\Pr(E)+o(\sqrt{1-\rho})
                    \qquad(\rho\uparrow1).
 \tag{TS4}
\]

Rotationally averaging halfspace disagreement gives
\(\mathbb E\vartheta=\arccos\rho\). For any fixed \(\eta>0\),
\(\Pr(\vartheta>\eta)=O_\eta(1-\rho)\). In fact a rotation
of the conditional mean \(\rho X\) by more than \(\eta\)
requires a displacement at least
\(\rho|X|\sin(\min\{\eta,\pi/2\})\). The independent noise
displacement is \(\sqrt{1-\rho^2}Z\); the small-ball
Rayleigh bound gives the estimate. As \(\Pr(E)=O(\sqrt{1-\rho})\),
splitting at \(\eta\), then taking \(\eta\) to zero, gives
\(\mathbb E[\vartheta1_E]=o(\sqrt{1-\rho})\) and justifies
the remainders in (TS4).

The Gaussian perimeter limit gives
\(\Pr(E)=(2/\sqrt\pi)P_\gamma(D)\sqrt{1-\rho}
+o(\sqrt{1-\rho})\). The simplex total interface area
is \(P_T=3/(2\sqrt{2\pi})\), and
\(P_\gamma(D)=\sqrt{\log(2)/2}\) because the median Rayleigh
radius is \(\sqrt{2\log2}\). Since
\(\arccos\rho\sim\sqrt{2(1-\rho)}\), (TS4) is

\[
 \frac2{\sqrt\pi}\big(maP_T-bP_\gamma(D)\big)
                   \sqrt{1-\rho}+o(\sqrt{1-\rho}).
\]

The ratio \(P_\gamma(D)/P_T=2\sqrt{\pi\log 2}/3\) exceeds
\(4/5\), and the coefficient in (TS1) is smaller than \(4/5\).
Thus the leading coefficient is negative for \(b>0\) and
nonpositive when \(b=0\). The comparison concerns this
symmetric two-type tangent class. Variable slopes on any finite
family of consecutive radial shells also admit a strict comparison.
A separate perturbation near the three-facet junction gives a gain
for the full PSD level-four relaxation near correlation one.
[The finite-shell variable-slope comparison](GAUSSIAN_AGREEMENT_VARIABLE_SLOPE_RADIAL_DUAL.md)
[The three-facet counterexample](GAUSSIAN_AGREEMENT_MULTI_FACET_COUNTEREXAMPLE.md)
[Gaussian perimeter and the two-sample limit](THREE_CELL_DEFICIT_AND_PERIMETER.md).

## Exact collapse for at most four states

Suppose \(F\) is constant on each product of two sets of a
measurable partition with at most four positive-mass states
\(B_i\). Let \(p_i=\gamma(B_i)>0\) and \(F_{ij}\) be the value
on \(B_i\times B_j\). Assume symmetry, \(0\le F_{ij}\le1\),
\(F_{ii}=1\), the row mean \(\sum_j p_jF_{ij}=1/3\), and
the four-point inequality \(\sum_{i<j}F(x_i,x_j)\ge1\).
Then \(F\) is exactly the agreement kernel of a balanced
three-cell partition, obtained by merging states if needed.
This holds on any probability space and uses neither triangles
nor positive semidefiniteness.

Every diagonal contribution gives \(p_i\le1/3\), excluding
fewer than three states. With three states their masses are
\(1/3\) and every off-diagonal entry is zero. With four states
choose \(q=p_4<1/3\), put \(u_{ij}=F_{ij}\) for \(i<j\le3\),
and put \(v_i=F_{i4}\) for \(i\le3\). Summing the first three
row equations and using \(\sum_{i\le3}p_i=1-q\) gives

\[
 q(1-\textstyle\sum_{i\le3}v_i)
   =\sum_{i<j\le3}(p_i+p_j)u_{ij}.
\]

Choosing one point from each state in the four-point inequality
gives \(1-\sum_i v_i\le\sum_{i<j}u_{ij}\). Each coefficient
on the right in the displayed equality satisfies
\(p_i+p_j\ge2/3-q>q\), so \(u_{ij}=0\) by nonnegativity.
It follows that \(\sum_i v_i=1\) and \(p_i=1/3-qv_i\)
for \(i\le3\). The remaining row gives

\[
 1/3=q+\sum_{i\le3}(1/3-qv_i)v_i
     =1/3+q(1-\textstyle\sum_{i\le3}v_i^2),
\]

so one \(v_i\) is one and the others are zero. Merging that
state with \(B_4\) gives the three equal masses and the claimed
agreement kernel.

## A common angular slope across finitely many radial shells

Partition the radii into any finite number of consecutive
positive-mass shells. Let \(F_\varepsilon\) be as in the
two-type tangent theorem above, with a separate continuous
angular function for each pair of shell types. In addition
to the range, balance and triangle conditions, impose the
four-point inequality. Write \(f_{ij}(\vartheta)\) for
the derivative at \(F_T\). Assume the difference quotients
are uniformly bounded, the derivatives exist away from
\(a\), and for some common \(A>0\),
\[
 f_{ii}(0)=0,\qquad
 f_{ii}(\vartheta)=A\vartheta/a+o(\vartheta),\qquad
 f_{ij}(\vartheta)=f_{ij}(0)+o(1)
             \quad(\vartheta\downarrow0).
\]
Let \(r_j\) be the radius of the boundary between consecutive
types \(j,j+1\), \(R=\sqrt{2\log2}\), and put
\(b_{ij}=-f_{ij}(0)\). Then
\[
 C:=\sum_j b_{j,j+1} r_j e^{-r_j^2/2}
                   \ \ge\ A\left(\frac R2+\frac{1}{16R}\right).
                                                               \tag{TS5}
\]
The leading near-one objective derivative is
\((2/\sqrt\pi)(A P_T-C)\sqrt{1-\rho}\), up to
\(o(\sqrt{1-\rho})\), with \(P_T=3/(2\sqrt{2\pi})\).
The leading coefficient is negative by (TS5).
The shell masses can be unequal and the median radius
need not be a shell boundary.

At angle zero the upper box and active triangles make
\(b_{ij}\) a pseudometric on the shell types. For
\(0\le t<a\), subdividing an active angular triangle
with one intermediate point of the same shell type
as its endpoint gives
\[
                   f_{ij}(t)\ge At/a-b_{ij}.
\]
Indeed the triangle implies
\(f_{ij}(t)\ge f_{ij}(0)+n f_{ii}(t/n)\), and the
right-hand side tends to the claimed bound. For
\(t>a\) the lower box gives \(f_{ij}(t)\ge0\).
The angular constraints can be evaluated on repeated
angles by continuity of the perturbed functions, even
if the original constraints hold almost everywhere.

Use the radial probability coordinate \(u=1-e^{-r^2/2}\)
and put \(R(u)=\sqrt{-2\log(1-u)}\). Regard \(b\) as a
step kernel of two radial probabilities. Put
\[
 \ell=\sum_j b_{j,j+1}\delta_{u_j},\qquad
 w(u)=(1-u)R(u),\qquad
 M=\int\min(u,1-u)\,d\ell(u),\qquad
 D=\int\min(u,1-u)|u-1/2|\,d\ell(u),
\]
where \(u_j\) is the probability coordinate of \(r_j\).
The pseudometric inequality gives \(b(u,v)\) bounded
above by the \(\ell\)-length between \(u,v\), and
\(C=\int w\,d\ell\). For fixed \(u\), the row equation
and the angular bounds give \(\int b(u,v)\,dv\ge A/2\).
Taking the median row, or the one-sided limit of rows
there, and bounding by path length gives \(M\ge A/2\).
An atom at the median has the same coefficient \(1/2\)
in either one-sided limit.

Let \(\overline b\) be the mean of \(b(u,v)\) for
independent uniform radial probabilities, and let
\(b_\times\) be its conditional mean for opposite
median halves. The path bound gives
\[
 \overline b\le M+2D,\qquad
 b_\times\ge 2\overline b-4D .                              \tag{TS6}
\]
Indeed a cut at \(u\) separates independent samples
with probability \(2u(1-u)\), and it separates samples
in the same median half with unconditional probability
\(2\min(u,1-u)|u-1/2|\). Integrate the two probabilities
over the length measure and use that the probability
of opposite halves is \(1/2\).

Subtract the angular baseline \(At/a-b_{ij}\) on
\(0\le t<a\) and zero on \(t>a\). The correction is
nonnegative and its total integral over angles and
radial pairs is \(a(\overline b-A/2)\) by balance.
Use the saturated four-point angular configuration
above, choosing the four radii independently in
alternating median halves. The baseline contribution
to the clique derivative is \(3A-4b_\times\).
Each neighboring angular density is at most \(2/a\)
and each diagonal density is at most \(4/a\).
There are four cross terms, whose radii are distributed
as conditional opposite halves, and two diagonal
terms, one in each same half. The integrated
correction is therefore at most \(16/a\) times the
total correction mass. Thus
\[
 5A\le16\overline b-4b_\times\le8M+32D.                     \tag{TS7}
\]

For \(u\le1/2\) the ratio \(w(u)/\min(u,1-u)\) is
\((1-u)R(u)/u\). Its derivative is
\((u-R(u)^2)/(R(u)u^2)\le-2/R\), since
\(R(u)^2\ge2u\) and \(R(u)u\le R/2\).
For \(u\ge1/2\) the ratio is \(R(u)\), whose derivative
is \(1/((1-u)R(u))\ge2/R\). Here \(w\) decreases
past the median because \(R>1\). Integration gives
\(C\ge RM+(2/R)D\). By (TS7) and \(M\ge A/2\),
\[
 5A\le (8-16R^2)M+16RC
       \le(4-8R^2)A+16RC .
\]
This proves (TS5).

Finally a change of radial type that crosses more than
one shell boundary requires a positive fixed lower
bound on \(|X-Y|\). Its probability is exponentially
small as \(\rho\uparrow1\). Crossing the boundary at
\(r_j\) has probability
\((2/\sqrt\pi)r_j e^{-r_j^2/2}\sqrt{1-\rho}
+o(\sqrt{1-\rho})\), by conditioning on the Rayleigh
radius. The angular small-gap estimate used in (TS4)
applies uniformly over the finitely many types and
gives the stated objective expansion. In particular
\(R/2+1/(16R)>P_T\), so the expansion excludes a
leading gain for this whole finite-shell class.
[Rayleigh conditioning at a radial boundary](GAUSSIAN_AGREEMENT_TRIANGLE_GRID_WITNESS.md#objective-slope-and-the-grid-threshold)

## Strict separation of finite-label consistency levels

On any nonatomic probability space partition into \(N>3\)
sets of mass \(1/N\). Give agreement one within a set and
agreement \(c_N=(N-3)/(3(N-1))\) between sets. The matrix
on these sets is positive semidefinite because it equals
\((1-c_N)I+c_N{\bf 1}{\bf 1}^{\mathsf T}\). The row mean
is \(1/3\). On the Gaussian plane these states may be
chosen as consecutive equal-mass radial shells.

For \(m<N\), the kernel satisfies level \(m\) of the
[finite-label consistency hierarchy](GAUSSIAN_AGREEMENT_TRIANGLE_DISCRETIZATION.md#a-convergent-hierarchy-using-finite-label-consistency)
whenever \(c_N\ge\beta_m\), where
\[
 \beta_m
  =\frac{\sum_{i=1}^3 n_i(n_i-1)}{m(m-1)}
  =\frac{m-3}{3(m-1)}
    +\frac{2\,1_{\,3\nmid m}}{3m(m-1)},\qquad m\ge2,               \tag{TS8}
\]
and \(n_1,n_2,n_3\) are the three integers as equal as
possible with sum \(m\). Indeed \(\beta_m\) is the minimum
pair agreement averaged over all pairs of \(m\) points
for three labels, since the agreement count
\(\sum_i\binom{n_i}{2}\) is minimized by these class sizes.
Permuting the minimizing label assignments uniformly
gives the exchangeable agreement matrix with off-diagonal
entry \(\beta_m\). Mixing with the all-one assignment
achieves every common entry from \(\beta_m\) to one.
Restrict this \(m\)-point distribution to the number of
distinct states represented in the test, and give repeated
states the same label. This proves finite-label membership
for every sample of at most \(m\) points.

For \(m\ge3\), one has \(c_{m+1}\ge\beta_m\), with equality
when \(3\nmid m\). Take \(N=m+1\) if \(3\nmid(m+1)\), and
\(N=m+2\) otherwise. Then \(c_N\ge\beta_m\) and \(c_N<\beta_N\).
A sample taking one point from each state violates the
level-\(N\) minimum agreement count. Such samples have
positive probability. In particular five equal-mass shells
with common cross agreement \(1/6\) pass level four
and fail level five, which requires cross agreement at
least \(1/5\) on exchangeable five-point tests. This is
a separation of feasible sets even with PSD, radial
symmetry and exact balance. The Gaussian objective
comparison is a separate question.

The obstruction persists for stronger finite systems
that require the local assignment laws to be consistent
under marginalization across all samples. Choose
\(N>3m\) with \(3\nmid N\). On any \(t\le m\) distinct
states give a specified label string with counts \(b_i\)
the probability
\[
                 \frac{\prod_{i=1}^3 (N/3)_{b_i}}{(N)_t},
                 \qquad (x)_j=x(x-1)\cdots(x-j+1).
\]
All factors are positive. Summing over the label of an
added state multiplies the numerator by
\(\sum_{i=1}^3(N/3-b_i)=N-t\), so the prescribed laws
marginalize to one another and have total probability
one. The two-state law has agreement \(c_N\) as above.
Repeated states inherit a single common label. There is
still no consistent law on all \(N\) states with these
pair marginals, by the level-\(N\) minimum count.
