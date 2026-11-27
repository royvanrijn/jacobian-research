# Exact angular pairing near equal circular masses

Fix a compact correlation interval \(J\subset(0,1)\).
There are \(\varepsilon_J,c_J>0\) and fixed
\(0<r_J<R_J<\infty\) such that, for \(\rho\in J\), any
ordinary balanced planar three-cell partition \(A\) with
circular mass vector \(p(r)\), \(b(r)=p(r)-(1/3)\mathbf1\),
and \(\operatorname*{ess\,sup}_{r_J\le r\le R_J}|b(r)|
\le\varepsilon_J\) satisfies
\[
 S_\rho(T)-S_\rho(A)
 \ge c_J\left(\int_0^{R_J} r|b(r)|^2\,d\mu(r)
      +\int_{R_J}^\infty |b(r)|^2\,d\mu(r)\right).
 \tag{JA1}
\]
Equality in the simplex comparison forces a common rotated
equal-arc partition at almost every radius. The comparison
allows arbitrary angular shapes. The circular masses below \(r_J\) and
beyond \(R_J\) are arbitrary. Use the notation of
[DG1](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md)
and \(p_*=(1/3)\mathbf1\).

## Exact pair envelope near balanced masses

For fixed \(K<\infty\), uniformly for \(0<\kappa\le K\) and
\(p,q\) sufficiently close to \(p_*\), the maximum of
\[
 \mathcal F_\kappa(u,v)=\sum_i\langle u_i,U_\kappa v_i\rangle_\sigma
 \tag{JA2}
\]
over circular partitions with prescribed masses \(p,q\) is attained
by two partitions into consecutive arcs with the same cyclic label
order. Every ordinary maximizing pair has this structure.

Allow measurable simplex-valued functions with prescribed masses.
The soft feasible sets are weak-star compact. Joint weak-star
continuity of the pairing follows by approximating the continuous
kernel uniformly by finite product sums. An ordinary maximizing
pair exists. First fix one argument, take an extreme point of the
maximizing face of its partner, then repeat with the other. An
extreme point of a mass slice is ordinary because otherwise two
components have simultaneous positive lower bounds on a set of
positive measure, supporting a nonzero bounded mean-zero transfer
with both signs feasible.

Use the normalized kernel \(\widehat k_\kappa=(k_\kappa-1)/\kappa\)
and \(\widehat k_0=\cos\). This only removes a constant pairing
\(p\cdot q\) and divides by a positive scalar when \(\kappa>0\).
The normalized kernels are continuous in each \(C^j\) on
\([0,K]\). They give positive semidefinite operators since
the nonconstant Fourier eigenvalues are
\(I_m(\kappa)/(\kappa I_0(\kappa))>0\) for \(\kappa>0\), and
at zero the two first real harmonics survive.

At \(p=q=p_*\), the normalized pairing is at most the mean
of the two self-pairings. For ordinary partitions each self
maximum is attained only by a partition into equal arcs.
For \(\kappa>0\), subtract the positive cosine component
\(a_\kappa\cos\vartheta\), where
\(a_\kappa=\kappa e^{-\kappa}/I_0(\kappa)\). In the variable
\(z=\cos\vartheta\), the derivative of the remainder is
\(\kappa(e^{\kappa z}-e^{-\kappa})/I_0(\kappa)\ge0\).
Add an irrelevant constant if needed and apply circular
rearrangement (DG0) to each cell for this remainder.
The cosine component is the squared first-moment norm and
has a strict fixed-mass cap at a centered arc, by projecting
in its first-moment direction and choosing the largest
cosine values. Thus equality for the full kernel forces
each cell to be an arc. At \(\kappa=0\) the normalized
kernel is cosine and the same argument applies. A fractional profile cannot maximize a self-pairing:
the feasible mean-zero transfer above strictly increases
the average of the two perturbed self-pairings. At positive
\(\kappa\) this uses strict positivity on mean-zero functions.
At zero choose the transfer with nonzero first angular
harmonic, possible because \(e^{i\theta}\) is not constant
almost everywhere on any set of positive measure. Equality
in the normalized Cauchy comparison therefore forces both
balanced extremizers to be the same arcs. At positive
\(\kappa\) this follows from positive definiteness on
mean-zero functions; at zero equal first moments fix all the
arc centers.

Otherwise take a contradicting sequence \(\kappa_j\to\kappa_*
\in[0,K]\), \(p_j,q_j\to p_*\) and ordinary maximizing pairs.
Their limiting normalized value is at least the equal-arc
value by testing close consecutive arcs of the given masses.
The preceding equality analysis and compactness show that
the pair converges weak-star, hence in \(L^1\), to common
equal arcs. The normalized convolved scores converge in
\(C^1\). Each slice separately maximizes the linear functional
given by the other slice's scores with its three prescribed
masses. Prices exist because the optimal linear assignment
value is finite and concave as a function of the mass
vector; a supergradient at an interior mass vector supplies
the prices and hence the unconstrained score-winner rule.
Price differences are uniformly bounded: a difference
exceeding twice the score supremum prevents a label from
winning despite its positive prescribed mass. Any centered
limit of the prices makes the reference arcs winners. At
their boundaries both incident scores agree by reflection
and are approached by winners from the two sides. Their
price differences must therefore be zero.

The equal-arc convolution score strictly decreases with
angular distance \(\theta\in(0,\pi)\) to its arc center.
For width \(\ell=2\pi/3\), the derivative is the difference
of kernel values at \(\theta+\ell/2\) and \(\theta-\ell/2\),
divided by \(2\pi\). Their cosines differ by
\(-2\sin\theta\sin(\ell/2)<0\), and the normalized kernel
is strictly increasing in the cosine. Thus the only
winning ties are the three arc endpoints, with simple
pair derivative of magnitude
\(2(k_\kappa(0)-k_\kappa(\ell))/(2\pi\kappa)\) and a strictly
lower third score. The limiting derivative at zero
parameter is also positive. Compactness away from these
points and monotonicity in neighborhoods of them show that
the nearby score-winner sets have exactly the three simple
switching points, in common cyclic order. This proves the
uniform pair-envelope assertion.

## Joint Hessian and remainder

Use equal-arc reference endpoints \(a_j=2\pi j/3\) and
perturb them by \(x_j,y_j\) in the two slices. Put
\(k_0=k_\kappa(0),k_\ell=k_\kappa(2\pi/3)\) and
\(\mathbf J=\mathbf1\mathbf1^T\). The raw pairing satisfies
\[
 \mathcal F_\kappa-\mathcal F_\kappa^*
 =\frac{1}{(2\pi)^2}\left[
 -(k_0-k_\ell)(|x|^2+|y|^2)
 +x^T((2k_0+k_\ell)I-k_\ell\mathbf J)y\right]
 +O_K\bigl(\kappa(|x|+|y|)^3\bigr).
 \tag{JA3}
\]
The first derivatives vanish by the endpoint ties. A
same-slice second endpoint derivative is
\(2(k_\ell-k_0)/(2\pi)^2\), from the derivative of the
two adjacent scores, and different same-slice endpoints
do not interact in this derivative. A mixed pair of
endpoint derivatives gives the kernel at those endpoints
times the dot product of the two label incidence
vectors, divided by \((2\pi)^2\); incidence products are
\(2\) on the diagonal and \(-1\) off it. For the remainder
subtract the constant kernel, whose contribution is
exactly \(p\cdot q\), quadratic in the endpoints.
Third derivatives of the remaining endpoint integrals
are bounded by \(C_K\kappa\) using the \(C^2\) kernel bound.
The derivative remainder is bounded in the corresponding
second-order scale.

The cyclic width map, divided by \(2\pi\), has squared
singular value \(3/(2\pi)^2\) on the mean-zero space.
The centered endpoint norms and dot product are therefore
\((2\pi)^2/3\) times \(|b|^2,|c|^2,b\cdot c\), where
\(b=p-p_*, c=q-p_*\). In (JA3) the mean rotations enter
only as \(-3(k_0-k_\ell)(\overline x-\overline y)^2/(2\pi)^2\).
Remove the common mean rotation. At a maximizing pair its
relative mean rotation tends uniformly to zero by the
compactness proof and its derivative vanishes. Dividing
that derivative in (JA3) by \(\kappa\) gives leading
coefficient bounded strictly away from zero on \([0,K]\)
and a remainder
\(O_K((|b|+|c|+|\overline x-\overline y|)^2)\). Absorbing
the square of the small relative rotation shows it is
\(O_K((|b|+|c|)^2)\). The joint pair envelope consequently
has expansion
\[
 \frac13\big[-(k_0-k_\ell)(|b|^2+|c|^2)
                 +(2k_0+k_\ell)b\cdot c\big]
         +O_K\big(\kappa(|b|+|c|)^3\big)
 \tag{JA4}
\]
relative to the common equal arcs. The expression gives the
Hessian of the joint pair envelope.

## Arc bounds at arbitrary tail masses

For any \(\kappa\) the two-set circular rearrangement upper
bound is \(\sum_i B_\kappa(p_i,q_i)\), using (DG9).
With
\(\Psi(v)=\int_0^{|v|}(|v|-w) k_\kappa(\pi w)\,dw\)
one has \(B_\kappa(p,q)=[\Psi(p+q)-\Psi(p-q)]/2\).
For all simplex profiles there are fixed \(K_0,C_0,c_0>0\)
such that the centered-arc upper bound minus the equal-arc value
is at most
\[
 \begin{cases}
 |b|\,|c|+C_0\kappa(|b|^2+|c|^2), & 0\le\kappa<K_0,\\
 C_0|b|\,|c|-c_0(|b|^2+|c|^2), & \kappa\ge K_0 .
 \end{cases}
 \tag{JA5}
\]
For bounded \(\kappa\), the constant kernel gives the exact
\(b\cdot c\). For the remainder kernel \(k_\kappa-1\),
all second derivatives in the bounded mass domain are
\(O_{K_0}(\kappa)\) by the formula for \(\Psi\), and the
linear terms at \(p=q=p_*\) sum to zero. This proves its
bound without a small-mass-deviation assumption.

For the large-\(\kappa\) bound first fix a small absolute
\(\varepsilon_0>0\), for instance \(1/12\). If
\(|b|+|c|\le\varepsilon_0\), the second derivative of
\(\Psi(2/3+b_i+c_i)\) is at most
\(C\sqrt{1+\kappa}e^{-c\kappa}\). This uses
\(I_0(\kappa)\asymp e^\kappa/\sqrt{1+\kappa}\), by comparing
\(1-\cos\theta\) with \(\theta^2\) on \([-\pi,\pi]\).
The linear terms sum to zero. On
\(0\le w\le c'/\sqrt{1+\kappa}\) the kernel
\(k_\kappa(\pi w)\) is at least \(c\sqrt{1+\kappa}\).
Restrict the integral for \(\Psi(d)\) to the smaller
of that interval and \([0,|d|/2]\) to get
\(\Psi(d)\ge c\min(\sqrt{1+\kappa}\,d^2,|d|)\ge c d^2\)
for \(|d|\le1\). This negative difference term gives
a fixed negative diagonal after expansion, absorbing
the exponentially small positive one.

If \(|b|+|c|>\varepsilon_0\), then
\(|b|^2+|c|^2\ge\varepsilon_0^2/2\). Always
\(\sum_i B_\kappa(p_i,q_i)\le\sum_i\min(p_i,q_i)
=1-\tfrac12\|b-c\|_1\), since the angular kernel is Markov.
For the reference pair, \(3B_\kappa(1/3,1/3)\ge
1-C/\sqrt{1+\kappa}\). To see this, condition on the signed
angular displacement \(\Theta\in[-\pi,\pi]\) with density
\(k_\kappa\) relative to \(\sigma\). The loss of intersection
under that displacement in the three reference arcs is at
most \(C|\Theta|\), and the Gaussian upper estimate for
\(k_\kappa\) just used gives
\(\mathbb E|\Theta|\le C/\sqrt{1+\kappa}\).
Since differences of simplex profiles have bounded norm,
\(\|b-c\|_1\ge c|b-c|^2\). Take \(K_0\) large enough that
\(C/\sqrt{1+\kappa}\le c'(|b|^2+|c|^2)\) for
\(\kappa\ge K_0\) in this case, with \(c'\) smaller than
the negative diagonal from expanding \(-c|b-c|^2\).
This proves the second line of (JA5) also for arbitrary
large-amplitude deviations. The choices of \(K_0\) and
all constants here are absolute.

## Angular comparison at small parameter, with arbitrary masses

Write \(\mathcal A_\kappa(p)=
\sum_i(B_\kappa(p_i,p_i)-p_i^2)\), the centered self-pairing
of consecutive arcs of widths \(p\).
For a sufficiently small fixed \(\kappa_0>0\) and every
\(0<\kappa\le\kappa_0\), one has uniformly for all \(p\)
\[
 \mathcal A_\kappa(p_*)-\mathcal A_\kappa(p)
                    \ge c\kappa|p-p_*|^2 .
 \tag{JA5a}
\]
The first cosine term has contribution
\(\kappa\sum_i\sin^2(\pi p_i)/\pi^2\).
For \(A_i=\pi p_i\), the identity
\(\sum_i\sin^2 A_i=2+2\prod_i\cos A_i\) holds since
\(\sum_i A_i=\pi\). If any angle is nonacute the product is
nonpositive and the deficit from \(9/4\) is at least \(1/4\).
If all angles are acute, \((\log\cos A)''=-\sec^2 A\le-1\)
gives
\(\prod_i\cos A_i\le(1/8)\exp(-\sum_i(A_i-\pi/3)^2/2)\).
In either case the cosine deficit is at least a constant
times \(|p-p_*|^2\).

The residual kernel \(k_\kappa-1-\kappa\cos\) is \(O(\kappa^2)\)
uniformly on the circle. The arc-pair formula via \(\Psi\)
shows that its self-pairing has Hessian \(O(\kappa^2)\)
as a function of the masses. Its linear terms about \(p_*\)
cancel on the sum-zero space by symmetry. Thus the residual
changes the cosine deficit by at most
\(C\kappa^2|p-p_*|^2\). Choosing \(\kappa_0\) small proves
(JA5a). Circular self-rearrangement bounds each centered
self-pairing by \(\mathcal A_\kappa\). Apply angular
Cauchy--Schwarz to the centered pairing of any two slices
to conclude that for \(0<\kappa\le\kappa_0\) their raw
pairing minus the simplex raw pairing is at most
\[
                  b\cdot c-c\kappa(|b|^2+|c|^2).
 \tag{JA5b}
\]
Both (JA5a) and (JA5b) allow arbitrary mass vectors.

## Quantitative radial comparison

Fix the compact interval \(J\) and split the radii at
\(0<r_0<1\) and \(R>2\), with \(r_0\) to be chosen after \(R\).
Put \(b_{\rm mid}=1_{\{r_0\le r\le R\}}b\),
\(b_{\rm low}=1_{\{r<r_0\}}b\), and \(b_{\rm out}=1_{\{r>R\}}b\).
Write
\[
 X=\int r|b_{\rm mid}|^2d\mu,\qquad
 Z=\int r|b_{\rm low}|^2d\mu,\qquad
 Y=\int |b_{\rm out}|^2d\mu .
\]
For \(r,u\in[r_0,R]\) use (JA4). Its integrated error
is at most \(C_{J,R}\varepsilon X\), where \(\varepsilon\)
bounds the middle-annulus deviation. Indeed
\(\kappa\le C_J r u\), \(u\le R\), and \(W_\rho\)
has marginal \(\mu\). The \(\kappa\) factor gives the
extra radial power even if the lower cutoff decreases.

For \(f(r)=r b_{\rm mid}(r)\), as a standard edge vector,
the quadratic of (JA4) integrated over all \(r,u\) with
\(b_{\rm mid}\) is exactly
\[
 -\frac{(2\pi)^2}{3}
       \left[\int a_\rho |f|^2d\nu-\langle f,K_+f\rangle_\nu\right].
 \tag{JA6}
\]
In fact \(k_\kappa(\vartheta)dW_\rho=(2\pi)^2 r u
M_{\cos\vartheta}(r,u)d\nu(r)d\nu(u)\) by (DG1) and (S14).
The diagonal in (JA6) uses
\(\int(M_1-M_{-1/2})(r,u)u\,d\nu(u)=a_\rho(r)r/2\)
by (S17). The mixed coefficient is \(K_+\).
On the actual middle--middle pairs only the negative diagonal
with \(u\notin[r_0,R]\) is missing. This costs at most
\((\epsilon_R+C_Jr_0^3)X\), with \(\epsilon_R\to0\)
uniformly in \(J\). Indeed \(k_\kappa(0)-k_\kappa(2\pi/3)
\le C\kappa\) for all \(\kappa\), and
\(\mathbb E[U1_{\{U>R\}}\mid r]\le C_Je^{-c_JR^2}\)
uniformly for \(r\le R\), because the conditionally Gaussian
radius has mean-vector length at most \((\sup J)R<R\)
and bounded variance. For \(u<r_0\), uniformly bounded
planar conditional density gives
\(\mathbb E[U1_{\{U<r_0\}}\mid r]\le C_Jr_0^3\).

The standard mean of \(f\) is controlled by balance:
\[
 \int f\,d\nu=(2\pi)^{-1}\int b_{\rm mid}\,d\mu
        =-(2\pi)^{-1}\int(b_{\rm low}+b_{\rm out})\,d\mu .
\]
Its squared norm is bounded by
\(C[r_0 Z+\mu((R,\infty))Y]\).
For the first term use Cauchy--Schwarz with weight \(r\)
and \(\int_0^{r_0}r^{-1}d\mu=O(r_0)\). For the second
use ordinary Cauchy--Schwarz. Subtract this mean and apply
the standard gap (S37), uniform on \(J\). Reinsert the
constant using Cauchy--Schwarz on the standard row
potential, which is uniformly bounded on \(J\). Thus the
bracket in (JA6) is at least \(c_JX-C_Jr_0 Z-\epsilon_R Y\).
The middle--middle contribution to the raw stability difference
is at most
\[
 (-c_J+\epsilon_R+C_J r_0^3+C_{J,R}\varepsilon)X
           +C_Jr_0 Z+\epsilon_R Y.
\]

For outer pairs \(\max(r,u)>R\), use (JA5) with no bound
on the outer mass deviations. Its positive mixed term is
at most \(\epsilon_R(X+Y+Z)\).
On inside--outside pairs apply Cauchy--Schwarz to
\(P_\rho(r,u)/\sqrt r\) and
\(\sqrt r\,|b(r)|\,|b_{\rm out}(u)|\) in the product measure.
On outside--outside pairs use no \(1/\sqrt r\) factor.
In both cases the squared kernel tails tend to zero
uniformly on \(J\), since \(I_0(\kappa)\le e^\kappa\) gives
\[
 P_\rho(r,u)^2\,d\mu(r)d\mu(u)
   \le C_J r u\,e^{-c_J(r^2+u^2)}\,dr\,du.
\]
The weighted \(1/r\) is harmless. The positive diagonal
allowed in (JA5) for \(\kappa<K_0\) is at most
\(\epsilon_R(X+Y+Z)\). For \(r>R\), that event implies
\(u\le C_J/r\), while the conditional Gaussian has
mean-vector length at least \((\inf J)r\), so its
probability is at most \(C_Je^{-c_JR^2}\).
For \(r\le R\) on an outer pair use
\(\kappa\le C_J r u\) and the uniform exponential bound
on \(\mathbb E[U1_{\{U>R\}}\mid r]\).
The negative term for \(\kappa\ge K_0\) supplies at least
\(c_JY\) by the same conditional probability estimate.
Thus outer pairs contribute at most
\(\epsilon_R(X+Z)+(-c_J+\epsilon_R)Y\).

For the remaining inner pairs \(\min(r,u)<r_0\) and
\(\max(r,u)\le R\), choose \(r_0\) so small that
\(C_Jr_0R\le\kappa_0\), and apply (JA5b).
The mixed integral is bounded by \(C_J\sqrt{r_0}(X+Z)\)
by Cauchy--Schwarz on \(P_\rho(r,u)/\sqrt{ru}\) and
\(\sqrt{ru}\,|b(r)|\,|b(u)|\).
The first squared integral is \(O_J(r_0)\) over these
pairs by the displayed Gaussian kernel bound after the
\(ru\) factor cancels. The negative diagonal supplies
at least \(c_JZ\): for \(r<r_0\) integrate over
\(u\in[1,2]\), which lies in this region and has conditional
probability bounded below uniformly on \(J\) and \(r<r_0\).
On it \(\kappa\ge c_Jr\). The contribution of inner pairs
is thus at most \(C_J\sqrt{r_0}(X+Z)-c_J Z\).

Choose \(R\) large depending only on \(J\), then \(r_0\)
small depending on \(J,R\), and then \(\varepsilon\)
small enough to apply (JA4) and absorb the integrated cubic
error. The three bounds prove (JA1) with \(r_J=r_0\),
\(R_J=R\), since \(X+Z=\int_0^R r|b(r)|^2d\mu\).
If stability equals the simplex value, (JA1) gives \(b=0\).
Then (DG3) and the equal-circle comparison make the
angular self-deficit and radial angular mismatch both
zero. Strict self rearrangement and positive angular
Fourier multipliers imply equal arcs with the same
orientation and labels at almost every radius.

An unrestricted global comparison still requires controlling
profiles with large annular mass deviation.

## Consecutive arcs can fail away from equal masses

Fix any $\kappa>0$. Take $p=(1/3,1/3,1/3)$ and
$q^0=(1/2,1/2,0)$. In the second slice let labels 1 and 2
be opposite halfcircles. In the first slice take the arcs of
length $1/3$ with those same two centers, giving their
two-component complement to label 3. This pair attains
$2B_\kappa(1/3,1/2)$ by simultaneous saturation of the two
individual circular rearrangements. It is therefore a global
pair maximizer at these masses.

No pair of consecutive-arc partitions realizes that value.
Equality in each of the two individual comparisons forces
both sets in that label to be arcs with the same center.
Indeed subtract the positive cosine component as in the proof
of (JA2) and use the strict first-moment cap and vector
Cauchy–Schwarz for the cosine pairing. This applies also
to different prescribed lengths. The two halfcircles have
opposite centers. The three consecutive arcs of equal length
in the first slice have pairwise angular center separation
$2\pi/3$, so the two required alignments are impossible.

The strict difference persists at positive mass vectors
$q^\delta=((1-\delta)/2,(1-\delta)/2,\delta)$ for all sufficiently
small $\delta>0$ at fixed $\kappa$. Perturb the second slice by
relabeling a total length $\delta/2$ from each halfcircle as
label 3. The unrestricted pairing tends to the boundary
maximum. On the consecutive side, the arc endpoints and
cyclic orders form a compact parameter space when empty arcs
are allowed. Any convergent sequence with masses tending to
$q^0$ thus has a consecutive limit at those masses. Continuity
of the pairing gives a strict upper bound below the boundary
maximum. Hence the exact consecutive-arc envelope near balanced
masses cannot extend to all interior mass vectors.

## Small squared radial deviation on the annulus suffices

The hypothesis in (JA1) can be weakened to
$\int_{r_J}^{R_J}|b(r)|^2d\mu(r)\le\delta_J$ for some
$\delta_J>0$. In particular any nonsimplex balanced
partition with stability at least the simplex has a positive
threshold of squared circular-mass discrepancy in this fixed
annulus, uniformly for $\rho\in J$.

Use the choices of $R,r_0$ in the radial proof, with the
exponential-in-$R^2$ and $O(\sqrt{r_0})$ errors also small
relative to $c_J/R$. Write $M=[r_0,R]$ and split $M$ into
$G=\{|b|\le\varepsilon\}\cap M$ and $B=M\setminus G$, where
$\varepsilon$ is the small fixed threshold for (JA4) at
$r,u\le R$. Chebyshev gives $\mu(B)\le\delta_J/\varepsilon^2$.
Put $X_G=\int_G r|b|^2d\mu$ and $X_B=\int_B r|b|^2d\mu$.

The $G\times G$ estimate is unchanged except for its missing
diagonal and mean. In the missing diagonal the additional
integral over $B$ costs $C_{J,R}\mu(B) X_G$, by the bounded
kernel and $k_0-k_\ell\le C\kappa$. The mean subtraction
costs at most $C\mu(B)X_B/r_0$ in addition to the inner
and outer terms, by Cauchy–Schwarz on $B$. Thus, for small
$\mu(B)$ at fixed $r_0,R$, the negative term on $X_G$
survives.

On the other pairs in $M\times M$ use Heilman's Lemma 6.1.
For the probability-normalized angular kernel $k_\kappa$,
it bounds the centered pairing minus the equal-arc centered
pairing by $-c\min(\kappa,1)(|b|^2+|c|^2)$ at every pair
of mass vectors. The prefactor in the lemma is a positive
weighted sum of functions $1-e^{-a\kappa}$, giving the
stated lower bound on the penalty. The raw pairing adds the
mixed mass term $b\cdot c$.
[Heilman, Lemma 6.1](https://arxiv.org/html/2306.03312v3)

For each $r\in B$, integrate that negative diagonal over
$u\in[1,R]$. The conditional probability of this interval is
at least $c_J$ uniformly for $r\le R$ when $R$ is large.
For the lower threshold use the independent one-dimensional
Gaussian coordinate perpendicular to the conditional mean,
whose variance is bounded below on $J$. For the upper
threshold use the exponential bound above, since the mean
has radius at most $(\sup J)R$. On this event
$\min(\kappa,1)\ge c_J r/R$. The negative term thus gives
at least $c_J X_B/R$. The possible positive mixed term on
$M^2\setminus G^2$ is at most
$C_{J,r_0,R}\sqrt{\mu(B)}(X_G+X_B)$, by Cauchy–Schwarz
with kernel weight $P_\rho(r,u)/\sqrt{ru}$, which is bounded
on this fixed annulus.

The outer and inner estimates depend only on total middle
energy, so replace $X$ there by $X_G+X_B$. First take $R$
large and then $r_0$ small as specified above. Then choose
$\varepsilon$ for (JA4) and $\delta_J$ sufficiently small
to absorb each term with a factor of $\mu(B)$ or
$\sqrt{\mu(B)}$. This leaves positive coefficients on
$X_G,X_B,Y,Z$ and proves (JA1) with the weaker hypothesis.

## Finite arc reduction for interior paired masses

Fix $0<k_-\le k_+<\infty$ and a lower bound $\delta>0$ on
each of the six masses in (JA2). Uniformly on these parameters,
there is a finite bound on the total number of arc components
of every ordinary maximizing pair. This gives a finite arc
endpoint description of the exact two-slice envelope on
compact interior parameter sets, allowing multiple arcs per label.

The key point is that a score contrast in either slice is
nonconstant even at a soft joint maximizer. Suppose, for example,
that $U_\kappa(v_i-v_j)$ is constant. The nonconstant Fourier
eigenvalues of $U_\kappa$ are positive, so $v_i-v_j$ is constant.
As both means are positive, the two components overlap on a
set of positive measure and are bounded below by a positive
constant on some positive-measure subset. Since $v$ maximizes
the linear functional with $u$ fixed, every bounded mean-zero
transfer between $v_i$ and $v_j$ in this subset has zero
first derivative. Thus $U_\kappa(u_i-u_j)$ is constant on the
subset and, by analyticity, everywhere. The same overlap
property then holds for $u_i,u_j$.

Choose bounded mean-zero transfers $h,k$ within those overlaps
of $u,v$. They can be chosen so that
$\langle h,U_\kappa k\rangle>0$. In fact start with any nonzero
$k$: its convolution is analytic and nonconstant. Its
difference from its mean on the other overlap, restricted
to that set, provides $h$ with strictly positive pairing.
Make both transfers with the same small parameter, preserving
the fixed masses and simplex-valued feasibility. Both linear
terms vanish by the constant-score identities. The mixed
quadratic term is $2\langle h,U_\kappa k\rangle$ times the
squared parameter, contradicting maximality.

To obtain a uniform bound, take a sequence of ordinary
maximizing pairs and extract a weak-star limit and convergent
parameters. This limit is a soft joint maximizer. The kernel
pairing is jointly weak-star continuous, and the feasible mass
vectors can be adjusted continuously in $L^1$ by transferring
the appropriate fractions of each surplus label to the deficit
labels. In each original slice use the bounded centered
assignment prices from the proof of (JA2), passing to a
convergent subsequence. The adjusted contrasts converge
uniformly on compact complex neighborhoods of the real circle,
since $e^{\kappa\cos(\theta-\varphi)}$ is entire in $\theta$
and the convolution kernels and their derivatives are uniformly
bounded on compact complex sets. The limiting contrasts are
nonzero analytic functions by the soft maximum argument. Each
has finitely many zeros on a period. In small complex disks
around those zeros the number of zeros in the approximants
is eventually constant by Rouché's theorem. Off the disks
there are no real zeros for large indices. Since a winner can
switch only at a contrast zero, this rules out unbounded arc
counts.

On the mass boundary there is a simple exact formula. If
$q_j=0$, the envelope is $\sum_i B_\kappa(p_i,q_i)$.
Take the two remaining arcs of the second slice complementary,
with opposite centers. Center the corresponding two arcs of
the first slice at those centers. They are disjoint because
their widths sum to at most one. The remaining points in the
first slice are assigned to $j$. This saturates both individual
rearrangement bounds. Changing a mass vector by total $\ell^1$
distance $d$ changes the envelope by at most $d/2$ in either
argument. Relabel a set of total circular mass $d/2$ from
surplus colors to deficit colors and use that the convolved
score difference at a relabeled point has absolute value at
most one. Apply this in each direction to compare optima.

### The circular kernel can increase sign changes

A uniform arc-count argument based on variation diminution of
$U_\kappa$ would require additional structure. For every
$\kappa>1/2$ the kernel has an explicit failure of this property.
Choose $\kappa<\lambda<\kappa^2+1/4$ and take the signed input
concentrated at $-t,0,t$ with respective weights
$1,-(2+\lambda t^2),1$. The convolved function divided by
$t^2 k_\kappa(\theta)$ converges as $t\to0$ to
\[
 \frac{k_\kappa''(\theta)}{k_\kappa(\theta)}-\lambda
       =-\kappa\cos\theta+\kappa^2\sin^2\theta-\lambda .
\]
It is negative at $0$ and $\pi$, and positive at the two
angles with $\cos\theta=-1/(2\kappa)$. For small nonzero $t$
the output has at least four cyclic sign changes. Replace
the three atoms by smooth disjoint bumps with the same signed
masses, narrow enough to preserve the four strict output
signs. The input then has two cyclic sign changes. This
example concerns signed convolution, independent of the
mass-constrained optimizer.

## Compatibility of exact paired maximizers across radii

Fix $0<\rho<1$. Suppose the angular slices of an ordinary planar
partition attain the exact paired envelope for almost every
pair of positive radii. Then almost every slice consists of a
single arc for each positive label, with empty or full arcs
allowed. In addition, the sets of radii with exactly two and
exactly three positive masses cannot both have positive radial
measure.

At almost every radius $r$, the slice is a Lebesgue point in
$L^1(S^1)^3$ on compact intervals of positive radii. For almost
every such $r$ the paired equality also holds for almost every
partner radius. Choose partner radii tending to $r$ along which
the slices converge in $L^1$ to that at $r$. The pairing is
continuous under this limit, since the kernels converge
uniformly. The envelope has the same continuity by its mass
Lipschitz bound and the uniform kernel convergence. Thus the
slice at $r$ attains its self-envelope. Consecutive arcs with
the prescribed masses attain the sum of the scalar self bounds
for any mass vector. Equality requires every nontrivial label
to be an arc by the positive cosine component in the proof of
(JA2).

For a two-label slice paired with a three-label slice, the
boundary envelope formula in the finite arc reduction shows
that equality must attain both scalar upper bounds. Again
subtract the positive cosine component. For each matched label,
equality of the cosine term requires that its arc centers agree
in the two slices, since both first moment norms are strictly
positive. In a two-label slice the two arc centers are
antipodal. In a three-label arc partition the distance between
the centers of any two arcs is $\pi(p_i+p_j)<\pi$, by measuring
along their common endpoint. Thus no pair of these two types
attains the envelope. Positive radial measure of both types
would contradict the assumed paired equality almost everywhere.

## Sharp circular mass penalty for the first moment

For every $p\in\Delta_3$, put $B(p)=\sum_i(p_i-1/3)^2$.
Then
\[
 \sum_i\sin^2(\pi p_i)\le \frac94-\frac32 B(p).
 \tag{JA7}
\]
Equality holds at the equal vector and at permutations of
$(0,1/2,1/2)$, so the coefficient of $B$ is optimal. To prove
it, use $\sum_i\sin^2(\pi p_i)=2+2\prod_i\cos(\pi p_i)$.
When all $p_i\le1/2$, put $q_i=1/2-p_i$, with $\sum q_i=1/2$.
The desired bound on the product becomes
$\prod_i\sin(\pi q_i)\le(3/2)\sum_{i<j}q_iq_j$.
The logarithm of $\sin(\pi q)/q^{2/3}$ is strictly concave
on $(0,1/2)$, since its second derivative is
$-\pi^2/\sin^2(\pi q)+2/(3q^2)<0$. Jensen's inequality gives
$\prod_i\sin(\pi q_i)\le(9/2)(q_1q_2q_3)^{2/3}$; apply the
arithmetic-geometric mean inequality to the three products
$q_iq_j$, and take boundary limits. When $p_1>1/2$, put
$s=1-p_1<1/2$. The cosine addition formula gives
$\cos(\pi p_2)\cos(\pi p_3)\ge\cos(\pi s)$, so the cosine
product is at most $-\cos^2(\pi s)\le-4(1/2-s)^2$. Also
$B(p)\le2/3-2s+2s^2$, which gives
$1/8-3B(p)/4\ge-(3/2)(1/2-s)^2$. This proves the bound also
in this case, with strict inequality.

The norm of the first Gaussian moment vector of an ordinary
planar partition satisfies
\[
 \sqrt{M_1(A)}\le\frac{3}{\sqrt{8\pi}}-
          \frac{1}{2\pi}\int r B(p(r))d\mu(r),\qquad
 M_1(A)=\sum_i\left|\int_{A_i}x\,d\gamma_2(x)\right|^2 .
 \tag{JA8}
\]
Indeed the per-circle moment caps are $\sin(\pi p_i)/\pi$.
The triangle inequality applied to the vector of the three
integrated moments bounds its norm by the integral of the
per-circle norm. Use (JA7) and the linear upper tangent to
the square root at $9/(4\pi^2)$. The reference moment is the
integral of the constant per-circle norm $3/(2\pi)$ times $r$.

For balanced partitions, (JA8) also gives an explicit exclusion
of some large radial deviations at intermediate correlations.
Let $B_0=\int B(p(r))d\mu(r)$, $W=\int r B(p(r))d\mu(r)$,
and put $a=\log(3/2)$, $b=\log 3$ and
$M=[4a^2+(2a-b)^2+b^2]/9$. Orthogonality of the radial and
angular-mean-zero parts of the indicators gives
\[
\begin{split}
 S_\rho(A)\le {}&\tfrac13+\rho^4 B_0+(\rho^2-\rho^4)\min(M,B_0)\\
 &+\min\left\{\rho(2/3-B_0),\ \rho^2(2/3-B_0)
       +(\rho-\rho^2)\left(\frac3{\sqrt{8\pi}}-\frac{W}{2\pi}\right)^2\right\} .
 \tag{JA9}
\end{split}
To bound the radial degree-two coefficient, use $t=r^2/2$
and $m_i=\int(t-1)p_i\,d\nu$. In the scalar product
$\sum_i m_i\int(t-1)p_i\,d\nu$, hold the weights $m_i$ fixed
and sort them. The optimal assignment allocates the three
ordered quantile intervals with endpoints $a,b$ to the
ordered weights. One proof writes the sorted weights as a
constant plus two nonnegative increments times upper
label sums, and maximizes each integral against the
increasing function $t-1$ by the upper tail of its prescribed
mass. The squared norm of the resulting moment vector is
$M$, so Cauchy–Schwarz gives $\sum_i m_i^2\le M$. It is also
at most $B_0$ by projection. All remaining nonconstant radial
Hermite degrees are at least four. The angular-mean-zero
part has total squared norm $2/3-B_0$, first-degree weight
$M_1(A)$, and all remaining degrees at least two. Apply the
corresponding eigenvalue bounds and (JA8) to obtain (JA9).
