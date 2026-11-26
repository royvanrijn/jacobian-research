# All-degree angular separation from opposite slots along matched tanh rays

Fix any finite ordinary two-monomial core. Additional rows can have
either orientation, positive integer gaps \(d_i\), and smaller bases
\(\ell_i(T)=a_i+s_iT\), with integers \(a_i\ge0,\ s_i>0\). At terminal
\(N=n+1\ge3\) let
\[
  \tau_N=\max_i-\frac{Na_i+1}{Ns_i}.
\]
The first-block cohort consists of the moving rows attaining this
maximum. Define the lower endpoint slot condition at this terminal
by requiring one Gauss slot in each orientation with base
\[
                       0<b_{n,i,h}(\tau_N)\le1.
\]
The marked upper slots always provide such a lower slot in their
own orientations. A private row of the fixed core provides one
in its orientation as well.

**Theorem.** Fix a terminal satisfying the lower endpoint slot
condition. For all but finitely many real \(T>\tau_N\), the exact
two arrays at orders \(n,N\) have disjoint root arguments. In
particular, if the condition holds at every terminal, all but
finitely many ordinary parameters \(T\in\mathbb Z_{\ge0}\) have
disjoint consecutive root arguments in every degree. The
exceptional set is effective given the
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).

A first-block cohort containing both orientations suffices,
regardless of the relative counts in those orientations. With a
private variable on each side in the fixed core, every affine ray
satisfies the condition.

For example, an arbitrary fixed core can be made angularly separated
in every degree by adjoining the two opposite shared rows
\[
   (T+1,T+1+d_+),\qquad (T+1+d_-,T+1),\qquad d_+,d_-\ge1
\]
for any sufficiently large integer \(T\). The gaps can differ and the
core requires no symmetry, curvature cap, or initial angular
separation. The appended curvature tends to zero. The choice
\(d_+=d_-=1\) appends just two matched slots.

## The exact adjacent endpoint

For a Gauss slot \(1\le h\le d_i\) the rising factorial base at
order \(m\) is \(b_{m,i,h}=(m\ell_i(T)+h)/d_i\). Its tanh coordinate
is \(T_{m,i,h}=(m-1)/(2b_{m,i,h}+m-1)\). At \(T=\tau_N\) a marked
first slot has \(T_{N,i,1}=1\) and \(b_{n,i,1}=1/(Nd_i)\).
The upper and lower values come from the same moving row.
[Exact matched tanh coordinates](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)

Write \(x=T-\tau_N\). For each \(m=n,N\), strip the smaller-end
factorials from the order-\(m\) array. The resulting polynomial
\(F_m(x,z)\) has coefficients in \(\mathbb Q[x]\); its factors at
index \(k\) are the rising factors \(m\ell_i(T)+h\),
\(1\le h\le d_i k\) or \(1\le h\le d_i(m-k)\), along with the
binomial coefficient. Let \(u,v\) count the marked increasing
and decreasing rows, respectively. At \(x=0\) the orders of
vanishing of the upper coefficients are \(v\) at index zero,
\(u\) at index \(N\), and \(u+v\) in the interior. Their first
nonzero coefficients are positive. For \(x>0\), rescale the
upper argument positively and divide out a scalar to write
\[
                   x^{-v}F_N\!\left(x,x^{(v-u)/N}z\right).
\]
Both endpoint coefficients now have positive limits. At an
interior index \(k\) the remaining power of \(x\) is
\(u(1-k/N)+vk/N>0\). Thus the upper endpoint array converges
to \(A+Bz^N\) with \(A,B>0\). Quotients of roots are invariant
under this scaling.

All lower factors at zero are positive. Reserve one slot of base
\(b\in(0,1]\) from each orientation, choosing a marked slot
wherever one is available. The ratios \((b)_j/j!\) are moments
of a \(\mathrm{Beta}(b,1-b)\) variable for \(b<1\) and of the
constant one for \(b=1\). They absorb the two binomial factorials.
The other rising factors are positive Gamma moments. After a
positive tilt and scalar normalization the lower endpoint has
coefficients \(\mathbb E X^k\), \(0\le k\le n\), with \(X>0\).
One reserved slot is marked, with base \(1/(Nd_i)<1\), so its
independent Beta factor makes \(X\) nonconstant. The same lower
specialization appears in the
[physical affine-ray theorem](MULTIVARIATE_AFFINE_BASE_RAY_ALL_DIRECTIONS.md#balanced-cohorts-have-a-positive-geometric-specialization).

In particular each lower upper-half-plane root has argument
\[
  \theta_j=\frac{2j\pi}{N}-\delta_j,\qquad
  0<\delta_j<\frac{(N-2j)\pi}{nN}<\frac{\pi}{N}.
  \tag{1}
\]
This is the strict
[geometric moment angle box](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry).
The quotient of an upper and a lower half-plane root of the
lower endpoint, in that order, has argument congruent to
\(\theta_j+\theta_k\) for some \(j,k\). By (1) the sum is strictly
between two consecutive multiples of \(2\pi/N\). Every quotient
of two roots of the upper endpoint is an \(N\)-th root of unity.
Thus the two endpoint arrays have disjoint **opposite-half-plane
root quotients**.

## Carrying the half-plane restriction through the ray

For real \(x>0\), all moment factors and both endpoint coefficients
are positive. The moment integral remains valid for the real
exponents on this entire interval. Indeed each smaller exponent
satisfies \(1+m\ell_i(T)>0\) for \(m=n,N\). For independent unit
exponentials, the two monomial values \(U,V>0\) therefore have all
mixed moments of total order \(m\) finite. On the real axis
\(\mathbb E(U+zV)^m>0\) when \(m\) is even, since \(U/V\) is
nonconstant. When \(m\) is odd the derivative
\(m\mathbb E[V(U+zV)^{m-1}]\) is strictly positive at every real
\(z\). The leading coefficient then gives exactly one real root,
negative and simple. These conclusions are unchanged by the
positive scalar normalizations of the two arrays.

For \(m=n,N\), let
\[
   D_m(x,u)=
   \prod_{\substack{F_m(x,\alpha)=F_m(x,\beta)=0\\
                   \Im\alpha>0,\;\Im\beta<0}}
           (u-\alpha/\beta),
   \qquad
   H(x)=\operatorname{Res}_u(D_n(x,u),D_N(x,u)).
   \tag{2}
\]
Count roots with multiplicity. Use the positively rescaled upper
array only for the limit as \(x\) decreases to zero. The degrees in \(u\) are fixed on \(x>0\).
The coefficients are real-analytic there, with complex values
allowed. At any given point, contours separate the two open
half-plane root sets from the real root, if present. The
logarithmic derivative of the original polynomial integrated on
those contours gives analytic power sums separately in the two
sets. This proves the asserted analyticity even at repeated
nonreal roots. The coefficients and \(H\) are algebraic over
\(\mathbb Q(x)\), since they are formed from algebraic roots by
finite rational operations.

Both limiting polynomials have simple roots, including the real
root in the odd degree. The nonreal roots of the rescaled arrays
converge within their respective half-planes as \(x\downarrow0\).
Thus \(H(x)\) has a nonzero limit. In
particular \(H\) is a nonzero analytic algebraic function on the
connected positive axis. A minimal polynomial relation over
\(\mathbb Q(x)\), cleared of denominators, has a nonzero
constant coefficient in its \(H\) variable. Every zero of
\(H\) on the ray must be a root of that fixed nonzero polynomial
in \(x\), so there are finitely many.

At a nonreal common argument the two roots give the same
quotient with their own conjugates, after conjugating both
roots if needed. This forces \(H(x)=0\). The positive real
axis is excluded by coefficient positivity, and the even
moment has no negative real root. This proves the
fixed-terminal assertion under independent positive rescalings.

There is also a degree bound for the number of exceptional points.
Write \(p=\lfloor n/2\rfloor,\ q=\lfloor N/2\rfloor,\ K=p^2q^2\)
and \(G=\sum_i d_i\) over the moving rows. For each assignment of
disjoint ordered pairs of subsets of sizes \(p,p\) and \(q,q\) to
the respective root multisets, form the numerator
\[
              R_\sigma=\prod
                         (\alpha_+\beta_- -\alpha_-\beta_+).
\]
The product runs through the four assigned subsets, with \(\alpha\)
the lower-order roots and \(\beta\) the upper-order roots.
The actual half-plane choice gives the numerator of (2), up to
sign. Let
\(M=\binom np\binom{n-p}p\binom Nq\binom{N-q}q\)
be the number of assignments. The coefficients of
\(\prod_{\sigma}(Y-R_\sigma)\) are symmetric separately in the
two root multisets. The coefficient of \(Y^{M-j}\) has degree
\(jK\) in each multiset. In the corresponding elementary
symmetric polynomials its weighted degree is \(jK\), so multiplying
by each leading coefficient to the power \(MK\) clears all
denominators. Every polynomial coefficient so obtained has
\(x\)-degree at most \(MK(n+N)G\), since each original coefficient
in order \(m\) has \(x\)-degree at most \(mG\).
Remove the largest common power of \(Y\) from the relation.
The actual analytic numerator still satisfies the reduced
relation because it is nonzero near zero. The reduced constant
coefficient is a nonzero polynomial and vanishes at every
exceptional point. Thus the number of exceptional real
parameters at this terminal is at most \(MK(n+N)G\).

For effectivity at a bounded terminal one may compute the true
exceptional set directly by real quantifier elimination. Ask for a
nonreal root \(u+iv\), \(v>0\), of the lower array and a root
\(\lambda(u+iv)\) of the upper array, with real \(\lambda>0\).
The coefficients of all four real and imaginary equations are
rational polynomials in \(x\). The resulting subset of the
positive axis is finite by (2). Terminal two separates by real
positivity of the quadratic. The absolute ordinary cutoff applies
to every higher terminal beyond its fixed computable bound
independently of \(T\). Thus only finitely many exceptional
ordinary parameters remain across all degrees.

This result uses the exact matched slots at a discrete-row
specialization. It leaves open angular separation at an arbitrary
fixed profile. The [physical affine-ray theorem](MULTIVARIATE_AFFINE_BASE_RAY_ALL_DIRECTIONS.md)
excludes physical common roots outside a finite set on every affine ray.
