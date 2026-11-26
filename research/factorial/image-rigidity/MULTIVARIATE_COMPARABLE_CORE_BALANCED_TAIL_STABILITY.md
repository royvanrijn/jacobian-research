# All-degree stability of comparable shared cores under balanced ordinary tails

Start with any ordinary two-monomial support with a private
variable on each side. Assume the exponent columns on its
shared variables are coordinatewise comparable. Append
arbitrarily many ordinary shared rows \((a_i,b_i)\), with
\(a_i,b_i\ge1\), satisfying
\[
           \sum_{\rm tail}(b_i-a_i)=0,\qquad
           \tau=\sum_{\rm tail}\frac{(b_i-a_i)^2}{\max(a_i,b_i)}.
\]

**Theorem.** There exists an absolute \(\tau_0>0\) such that
\(\tau<\tau_0\) gives consecutive coprimality in every degree
over every characteristic-zero field. The tolerance is
independent of the number and exponents of the core and tail
rows, and of the total degree difference of the core.
The appended rows may have either orientation. Thus the full
shared support can have genuinely mixed orientations.

The positive tolerance here is unevaluated. The proof uses a
computable absolute eventual cutoff and a compactness argument
at the remaining orders. It does not evaluate the compactness
margin or test individual terminal polynomials.

## Uniform normalized endpoint data

Suppose the assertion fails along a sequence of tails with
\(\tau\to0\), and choose a shared physical root in each example.
The [absolute eventual theorem](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
bounds the terminal \(N=n+1\) in these examples by one absolute
integer. Pass to a subsequence with \(n\) fixed. The terminal
\(N=2\) already separates by positivity of the real quadratic
moment, so \(n\ge2\). The
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md)
also lets us assume that the total curvature of each full
support is below \(27\), hence the same bound holds for its core.

For each core reserve one private factorial on each side,
as in the
[both-private geometric representation](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry).
Let \(U,V\) be the resulting positive bases. For every private
row one may group its \(d\) factors with Gamma shapes
\(h/d\), \(1\le h\le d\), deleting the shape-one factor at
each of the two reserved rows. Thus the residual private
factors have exponent one on their respective sides.
Their number is below \(27\) and their shapes are at
least \(1/27\). Deterministic positive scales can be
included in \(U,V\). Put \(T=V/U\) and let
\(\mathbb E_j\) denote the probability tilt by \(U^j\).
Write
\[
                   g_j=(\mathbb E_j T^j)^{1/j}
\]
for the core normalization, so the core normalized polynomial
has coefficients \(\mathbb E_j(T/g_j)^k\), \(0\le k\le j\).
Orient the core so its shared gaps \(d_i=b_i-a_i\) are
nonnegative.

We will use uniform moments in a fixed neighborhood of
the two grids. Set \(\sigma=1/128\). For
\(-\sigma\le s\le j+\sigma\), write
\[
                  F_j(s)=\log\mathbb E_j T^s .
\]
All arguments of its Gamma factors are positive. For each
shared unequal row \(a_i\ge1\), the bound on total curvature
gives \(d_i/a_i\le1+d_i^2/b_i\le28\) and hence \(b_i/a_i\le29\).
Thus \(j a_i+s d_i\ge j a_i/2\), and its contribution to
\(F_j''(s)\) is at most \(2d_i^2/(j a_i)\), by
\(\psi'(1+x)\le1/x\) for \(x>0\).
The sum of these bounds is at most \(2\cdot29\cdot27/j\).
Each residual private argument is at least
\(1/27-\sigma\); the estimate
\(\psi'(x)\le x^{-1}+x^{-2}\) bounds its contribution
by \(1300\). The two trigamma estimates follow by
comparison of the series with an integral, keeping at
most the first term separately. Consequently
\[
          0\le F_j''(s)\le C:=40000
                         \qquad(-\sigma\le s\le j+\sigma).
\]
The centered log moment function agrees with zero at
\(s=0,j\). The interpolation remainder at these two
points gives
\[
 \left|\log\mathbb E_j(T/g_j)^s\right|
              \le\frac C2\,|s(s-j)|
                         \qquad(-\sigma\le s\le j+\sigma).     \tag{1}
\]
This bound is independent of the private exponent choices.

## The balanced tail disappears at both orders

Write \(D_j\) for the tail normalization, so
\[
                 \log D_j=\sum_{\rm tail} d_i
                      \int_0^1\psi(1+j(a_i+u d_i))\,du ,
                  \qquad d_i=b_i-a_i .
\]
All normalized tail coefficients tend to one for the fixed
orders \(n,N\). Indeed on the continuous index interval
\([0,j]\) the logarithmic second derivative of the product
of tail Gamma factors is at most
\(\sum_{\rm tail} d_i^2/(j\min(a_i,b_i))
\le\tau/[j(1-\tau)]\) for \(\tau<1\).
The last step uses \(|d_i|/\max(a_i,b_i)\le\tau\).
The endpoint chord bound then gives convergence to the
constant normalized log array.

The two tail normalizations also match asymptotically.
Differentiation in a real order \(j\) and the comparison
\[
                 0\le \frac1j-x\psi'(1+jx)
                           \le\frac{1}{j^2x}\qquad(x>0)
\]
give
\[
        \left|\log\frac{D_N}{D_n}\right|
                  \le\frac{\tau}{nN(1-\tau)}.                 \tag{2}
\]
Indeed the sum of the \(d_i/j\) terms vanishes by balance,
and \(\sum_{\rm tail}|d_i|/\min(a_i,b_i)\le\tau/(1-\tau)\)
by integrality of the unequal gaps.

The full normalized polynomials have positive
coefficients at most one and both endpoints one, by
Hölder after private cancellation. Hence their roots have
radii in \((1/2,2)\), by comparing the leading term
against the other terms for \(|z|\ge2\) and reversing.
For the alleged physical common root \(z\) put
\(w=g_nD_nz\), \(w'=g_ND_Nz\).
Pass to a subsequence with \(w\to w_*\ne0\) and
\(g_N/g_n\to c_*>0\); the ratios are bounded above and
below because \(w'/w=g_ND_N/(g_nD_n)\) and (2) holds.
Tail coefficient convergence shows that the normalized
lower core polynomial at \(w\) tends to zero. The normalized
upper core polynomial at \((g_N/g_n)w\) also tends to zero:
it differs from evaluating at \(w'\) by an argument change
tending to zero, and its coefficients are uniformly bounded.

## The limiting pair retains the conditional order

For each core along the sequence, let \(\mu\) be the
law of \(X=T/g_n\) under \(U^n\), and let \(\nu\) be its
law under \(U^nV\), with both tilts normalized.
The
[conditional monotonicity for comparable cores](MULTIVARIATE_BOTH_PRIVATE_COMPARABLE_SHARED_INTERLACING.md#why-comparable-shared-ordinary-rows-satisfy-the-test)
makes \(d\nu/d\mu\) nondecreasing.
Equivalently, \(\nu\otimes\mu-\mu\otimes\nu\) is a
nonnegative measure on the ordered region \(x>y\).
This formulation of the likelihood-ratio comparison
survives weak convergence, by testing against nonnegative
continuous functions supported in that open region.

Under \(U^N\) put \(b=\mathbb E_N X\).
By (1) and the bounds on \(g_N/g_n\), these numbers
are bounded above and bounded away from zero, and
the \(N+\sigma\) moments of \(X\) under \(U^N\)
are uniformly bounded. The measure \(\nu\) is also
the \(X\)-size-biased law under \(U^N\), since
\(U^n V=U^N T\). Therefore both \(\mu,\nu\) have
uniformly bounded \(n+\sigma\) moments.
The law \(\mu\) has uniformly bounded negative
\(\sigma\)-moment by (1). Pass to a subsequence with
weak limits \(\mu_*,\nu_*\) and \(b\to b_*>0\).
The negative moment bound excludes an atom at zero
for \(\mu_*\). The positive moment bounds give
uniform integrability for all powers through order \(n\).
Consequently the two vanishing limits in the previous
section give
\[
                    \mathbb E_{\mu_*}G_n(w_*X)=0,\qquad
                 1+w_* b_*\,\mathbb E_{\nu_*}G_n(w_*X)=0,      \tag{3}
\]
where \(G_n(s)=1+s+\cdots+s^n\).
The upper equality before passing to the limit
uses \(G_N(s)=1+sG_n(s)\) and the size-biased law.

Conjugate the roots to the upper half-plane. By the
[angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
a shared angle must have the same positive label \(k\)
in the two orders. Pass to a further subsequence with
fixed \(k\). Its limit \(\theta=\arg w_*\) satisfies
\[
        0<\frac{(2k-1)\pi}{n}\le\theta
                         \le\frac{2k\pi}{N+1}\le\frac\pi2 .
\]
In particular \(A=-\sin(N\theta)>0\) and
\(B=-\sin(n\theta)\ge0\). With \(q=\sin\theta>0\),
the imaginary part of the lower geometric sum is
\[
  I(s)=\operatorname{Im}G_n(se^{i\theta})
       =\frac{s[q+s^n(A-Bs)]}{1-2s\cos\theta+s^2}.
\]
If \(B=0\), then \(I>0\) for \(s>0\), contradicting
the first equality of (3) and the absence of an atom
at zero. Otherwise \(I\) is positive up to one
positive crossing and negative thereafter.
Write \(I_\pm\) for its positive and negative parts,
applied to \(s=|w_*|X\). The first equality of (3)
gives \(\int I_+\,d\mu_*=\int I_-\,d\mu_*=a_*\).
Here \(a_*>0\): otherwise \(\mu_*\) would be supported
on the single crossing, with the atom at zero already
excluded. The first equality would then make the angle
a lower geometric-sum root angle, and the lower box
permits only \(2k\pi/N\), above the displayed upper bound.

On the ordered region \(x>y\), test the limiting
likelihood-ratio inequality with \(I_-(x)I_+(y)\).
One may first use bounded continuous cutoffs supported
there and then pass by monotone convergence. This gives
\[
                a_*\int I_-\,d\nu_*\ge a_*\int I_+\,d\nu_* .
\]
Thus \(\operatorname{Im}\mathbb E_{\nu_*}G_n(w_*X)\le0\).
The second equality in (3) makes this imaginary part
\(\operatorname{Im}(-1/(w_* b_*))>0\), a contradiction.

The contradiction rules out sequences of failures with
\(\tau\to0\) over all core and tail profiles. It proves
the stated absolute positive tolerance. The original
coefficient arrays are rational and the positive
normalizations were used only to test complex roots,
so the rational Bézout identity gives every
characteristic-zero field.
