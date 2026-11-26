# The universal Gaussian clock for adjacent circular interlacing

For a polynomial \(P_m(z)=\sum_{k=0}^m c_k z^k\)
write its angular heat flow at time \(s\ge0\) as
\[
          {\cal G}_{m,s}P_m=
             \sum_{k=0}^m c_k e^{-s k(m-k)/2}z^k .
\]
The physical degree normalization uses
\({\cal H}_{m,t}={\cal G}_{m,t/m}\), so
\({\cal H}\) multiplies the normalized coefficients by
the common symbol \(e^{mL_t(k/m)}\) with
\(L_t(q)=-t q(1-q)/2\).

**Theorem.** For every \(n\ge1\), initially strictly
interlacing positive palindromic polynomials with simple
unit-circle roots at degrees \(n,N=n+1\) keep strict
interlacing under \({\cal G}_{n,s}\) and \({\cal G}_{N,s}\)
at equal angular times \(s\).
For every \(n\ge4\) and distinct positive rates \(c_n,c_N\),
there exist such polynomials, also with strictly
log-convex binomial-normalized coefficients, for which
\({\cal G}_{n,c_nt}\) and \({\cal G}_{N,c_Nt}\) give a
common root at some \(t>0\).
In particular the universal preservation fails for
\({\cal H}_{n,t}\) and \({\cal H}_{N,t}\).
The counterexamples concern abstract coefficient arrays.

## The angular root flow

Near a simple unit-circle root, a positive real
palindromic input continues to have its root on
the circle under sufficiently small real smoothing.
Indeed reciprocal conjugation preserves the root
set and the nearby root is unique.
For each continued root angle \(\theta_j(t)\),
the backward heat equation at angular time \(s\) gives
\[
              \frac{d\theta_j}{ds}=\frac12
                         \sum_{\ell\ne j}
                         \cot\frac{\theta_j-\theta_\ell}{2}.
\]
To see this, center the circle value by
\(e^{-im\theta/2}\) and multiply by \(e^{s m^2/8}\).
The resulting trigonometric polynomial satisfies
\(\partial_s f=-f''/2\).
At a simple root \(f''/f'\) is the sum of the
indicated cotangents by factoring the centered
polynomial into half-angle sine factors.

In the upper angular half put \(x_j=\cos\theta_j\),
excluding the odd structural root at \(-1\).
Pairing the terms from \(\theta_\ell\) and
\(-\theta_\ell\) gives
\[
 \frac{d x_j}{ds}=(1-x_j^2)\sum_{\ell\ne j}\frac1{x_j-x_\ell}
                  -\frac{x_j-\mathbf{1}_{m\ {\rm odd}}(1-x_j)}2 .
                                                               \tag{1}
\]
The single reflected root contributes \(x_j\)
to the last numerator after multiplication by
\(\sin\theta_j\). The odd root contributes
\(-(1-x_j)\).

The flow stays simple and on the circle for
every finite \(s\). On each interval of simple
roots the angular equation is the gradient of
\[
              \sum_{j<\ell}\log\left|
                  \sin\frac{\theta_j-\theta_\ell}{2}\right|.
\]
The potential increases and cannot tend to
minus infinity, whereas a collision would
force that limit. Thus local simple reciprocal
root continuation extends for finite times.
Positive coefficients exclude the positive cut
throughout. The odd root stays at \(-1\), and an
even polynomial cannot have a simple root there
by reciprocity.

## Equal angular times

Use decreasing cosine coordinates \(A_i\) for
the lower positive-half roots and \(B_i\) for
the upper positive-half roots, excluding the
structural odd root. At a hypothetical first
contact \(x\), the function \(q(y)=1/(x-y)\)
is strictly increasing on either side of
\(x\). In (1) compare the two sums over roots
other than the contacted pair.

If \(n=2p\), the ordering is
\(B_1>A_1>\cdots>B_p>A_p>-1\), with the upper
structural root at \(-1\).
At a left contact \(B_k=A_k=x\), each uncontacted
\(B_i\ge A_i\) lies on the same side of \(x\),
so the upper sum is at least the lower one.
Its additional endpoint drift \((1-x)/2\)
makes \(d(B_k-A_k)/ds>0\).
At a right contact \(A_k=B_{k+1}=x\), pair
each uncontacted \(A_i\), \(i<p\), with
\(B_{i+1}\). Each lower summand is at least
its upper partner. The extra pair \(A_p,B_1\)
contributes to the lower-minus-upper sum at least
\[
                      q(A_p)-q(B_1)
                        \ge\frac{1}{x+1}+\frac{1}{1-x}
                        =\frac{2}{1-x^2}.
\]
After multiplying by \(1-x^2\) this dominates
the adverse endpoint drift \((1-x)/2\).
Thus \(d(A_k-B_{k+1})/ds>0\).

If \(n=2p+1\), the upper has \(p+1\) positive-half
cosines, alternating with the \(p\) lower ones
before the lower structural root.
At a left contact pair equal indices and use
the extra upper \(B_{p+1}<x\).
It contributes at least \(1/(x+1)\)
to the upper-minus-lower sum.
Multiplication by \(1-x^2\) gives at least
\(1-x\), greater than the adverse lower
endpoint drift \((1-x)/2\).
At a right contact pair \(A_i\) with \(B_{i+1}\)
and use the extra upper \(B_1>x\), whose
summand is negative. The lower sum is
larger and the endpoint drift favors
the lower root as well.
The case \(p=0\) has just the lower
structural root and is immediate.
In every possible first-contact case the
vanishing gap has strictly positive derivative,
which is incompatible with the gap approaching
zero from above. This proves equal-time preservation.

## Unequal angular times

Write \(\chi=c_N/c_n>0\).
For any fixed \(n\ge4\) there are at least two
lower positive-half roots. Build a left contact
\(B_k=A_k=x\), keeping the other required
inequalities strict. When \(\chi<1\), put the
next lower and upper roots at
\(x-e\) and \(x-(1-\nu)e\), with
\(0<\nu<1-\chi\) and \(e>0\) small.
Take \(k=1\). In the upper-minus-lower
velocity, divided by \(c_n\), these two
force contributions from (1) give
\[
                \frac{1-x^2}{e}
                       \left(\frac{\chi}{1-\nu}-1\right)<0.
\]
When \(\chi>1\), use a preceding lower
root at \(x+e\) and its upper partner at
\(x+(1+\nu)e\), with \(0<\nu<\chi-1\).
Take \(k\) to be the last lower index.
The contribution is now
\[
                \frac{1-x^2}{e}
                       \left(1-\frac{\chi}{1+\nu}\right)<0 .
\]
In each case keep all other roots in
the prescribed alternating order a fixed
distance from \(x\), as \(e\) decreases.
Their contributions and the endpoint
drifts stay bounded. Thus the left
gap has strictly negative derivative
at the chosen contact.

The construction can have all its
cosines close to a fixed negative
\(x\) close to \(-1\).
The repeated-cosine input is a power of
\(z^2-2xz+1=(1+z)^2-u z\),
\(u=2(1+x)>0\), times the odd factor
if required. Its binomial-normalized
logarithms have strictly positive
second differences at small \(u\),
by the
[repeated-root coefficient expansion](MULTIVARIATE_LINEAR_CAYLEY_COMPRESSION_OBSTRUCTION.md#compression-on-the-parity-spaces).
Positivity and these strict inequalities
persist in a neighborhood of the
repeated inputs. Choose the indicated
distinct cosine locations there.
Finally move only \(B_k\) a sufficiently
small positive amount to make the initial
pair strictly interlacing. By the analytic
implicit function theorem at the transverse
contact, this pair has a common cosine
at a positive time. This proves the rate
obstruction in each adjacent degree.

## An exact instance at the physical normalization

Let \(x=-3/4\) and \(Q_a(z)=z^2-2az+1\).
Use the degree-four polynomial \(Q_xQ_{-19/25}\).
At the boundary collision use the degree-five
polynomial \((1+z)Q_xQ_{-759/1000}\).
Both have simple unit-circle roots.
Their binomial-normalized sequences are
\[
 \begin{split}
 (1,151/200,107/150,151/200,1),\\
 (1,2009/2500,1459/2000,1459/2000,2009/2500,1).
 \end{split}
\]
Strict log convexity follows by taking successive
cross-products minus the squared middle term.
For the first sequence the distinct internal
margins are \(17197/120000\) and \(881/14400\).
For the second they are \(261647/3125000\)
and \(1081119/20000000\).

At the common value \(x=-3/4\), formula (1)
at the physical rates \(c_4=1/4,c_5=1/5\) gives
\[
 \dot x_5-\dot x_4
  =\frac7{16}\left(\frac{200}{9}-25\right)
       -\frac{2x-1}{10}+\frac{x}{8}
  =-\frac{305}{288}.                                         \tag{2}
\]
The two positive-angle roots involved are simple
in their respective polynomials.

## An initially strict pair with an exact later contact

Replace \(Q_x\) in the degree-five polynomial
by \(Q_{x+\epsilon}\), with \(\epsilon>0\) sufficiently
small and rational. The upper-half cosine roots
now satisfy
\[
                 x+\epsilon>x>-759/1000>-19/25>-1,
\]
where the last endpoint is the structural
root of degree five. Reflection gives strict
alternation around the whole circle cut at \(1\).
Positivity and the strict logarithmic margins
persist by continuity.

The difference between the continued cosine
roots from \(x+\epsilon\) and \(x\) is analytic
near \(t=\epsilon=0\), by simplicity.
At \(t=0\) the difference is \(\epsilon\).
At the boundary collision its \(t\)-derivative
is (2). The implicit function theorem thus
gives an exact common cosine, hence an exact
common circle root, at
\[
                        t(\epsilon)=\frac{288}{305}\epsilon
                                               +O(\epsilon^2)>0 .
\]
Individual root simplicity and their unit-circle
location persist on this small parameter
neighborhood. The construction applies to every
sufficiently small positive rational \(\epsilon\).
