# Thin-strip bulk extension to the gap-two endpoint neighborhood

This proof narrows the remaining large-a gap-two problem to the two compact endpoint windows. It does not assert the needed CFU cycle/remainder theorem in those windows.

Assume integer a>=4*10^16, delta=2/sqrt(a)<=10^(-8), and N>a/4+3/2. Let n=N-1 and use the exact actual moment integrals and z normalization in [the complete bulk contour proof](COALESCING_HERMITE_BULK_CONTOUR_EXCLUSION.md).

**Lemma ThinStripExtendedBulkExclusion.** The actual pair cannot have a common zero when |Re z|<=10^(-7) and |Im z|<=19/10.

For theta in[0,1], put kappa=1-theta/N. The exact radial shift is

z_theta=kappa^(-2)*z+(1-kappa^(-2))/delta.

Since 1/N<delta^2, kappa^(-2)-1<3*delta^2. Thus every z_theta lies in the intermediate rectangle

|Re z_theta|<=2*10^(-7), |Im z_theta|<=191/100.

Write z=x+iy in this rectangle. At x=0, the Hermite roots are v_0,+/-=+/-sqrt(1-y^2/4)-iy/2, of modulus1. Their real parts have magnitude greater than37/125=.296, since 1-(191/200)^2>(37/125)^2. On circles of radius r=10^(-6) around these roots, |P_0|>r*(74/125-r). The exact stationary polynomial is

P=v^2+z*v-1+delta*(v^3/4-v)-delta^2*v^2/4.

The perturbation on either circle has modulus below

2*10^(-7)*(1+r)+delta*[(1+r)^3/4+1+r]+delta^2*(1+r)^2/4,

which is smaller than that boundary lower bound. Each circle contains exactly one simple stationary root V_+/-; |V-v_0|<r. The selected branches are analytic by the same nonvanishing derivative bound. Their real-part separation exceeds59/100.

At either actual saddle, compare its curvature A=-F_vv with A_0=1+v_0^2. The actual logarithmic denominator differs from 1/v_0 by at most

t=2*10^(-7)+r+delta*(1+r)^2/4.

The inverse-square difference is bounded by (2t+t^2)/(1-t)^2. Adding the change in h'^2, the h'' term, and the U'' correction bounds |A-A_0| by3*10^(-6). Since Re A_0=2-y^2/2>=2-(191/100)^2/2, this yields Re A>17/100 and |A|<3. On the complex radius1/20 around V, the logarithmic denominator exceeds47/50 and the third derivative is below5/2. Therefore the local horizontal decay is at least3*u^2/50.

For the complete horizontal profile use the symmetric Hermite line Im v=-y/2. Its real phase is

H(t)=-t^2/2+y^2/8+(1/2)*log(t^2+y^2/4).

On either half-line it increases up to the corresponding saddle and then decreases. Its peak S_0=-1/2+y^2/4. The exact centre barrier is

G(u)=-u/2-(1/2)*log(1-u), u=1-y^2/4.

It increases with u, and at u_*=1-(191/200)^2 satisfies G(u_*)>1/600. This follows from (1-u_*)*exp(u_*+1/300)<1. A zero denominator at y=0 has phase minus infinity.

On |v|<=16/5, the actual horizontal heights and the Hermite height differ by less than r. The first log term changes by at most 6*delta+(16/5)*r+r^2/2. The second denominator changes by less than

r+2*10^(-7)+6*delta.

If the Hermite denominator is at least1/10, the total real-phase change is less than1/50000. If it is smaller, the actual phase is below -3/2, while every actual peak exceeds -1. The same comparison holds for the saddle values.

Outside actual horizontal distance1/20, the corresponding Hermite distance exceeds49/1000. The Hermite local decay bound3/50 and its monotone profile yield loss at least (3/50)*(49/1000)^2. Hence the actual nonlocal horizontal loss exceeds

(3/50)*(49/1000)^2-2/50000>1/10000.

Both actual heights lie within r of the SAME Hermite height. The short centre connector therefore has phase at least1/600-2/50000>1/1000 below the larger actual peak. Assign it to that stronger packet, after cancelling the common central vertical exactly.

The original outer real tails have Re F_v<-2 on the right and >2 on the left: the small real part |x|<=2*10^(-7) permits the same bound using Re[h'/(z+h)]<=h'/(h-|x|). At v=+/-3 their phases are below -2. The two far verticals also have phase below -2, since |Im v|<191/200+r, |z+h(v)|<6, log6<2, and the Gaussian real part is at most -(9-(191/200+r)^2)/2+6*delta. Each peak exceeds -1.

Use the entire cancelled contour of the central-bulk proof. The complete local Gaussian relative error is bounded by

[120+240*exp(-m/5000)]/sqrt(m),

and the remaining horizontal, far-tail, and stronger-connector errors are bounded respectively by

3*sqrt(m/2)*exp(-m/10000),
4*sqrt(m/2)*exp(-m),
2*sqrt(m/2)*exp(-m/1000).

All are below1/1000 for m>10^12. Thus the complete selected packets have relative error below1/100 and curvature prefactor ratio below sqrt(300/17)<5.

Finally the exact action derivative from the retained local bulk account, now using Re(V_+-V_-)>59/100, gives

Delta_1<(n/N)*Delta_0-(2/5)*sqrt(a).

At an upper-order zero, the prefactor ratio bound5 and relative errors below1/100 give |Delta_0|<log6. The lower-order packet ratio is then below40*exp(-(2/5)*sqrt(a))<1. The actual lower order cannot vanish. This proves the thin-strip bulk lemma; the displayed elementary bounds are independently replayed by the accompanying exact guards.

Together with the independently proved global Grace rectangle and the outside-endpoint one-packet lemma, any remaining common zero in the present large-a hard regime must satisfy

|Re z|<=10^(-7), 19/10<|Im z|<21/10.

This reduction alone proves no conclusion about that compact turning-point window. Explicit endpoint control for the fixed a>=4*10^16 range is supplied by [the effective gap-two theorem](GAP_TWO_EFFECTIVE_LARGE_LOWER_EXPONENT_THEOREM.md). The broader a>=40000 assertion remains separate.


## A genuine equivalence on the stated hard-order domain

**Corollary HardOrderEndpointWindowEquivalence.**

Let a>=4*10^16 and N>a/4+3/2 be integers, n=N-1, and define the
ACTUAL ordinary factorial polynomials
\[
 A_{m,a}(x)=\sum_{k=0}^{m}\binom mk(am+1)_{2k}x^{m-k}.
\]
Use the COMMON upper-order affine argument
\[
 s_N(z)=N^2[-a^2+2a^{3/2}z],\quad
 U(z)=A_{N,a}(s_N(z)),\quad L(z)=A_{n,a}(s_N(z)).
\]
The lower polynomial retains its own degree n and parameter an+1.
The exponent gap is2 and the moment orders are consecutive.

Then
\[
 \gcd_{\mathbb Q[x]}(A_{N,a},A_{n,a})=1
 \quad\Longleftrightarrow\quad
 \text{there is no common complex zero of }U,L\text{ in }
 \{|\Re z|\le10^{-7},\ 19/10<|\Im z|<21/10\}.                \tag{E1}
\]

To prove the nontrivial direction, apply
[the all-root Grace theorem](GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md).
Its real bound is
8/sqrt(a)+18/a^(3/2)<10^(-7) on the present range, and every upper
root has |Im z|<4. By the
[outside-endpoint one-packet lemma](GAP_TWO_OUTSIDE_ENDPOINT_ONE_PACKET.md),
the upper polynomial has no zero with21/10<=|Im z|<=4.
The thin-strip lemma proved above excludes a COMMON zero when
|Im z|<=19/10. Every remaining common zero must therefore be inside
the set in(E1). Coprimality is equivalent to no common complex root,
and s_N is an invertible complex affine coordinate. This proves both
directions.

Quantifying(E1) over ALL integers a>=4*10^16 and HARD orders
N>a/4+3/2 gives an equivalence between two assertions:
global consecutive coprimality on that domain and noncollision in
the two compact endpoint windows. This does not cover the
complementary soft-order range or prove either assertion by itself.
The effective gap-two theorem separately establishes this fixed4*10^16
range; it does not imply the all-a>=40000 goal. No unproved
Airy-cycle identification or remainder comparison is imported.

The rational coprimality statement is equivalently common-root
exclusion over C or over EVERY characteristic-zero extension field.
The windows themselves are complex analytic regions. Failure of
a rational gcd need not yield a root in a specified nonclosed field K.

## Extending the equivalent goal to all moment orders

**Corollary AllTerminalHugeLowerEndpointEquivalence.**

For the same ordinary exponent-gap-two family, the following universal
assertions are equivalent:

1. For every integer a>=4*10^16 and EVERY N>=2,
   gcd_Q(A_(N,a),A_(N-1,a))=1.
2. For every integer a>=4*10^16 and HARD N>a/4+3/2,
   the actual pair U,L defined above has no common complex zero in
   |Re z|<=10^(-7),19/10<|Im z|<21/10.

The forward implication is restriction followed by(E1). For the reverse,
the hard orders are covered by(E1). For the complementary SOFT orders
N<=a/4+3/2, one has a>=2(2N-3). Apply
[the retained linear lower-exponent criterion](GAMMA_PRODUCT_ROOT_MODULUS.tex),
Corollary cor:linear-lower-exponent, with exponent gap d=2 and
Gamma shape beta=1. It proves coprimality of all normalized ordinary
moment polynomials through order N, including the equality boundary.

In that corollary P_j(t)=T_(j,aj+1,2)(t), whereas our reversed
polynomial is exactly A_(j,a)(x)=x^j P_j(1/x). Both constant terms
of A_j and P_j are nonzero, so reciprocal change preserves common-root
exclusion. This proves soft-order coprimality in the precise convention
of(E1), and completes the reverse implication.

The all-terminal huge-a assertion and its hard-order endpoint-window
formulation are not proved by this reduction; both are established in
[the effective gap-two theorem](GAP_TWO_EFFECTIVE_LARGE_LOWER_EXPONENT_THEOREM.md).
The equivalence has an explicit separate soft-order premise, and no
solution is inferred from the contour estimates alone. The broader
all-a>=40000 gap-two assertion remains open.

## Exact reproduction

    python3 research/factorial/image-rigidity/verify_gap_two_endpoint_reduction_guards.py \
      --repo-root . --output /tmp/gap_two_endpoint_reduction_review.json

The ledger explicitly pins this proof, the outside one-packet proof,
the complete central-contour proof and its checker, the retained
local Hermite/action proof and its checker, and the global Grace
proof/checker, together with the exact retained linear lower-exponent
criterion used only for soft orders. The finite exact guards support the written contour
proof and the elementary endpoint reduction; they do not certify
the remaining endpoint-window noncollision assertion.

Updated 20 November 2026.
