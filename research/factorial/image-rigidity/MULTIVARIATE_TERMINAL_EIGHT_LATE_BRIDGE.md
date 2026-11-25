# Shared terminal-eight separation from outer mean .72 to .88

Let \(a,b\) be distinct ordinary exponent columns whose unequal
rows are shared. Let \(u\) be the order-seven outer geometric mean
after normalizing the two endpoints to one and removing the binomial
factors. If \(18/25\le u\le22/25\), then the order-seven and order-eight
polynomial root arguments are disjoint under independent positive
rescalings. The dimension and shared exponents are unrestricted.

The [exact checker](verify_multivariate_terminal_eight_late_bridge.py)
separates the possible row \((1,2)\), and separates the possible
row \((2,3)\) on the interval \([.80,.85]\).

## Profiles without a marked row

With lower exponent \(\ell\) and gap \(d\) on each row, use the
[upper-bridge reciprocal cap](MULTIVARIATE_TERMINAL_EIGHT_UPPER_BRIDGE.md#ordinary-row-bounds)
on each outer-mean interval. A unit gap in the present branch has
\(\ell\ge2\), so its order-seven reciprocal is at most \(1/15\).
For \(d\ge2\) the upper-bridge proof gives \(d/(7\ell+h)\le(1-u)/(7u-1)\).
Thus intersect the upper-bridge cap with
\(\max(1/15,(1-u_-)/(7u_--1))\), where \(u_-\) is the
interval's lower endpoint. On \([4/5,17/20]\) also put the
row \((2,3)\) in a marked case, so that in the unmarked case
the unit-gap bound is \(1/22\) instead of \(1/15\).
Use this reciprocal cap in the same
[quadratic Cayley coefficient model](MULTIVARIATE_TERMINAL_EIGHT_UPPER_BRIDGE.md#quadratic-cayley-families).

Every actual normalized interior multiplier is at most one by
log convexity of the factorial coefficients. This bound also
holds for the symmetric reference arrays and the straight
coefficient segments used in the deformation. Since
\((13/12)^8<2\), dominance by the constant and leading terms
bounds all the physical roots within \(1/12<|t|<12\).

Put \(p=1-u=\lambda^2\). The positive ray parameters \(b\) are
piecewise affine at the checker's rational knots, with physical
argument \(\pi-2\arctan(\lambda b)\). The checker verifies
\(p b^2<1\) on every knot interval, so these are negative-cosine
rays throughout. Because coefficient reversal gives
\(t\mapsto1/\overline t\) on the upper ray and preserves the
model, it suffices to check radii \(1/12\le r\le1\).
For \(X=(1+t)/(\lambda(1-t))=x+iY\), the horizontal coordinate
\(\xi=\lambda x\) satisfies
\[
 0\le\xi=\frac{1-r^2}{1-2r\cos\theta+r^2}
          \le\frac{1-r^2}{1+r^2}\le\frac{143}{145}.
\]
The exact ray equation is
\[
 (Y-b)(1+pbY)=-b\xi^2,\qquad 0<Y\le b .
\]
The checker bounds \(Y\) with three interval updates from
\(0\le Y\le b\) and uses \(x=\xi/\lambda\).
This avoids needing a scaled-root disk bound.

The checker bisects rectangles in \(u,\xi\), testing the strict
sign of \(\Im(\widetilde P_j(X)\overline X^{\,h})\) for at least
one \(0\le h\le j\). It first tries a combined coefficient interval
and then uses the individually projected contribution of each
shared scalar and product. A cell is accepted only with a strict
sign, and reaching the subdivision limit is an error.

At fixed \(u\), deform the actual array to the symmetric
reference at ratio \(49/48\), along the straight coefficient
segment in the convex coefficient enclosure. Then increase \(u\)
to \(22/25\) along the symmetric arrays. The ray parameters
match the [upper bridge anchor](MULTIVARIATE_TERMINAL_EIGHT_UPPER_BRIDGE.md#ray-exclusions-and-root-counts)
at that endpoint. The verified ray exclusions preserve its
seven upper sector counts.

## Profiles with a marked unit-gap row

At order seven a row \((1,2)\) contributes \(4/7\) to \(u^2\),
so two such rows would imply \(u\le4/7<18/25\). Its residual
outer mean \(v\) satisfies
\(v^2\ge(18/25)^2/(4/7)>(19/20)^2\).

A row \((1,2)\) is impossible for \(u\ge4/5\). In the interval
\([4/5,17/20]\), a row \((2,3)\) contributes \(5/7\) to \(u^2\)
and is also unique, since two such rows would give \(u\le5/7\).
Its residual outer mean satisfies
\(v^2\ge(4/5)^2/(5/7)>(47/50)^2\).

Orient the marked row as \((a,a+1)\), with \(a=1\) or \(2\).
It contributes the exact normalized binomial array
\[
 B_{j,k}=\binom jk\frac{(ja+k)!}{((ja)!)^{1-k/j}((j(a+1))!)^{k/j}} .
\]
On each interval with lower endpoint \(v_-\), every residual
order-seven reciprocal is at most \(c=(v_-^{-2}-1)/6\),
by the outer mean product. The order-eight cap is \(7c/(8-c)\).

Use the
[paired marked coefficient method](MULTIVARIATE_TERMINAL_EIGHT_MARKED_SINGLETONS.md)
with this retained array and these caps. On each interval of \(v\)
use rational disk centers proposed by affine interpolation at its
midpoint. The fixed radii are \(R=3/500\) for \((1,2)\) and
\(R=1/100\) for \((2,3)\). For each disk the checker
proves the strict estimate
\[
 |P'(z_0)|R > |P(z_0)|+
       \sum_{h=2}^j \frac{|P^{(h)}(z_0)|}{h!}R^h
\]
for every coefficient array in the model. It combines each shared
scalar's contribution before applying its interval in each Taylor
coefficient. Rouché's theorem gives one root in each disk. The
upper disks, their conjugates, and the additional negative real
septic disk are disjoint and account for all the roots.

The six marked ray cosines are affine on their own rational knots.
For a disk center \(x+iY\) the checker bounds its signed distance
\(\cos\theta\,Y-\sin\theta\,x\) to the line of each ray over the
entire \(v\) interval. A bound with magnitude greater than \(R\)
places the whole disk on the indicated side. The strict bounds
give the required sector counts directly at every \(v\), using
26 root-disk intervals for \((1,2)\) and 50 for \((2,3)\).

## Sector conclusion

The unmarked anchor and the marked disks give successive upper sectors of orders
\(8,7,8,7,8,7,8\), and the septic has the additional negative
real root. Join the first and last sectors to their lower
reflections during continuation. The only septic root in the last
joined sector remains real by conjugation. The two octic roots
there are nonreal since the actual eighth factorial moment is a
strictly positive Gamma integral on the real line. Hence the
orders have disjoint root arguments, and rational Bézout gives
coprimality over characteristic zero.

The replay verifies 327,414 unmarked Cayley projection cells and
76 marked root-disk intervals.

Replay from the repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eight_late_bridge.py
~~~
