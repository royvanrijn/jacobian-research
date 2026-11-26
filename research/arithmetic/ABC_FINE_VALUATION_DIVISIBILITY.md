# A count saving for high divisibility of the primitive valuations

The coefficient-uniform generalized Fermat estimate gives a small-budget
saving on a sector whose defining restriction is invisible to the weak
scaled-degree measures. The large common exponent may have scaled width
tending to zero, so this is compatible with formal diffuse measures.

## Exact sector bound

Let \(K\ge2\) be an integer and \(f\ge0\). Count positive primitive
\(a+b=c\le X\) admitting factorizations
\[
 a=s_1x_1^{k_1},\qquad b=s_2x_2^{k_2},\qquad
 c=s_3x_3^{k_3},\qquad k_1,k_2,k_3\ge K,
 \qquad \prod_{i=1}^3\operatorname{rad}(s_i)\le X^f .
\]
All factors here are positive integers. The orders can depend on the
triple and need not be bounded independently of \(X\). Then
\[
 \#\{(a,b,c)\}\ll_{K,\varepsilon}
        X^{f+3K^{-3/2}+\varepsilon}
\]
uniformly in \(f\). The implied constant can be taken uniform in \(f\)
by replacing \(X^f\) by its minimum with \(X^3\) when counting the
coefficient choices.

To see this, choose the joint coefficient radical and assign its
primes to the coefficient letters. Primitivity makes their supports
disjoint, so there are \(3^{\omega(R)}=X^{o(1)}\) assignments uniformly
for \(R\le X^3\). The fixed-support estimate gives \(X^{o(1)}\) triples
of exact coefficients per allocation. There are therefore at most
\(X^{f+o(1)}\) coefficient choices. Fix one such combination. The
variables \(x_i\) are pairwise coprime.

Put \(H=\lceil K^{3/2}\rceil\). If all \(k_i\le H\), there are only
finitely many choices of the exponents, depending on \(K\). For each
choice use Browning–Verzobio, Theorem 1.2. Its two ordinary input
boxes have heights at most \(X^{1/K}\) and the distinguished exponent is
at least \(K\). The auxiliary \(W\) of that theorem is at most
\(X^{K^{-3/2}}\). The three terms in the estimate are bounded,
respectively, by
\[
 X^{2K^{-3/2}},\qquad X^{3K^{-3/2}},\qquad
 X^{K^{-3/2}+K^{-2}} .
\]
The theorem is uniform in the fixed coefficients. If a root box has
trivial size, one can enlarge it to \(X^{1/K}\) since the theorem is
monotone and applies when box sides are at least two.
[Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)

If some \(k_i>H\), enumerate its actual term
\(s_i x_i^{k_i}\) in at most \(O(X^{1/H}\log X)\) ways,
including the term \(s_i\) when \(x_i=1\). For \(x_i>1\) the
order is at most \(\log_2 X\).
If two or three orders exceed \(H\), enumerate two such terms and
the third is determined. The count is \(X^{2/H+o(1)}\).

If just one order exceeds \(H\), after fixing its term the other two
roots satisfy
\[
 \alpha y^{k}+\beta z^{j}=\gamma,\qquad K\le k,j\le H,
 \qquad \alpha\beta\gamma\ne0 .
\]
This plane curve is absolutely irreducible. Over
\(\mathbf C(y)\), the function
\((\gamma-\alpha y^k)/\beta\) has simple zeros, so adjoining its
\(j\)-th root has degree \(j\). Equivalently the binomial is
Eisenstein at one of those simple zeros. The plane equation follows
by Gauss's lemma. The coefficient-uniform integral plane-curve estimate
in the \(X^{1/K}\)-square gives at most
\(X^{1/(K\max(k,j))+\varepsilon}\le X^{K^{-2}+\varepsilon}\)
points. Thus this case costs at most \(X^{1/H+K^{-2}+\varepsilon}\).
Every branch fits \(3K^{-3/2}\).
[Browning–Verzobio, Lemma 2.2](https://arxiv.org/html/2608.24512v1#S2.SS1)

## Consequence for a genuinely fine valuation grid

At fixed \(0<\lambda<1\), take
\(K=\lceil(72/\lambda)^{2/3}\rceil\). On any family where the primes
whose valuations are *not* multiples of their own letter's order
\(k_i\ge K\) have joint radical at most \(X^{\lambda/4}\), place
all these prime powers in the \(s_i\). The displayed estimate gives
\[
 \#\{(a,b,c)\}\ll_{\lambda,\varepsilon}
 X^{(7/24)\lambda+\varepsilon}.
\]
It does not require the radical cap on the other primes. In particular
one may impose this divisibility on at least three quarters of the
normalized radical budget of an abc family. The exponent \(7/24\)
then leaves a fixed gap below \(1/3\).

The condition can select a grid with scaled mesh
\(\lambda K=O(\lambda^{1/3})\to0\). For example, the measures
obtained by assigning mass \((1/3)/J_\lambda\) to each of the
\(J_\lambda\) degree positions in \([2,4]\) that belong to
\(\lambda K\mathbf N\) converge weakly to
\((1/3)\mathrm{Unif}[2,4]\). Thus diffuse limiting measures
do not rule out the grid condition at the profile level. The estimate
excludes arithmetic realization at the endpoint for this sector.

## A bounded number of adaptive valuation grids

Let \(m\ge3\) and \(K\ge2\) be integers. At each positive
primitive triple \(a_1+a_2=a_3\le X\), allow at most \(m\)
grid orders in total, distributed among the three letters.
Every chosen order is an integer at least \(K\), and the
orders and their distribution may depend on the triple.
A prime dividing \(a_i\) is exceptional if its valuation
is divisible by none of the orders chosen for letter \(i\).
In particular, every prime on a letter with no assigned
order is exceptional. Impose a cap \(X^f\) on the product
of all exceptional primes, for \(f\ge0\). Then
\[
 \#\{(a_1,a_2,a_3)\}\ll_{m,K,\varepsilon}
    X^{\,f+C'_m K^{-m^2/(m^2-1)}+\varepsilon}.
 \tag{FV1}
\]
Here \(C'_m\) depends only on \(m\). The implied constant is
uniform in \(f\) and the choices of grid orders.

Assign each nonexceptional prime to one order on its own
letter dividing the valuation. Write
\[
       a_i=s_i\prod_{j\in I_i}x_j^{e_j},
       \qquad e_j\ge K,\qquad \sum_i |I_i|\le m.
\]
Put the entire original power of an exceptional prime into
\(s_i\). Each moving root has only the primes assigned to
its order, at exponents equal to the original valuation
divided by that order. The roots are pairwise coprime,
including roots on the same letter. The coefficients
\(s_i\) have joint radical at most \(X^f\) and can be fixed
at cost \(X^{f+o_X(1)}\), by the same radical enumeration
as for the single-grid count. This is uniform in \(f\)
because \(s_1s_2s_3\le X^3\).

Delete unit roots and put
\[
              T=\left\lceil K^{m^2/(m^2-1)}\right\rceil.
\]
If \(e_j>T\), enumerate the actual prime-power product
\(x_j^{e_j}\), then absorb it into its letter's coefficient.
There are \(O(X^{1/T}\log X)\) candidates for each such
product, since \(x_j\ge2\) and \(e_j\le\log_2 X\).
The bounded number of choices for the allocation of slots
to letters and the set of frozen slots depends only on
\(m\). The remaining orders all belong to \([K,T]\), so
they can be fixed at a cost depending only on \(m,K\).
Every remaining root is at most \(X^{1/K}\).

The residual equation has three nonzero signed coefficients
multiplying the three monomials on their disjoint sets
of roots. A letter with no remaining roots contributes
its fixed nonzero coefficient as a constant. The
bounded-root release (30) allows empty root sets. If no
monomial varies, the selected coefficients give at most
one triple. If only one monomial varies, the equation
fixes the value of that product. For any solution it
is a positive integer at most \(X\), so the finitely many
roots can be read from its divisors at subpower cost.
In all other cases apply (30) with root bound \(X^{1/K}\)
and lower degree \(K\), uniformly in the frozen
coefficients. Its loss in the \(X\)-exponent is at most
\(C_m K^{-1-1/(m^2-1)}\), with the arbitrary subpower
allowance made smaller if necessary. Adding the freezing
cost \(m/T\) proves (FV1).
[Bounded-root release (30)](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#releasing-a-bounded-number-of-high-degree-coordinates)

For fixed \(m\), any exponent
\((m^2-1)/m^2<\alpha<1\) permits
\(K=\lceil\lambda^{-\alpha}\rceil\).
The additional exponent in (FV1) is \(o(\lambda)\) as
\(\lambda\to0\), and the scaled step of the base order
\(K\) is \(\lambda K=o(1)\). Consequently, for any fixed
\(\rho>0\), exceptional prime product at most
\(X^{(1/3-\rho)\lambda}\) gives the count
\[
 O_{m,\alpha,\rho,\lambda,\varepsilon}
      \left(X^{(1/3-\rho/2)\lambda+\varepsilon}\right)
\]
for all sufficiently small fixed \(\lambda\). This permits
any adaptive choice of at most \(m\) orders at each row.
As in the single-grid result, no radical cap is used on
the nonexceptional primes.

One can let \(m=m(\lambda)\to\infty\) sufficiently slowly.
For each \(m\) take
\(\alpha_m=1-1/(2m^2)\) and choose a decreasing sequence
\(\lambda_m\to0\) so that for \(\lambda\le\lambda_m\),
\[
 C'_m\lambda^{1/(2(m^2-1))}\le 1/m,
 \qquad \lambda^{1/(2m^2)}\le 1/m.
\]
Use \(m(\lambda)=m\) on
\((\lambda_{m+1},\lambda_m]\) and
\(K=\lceil\lambda^{-\alpha_m}\rceil\).
Both the \(o(\lambda)\) count loss and the vanishing
base mesh persist. The chosen \(m,K\) are fixed before
the height limit at each budget.

### The same limiting measure can avoid all small grid collections

Set
\(J_\lambda=\#\{r\text{ prime}:2/\lambda\le r\le4/\lambda\}\).
On each letter give mass \(1/(3J_\lambda)\) to each scaled
degree \(\lambda r\) in this set. By the prime number theorem,
these measures tend weakly to \((1/3)\mathrm{Unif}[2,4]\),
with first moments tending to one. This is the same limit
as in the single-grid example. An integer order at least
two can divide at most one of the prime degree indices.
Thus any \(o(J_\lambda)\) grid orders, however chosen,
capture only \(o(1)\) of the normalized mass on these
measures.

These finite degree data can be realized by three pairwise
coprime integers without claiming an additive relation.
At fixed \(\lambda\), choose \(3J_\lambda\) distinct prime
bases in a dyadic interval of size
\(X^{\lambda/(3J_\lambda)}\). Give each letter one prime
base to each of the indicated exponents. The resulting
root masses and heights tend to the stated data in the
height limit, before taking \(\lambda\to0\). The obstruction
is to recovering a genuine exact grid from weak degree
measures alone.


### Prime degree indices also evade sparse signed grids

The same example is an explicit obstruction to bounded-rank
signed digit assignments with a small product of widths.
Write \(J=1/\lambda\) and let \(\mathcal P_J\) be the set
of primes in \([2J,4J]\). If \(U,V\ge1\), \(UV\le J\),
and the integer steps \(P,Q\) are arbitrary, then for any
fixed \(0<\theta<1\),
\[
 \left|\mathcal P_J\cap
   \{Pu+Qv:u,v\in\mathbb Z,\ |u|\le U,\ |v|\le V\}\right|
 \ll_\theta J^{3/4}+
           \frac{J}{\log J}\left(\frac{UV}{J}\right)^{1-\theta}.
 \tag{FV2}
\]
The estimate is uniform in the steps and retains its form
when one digit width stays bounded.

Order the widths so that \(U\ge V\). If \(U<J^{1/3}\),
there are \(O(J^{2/3})\) pairs. If \(P=0\), there are
\(O(V)\) distinct outputs. If \(\gcd(P,Q)>1\), every
output is divisible by this gcd, so at most one output
is prime. Orient \(P>0\) for the other cases. If
\(P>J^{3/4}\), each \(v\) supplies \(O(1+J/P)\) possible
values in the interval, and \(V\le J^{1/2}\) gives the
first term of (FV2).

It remains to take \(P\le J^{3/4}\) and \(U\ge J^{1/3}\),
so \(D=\min(U,J/P)\ge J^{1/4}\). For fixed \(v\) with
\(\gcd(v,P)=1\) the eligible output values are in a
single reduced residue class modulo \(P\), in an interval
of length \(O(PD)\). The interval may start at a
positive real and be enlarged to length at least \(PD\).
Brun–Titchmarsh gives \(O(PD/(\varphi(P)\log J))\)
primes in that class. If \(\gcd(v,P)>1\), there are
no eligible primes since \(P<2J\).
[Montgomery–Vaughan, Theorem 2, equation (1.10)](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf#page=3)

If \(P\le V^4\), inclusion–exclusion gives
\(O(V\varphi(P)/P+2^{\omega(P)})\) eligible \(v\).
For instance, using the elementary inequalities
\(2^{\omega(P)}\ll P^{1/16}\) and
\(\varphi(P)/P\gg P^{-1/16}\), this is \(O(V\varphi(P)/P)\).
Multiplying the counts gives \(O(UV/\log J)\).
If \(P>V^4\), bound the number of \(v\)'s by \(O(V)\).
For \(\nu=3\theta/4\), the bound
\(P/\varphi(P)\ll_\nu P^\nu\) then gives, after division
by \(J/\log J\), at most
\[
  O_\theta\!\left(
      P^\nu\min\!\left\{\frac{UV}{J},\frac{V}{P}\right\}
  \right)
  \le O_\theta\!\left(
      P^\nu\min\!\left\{\frac{UV}{J},P^{-3/4}\right\}
  \right)
  \le O_\theta\!\left((UV/J)^{1-\theta}\right).
\]
The last estimate follows by splitting at
\(P=(UV/J)^{-4/3}\). This proves (FV2).

For a bounded-rank box in signed integer steps of arbitrary size,
apply (VG1) to its output clipped to \([-4J,4J]\).
If its coefficient-width product is \(o(J)\), the
two resulting widths still have product \(o(J)\).
A fixed translate can be included by adding one signed
step of width one. A bounded collection of unrelated
grids can be concatenated when the product of their
widths is \(o(J)\). By (FV2) and the prime number
theorem, any such adaptive box meets only
\(o(|\mathcal P_J|)\) prime positions. In particular the
valuation law of the example has entropy-effective
support size \(\lvert\mathcal P_J\rvert\asymp J/\log J=o(J)\)
and a diffuse weak limit, yet every admissible sparse
signed digit assignment misses asymptotically all
of its radical mass. The example concerns actual
degree data on individual integers. The additive
condition at the endpoint requires a separate estimate.
[Rank compression (VG1)](ABC_SPARSE_VALUATION_GRID_COMPRESSION.md#rank-reduction-inside-a-bounded-valuation-strip)

## Effect of one-target minimal repair

The least multiplier \(Z_M(c)\) is coprime to \(ab\). At a prime of a
primitive summand the valuation in \(Z_M(c)a\) or \(Z_M(c)b\)
is therefore unchanged. In particular the repair cannot make a
primitive summand into a high exact power unless its original
valuations already satisfy the corresponding divisibility.

When all original target valuations are at most \(M\), the repaired
target is \(\operatorname{rad}(c)^M\). If \(k\mid M\) and either
repaired summand is also an exact \(k\)-th power, then at every prime
\(p\mid c\) its valuation is \(M-v_p(c)\). Consequently \(k\mid
v_p(c)\), so the target was already an exact \(k\)-th power before
repair. One can instead put \(Z_M(c)\) in the nonpower coefficients,
but its varying support then enters the coefficient radical cost in
the sector bound. These facts isolate the genuine valuation-grid
input needed for this use of high-power arithmetic.

## Unequal grid spacings

A precise anisotropic statement follows by keeping the input boxes of
Browning–Verzobio separate. Suppose the orders \(p,q,r\ge2\) are fixed
before the height limit, and the three corresponding complete
height caps in base \(X\) are \(h_1,h_2,h_3\). Put
\(u=h_1/p\), \(v=h_2/q\), and \(w=\sqrt{uv/r}\).
The same coefficient enumeration gives exponent at most \(f+E\), where
\[
 E=\max\left(2w,\quad
             w+\frac{2\max(u,v)}{\sqrt{\max\{p,q,36\}}},\quad
             w+\frac{\max(u,v)}{r}\right).
\]
One may minimize over the three choices of distinguished letter.
When an indicated height cap is zero, take arbitrarily small positive
enlargements and pass to the limiting inequality. In particular with
all complete heights at most one, two orders of size at least
\(\lambda^{-\alpha}\) and a third at least
\(\lambda^{-\beta}\) incur \(E=o(\lambda)\) provided
\[
 \alpha>2/3,\qquad \beta>0,\qquad \alpha+\beta/2>1 .
\]
Indeed \(w\le\lambda^{\alpha+\beta/2}\), the second increment is at
most \(2\lambda^{3\alpha/2}\), and the third at most
\(\lambda^{\alpha+\beta}\).
This permits a substantially finer valuation grid on the third
letter while retaining the same coefficient-budget saving.


## Relation to the external count estimates

The input from Browning and Verzobio is the coefficient-uniform
fixed-exponent bound. The radical enumeration and the split at
\(K^{3/2}\) give the estimate for triples with variable grid orders.
Bernert, Browning, Lichtman and Teräväinen give an unrestricted
\(X^{3/5+\varepsilon}\) estimate for \(0<\lambda<1\).
Zeng's complete-modulus-locking preprint reports an unrestricted
exponent \(11680007/20000012\) independent of \(\lambda\). These
constant exponents do not address a saving proportional to
\(\lambda\) as \(\lambda\downarrow0\). Zeng's separate Pell preprint
lists a \(\lambda/3-C\lambda^3\) exponent on its three fixed linear
forms \(n,n+1,2n+1\). Its input uses unit differences.
[Generalized Fermat surface estimate](https://arxiv.org/html/2608.24512v1#S1)
· [Bernert–Browning–Lichtman–Teräväinen](https://arxiv.org/html/2410.12234v2#S1)
· [Complete modulus locking](https://hexagonmath.org/pdf/2611.00621v1)
· [Small Pell discriminants](https://hexagonmath.org/pdf/2611.00611v1)

The fine-grid saving concerns unrestricted signed coefficients with
small joint radical and variable exact power orders. Arbitrary
nonatomic limiting spectra do not give those divisibility
conditions. In particular a model with weights spread evenly across
all integer degree positions in the compact band has asymptotically
negligible mass on every grid of mesh \(k\ge K\), while its weak
limit is the same as that of the grid-supported example above.

Two exact high orders can also be handled with an arbitrary
third-letter spectrum. A rectangular finite-place determinant
uses a low-cost full factor of that letter, or the generalized
Fermat estimate uses a single macroscopic prime on it. In the
symmetric case fixed large multiples of \(\lambda^{-1/2}\) suffice.
[Two valuation grids](ABC_TWO_VALUATION_GRIDS.md)
