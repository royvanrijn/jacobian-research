# Cyclotomic splitting in the primitive quotient of target repair

Write \(N_{\lambda,q}^{\rm opp}(X)\) for the primitive positive
exceptional triples \(a+b=c\le X\) with \(c\) and at least one
summand exact \(q\)-th powers. For every fixed integer \(q\ge2\)
and every \(\eta>0\), for all sufficiently small fixed \(\lambda>0\),
\[
 N_{\lambda,q}^{\rm opp}(X)
 \ll_{q,\lambda,\eta,\varepsilon}
 X^{(1/(\tau(q)+2)+\eta)\lambda+\varepsilon}.
 \tag{CS1}
\]
Here \(\tau(q)\) counts the positive divisors of \(q\). In
particular the coefficient is \(1/4\) for opposite primitive
squares and \(1/6\) for opposite primitive sixth powers.
The proof keeps the full radical constraints on all three inputs.

## Primitive cyclotomic factors and their fixed-value counts

Write \(a=x^q\), \(c=z^q\), with \(0<x<z\) and
\(\gcd(x,z)=1\). The third letter is
\[
             b=z^q-x^q=\prod_{e\mid q}\Phi_e(z,x),
\]
where \(\Phi_e(Z,X)\) is the homogeneous cyclotomic polynomial.
Every factor is positive. Distinct factors can have common
prime support only among primes dividing \(q\).
Indeed, if a prime \(p\nmid q\) divided two factors, then
\(p\nmid xz\) by primitivity and the monic leading and constant
coefficients of the cyclotomic factors. The point \(z/x\bmod p\)
would be a multiple root of \(Z^q-1\bmod p\), which is
separable for \(p\nmid q\). Consequently the sum of the logs of
the factor radicals is
\[
                   \log\operatorname{rad}(b)+O_q(1).
 \tag{CS2}
\]

For \(e\ge3\), fixing \(F=\Phi_e(z,x)\) leaves
\(O_{q,\varepsilon}(X^\varepsilon)\) pairs with \(0<x<z\le X\)
and \(1\le F\le X^q\). To see this use the integer
\(z-x\zeta_e\) in the fixed cyclotomic field, whose norm is
\(F\), where \(\zeta_e\) is a primitive \(e\)-th root of
unity. Its principal ideal can be counted among the integral
ideals of norm \(F\). Above a rational prime there are at most
\(\varphi(e)\) prime ideals, so the number is at most
\(\tau(F)^{\varphi(e)}\). All embeddings of a possible
generator have absolute value at most \(2X\). Their product
has absolute value \(F\ge1\), so every absolute value is at
least \((2X)^{-(\varphi(e)-1)}\). Ratios of generators of the
same ideal are units, with all logarithmic embedding
coordinates \(O_q(\log X)\). These vectors form a discrete
group with finite kernel. For this counting assertion it
suffices that algebraic integers in the fixed field with all
embedding absolute values bounded form a finite set.
This gives a minimum spacing for nonzero logarithmic
differences in any fixed box about zero and a finite kernel.
Thus there are only \(O_q((1+\log X)^{\varphi(e)})\) generator
choices per ideal. Since \(\zeta_e\notin\mathbb Q\), each
generator determines \(x,z\). The bound follows from the
ordinary divisor estimate. This is the fixed-field norm
version of the recovery mechanism used for mixed powers.
[Elementary quadratic norm recovery in Erdős #943](ERDOS_943_PELL_RECOVERED_KERNEL_GCD_WEIGHT.md#1-a-uniform-elementary-norm-lemma)

## Paying for every nonlinear factor

Use dyadic \(c\)-height boxes and refine the logarithmic
radical caps of \(a,c\) and all factors of \(b\), at fixed \(q\).
Write \(L\le (1+o(1))\lambda\) for the resulting total
radical exponent in base \(X\) and \(\beta_e\) for the
factor radical exponents in the same base. By (CS2) the
coordinate caps may be taken so that the sum of all the
\(\beta_e\) and the two other radical exponents is
\(L+o(\lambda)\). For a precise limiting argument, use cap
mesh \(o(\lambda)\) chosen before the inner height limit
and sum over the finitely many bins at each fixed mesh.
An auxiliary \(\varepsilon\) in a count can be chosen
smaller than any prospective failure margin at that
fixed budget.

Let \(\Delta\) be the count exponent in such a bin after
arbitrarily small fixed truncation losses, in base \(X\).
Numbering the possible \(F\)'s by their radicals and then
using the norm bound gives \(\Delta\le\beta_e+\varepsilon\)
for every \(e\ge3\) dividing \(q\). In particular, with
\(S=\sum_{e\ge3,e\mid q}\beta_e\), one has
\(S\ge(\tau(q)-1)\Delta-o(\lambda)\) for odd \(q\), and
\(S\ge(\tau(q)-2)\Delta-o(\lambda)\) for even \(q\).
The fixed-support factor count is uniform since the factor
values have polynomial height.

For odd \(q\), the derived primitive triple
\(x+(z-x)=z\) has logarithmic radical exponent at most
\(L-S+o(\lambda)\) in base \(X\), because
\(\operatorname{rad}(x)=\operatorname{rad}(a)\) and
\(\operatorname{rad}(z)=\operatorname{rad}(c)\).
Its height is comparable to \(X^{1/q}\) in the dyadic
original box, and it determines the original triple.
The unrestricted one-third limiting count applied to the
adjusted budget gives \(3\Delta\le L-S+o(\lambda)\).
To pass between absolute caps in base \(X\) and the
height-relative budget, enlarge the derived budget by
an arbitrarily small fixed multiple of \(\lambda\).
The factor \(q\) in the budget cancels the \(1/q\) in
height, and the derived budget tends to zero as \(q\) is
fixed before \(\lambda\downarrow0\).

For even \(q\), select \(z-x\) or \(z+x\), whichever has
smaller radical. Use the respective ordered primitive
triple \(x+(z-x)=z\) or \(x+z=z+x\), each injectively
within its choice. If \(s\) denotes the selected linear
factor's exponent, then \(2s\le\beta_1+\beta_2+o(\lambda)\).
The derived triple has total logarithmic radical exponent
at most \(L-S-\beta_1-\beta_2+s+o(\lambda)\) in base \(X\).
One of its letters has cap \(s\). Apply the global
one-letter quarter inequality to the derived triple.
It gives
\[
               4\Delta\le L-S-\beta_1-\beta_2+2s+o(\lambda)
                         \le L-S+o(\lambda).
\]
This use of the quarter inequality is in base \(X\);
the derived positive height is comparable to \(X^{1/q}\)
and its budget is again \(O_q(\lambda)\).
Combining the odd or even inequality with the count from
every nonlinear factor gives \((\tau(q)+2)\Delta\le
L+o(\lambda)\). This proves (CS1), with any fixed \(\eta\)
absorbing the vanishing normalized losses.
[Unrestricted one-third limiting count](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)
· [One-letter quarter inequality](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-quarter-inequality-for-a-letter-of-small-radical-mass)

## Additive identities for opposite even powers

Fix an even integer \(q\ge2\) and nonnegative caps \(\ell_a,\ell_b,\ell_c\). Consider positive primitive \(a+b=c\le X\), with \(a=x^q\), \(c=z^q\) and \(\operatorname{rad}(j)\le c^{\ell_j\lambda}\) for \(j=a,b,c\). For every \(\eta>0\) and all sufficiently small fixed \(\lambda>0\), the count is
\[
 \ll_{q,\boldsymbol\ell,\eta,\lambda,\varepsilon}
 X^{\left(
 \frac{2\min(\ell_a,\ell_c)+3\ell_b}{3\tau(q)+2}+\eta
 \right)\lambda+\varepsilon}.
 \tag{CS5}
\]
Under the joint cap \(\operatorname{rad}(ab)\le c^{\ell_{ab}\lambda}\)
alone, the count also satisfies
\[
 \ll_{q,\ell_{ab},\eta,\lambda,\varepsilon}
 X^{(\ell_{ab}/(\tau(q)+1)+\eta)\lambda+\varepsilon}.
 \tag{CS7}
\]
If \(3\mid q\) or \(4\mid q\), there is a further bound
\[
 \ll_{q,\ell_b,\eta,\lambda,\varepsilon}
 X^{(\ell_b/\tau(q)+\eta)\lambda+\varepsilon}.
 \tag{CS6}
\]
For (CS6) only the cap on \(\operatorname{rad}(b)\) is needed. At balanced caps \(\ell_j=1/3\), (CS5) gives \(5/(3(3\tau(q)+2))\) up to the allowance \(\eta\). The coefficient is \(5/24\) for opposite squares. (CS6) gives \(1/9\) for opposite fourth powers and \(1/12\) for opposite sixth powers.

Set \(u=z-x\), \(v=z+x\) and \(F_e=\Phi_e(z,x)\) for \(e\mid q\). In particular \(F_1=u\), \(F_2=v\), and \(\gcd(u,v)\mid2\) because \(\gcd(x,z)=1\). Write \(\beta_e\) for the logarithmic radical masses of the \(F_e\) in base \(X^\lambda\), in a dyadic \(c\)-shell with endpoint \(X\). Use upper caps on a mesh of fixed width \(\delta>0\). At fixed \(\lambda\), for the upper limiting normalized count exponent \(d\) in one bin,
\[
 \sum_{e\mid q}\beta_e\le\ell_b+O_q(\delta),\qquad
 d\le\beta_e\quad(e\ge3,\ e\mid q).
\]
The first inequality holds because only primes dividing the fixed \(q\) can occur in two factor radicals. For the second, enumerate the factor radical and then the complete value at subpower cost on its fixed support. A fixed \(F_e\), \(e\ge3\), gives the norm of \(z-x\zeta_e\) in a fixed cyclotomic field. The divisor count for ideals of that norm and the logarithmic count of bounded generators in each unit orbit give subpower many \(x,z\), as in the fixed-value argument above.

The identities \(u+2x=v\) and \(u+v=2z\) give positive primitive triples after dividing their common gcd, which divides \(2\) in both cases. Each derived triple has at most two original lifts. Choose the identity that uses the opposite root with cap \(\ell_*=\min(\ell_a,\ell_c)\). Its radical product is at most \(2\operatorname{rad}(uv)\operatorname{rad}(x)\) or \(2\operatorname{rad}(uv)\operatorname{rad}(z)\), respectively. The derived height is comparable to \(X^{1/q}\). One of its letters has radical cap \(\min(\beta_1,\beta_2)\) in base \(X^\lambda\). The quarter inequality therefore gives
\[
 4d\le \ell_*+\beta_1+\beta_2+\min(\beta_1,\beta_2)+O_{q,\boldsymbol\ell}(\delta)
             \le \ell_*+\tfrac32(\beta_1+\beta_2)+O_{q,\boldsymbol\ell}(\delta).
\]
There are \(\tau(q)-2\) nonlinear factors, so \(\beta_1+\beta_2\le\ell_b-(\tau(q)-2)d+O_q(\delta)\). Solving for \(d\) proves (CS5). [One-letter quarter bound (LQ7)](ABC_WIDE_COFACTOR_CONGRUENCES.md#small-radical-on-any-one-letter)

For (CS7), always use \(u+2x=v\). Its total radical cap is the cap on \(\operatorname{rad}(a)\) plus \(\beta_1+\beta_2\), up to the fixed factor \(2\). In a bin let \(\beta_a\) be the logarithmic radical mass of \(a\). Primitivity and the joint cap give \(\beta_a+\sum_{e\mid q}\beta_e\le\ell_{ab}\) after bin losses. The one-third bound gives
\[
 3d\le\beta_a+\beta_1+\beta_2+O_{q,\ell_{ab}}(\delta)
                  \le\ell_{ab}-(\tau(q)-2)d+O_{q,\ell_{ab}}(\delta).
\]
This holds with the target radical unrestricted. [Unrestricted one-third limiting count](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

If \(4\mid q\), use \(u^2+v^2=2F_4\) with \(F_4=z^2+x^2\). If \(3\mid q\), use \(u^2+3v^2=4F_3\) with \(F_3=z^2+xz+x^2\). These equations give positive primitive triples after dividing a gcd at most \(4\) in the first case and at most \(12\) in the second. Given the derived triple and the gcd, recover \(u,v\) by positive square roots, hence \(x,z\). The derived height is comparable to \(X^{2/q}\). The radical product is supported on the three factor radicals of \(F_1,F_2,F_e\), with \(e=4\) or \(3\), and the fixed primes dividing \(6\). The one-third bound gives
\[
 3d\le\beta_1+\beta_2+\beta_e+O_{q,\ell_b}(\delta).
\]
Each of the other \(\tau(q)-3\) factors is nonlinear and costs at least \(d\). Summing gives \(\tau(q)d\le\ell_b+O_{q,\ell_b}(\delta)\), which proves (CS6). [Unrestricted one-third limiting count](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

Enlarge each derived cap by \(O(\delta)\) before the height limit to absorb fixed factors and comparisons between shell heights. If the derived height is comparable to \(X^\alpha\), with \(\alpha=1/q\) or \(2/q\), the enlarged total cap \(T\) gives the derived budget \(T\lambda/\alpha\). The factor \(\alpha\) cancels when expressing the count in base \(X\). Choose the mesh width and the one-third or quarter allowance small enough in terms of \(\eta\). There are finitely many resulting cap profiles. Choose the original budget small enough for all their derived budgets and then take the height limit at that fixed budget. Summing the dyadic shells and bins proves the stated bounds.

The angular estimate below applies to genuine near-power targets when a summand is short. It can give a smaller exponent when the parity kernels and target coefficient have very small radicals. The identities above allow comparable summand heights and use the radical of the actual difference \(b\).

## Count after minimal target repair

Keep the sector \(\mathcal S_{\lambda,\xi}\) of the minimal
repair equivalence, and require the two opposite primitive
\(q\)-th powers before repair. In particular the original
three-letter radical is at most \(c^{(1+2\xi)\lambda}\).
For every fixed \(M,\lambda\) in the admissible small-budget regime,
\[
 H_M(c)\le Y\quad\Longrightarrow\quad
 c\le Y^{1/(M(1/3-\xi)\lambda)} .
\]
Use (CS1) with an enlarged budget
\((1+2\xi+\eta)\lambda\). The repaired exponent, up to
arbitrarily small height losses, is at most
\[
 \frac{(1/(\tau(q)+2)+\eta)(1+2\xi+\eta)}
 {M(1/3-\xi)} .
 \tag{CS3}
\]
For \(\xi,\eta\) sufficiently small this is strictly less
than \(1/M\), and the coefficient approaches
\(3/(\tau(q)+2)\) as \(\xi,\eta\downarrow0\).
For opposite squares this is \(3/4\).

For every even \(q\), (CS7) improves this repaired exponent on the entire \(\mathcal S_{\lambda,\xi}\) sector with opposite primitive \(q\)-th powers, since the sector supplies \(\ell_{ab}=2/3+\xi\). The exponent is at most
\[
 \frac{(2/3+\xi)/(\tau(q)+1)+\eta}{M(1/3-\xi)}+\varepsilon.
 \tag{CS8}
\]
As \(\xi,\eta\) decrease, its coefficient relative to \(1/M\) approaches \(2/(\tau(q)+1)\), giving \(2/3\) for squares. This uses the joint cap on \(ab\) and the lower cap on the target radical.

If each original radical also satisfies the upper cap \(c^{(1/3+\xi)\lambda}\), (CS5) gives repaired exponent
\[
 \frac{5(1/3+\xi)/(3\tau(q)+2)+\eta}{M(1/3-\xi)}+\varepsilon .
\]
The corresponding (CS6) exponent is
\[
 \frac{(1/3+\xi)/\tau(q)+\eta}{M(1/3-\xi)}+\varepsilon
 \qquad (3\mid q\ \text{or}\ 4\mid q).
\]
These coefficients relative to \(1/M\) approach \(5/8\) for squares, \(1/3\) for fourth powers and \(1/4\) for sixth powers. All of the counts use the unique mandatory repair of each primitive quotient. [Minimal target-repair equivalence](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

## What remains of a quadratic norm under repair

Let \(a=A x^2\), \(b=B y^2\), with \(A,B\) their positive
squarefree kernels. Primitivity makes their supports
disjoint and coprime to \(c\). \(A x+y\sqrt{-AB}\) belongs to
\(\mathbb Q(\sqrt{-AB})\) and has norm \(A c\).
Let \(Z=Z_M(c)\) and write \(Z=Z_0 u^2\) with \(Z_0\)
squarefree. When the repaired summands are expressed
in squarefree-kernel form, their kernels are \(Z_0 A\)
and \(Z_0 B\). The quadratic field is unchanged and the
new norm element in this same field is
\[
                  Z_0 u\,(A x+y\sqrt{-AB}).
\]
Thus at a split target prime the two valuations each
increase by \(v_p(Z_0u)=\lceil v_p(Z)/2\rceil\).
Their difference in the repaired norm is the original
difference, even if the repaired target has valuation \(M\).
This explains why an argument based on the primitive
imbalance of the two localized norm factors uses the
original orders. Common multiplication also leaves
the squareclass of the ratio of the opposite letters
unchanged. In particular, the real split
\(z^2-x^2=(z-x)(z+x)\) is a primitive square condition.

### Actual norm-family saturation when only kernels are controlled

The cap on \(\operatorname{rad}(ab)\) cannot be replaced in
(MRG) by a bound on the two squarefree kernels and by
local norm data, even with a wholly diffuse target degree
limit and comparable heights.

Fix \(K>4\), take \(M=\lceil K/\lambda\rceil\), and choose
\(J=J(\lambda)\to\infty\) as \(\lambda\downarrow0\), fixed at
each \(\lambda\). Take positive integer exponents \(e_j\le M\)
with \(\lambda e_j\) within \(O(\lambda)\) of the \(J\) equally
spaced midpoints of \([2,4]\). Fix a large \(Y\), and choose
\(J\) distinct primes \(p_j\equiv1\bmod4\) other than 5
in \([B_0,2B_0]\), where
\(B_0=\tfrac12 (Y^{1/M}/5)^{1/J}\).
Order the primes increasingly before assigning the \(e_j\).
Choose one Gaussian prime \(\pi_j\) over each \(p_j\) and put
\[
              x+iy=(2\pm i)\prod_{j=1}^J\pi_j^{e_j}.
\]
One of the two signs makes both coordinates have absolute
value bounded below by a fixed positive multiple of the
complex modulus. Indeed the two possible arguments differ
by \(2\arctan(1/2)\), which has positive distance from
\((\pi/2)\mathbb Z\); choose the argument that stays at
least a third of that distance from the coordinate axes.
Then \(a=x^2\), \(b=y^2\), \(c=5\prod p_j^{e_j}\) form a
positive comparable-height triple. They are primitive
because the Gaussian element includes exactly one of the
two prime factors over each selected split rational prime.
A rational common divisor of \(x,y\) would give both.

Here \(R=\operatorname{rad}(c)=5\prod p_j\),
\(H_M(c)=R^M\le Y\), and both summands have squarefree
kernel one. The radical \(R\) identifies the set of primes,
so choosing one triple per prime set makes distinct
primitive quotients. The prime number theorem in the
fixed progression \(1\bmod4\) gives
\[
      \gg_{\lambda,M,J}Y^{1/M}(\log Y)^{-J}
\]
such quotients, solely from those prime choices.
Their \(\log R/(\lambda\log c)\) tends to
\(1/(\lambda\overline e)=1/3+O(\lambda)\) in the inner
height limit, where \(\overline e=\sum e_j/J\). At fixed
\(\lambda\) the primes have asymptotically equal root log
masses, so the \(\lambda\)-normalized scaled target
distributions tend in the outer limit to mass \(1/3\)
uniform on \([2,4]\). The extra prime 5 has vanishing
mass in the inner limit. [Prime number theorem in progressions](https://academic.csuohio.edu/soprunov-ivan/wp-content/uploads/sites/93/2023/02/primes.pdf) Thus the relaxed exponent
is at least \(1/M\) even for this actual primitive
norm family. The full original summand radical cap is
the condition left out in this construction.


## Direct target bounds for reducible weighted powers

Suppose \(a=A x^q\), \(b=B y^q\) with \(q\ge2\) fixed,
\(A,B,x,y\) positive integers and \(a,b\) coprime.
The coefficients may vary. Suppose the homogeneous
binomial \(A X^q+B Y^q\) has at least \(m\ge2\) irreducible
factors of degree at least two over \(\mathbb Q\).
For any fixed \(M\), count such primitive quotients with
\(H_M(c)\le Y\) and \(\operatorname{rad}(AB)\le Y^{f/M}\).
Then
\[
           \#\{(a,b,c)\}\ll_{q,M,\varepsilon}
                         Y^{(f+1/m)/M+\varepsilon}.
 \tag{CS4}
\]
In particular a fixed margin \(f<1-1/m\) gives a genuine
saving in the minimally repaired count. The full original
summand radical caps may be imposed in addition.

Choose \(A,B\) first. They are coprime and at most \(Y\),
so the joint coefficient radical allows
\(O(Y^{f/M+\varepsilon})\) choices by the uniform supported
factor count. The homogeneous binomial is primitive.
Its primitive integral factors have polynomial-height
coefficients at fixed \(q\). This follows by normalizing
coefficient sup norms and using compactness in the fixed
degree spaces, since two nonzero polynomials have nonzero
product. If \(p\mid c\) then \(A,B,x,y\) are units modulo
\(p\). If also \(p\nmid q\), the binomial is separable
after setting \(Y=1\), so at most one factor vanishes at
\(x/y\bmod p\). Thus the radicals of the \(m\) nonlinear
factor values multiply to at most \(O_q(\operatorname{rad}(c))\).
One of the values has radical at most
\(O_q(\operatorname{rad}(c)^{1/m})\le O_q(Y^{1/(mM)})\).

The fixed-value bound for each nonlinear factor can be made
uniform in the coefficients. If \(F\) is such a factor, of
degree \(d\ge2\) and leading coefficient \(A_0\ne0\), take
a root \(\theta\) of \(F(T,1)\). Then
\(A_0 x-A_0\theta y\) is an algebraic integer of norm
\(A_0^{d-1}F(x,y)\) in degree \(d\). Its conjugates have
polynomial height in \(Y\), since both the coefficients
and \(x,y\) do. At a fixed nonzero value \(F(x,y)\),
the possible principal ideals number at most a fixed
power of the corresponding norm's divisor count.

The unit-spacing argument is uniform in the field here.
A unit of degree at most \(d\) whose logarithmic embedding
coordinates all have absolute value at most one has
all conjugates bounded by \(e\). It therefore belongs to
a fixed finite set by the minimal polynomial coefficient
bounds. Among those with a nonzero logarithmic vector
there is an absolute positive minimum spacing depending
only on \(d\). The kernel is uniformly finite by the same
argument. Comparing two generators of the same ideal
now gives at most \(O_q((1+\log Y)^d)\) choices. A norm
element determines \(x,y\), since \(\theta\) is irrational.
Thus the chosen factor value costs at most
\(Y^{1/(mM)+\varepsilon}\) after counting its radical and
then its exact supported value. Recovering the primitive
quotient determines the mandatory multiplier uniquely.
This proves (CS4).

For example \(a=x^4\), \(b=4y^4\) has the two factors
\((x-y)^2+y^2\) and \((x+y)^2+y^2\), whose fixed values are
quadratic norms. The coefficient radical is fixed, giving
exponent \(1/(2M)+\varepsilon\) without a budget loss.
For the exact even order \(q=2^j s\) with \(s\) odd,
\(a+b=x^q+y^q\) has \(\tau(s)\) nonlinear cyclotomic factors.
Thus the direct repaired exponent is
\(1/(\tau(s)M)+\varepsilon\) when \(s>1\). The pure
power-of-two case has only one factor.

For completeness, if the two summands are exact odd
\(q\)-th powers, \(x^q+y^q\) has \(\tau(q)\) cyclotomic
factors, with exactly one linear factor \(x+y\).
The argument used for the odd opposite-power case
applies to \(x+y=z'\), where \(z'=x+y\) is this factor.
It gives the coefficient \(1/(\tau(q)+2)\) on that
same-side exceptional sector as well.


## Quadratic kernel recovery inside a sparse affine record

The norm bound also gives a local count that does not need
an exact square on a complete letter. In the setting of the
[affine interval count](ABC_SPARSE_AFFINE_CONTENT.md), let
a fixed pair record and an optional target lock \(M_0\)
force the primitive cofactors onto a line \(a_0u+b_0v=c_0\).
As there, assume the sample is a matching, or that
\(\log_X M_0+\min(w_F,w_G)>h\) so the axial lines are
singletons. Suppose the remaining two cofactor squarefree
kernels \(A_0,B_0\) have joint radical cap \(n^p\).
If there are at most \(n^{\rho+o_X(1)}\) choices of
\(M_0\), the count per fixed pair record and interval
is at most \(n^{\rho+p+o_X(1)}\).

The origin line gives one primitive pair. On any nonaxial
line not through the origin, write \(u=A_0 x^2\), \(v=B_0 y^2\)
and fix both kernels. Their values cost at most \(n^{p+o_X(1)}\),
since the kernels are squarefree and coprime.
The equation becomes the fixed norm
\[
            (a_0A_0x)^2-(-a_0b_0 A_0B_0)y^2
                            =a_0A_0 c_0 .
\]
Its right side is nonzero and all data have polynomial height.
If the discriminant is a square, factor the difference of
squares and use the divisor bound. Otherwise use the uniform
quadratic norm count, real or imaginary, proved by the ideal
and unit argument above. Each solution determines the original
cofactors. Thus the parity-kernel cost replaces the smaller
nonatomic cap when the kernels are cheaper. This local
bound still requires a count of the occupied intervals to
cover arbitrary affordable records.

## Angular packing of primitive square kernels and near-power targets

Fix \(0\le h<1\) and an integer \(k>2/(1-h)\). Consider primitive \(a+b=c\le X\) with
\[
 a=A x^2\le c^h,\quad b=B y^2,\quad c=C R^k,\quad
 AB\le X^{f\lambda},\quad \operatorname{rad}(C)\le X^{g\lambda},
\]
where \(A,B\) are positive squarefree. This is an exact parity and defect condition before common repair. The count is \(O(X^{(3f/2+g)\lambda+\varepsilon})\) for general fixed \(k\) and \(O(X^{(f+g)\lambda+\varepsilon})\) for fixed dyadic \(k\). The all-three original radical caps can also be imposed. Either summand can be the short one.

Fix a dyadic \(c\asymp X\) shell. Fix \(A,B\) at cost \(X^{f\lambda+o(1)}\) by allocating primes of their squarefree product, and fix \(C\) at cost \(X^{g\lambda+o(1)}\) by enumerating its radical and the fixed-support count. By primitivity \(D=AB\) is squarefree. Let \(K=\mathbb Q(\sqrt{-D})\), \(\alpha=A x+y\sqrt{-D}\), and work in its full ring of integers. Then \(N\alpha=A c\) and its argument differs from \(\pi/2\) by \(O(X^{-(1-h)/2})\), because \(\Re\alpha/|\alpha|=\sqrt{a/c}\).

Away from primes above \(2AB\), the ideals of \(\alpha,\bar\alpha\) are coprime. A common prime ideal would divide \(2Ax\) and \(2y\sqrt{-D}\), forcing its rational prime to divide the coprime \(x,y\). If \(p\mid AB\), primitivity makes \(p\nmid c\). The factor of \(A\) and the exponents on bad support above \(2AB\) can therefore be fixed at subpower cost, including the varying contribution at \(2\). At good primes dividing \(C\), record which prime ideal carries its contribution (at most two choices per rational prime). All remaining exponents are divisible by \(k\), since the full rational norm exponent at a good prime is carried by exactly one ideal. Hence
\[
             (\alpha)=\mathfrak a\mathfrak b^k,\qquad
             N\mathfrak b\le X^{1/k},
\]
where \(\mathfrak a\) is fixed after subpower many choices. The norm bound uses that the coefficient \(A\) was absorbed on bad support.

For each nonempty ideal class of \(\mathfrak b\), choose an actual representative solution \(\alpha_0,\mathfrak b_0\). Then \(\mathfrak b/\mathfrak b_0=(\gamma)\), with \(\gamma\in\mathfrak b_0^{-1}\), and
\[
       \alpha/\alpha_0=u\gamma^k,\quad
       |\gamma|\le X^{1/(2k)}/\sqrt{N\mathfrak b_0},
\]
where \(u\) is one of at most six imaginary quadratic units. The angle of the quotient \(\alpha/\alpha_0\) is within \(O(X^{-(1-h)/2})\) of zero. For each \(u\), the possible \(\arg\gamma\) lie in \(O(k)\) arcs of that length divided by \(k\). Any two nonparallel vectors in the lattice \(\mathfrak b_0^{-1}\) have sine separation at least their covolume divided by their lengths' product, which here is \(\gg\sqrt{|D_K|}\,X^{-1/k}\) using covolume \(\sqrt{|D_K|}/(2N\mathfrak b_0)\). Hence each arc contains at most one ray for large height. On one ray the ratio of two \(\gamma\)'s is positive rational, because \(K\cap\mathbb R=\mathbb Q\). For fixed \(u\), this makes their \(\alpha\)'s proportional, and coprimality of \(x,y\) determines them. There are \(O(k)\) completions per ideal class. This is the uniform angular-spacing step analogous to the quadratic norm completion for Erdős #943.

For general fixed \(k\), the elementary Minkowski bound gives \(O(D^{1/2}X^{o(1)})\) ideal classes: choose an integral ideal of norm \(O(\sqrt{|D_K|})\) per class and use the uniform divisor bound on ideals of specified norm. This gives the extra \(f/2\). For \(k=2^j\), the equation \([\mathfrak b]^k=[\mathfrak a]^{-1}\) is a coset of \(2^j\)-torsion if nonempty. The bound \(|\mathrm{Cl}(K)[2]|\le2^{\omega(D_K)}\) has an elementary proof. A representative \(\mathfrak d\) of a 2-torsion class satisfies \(\mathfrak d=(v)\bar{\mathfrak d}\) and \(N(v)=1\). Write \(v=w/\bar w\), taking \(w=1+v\) unless \(v=-1\), when take \(w=\sqrt{-D}\). The ideal \(\mathfrak d/(w)\) is fixed by conjugation. Its class is represented by a product of ramified prime ideals, each used at most once, because inert primes and paired split primes contribute rational ideals. Hence \(|\mathrm{Cl}(K)[2^j]|\le2^{j\omega(D_K)}=X^{o(1)}\) for fixed \(j\). Sum over the shells to obtain the stated exponents.

Under the lower repaired target cap \(\operatorname{rad}(c)\ge c^{(1/3-\xi)\lambda}\), the inequality \(H_M(c)\ge\operatorname{rad}(c)^M\) gives dyadic-\(k\) repaired exponent at most \((f+g)/((1/3-\xi)M)+\varepsilon\). For example \(h=3/5,\ k=8,\ f+g\le1/4,\ \xi=1/100\) gives at most \(75/(97M)+\varepsilon\). The condition depends on primitive target power defects. A common repair makes the sum a high power, but introduces shared support in the two norm coordinates, so the coprime ideal decomposition above does not follow from repaired fullness alone.
