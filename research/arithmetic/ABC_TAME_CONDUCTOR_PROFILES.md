# Tame conductor loss on dense valuation profiles

Let \(\pi:Y\to\mathbb P^1\) be a fixed geometrically connected
regular cover over a number field \(K\), with all deck
transformations defined over \(K\), genus \(g(Y)\ge2\),
and branch points among \(0,1,\infty\). Denote the
respective inertia orders by \(m_1,m_2,m_3\), using \(1\)
if a point is unramified. Riemann–Hurwitz gives
\(\sum_i1/m_i<1\).

For a primitive \(a+b=c\) at height \(X\), put \(t=a/c\)
and \(T_1=a,T_2=b,T_3=c\). Let \(E_t/K\) be the splitting
field of the geometric fibre at \(t\). For the fixed
cover exclude a fixed finite set of rational primes,
including the ramification of \(K\) and the bad and
wild primes of specialization. Define the normalized
weights
\[
  w_{i,e}=\frac1{\log X}
                    \sum_{\substack{q\ \mathrm{good}\\v_q(T_i)=e}}
                        \log q,\qquad
  L_i=\sum_{e\ge1} w_{i,e},\qquad L=\sum_i L_i .
\]
There is an exact equality for the good-prime contribution
to the logarithmic root discriminant in units of \(\log X\):
\[
     \mathscr D_{\rm good}
           = \sum_i\sum_{e\ge1} w_{i,e}
                             \left(1-\frac{\gcd(m_i,e)}{m_i}\right).
\]
The corresponding ramified radical weight is
\[
     \mathscr R_{\rm good}
                    =\sum_i\sum_{m_i\nmid e} w_{i,e}.
\]
Indeed at a good prime of support colour \(i\), the
intersection multiplicity with the corresponding branch
section is exactly \(e\), by primitivity and the
unramified base. Beckmann's specialization inertia theorem
makes the tame inertia a cyclic group of order
\(m_i/\gcd(m_i,e)\), which gives the root discriminant
exponent \(1-\gcd(m_i,e)/m_i\). This statement uses
the actual arithmetic inertia in the Galois fibre field,
including cases where the fibre field has degree smaller
than the regular cover.
[Conrad, Inertia groups and fibers, Theorem 0.1 reporting Beckmann Proposition 4.2](https://math.stanford.edu/~conrad/papers/inertiapaper.pdf)

## Prime degrees and dense consecutive degrees

For a balanced profile \(L_i=L/3+o(L)\) supported on
prime valuation degrees, suppose that
\(\max_e w_{i,e}=o(L)\) on each letter. Then, uniformly
over hyperbolic inertia triples,
\[
               \mathscr R_{\rm good}\ge (1-o(1))L,\qquad
               \mathscr D_{\rm good}\ge (85/126-o(1))L .
\]
Every \(m_i>1\) divides at most one prime degree, so the
first assertion is immediate. For the second, use
\[
  \sum_{e\ {\rm prime}} w_{i,e}\frac{\gcd(m_i,e)}{m_i}
        \le \frac{L_i}{m_i}+\max_e w_{i,e},
\]
which follows from
\(\sum_{q\mid m,\ q\ {\rm prime}}(q-1)\le m-1\) by
expanding the product of \(1+(q-1)\).

The largest possible reciprocal sum of a hyperbolic
integer triple is \(41/42\), attained at \((2,3,7)\):
after sorting, separate \((2,3,m\ge7)\), \((2,4,m\ge5)\),
\((2,m_2\ge5,m_3)\), and \(m_1\ge3\) with
\((3,3,3)\) excluded. The balanced saved charge is at most
\((41/42)L/3+o(L)\), giving the constant \(85/126\).

Even the dense profile which is uniform on \(J\)
consecutive integer valuation degrees on each letter up
to total variation \(o(L)\), with \(J\to\infty\), forces
\[
                \mathscr D_{\rm good}\ge (4/9-o(1))L .
\]
The estimate is uniform in the inertia orders, including
orders varying with \(J\). For each \(m\ge2\) and any
interval \(I\) of \(J\) integers, the identity
\(\gcd(m,e)=\sum_{d\mid(m,e)}\varphi(d)\) gives
\[
    \frac1{Jm}\sum_{e\in I}\gcd(m,e)
       \le A(m)+\frac1J,\qquad
    A(m)=\frac1m\sum_{d\mid m}\frac{\varphi(d)}d
        =\prod_{q^v\Vert m}\frac{1+v(1-1/q)}{q^v}.
\]
The error bound uses \(\sum_{d\mid m}\varphi(d)=m\).
In fact the last expectation is \(A(m)+O(1/J)\),
uniformly in \(m\). For weights that are within \(o(L)\)
in summed absolute mass of uniform interval weights, the
same estimate for each weighted sum costs \(o(L)\).

We obtain the exact best good-prime charge for arbitrary
letter masses in this consecutive model. Write the three
masses in decreasing order \(x\ge y\ge z\). Then the
least possible conductor charge over regular hyperbolic
three-point covers, on the formal profile, equals
\[
 L-\max\!\left\{\frac34x+\frac59y+\frac5{16}z,\
                       \frac34x+\frac12y+\frac5{12}z\right\}
             +o(L).                                         \tag{TC1}
\]
The two expressions use \((2,3,8)\) and \((2,4,6)\)
respectively, assigning the largest saved coefficients
to the largest masses. Both signatures are hyperbolic and
realized by the three-point cover theorem, so the
nominal minimum is attained in the formal calculation.

To prove the upper bounds on saved coefficients, note
\(A(2)=3/4,\ A(3)=5/9,\ A(4)=1/2,\ A(5)=9/25,\ A(6)=5/12\).
Also \(A(m)\le5/16\) for \(m\ge7\), with equality at \(8\).
Indeed for each odd prime \(p\) and exponent \(a\ge1\)
one has \(1+a(1-1/p)\le p^{a/2}\), by checking \(a=1\)
and comparing the successive ratios. For \(p=2\), the
ratio of the left side to \(2^{a/2}\) is at most
\(3/(2\sqrt2)\). Thus \(A(m)\le 3/(2\sqrt{2m})<5/16\)
for \(m\ge12\), and the values at \(7,8,9,10,11\) give
the assertion.

Sort the inertia orders by their size for the case split,
and their \(A\)-values by size for pairing with the masses.
If the two smallest orders are \(2,3\), the third is at
least seven, and the first expression in (TC1) bounds the
saved charge. If they are \(2,4\), the third is at least
five and its maximal \(A\)-value is \(A(6)\). If the
smallest is \(2\) and the other two are at least five,
their \(A\)-values are both at most \(5/12\). These two
cases are bounded by the second expression. If all three
orders are at least three, the sorted \(A\)-values are
bounded coordinatewise by \((5/9,5/9,1/2)\), since
\((3,3,3)\) is excluded. Pairing this vector with
\(x,y,z\) is bounded by the second expression too:
its difference is \((7x-2y-3z)/36\ge0\).
This proves (TC1). The second expression wins exactly
when \(z\ge(8/15)y\). For \(x=y=z=L/3\) it equals
\(5L/9\), giving the sharp formal lower charge \(4L/9\).
[Darmon–Granville, Proposition 3.1](https://www.math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf)

For example, consecutive degrees in \([J,2J)\) with
equal radical mass and \(L\) chosen near \(2/J\) are
compatible with the numerical height moments. The
estimates apply to profile data and do not supply
families of primitive triples with such distributions.

The rotation among covers may be taken after the
profile parameters are fixed. For a fixed budget and
finite list of covers chosen before the inner height
limit, their finite bad-prime sets have subpower
weight. Allowing arbitrary covers or base fields varying
with the actual triple requires separate control of
their bad primes.

## Charging the same tame conductor to a moving twist

Keep the original regular cover over \(K\) with
constant deck group \(G\). At a fixed rational parameter
\(t\), choose an element of its geometric fibre and
write \(\rho_t:G_K\to G\) for the specialization
homomorphism. For a twist coming from a \(G\)-torsor
cocycle \(\sigma:G_K\to G\), the twisted fibre has
action \(x\mapsto\rho_t(h)x\sigma(h)^{-1}\), up to the
choice of left and right actions. At a good tame prime,
if the twisted fibre has even one point over the
maximal unramified local extension, inertia fixes some
\(x\). Hence
\[
                \rho_t(i)=x\sigma(i)x^{-1}
                  \qquad\text{for every inertia element }i.
\]
In particular the twist's torsor splitting field has the
same inertia order as the original fibre field at that
prime, and the same tame root-discriminant exponent.
A twist unramified at that prime cannot cancel any
nontrivial specialization inertia. If one charges the
twisting field in addition to the lifted points' field,
the local tame conductor does not disappear. The calculation
concerns twists of the specified cover. It leaves open
covers with new branch data tied to the arithmetic family.
