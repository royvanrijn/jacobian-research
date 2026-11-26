# Character moments for correlated complete-value dictionaries

A dictionary can occupy almost the full exponent of its root box
while its root coordinates are strongly dependent. A character
frame on the full box gives high-order moments on the actual
dictionary without replacing its law by a Cartesian product.
This supplies nearly uniform products of actual complete values
at most auxiliary primes. A positive multiplier from \(\{1,2\}\)
can remove the quadratic obstruction when the test primes are
chosen in an appropriate progression.

[Independent products of radical labels (MSU1–MSU2)](ABC_MULTIPLICATIVE_SPECTATOR_UNIT_SUPPLY.md)
· [Signed auxiliary tree comparison](ABC_SIGNED_AUXILIARY_TREE_COMPARISON.md)

## A large sieve for weighted tuples inside a root box

**CFM1.** Fix \(L>0\), set \(n=X^L\), and let
\[
 \Omega=\prod_{i=1}^{J}[1,B_i]_{\mathbb Z},\qquad
 B_i=\lfloor n^{b_i}\rfloor,\quad b_i>0,\qquad
 t=\sum_i b_i,\quad b_*=\max_i b_i .
\]
The list of caps and the nonzero integral orders \(e_i\)
are fixed before the height limit. Let \(a:\Omega\to\mathbb C\)
be any weight, and fix \(\beta>b_*\) and \(a_0>0\).
For a prime \(q\in[Q,2Q]\), \(Q=n^\beta\), and a multiplicative
character \(\chi\bmod q\), put
\[
  \mathcal A_{q,\chi}
       =\sum_{i:\,\chi^{e_i}\ne1} b_i,\qquad
  \widehat a(q,\chi)
       =\sum_{\boldsymbol x\in\Omega}
                         a(\boldsymbol x)\prod_i\chi(x_i)^{e_i}.
\]
The first condition compares each powered character with
the principal character on \(\mathbb F_q^\times\).
All the arguments \(x_i\) are units for sufficiently large \(X\).

For any fixed integer \(h\ge1\) with
\[
                     h a_0\ge 2(2\beta+b_*),
 \tag{CFM1a}
\]
one has, uniformly in \(a\),
\[
 \sum_{\substack{q\in[Q,2Q]\ {\rm prime},\ \chi\bmod q\\
                 \mathcal A_{q,\chi}\ge a_0}}
               |\widehat a(q,\chi)|^{2h}
   \ \ll\ n^{o_X(1)}
     \left(|\Omega|
                  \sum_{\boldsymbol x\in\Omega}|a(\boldsymbol x)|^2\right)^h.
 \tag{CFM1}
\]
The inner subpower allowance is at fixed \(L\) and fixed
root lists. The same result holds with a subset of the test
primes or of the eligible characters. Terms with \(B_i=1\)
from zero-cap coordinates can be omitted.

### A uniform twisted moment on one slot

Index the eligible prime–character pairs by \(\omega\).
For two such indices \(\omega=(q,\chi)\), \(\omega'=(r,\psi)\)
write
\[
   \rho_i(\omega,\omega')
       =B_i^{-1}\sum_{x\le B_i}
                     \chi(x)^{e_i}\overline{\psi(x)^{e_i}}.
\]
With \(k_i=\lceil 2\beta/b_i\rceil\), for any fixed \(\omega\)
one has
\[
     \sum_{\omega':\,\psi^{e_i}\ne1}
                  |\rho_i(\omega,\omega')|^{2k_i}
        \ll n^{o_X(1)}.
 \tag{CFM2}
\]
Indeed let \(d_{k_i,B_i}(m)\) be the number of ordered
products of \(k_i\) integers in \([1,B_i]\) equal to \(m\).
Then \(B_i^{k_i}\rho_i(\omega,\omega')^{k_i}\) is a
Dirichlet polynomial in \(\overline{\psi^{e_i}}\) with
coefficients \(d_{k_i,B_i}(m)\chi(m)^{e_i}\).
They are supported on \(m\le B_i^{k_i}\) and their squared
sum is at most \(B_i^{k_i}n^{o_X(1)}\) by the fixed-order
divisor bound. The anchor \(\omega\) contributes only
a fixed bounded coefficient twist while \(r,\psi\) vary.
Every surviving \(\psi^{e_i}\) is nonprincipal and
therefore primitive modulo its prime \(r\). The power map
on characters has fibres of size at most \(|e_i|\).
The character large sieve gives
\[
 \sum_{\omega':\,\psi^{e_i}\ne1} |\rho_i(\omega,\omega')|^{2k_i}
   \ll_{e_i} n^{o_X(1)}
                   \frac{Q^2+B_i^{k_i}}{B_i^{k_i}} .
\]
This is (CFM2) because \(k_i b_i\ge2\beta\).
The argument includes \(r=q\), since primitivity is needed for
\(\psi^{e_i}\) itself even when its twist against the fixed
anchor is principal. It also includes the one diagonal term
\(\omega'=\omega\) when that term is present.
[Montgomery–Vaughan, Theorem 1, and the Gauss-sum reduction to prime characters](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf)
· [Explicit prime-character reduction in MF5](ABC_MODULAR_JACOBI_INCIDENCE.md#higher-moments-on-an-actual-factorized-letter)

### The character frame on the full box

Group indices \(\omega=(q,\chi)\) by their exact active
slot set \(I=\{i:\chi^{e_i}\ne1\}\). Each group has
\(\sum_{i\in I}b_i\ge a_0\). There are at most \(2^J\)
groups, a constant in the height limit. In the probability
inner product on \(\Omega\), regard the product
\(\xi_\omega(\boldsymbol x)=\prod_i\chi(x_i)^{e_i}\) as a
unit vector, and use the vectors \(\xi_\omega^{\otimes h}\)
on \(\Omega^h\). Their Gram entry is
\(\prod_i \rho_i(\omega,\omega')^h\), up to complex conjugation.
For a fixed \(\omega\), discard the inactive slots from the
absolute row sum. Since
\[
        \sum_{i\in I}\frac1{2k_i}
                  \ge \frac{\sum_{i\in I} b_i}{2(2\beta+b_*)}
                  \ge 1/h ,
\]
choose positive weights \(w_i\) with sum one and
\(2k_iw_i\le h\). Every \(|\rho_i|\) is at most one.
Hölder and (CFM2) now give
\[
 \sum_{\omega'\ {\rm in\ the\ group}}
           \prod_{i\in I}|\rho_i(\omega,\omega')|^h
   \le \prod_{i\in I}
            \left(\sum_{\omega'\ {\rm in\ the\ group}}
                       |\rho_i(\omega,\omega')|^{2k_i}\right)^{w_i}
   \ll n^{o_X(1)} .
\]
The powers in (CFM2) and the weights are fixed at the
height limit. The Hermitian Gram matrix therefore has
operator norm at most \(n^{o_X(1)}\), by the bound on its
absolute row sums.

For the test vector use the \(h\)-fold tensor of the function
\(|\Omega|a\) on \(\Omega\). Its squared norm in the
probability inner product is
\((|\Omega|\sum_{\boldsymbol x}|a(\boldsymbol x)|^2)^h\).
Use the conjugate frame vectors \(\overline{\xi_\omega}^{\otimes h}\),
whose absolute Gram entries are the same. The squared matrix
coefficient of the test vector against them is
\(|\widehat a(q,\chi)|^{2h}\) in the inner-product convention
linear in the first entry. The analysis-operator squared
norm equals the Gram matrix norm. Apply the bound on each
group and sum over the fixed number of groups to obtain
(CFM1). The arbitrary tuple weight is used in this test
vector and need not factor into slot weights.

## Actual complete values at the occurrence scale

Let \(\mathcal D\) be a set of distinct complete values,
each represented by
\[
                  F\prod_i x_i^{e_i},\qquad \boldsymbol x\in\Omega ,
\]
with positive integer orders and a common positive coefficient
\(F\) of polynomial height. Choose one root representation
per complete value. Suppose \(|\mathcal D|\ge n^d\).
Use the uniform probability weight on those root tuples in
(CFM1). Discard the primes that divide \(F\).
The extra character value \(\chi(F)\) has modulus one,
and the discarded primes number \(O_{L,\beta}(1)\) under
a fixed polynomial height cap. For the law \(\mu_q\) of
the actual complete values modulo \(q\) this gives
\[
 \sum_{q}\sum_{\chi:\mathcal A_{q,\chi}\ge a_0}
                |\widehat\mu_q(\chi)|^{2h}
       \ll n^{h(t-d)+o_X(1)}.
 \tag{CFM3}
\]
If the test family has cardinality \(n^{\beta-o_X(1)}\)
and \(\beta>h(t-d)\), then for every fixed
\(0<\gamma<\beta-h(t-d)\), outside a proportion at most
\(n^{-\beta+h(t-d)+\gamma+o_X(1)}\) of its primes,
\[
   \sum_{\chi:\mathcal A_{q,\chi}\ge a_0}
                    |\widehat\mu_q(\chi)|^{2h}\le n^{-\gamma}.
 \tag{CFM4}
\]
This is also a uniform bound on the absolute contribution
of these characters to the relative density of
\(\mu_q^{*2h}\), the law of the product of \(2h\) independent
draws of a complete value. The root coordinates *within*
each draw can remain correlated.

At the balanced complete occurrence scale, each of the
three matching dictionaries has exponent tending to \(1/3\).
On the common anatomic box its sum of root caps also tends
to \(1/3\), with the exceptional coefficient on each role
fixed before using the frame. Passing to a box and deleting
rows costs only the prescribed small normalized allowance,
and no original radical or repair condition is removed.
For any fixed \(\beta>1/3\) one has \(\beta>b_*\) at
sufficiently small balance slack. Taking \(a_0\) a fixed
positive number below \(1/3\) and \(h\) fixed sufficiently
large, (CFM4) gives a polynomial saving on all modes with
that much active raw mass, simultaneously in the three
dictionaries outside a power-small prime set. A finite mesh
of root height caps, chosen before the height limit, absorbs
rounding and ensures that the lists used in (CFM1) are fixed.
[The anatomic decomposition in BBLT, Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)
· [Complete matching and preserved caps (IM1)](ABC_OCCURRENCE_INDUCED_ADDITIVE_MATCHINGS.md#extracting-an-exact-induced-matching-under-all-three-caps)

### Full unit slots and a positive quadratic correction

Choose test primes in the fixed reduced progression
\[
       q\equiv3\pmod8,\qquad
       q\equiv2\pmod p
          \quad\text{for every odd }p\mid\prod_i |e_i|.
\]
The conditions are compatible by CRT and the resulting
dyadic prime set has size \(n^{\beta-o_X(1)}\) by the
prime number theorem in arithmetic progressions.
Its modulus may depend on the original fixed budget
and root lists. For each such prime, \(\gcd(e_i,q-1)\)
divides \(2\). The only nonprincipal character that can
be inactive anywhere is the quadratic \(\chi_2\).
If the odd-order slots have total cap at least \(a_0\),
it also lies in (CFM4), so Fourier inversion gives
\[
   \sup_{u\in\mathbb F_q^\times}
          |(q-1)\mu_q^{*2h}(u)-1|\le n^{-\gamma}
\]
on the good primes. In particular this holds if all the
orders are odd and \(t\ge a_0\).

Without the odd-mass condition, (CFM4) still bounds all
modes except the principal and possibly \(\chi_2\), provided
\(t\ge a_0\). The relative density is
\[
             1+\widehat\mu_q(\chi_2)^{2h}\chi_2(u)
                     +O(n^{-\gamma})
\]
uniformly on the units. Its first two terms give a
mixture of the two quadratic cosets, since
\(\widehat\mu_q(\chi_2)\) is real and its modulus is at most
one. Alternatively multiply the product by an independent
uniform choice from \(\{1,2\}\). Since
\((2/q)=-1\) for the chosen primes, the additional variable
kills the quadratic coefficient exactly. Its other
Fourier multipliers have modulus at most one. Thus the
corrected product has full-unit relative density
\(1+O(n^{-\gamma})\). The correction uses a positive
integer and changes the product height and radical
by at most a factor two.

## Real windows on the same auxiliary products

**CFM5.** Keep the hypotheses and margin of (CFM3)–(CFM4).
Let \(W\) have the uniform law on the actual positive dictionary,
and let \(Z=W_1\cdots W_{2h}\) for \(2h\) independent copies of \(W\).
For a scale \(H>0\), define the triangular window
\[
                    K_H(s)=(1-H|s|)_+ .
\]
The scale can depend on \(X\). Fix \(0<\gamma<\beta-h(t-d)\).
Outside a proportion at most
\(n^{-\beta+h(t-d)+\gamma+o_X(1)}\) of the test primes,
simultaneously for every real centre \(c\) and every
\(u\in\mathbb F_q^\times\) one has
\[
\begin{aligned}
 &(q-1)\mathbb E\!\left[
       {\bf1}_{Z\equiv u\bmod q}K_H(\log Z-c)\right]\\
 &\quad=\mathbb E[K_H(\log Z-c)]
   +\sum_{\substack{\chi\ne 1\\\mathcal A_{q,\chi}<a_0}}
       \overline{\chi(u)}\,
         \mathbb E\!\left[\chi(Z)K_H(\log Z-c)\right]
   +O(n^{-\gamma}).
\end{aligned}
\tag{CFM5}
\]
The bound also holds simultaneously for any prescribed
\(n^{o_X(1)}\) scales, with the same exponent and possibly a
different inner subpower factor in the exceptional proportion.
For example, this includes all dyadic \(H\) in
\([n^{-C},n^C]\) for fixed \(C>0\).

In fact, set
\(\phi_{q,\chi}(\tau)=\mathbb E[\chi(W)e^{2\pi i\tau\log W}]\).
Twist the tuple probability weight in (CFM1) by
\(e^{2\pi i\tau\log W}\). Its squared norm is unchanged,
so the \(h\)-tensor Bessel bound gives
\[
   \sum_q\sum_{\mathcal A_{q,\chi}\ge a_0}
                  |\phi_{q,\chi}(\tau)|^{2h}
       \ll n^{h(t-d)+o_X(1)}
       \qquad(\tau\in\mathbb R).
\]
The constant and subpower bound are uniform in \(\tau\).
The Fourier representation of \(K_H\) is
\[
     K_H(s)=\int_{\mathbb R} e^{2\pi i\tau s}\,d\nu_H(\tau),
     \qquad
     d\nu_H(\tau)=\frac{1}{H}
              \left(\frac{\sin(\pi\tau/H)}{\pi\tau/H}\right)^2 d\tau .
\]
This is the transform of the convolution of two interval indicators
with the normalization \(K_H(0)=1\), so \(\nu_H\) is a positive
measure of total mass one for every scale. Integrating the Bessel
bound proves
\[
   \sum_q\sum_{\mathcal A_{q,\chi}\ge a_0}
             \int |\phi_{q,\chi}(\tau)|^{2h}\,d\nu_H(\tau)
          \ll n^{h(t-d)+o_X(1)}.
\]
Markov's inequality gives the asserted exceptional proportion, and a
union bound handles the prescribed family of scales. Independence of
the \(2h\) draws gives
\[
  \mathbb E[\chi(Z)K_H(\log Z-c)]
     =\int e^{-2\pi i\tau c}
                      \phi_{q,\chi}(\tau)^{2h}\,d\nu_H(\tau).
\]
The magnitude of the phase from the centre is one. Summing its
absolute upper bound over high-active modes at a good prime
therefore costs at most \(n^{-\gamma}\), for every centre.
Character inversion then gives (CFM5).

If all nonprincipal characters are high-active, divide (CFM5) by
the window weight \(w_c=\mathbb E[K_H(\log Z-c)]>0\) to get
relative error \(O(n^{-\gamma}/w_c)\) in the residue law under
window weighting. In particular, \(w_c\ge n^{-\sigma}\) gives
a power saving when \(0\le\sigma<\gamma\).
When low-active modes are present their terms require separate
control, regardless of the window weight. In the progression of
the previous section the only possible low nonprincipal mode is
quadratic. An independent uniform sign multiplying \(Z\) kills that
term at every \(\tau\), since \(\chi_2(-1)=-1\) and
\(\log|\pm Z|=\log Z\). This is an auxiliary signed law with
the same absolute height, applicable when either sign is permitted.
An independent uniform multiplier in \(\{1,2\}\) instead kills
the quadratic term when the window is on the uncorrected
\(\log Z\). For a window on the corrected logarithm the
quadratic frequency multiplier is
\((1-e^{2\pi i\tau\log2})/2\), which vanishes at \(\tau=0\)
and can survive the window integral.

In (CFM5) every dictionary value is a unit for the test primes
because \(\beta>b_*\) and primes dividing \(F\) were discarded.
There is also a unit-conditioned version for any fixed \(\beta>0\).
Extend each powered character in (CFM1) by zero on nonunits,
including the powered principal character. Each frame vector then
has norm at most one, and the slot moments and their coefficient twists
remain valid. The same proof of (CFM1) follows. Apply the resulting
bound to \(\phi_{q,\chi}\) computed with extended characters.
Let \(p_q=\mathbb P(q\nmid W)\). Every polynomial-height value is
divisible by \(O(1)\) primes in \([Q,2Q]\), so the mean of \(1-p_q\)
over the test family is at most \(n^{-\beta+o_X(1)}\).
Primes with \(p_q<1/2\) form at most this proportion.
At each other prime use independent draws from the law of \(W\)
conditioned on \(q\nmid W\). Its twisted transform is
\(\phi_{q,\chi}(\tau)/p_q\), and the integrated bounds change by
at most \(2^{2h}\). This preserves the asserted exponent for the
exceptional proportion. Conditioning supplies units and still
leaves every low-active character in (CFM5); it does not itself
remove a proper-subgroup obstruction.

### Retaining a conditioning record

**CFM6.** Keep \(\beta>b_*\) as in (CFM5), with the same fixed
box and units. Let \(Y\) be a record with any probability law, and
let \(p_y\) be a measurable family of conditional probability
weights on the root tuples of the actual complete values. Put
\(C_y=\sum_{\boldsymbol x}p_y(\boldsymbol x)^2\). Suppose
\[
                   \mathbb E_Y C_Y^h
                                  \le n^{-h(d-\kappa)} .
\]
For each \(y\), form \(Z_y\) from \(2h\) independent draws from
the *same* law \(p_y\), using the common complete-value map.
Let \(E_{q,y}(u,c;H)\) denote the left side of (CFM5) minus its
principal and low-active terms, all three computed with \(Z_y\).
For any fixed \(0<\gamma<\beta-h(t-d+\kappa)\), outside a
proportion at most \(n^{-\beta+h(t-d+\kappa)+\gamma+o_X(1)}\)
of test primes one has
\[
              \mathbb E_Y\sup_{u,c}|E_{q,Y}(u,c;H)|
                      \le n^{-\gamma}.
\tag{CFM6}
\]
This also holds for each scale in a prescribed subpower family
on a common good prime set. In particular the window centre
may depend measurably on \(Y\).

For the proof, apply (CFM1) to each logarithmically twisted \(p_y\),
with squared counting norm \(C_y\). Integrate first over the
same positive Fourier probability measure \(\nu_H\) as in (CFM5)
and then over \(Y\). The total summed moment is at most
\(n^{o_X(1)}|\Omega|^h\mathbb E C_Y^h\). Markov's inequality
at the prime variable and the uniform absolute bound on the
centre phase give (CFM6). The displayed collision condition
holds if \(C_y\le n^{-d+\kappa}\) for every retained history.
If it holds after restricting to a chosen subset of histories,
the conclusion is for the law on that subset. The estimate
keeps \(Y\) in averaged error and includes no conclusion for
discarded histories. The nonprincipal low-active terms remain
explicit. If they are absent or killed without changing the
window, histories of window weight at least \(n^{-\sigma}\)
have averaged relative residue error at most
\(n^{-\gamma+\sigma}\) when averaged with their original
history probabilities. This follows by division on those
histories, for any \(0\le\sigma<\gamma\).

For arbitrary coupled weights \(b\) on \(\Omega^h\) the Gram
argument also gives
\[
 \sum_{\substack{q,\chi\\\mathcal A_{q,\chi}\ge a_0}}
 \left|\sum_{\boldsymbol x^{(1)},\ldots,\boldsymbol x^{(h)}}
             b(\boldsymbol x^{(1)},\ldots,\boldsymbol x^{(h)})
                      \prod_{j=1}^h \xi_{q,\chi}(\boldsymbol x^{(j)})
 \right|^2
 \ll n^{o_X(1)}|\Omega|^h
       \sum_{\boldsymbol x^{(1)},\ldots,\boldsymbol x^{(h)}}
                    |b(\boldsymbol x^{(1)},\ldots,\boldsymbol x^{(h)})|^2 .
\]
Here the test vector on the probability space \(\Omega^h\)
is \(|\Omega|^h b\), without a tensor product assumption.
This is an aggregate squared bound on the coupled weights.
The relative-density estimates above use the higher power
arising from independent repeated conditional draws.

## Comparison and the arithmetic transfer still needed

(MF5) bounds the unweighted second moment on a Cartesian
root box by Hölder over prime–character pairs. In (CFM1)
the Hölder step instead bounds a row of a Gram matrix,
with one prime–character pair held fixed as a twist.
This allows arbitrary weights on the actual correlated
root tuples, at the cost of taking higher moments.
The ray class norm argument uses a fixed modulus and
diagonal moment estimates of ideal products. Here the
average is over prime moduli, and the powers of the
complete values can have large numerical height because
the large sieve is applied to products of the individual
small roots.
[Factorized second moment (MF5)](ABC_MODULAR_JACOBI_INCIDENCE.md#higher-moments-on-an-actual-factorized-letter)
· [Weighted ray class moments](ABC_PROJECTIVE_RAY_CLASS_NORM_MOMENTS.md#the-compact-group-and-the-diagonal-moment)

Independent products of actual complete values now supply
nearly uniform auxiliary units on typical primes when
the corresponding slots can be sampled independently.
A positive correction is available when it is allowed
by the slot's arithmetic construction. A signed tree
application still needs histories in which the original
true incidence gives a quantitative lower bound and the
later conditioning preserves those auxiliary unit laws.
The moment theorem concerns complete-letter marginals
or any other prescribed tuple law with enough mass in its
root box. A law supported on all three letters of a
matching can instead have exponent near \(1/3\) in a
joint box of exponent near one. Estimate (CFM3) does not
give its joint decorrelation at that density.
[Signed-tree arithmetic obligations](ABC_SIGNED_AUXILIARY_TREE_COMPARISON.md#what-the-positive-input-would-need-to-retain)
