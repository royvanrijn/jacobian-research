# All-order factorial windows from uncontaminated residue sections

Fix a prime \(p\) and a finite support of nonzero ordinary
exponent columns. Work at a place over \(p\) where the
natural coefficients \(u_\nu\) in
\[
                 f=\sum_\nu (u_\nu/\nu!)X^\nu,\qquad
                 E(T)=\sum_{n\ge0}\frac{\Lambda(f^n)}{n!}T^n
\]
are integral. Let \(S\) be the unit columns of the form
\(p^b e_i\). Assume \(S\) is nonempty. Let \(H\) be all
the other unit columns. Coefficients that are nonunits
can be discarded when computing the residue of \(E\).

For each row carrying columns of \(S\), let \(b_i,h_i\)
be their smallest and largest digit positions and put
\[
                           d=\sum_i p^{h_i-b_i}.
\]
There may be several columns in \(S\) in one row.

**Theorem.** If every column in \(H\) has digit \(p-1\)
at the own slot \((i,b)\) of every \(p^b e_i\in S\),
then every \(pd\) consecutive positive orders contain
a nonzero ordinary factorial moment. In particular
this gives the actual-support window if the support
size is at least \(pd\).

More generally, the count test below identifies each
residue class of orders whose moments modulo the
maximal ideal equal those from the \(S\)-part alone.
Each such class gives a finite rational-function
certificate for a window valid at every start.

For example at \(p=2\), allow unit singleton columns
\(X,X^2\) and at least two other nonzero monomials.
Every unit nonsingleton must have \(X\)-exponent
congruent to \(3\) modulo \(4\), but its higher
exponent digits and all other coordinates are
unrestricted. Then \(d=2\) and every four-window
is nonzero in the stated sense. Thus the theorem
covers two different singleton positions in the
same row. With one unit singleton in each of \(q\)
distinct rows, the dyadic bound is \(2q\) even
when the nonsingleton heights above the specified
digits are unbounded.

There is a sharp limitation of this residue method.
For pure two-term supports \(\{1,Q=p^h\}\) in one
variable and natural coefficient residues \(1,1\),
the largest residue zero-window length is exactly
\(Q-1\). Thus no bound depending only on the support
size follows from the first residue, even for pure
prime-power columns.

More generally, any fixed finite list of primes and
precisions can miss an arbitrarily long window on
one rational two-term phase, as shown below.

## A count test for a residue section

Use the
[native digit formula](FACTORIAL_MOMENT_SERIES_PRIME_POWER_RATIONALITY.md#degree-bounds-in-characteristic-p).
Every count vector contributing a nonzero residue
has counts \(m_\nu<p\) on \(H\), and the output
digit sums are carry-free. For a fixed count
vector on \(H\), write \(h=\sum_{\nu\in H}m_\nu\)
and let \(\kappa_{i,b}\) be its static output load
at digit \(b\) in coordinate \(i\).
Only vectors with every static load at most \(p-1\)
need be considered.

If \(p^b e_i\in S\) has count \(k_{i,b}\), the lowest
digit of that count appears at its own slot
\((i,b)\). Carries are forbidden there, so
\[
 0\le k_{i,b}\bmod p\le p-1-\kappa_{i,b}.
\]
Other singleton counts, whose higher digits can also
land in that slot, can only reduce the right-hand
capacity available to this digit. Hence a count
vector extending the fixed vector on \(H\) can
contribute only at orders \(n\) with
\[
       n\bmod p\ \in\
       \left\{h,h+1,\ldots,
          h+\sum_{p^b e_i\in S}(p-1-\kappa_{i,b})\right\}\bmod p .
                                                               \tag{1}
\]
Call a residue class available if it is outside
these intervals for every nonzero, admissible
count vector on \(H\). At those orders the
residue series agrees with the series from the
\(S\)-part, term by term in the moment index.
The test is finite in the support and the prime
and does not run through starting orders.

If every column in \(H\) saturates every own
singleton slot, the only admissible nonzero
vector on \(H\) consists of one count of one
column. Indeed \(S\ne\varnothing\) and a shared
saturated slot forbids any further insertion.
Here \(h=1\) and all capacities in (1) are zero.
Every class except \(1\bmod p\) is therefore
available. If \(H\) is empty every class is
available without a test.

## Windows on the available sections

Work over the residue field. In a row \(i\) used
by \(S\) define
\[
 U_{i,j}(T)=1-\sum_{\substack{p^b e_i\in S\\b\le j}}
                                (\bar u_{p^b e_i}T)^{p^{j-b}} .
\]
Independently summing the pure column counts in
each row, or setting the nonsingleton counts
to zero in the native digit formula, gives
\[
       E_S(T)=\frac{P(T)}{D(T)},\qquad
       P=\prod_i\prod_{j=b_i}^{h_i-1}U_{i,j}^{p-1},\qquad
       D=\prod_i U_{i,h_i}.
                                                               \tag{2}
\]
Indeed the factor at each digit is the sum of
carry-free multinomials there, equal to
\(U_{i,j}^{p-1}\). Above position \(h_i\) the
factors are successive Frobenius powers and
their infinite product from \(h_i\) is
\(1/U_{i,h_i}\). The lowest singleton position
in a row gives the unique highest-degree term
of \(U_{i,j}\), so \(\deg D=d\) and
\(\deg P=d-\#\{i:\text{row }i\text{ carries }S\}<d\).
The numerator and denominator have constant
coefficient one.

For \(0\le a<p\), let
\(\mathcal S_a(\sum_j c_jT^j)=\sum_j c_{pj+a}Z^j\).
Write \(D(T)=\sum_j d_jT^j\) and
\(D^{[p]}(Z)=\sum_j d_j^pZ^j\).
Using \(D(T)^p=D^{[p]}(T^p)\) gives
\[
           \mathcal S_a(E_S)=
      \frac{\mathcal S_a(PD^{p-1})}{D^{[p]}(Z)} .
                                                               \tag{3}
\]
The numerator has degree less than \(d\), since
\(\deg(PD^{p-1})<pd\). This is the elementary
section identity for rational series. See
[Bostan, Caruso, Christol and Dumas, §2](https://msp.org/obs/2019/2-1/obs-v2-n1-p08-s.pdf).

A nonzero section is proper and its reduced
denominator has degree \(d_a\) with \(1\le d_a\le d\).
It has no block of \(d_a\) consecutive zero
coefficients. Indeed its recurrence with nonzero
leading and constant denominator coefficients
propagates such a block both backwards to the
first coefficients and forwards to all indices.
Thus any ordinary window containing \(d_a\) indices
in an available class with nonzero section
contains a nonzero ordinary moment.

This supplies a numerical certificate. For a
beginning residue \(r\) let \(\delta(a,r)\) be
the least nonnegative integer congruent to
\(a-r\bmod p\). If \(J\) is the set of available
classes with a nonzero section in (3), a
sufficient length is
\[
             1+\max_{0\le r<p}\min_{a\in J}
                                 \bigl(\delta(a,r)+p(d_a-1)\bigr).
                                                               \tag{4}
\]
In the saturated situation, \(0\in J\) because
the constant term of \(\mathcal S_0(E_S)\) is
one. Taking only class zero in (4) gives at
most \(pd\). The proof works over the residue
field of any finite extension, with the usual
lifting of a nonzero residue to a nonzero
characteristic-zero moment.

## Arbitrarily long residue zero blocks on two monomials

For \(S=\{1,Q\}\) in one row, \(Q=p^h\) with
\(h\ge1\), set both natural coefficient residues
to one. Formula (2) becomes
\[
                   E_S(T)=\frac{(1-T)^{Q-1}}{1-T-T^Q}.
                                                               \tag{5}
\]
In fact \(U_j=(1-T)^{p^j}\) for \(j<h\), and
\(U_h=1-T-T^Q\). Denote its coefficients by \(a_n\).
For \(0\le r<Q\), reduction of (5) through degree
\(Q-1\) gives \(a_r=r+1\), since
\((1-T)^{Q-1}=(1-T^Q)/(1-T)\) in characteristic
\(p\). For \(n\ge Q\) the recurrence is
\(a_n=a_{n-1}+a_{n-Q}\).

Induction by blocks, using Pascal's identity,
therefore gives
\[
        a_{kQ+r}=\binom{r+k+1}{k+1}\pmod p,
             \qquad 0\le k\le Q-2,\quad 0\le r<Q .
                                                               \tag{6}
\]
The induction has zero initial value in each new
block, since the preceding terminal value has
the form \(\binom{Q+j-1}{j}\) with \(1\le j<Q\),
which vanishes by the carry in adding \(j\) to
\(Q-1\). For \(k=Q-2\), (6) is one at \(r=0\)
and zero at every \(1\le r<Q\) by the same
carry test. This is a block of \(Q-1\) zero
residues. The numerator in (5) has no common
factor with the denominator, which does not
vanish at \(T=1\). The degree-\(Q\) recurrence
excludes any block of length \(Q\), proving
sharpness. Its state evolution is invertible
over the finite field, so the exhibited block
also recurs infinitely often.

For example the positive rational phase
\(f=X+X^Q/Q!\) realizes (5). Its ordinary
factorial moments are all positive. The long
blocks here are failures of the first residue
certificate, not counterexamples to the
[ordinary strong factorial conjecture](https://arxiv.org/abs/1304.3956).

## Simultaneous long zeros at prescribed primes and precisions

Let \(L\ge1\) be arbitrary. Fix any finite list of
distinct primes \(p\) and positive precision integers
\(k_p\). There is an ordinary phase \(f=X-X^D/D!\)
with \(D>1\) whose normalized moments satisfy
\[
                  \Lambda(f^n)/n!\equiv0\pmod{p^{k_p}}
                    \qquad(1\le n\le L)
                                                               \tag{7}
\]
simultaneously over the entire list. Its natural
coefficients are \(1,-1\), so both are units at
all the chosen places. The phase satisfies the
ordinary two-window because its even powers
have strictly positive factorial integrals.

For each chosen \(p\), take
\[
 b_p\ge k_p+\left\lfloor\log_p\max(1,L-1)\right\rfloor .
\]
Use the Chinese remainder theorem to choose \(D>1\)
with \(b_p\) trailing zero digits followed by
\(k_p+1\) digits all equal to \(p-1\) in each base.
The required congruence is
\[
         D\equiv p^{b_p}(p^{k_p+1}-1)
                            \pmod {p^{b_p+k_p+1}} .
\]

The exact insertion at the unit linear column,
or direct binomial expansion, gives
\[
              E(T)=\sum_{m\ge0}
              \frac{(Dm)!}{m!\,(D!)^m}
                           \frac{(-T)^m}{(1-T)^{Dm+1}} .
                                                               \tag{8}
\]
See the [exact linear insertion identity](FACTORIAL_SINGLETON_SLOT_WINDOWS_AND_EXACT_INTERPOLATION.md#3-exact-interpolation-with-unit-linears).
At a chosen \(p\), the digit-weight lower bound
for the scalar in (8) is
\((s_p(D)-1)v_p(m!)\). This is at least \(k_p\)
for \(m\ge p\), since \(s_p(D)\ge(k_p+1)(p-1)\).
When \(2\le m<p\), its valuation is the number
of base-\(p\) carry operations in multiplying
\(D\) by \(m\), with their multiplicities.
Each of the \(k_p+1\) prescribed saturated
digits already exceeds the digit capacity after
multiplication by \(m\), so there are at least
\(k_p+1\) carries. Therefore only \(m=0,1\)
remain in (8) modulo \(p^{k_p}\).

At a positive order \(n\), the two surviving
summands yield \(1-\binom{D+n-1}{n-1}\).
For \(1\le n\le L\) all factors in
\[
                \binom{D+n-1}{n-1}
                            =\prod_{j=1}^{n-1}(1+D/j)
\]
are one modulo \(p^{k_p}\) by the choice of
the trailing digits. This gives (7) at every
chosen prime. In particular none of these
fixed-precision tests can certify a nonzero
moment in the specified window. The construction
does not preclude selecting primes or precisions
depending on the exponent profile.
