# Consecutive reciprocal powers and moving targets of small radical

The additive relation count can retain both a moving target and
an even number of distinct active rows. The useful coordinates
are two consecutive powers of the same signed reciprocal. At a
prime shared by exactly two active rows, the two equations give
an amplified congruence between actual complete letters. At a
prime whose last occurrence is on the current row, the earlier
rows determine the current signed complete value at the same
precision. The arguments use both congruences before allowing
any numerical weight to escape to high orders.

## Statements and order of limits

For positive primitive \(x=(a,b,c)\) with \(a+b=c\) and
\(a\le b\), put
\[
       A_s(x)=a^{-s}+b^{-s}+(-c)^{-s}\qquad(s\ge1).
\]
Fix a polynomial-height allowance \(K>0\) for nonzero rational
weights, and put \(n=X^L\), with \(L>0\) fixed before taking
the height limit. All rows lie in a shell \(X/2<c\le X\).
For fixed \(k\ge2\) weights \(w_1,\ldots,w_k\) of reduced
height at most \(X^K\) and distinct rows \(x_1,\ldots,x_k\)
write
\[
 t_s=\sum_{i=1}^k w_i A_s(x_i).
\]
Let \(P_0\) be the union of the prime supports of all weight
numerators and denominators. In the fixed-letter assertion add
the prime support of the fixed complete letter \(T\). The
target support is the union of the prime divisors of the
reduced denominators of \(t_q,t_{q+1}\) outside \(P_0\).
The targets need not be fixed. The supports and weights can
depend on height. All implied constants are uniform in the
specific choices subject to the stated bounds.

Under the total cap \(\operatorname{rad}(abc)\le n\), there
is an absolute integer \(q\ge1\) such that for every fixed
\(k,K\) there are \(\delta,\mu_0,L_0>0\) for which the number
of tuples whose target support has product at most \(n^{\mu_0}\)
is
\[
           O_{k,K,L,\varepsilon}
                \bigl(X^{(k/6-\delta)L+\varepsilon}\bigr)
                   \qquad (0<L<L_0).                          \tag{VP1}
\]
The estimates apply to arbitrary subfamilies of the shell.

There is a fixed-letter version. Prescribe any complete \(T\)
with radical at most \(n^C\), and let \(U,V\) be the two
moving roles. Assume raw caps \(\ell_U,\ell_V\) in base \(n\)
with \(B=\ell_U+\ell_V>0\), where the caps and \(B,C\) are
fixed independently of \(L\). For an integer \(q\) depending
only on \(B,C\), the analogous count is
\[
         O_{B,C,k,K,L,\varepsilon}
            \bigl(X^{(kB/8-\delta)L+\varepsilon}\bigr)
                   \qquad (0<L<L_0),                          \tag{VP2}
\]
with fixed \(\delta,\mu_0,L_0>0\) depending on the displayed
fixed parameters. In this version one may omit the reciprocal
of the fixed signed value in each \(A_s\), since shifting the
target by its known weighted contribution does not change the
denominator orders outside \(P_0\).
For a fixed total cap \(n^B\) on a global family use (VP1)
with budget \(BL\), and rescale the target support allowance.

The small target support can include any prescribed pool
supported on an integer of polynomial height, by adding that
pool to \(P_0\). In particular the sum targets may have
arbitrary numerators and denominators supported in such a
pool. A relation with repeated rows can first merge equal
rows and add their weights. For each resulting equality
pattern with \(k\) distinct active rows, these estimates apply
to those \(k\) and the nonzero merged weights. Each additional
cancelled group can be counted by the one-row bound.

## Compact positive valuation mass

The required regularity is a common finite band of scaled
orders. In (VP1) there are absolute \(r_0,R_0,\gamma,\sigma>0\)
such that, outside a class of size
\(O_{L,\varepsilon}(X^{(1/3-\sigma)L+\varepsilon})\),
each role \(W\) has
\[
       \sum_{\substack{p\mid W\\ r_0\le L v_p(W)\le R_0}}
                          \log_n p\ \ge\ \gamma.               \tag{VP3}
\]
In (VP2), the same statement on the two moving roles holds
with exponent \(B/4-\sigma\), with constants depending on
\(B,C\). These assertions are consequences of the height-first
profile lemmas. For the global assertion use the uniform
positive-band result at the balanced one-third limit. For the
fibre use (FT1), which gives positive nonatomic mass on each
moving role in any fixed positive-exponent limiting fibre.
The raw caps in (VP2) are bounded, so the mass measures are
tight by their bounded numerical first moments. To deduce
the uniform finite statements, suppose they fail at every
finite band and tolerance. Choose bands expanding to
\((0,\infty)\) and band mass thresholds tending to zero.
Fix the offending role on a subsequence. For each fixed budget
choose arbitrarily large heights violating the count and make
the exact-layer box refinement before choosing a diagonal
across budgets. The offending role's weak profile is then
supported at zero. At the global coefficient it contradicts
the positive-band result. On the fibre the limiting exponent
can be taken arbitrarily close to \(B/4>0\), so the zero
nonatomic mass contradicts (FT1). This proves (VP3) at all
sufficiently small fixed budgets. One can choose \(\sigma<B/4\)
in the fibre version.
[Positive scales in every endpoint letter](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#positive-scales-in-every-endpoint-letter)
· [Nonatomic mass on a complete fibre](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#nonatomic-mass-outside-prescribed-prime-pools)

## Two exact amplifications

At a prime outside the weight and target supports, exactly
two appearances of \(p\) in the row radicals mean that,
writing \(U\) and \(V\) for the two signed complete values
it divides, both
\[
           w_i/U^q+w_j/V^q,\qquad
           w_i/U^{q+1}+w_j/V^{q+1}
\]
are \(p\)-integral. Subtract the first expression from the
second multiplied by \(V\). Since \(w_i\) is a \(p\)-unit,
this gives the exact divisibility
\[
                  p^{(q+1)v_p(U)}\mid V-U .                  \tag{VP4}
\]
In particular the two prime orders agree. The negative sign
on the target is part of \(U,V\). This allows one to recover
signed value equality across roles. Two distinct signed
values determine the whole sorted primitive triple, because
the third is minus their sum. If one complete letter is
fixed, one of the two moving signed values determines the
sorted triple by the same rule. In a shell at large height,
the signed values in a row are distinct. Equality of the
positive summands would force \((a,b,c)=(1,1,2)\).

For the other amplification, fix an order of exposing rows,
with \(x_i\) the current row, and put
\[
            S_s=-\sum_{\ell\ {\rm earlier}} w_\ell A_s(x_\ell).
\]
For a prime \(p\) outside the weight and target supports that
occurs on the current row and is absent from every later
row, the unknown targets are \(p\)-integral. If \(p^e\)
is its full power on current signed role \(U\), then
\(S_s\equiv w_i/U^s\bmod\mathbb Z_p\), for \(s=q,q+1\).
Consequently
\[
       v_p(S_{q+1})=-(q+1)e,\qquad
       U\equiv S_q/S_{q+1}\pmod{p^{(q+1)e}}.                 \tag{VP5}
\]
Indeed \(U S_{q+1}-S_q\) is \(p\)-integral, and the displayed
ratio is \(p\)-integral as well. Thus the precision and
residue come from the actual earlier rows, even though the
targets can move. The denominator of the known \(S_{q+1}\)
has polynomial height at fixed \(k,q,K\). Choosing any
subset of its primes and assigning their roles costs
\(X^{o_X(1)}\) uniformly. If \(S_{q+1}=0\), there are no such
marked primes. The CRT applied to any one role fixes its
complete signed value as soon as its product modulus exceeds
\(2X\).
[Single projective coordinate and the squared last-prime modulus](ABC_THREE_CAP_RATIO_SUMSETS.md#full-congruences-from-the-last-occurrence-of-a-denominator-prime)


## Paying for a reconstructed intermediate row

Use \(B=1\) in (VP1) and \(B=\ell_U+\ell_V\) in (VP2).
Write \(\kappa=3\) for the global count and \(\kappa=4\)
for the fixed-letter count. All radical prices in this
argument are for the three complete letters in (VP1) and
only the two free letters in (VP2). If the earlier rows
and any polynomial-height prime pool are fixed, the number
of next rows with outside radical mass at most \(v\) is
\[
                   O(X^{(v/\kappa+\eta)L+\varepsilon})
                                                               \tag{VP6}
\]
for any fixed \(\eta,\varepsilon>0\) at sufficiently small
fixed \(L\). This is (EB4) globally. On the fibre, if the
two outside masses are \(u,v'\), (FL1) gives
\(\min(u,v')/2\le(u+v')/4\). Bin the separately capped
masses to deduce (VP6). Both bounds allow the earlier row
primes and \(P_0\) in the prescribed pool.
[Effective price on a fixed prime pool](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Anisotropic free-letter price](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

For each row let \(b_i\le B\) be its tracked raw radical
mass. At a prime let \(m_p\) be its multiplicity across
these rows. Let \(D=kB-\sum_i b_i\), and let \(S\) be the
incidence mass in the fixed pool plus the sum of
\((m_p-2)\log_n p\) over the outside primes with \(m_p\ge3\).
Let \(U_1\) be the mass of outside primes with \(m_p=1\).
A singleton survives in the denominator of either sum,
by the unique negative term of its \(p\)-adic expansion.
Thus \(U_1\le\mu_0\). In any order let \(v_i\) be the new
mass on row \(i\), excluding \(P_0\) and all earlier row
primes. Then
\[
       2\sum_i v_i=kB-D-S+U_1,\qquad
       \sum_i v_i/\kappa
         \le \frac{kB-D-S+\mu_0}{2\kappa}.                    \tag{VP7}
\]
The pricing in (VP6) counts this succession of rows at
the right-hand exponent with an arbitrary fixed loss,
after binning all the new masses on a fixed finite mesh.
The target support does not have to be enumerated.

Choose \(\rho=B/(4k)\),
\(0<\tau<\min(B/16,\gamma/24,\kappa\sigma/8)\) and
\(\mu_0<\min(\tau/4,\rho/4,\gamma/8,\kappa\sigma/8)\).
Discard the classes with \(D>\tau\) or \(S>\tau\) using
(VP7). In the remaining class, any row failing (VP3)
can be placed first. Its charged cost is
\(v_1/\kappa\ge(B-D-S)/\kappa\), whereas its actual cost
is at most \(B/\kappa-\sigma\) by (VP3). This saves at
least \(\sigma-2\tau/\kappa\) in that order and absorbs
the \(\mu_0/(2\kappa)\) target allowance. The number of
orders considered is bounded at fixed \(k\).

All rows now have the compact band on all their tracked
roles. Fix an active row \(i\). The union of primes
outside \(P_0\) occurring exactly on \(i\) and one other
row has mass at least \(B-D-3S-\mu_0\). Indeed the total
incidence in \(P_0\) and on outside primes of multiplicity
at least three is at most \(3S\). There is a partner \(j\)
whose exclusively repeated packet carries mass at least
\(\rho\). Choose \(i\) penultimate and \(j\) last. Then
\(v_i\ge\rho\), because the exclusively repeated packet
has no earlier occurrence.

In any tracked role on row \(i\), call the primes in
its band eligible when they are outside the weight pool
and both target supports and occur on exactly two active
rows. The eligibility deletion removes at most
\(3S+\mu_0\) from the band, leaving at least \(\gamma/2\).
Say that \(j\) occupies a role if at least half of this
eligible band mass on \(i\) belongs to its exclusive
packet with \(j\). On an occupied role, at least
\(\gamma/4\) of the mass is available. Partition by its
signed role on \(j\). One part has mass at least
\(\gamma/12\) and numerical weight at least
\(r_0\gamma/12\). Fix \(q\) so that
\[
                           (q+1)r_0\gamma/12>2.
\]
By (VP4) this forces equality of the corresponding signed
complete values. If \(j\) occupies two roles in the global
case, the two rows are identical. If \(j\) occupies a
single free role in the fixed-letter case the same follows
using \(T\). Both outcomes contradict distinctness.

In the global case at least two roles are therefore
unoccupied. In the fixed-letter case both free roles
are unoccupied. On each unoccupied role, eligible compact
mass at least \(\gamma/4\) is absent from the last row.
Its primes are last occurrences at step \(i\), so (VP5)
recovers the known signed value modulo a product larger
than \(2X\). Thus two complete values are fixed on the
global row and at least one free value is fixed on the
fibre row. Either fixes \(x_i\). The only extra choices
at this step are the marked primes and roles from known
earlier denominators, at subpower cost. One can drop the
eligibility and absence restrictions after fixing this
record in an upper bound. Hence the ordinary cost
\(v_i/\kappa\ge\rho/\kappa\) can be deleted. Choose the
fixed mesh, the small-budget threshold and the inner
height allowances below this gap and the discard margins.
This proves (VP1) and (VP2).

When a weighted equality pattern leaves only one active
row, outside \(P_0\) the denominator support directly
records its outside radical, so (VP6) bounds that
row at exponent \(\mu_0/\kappa\) with arbitrarily small
fixed allowances. The cost of zero active rows is only
the one-row bounds for the cancelled groups.

## Powerful triples and the remaining incidence question

For prescribed powerful orders \(p_a,p_b,p_c\), the primitive
triple has radical at most \(X^{1/p_a+1/p_b+1/p_c}\).
On a fixed complete target the sum of the two free
reciprocal order bounds replaces this budget in (VP2).
Sort the original summands when using these assertions.
If their orders differ, split by the orientation, so
order sorting does not alter any exponent bound.
These are estimates for relations among genuine integer
representations, with a variable output. They provide a
saving on each complete-target graph from the powerful
problem whenever the available sums have small denominator
radical outside a prescribed pool.

They also apply to the uniquely minimally repaired
primitive quotients by evaluating \(A_s\) on the quotients.
This requirement matters for nonhomogeneous reciprocal
coordinates. Multiplying a row by an arbitrary scalar
multiplies \(A_s\) by the inverse \(s\)-th power, which
varies between the two coordinates. A direct use on
repaired complete values would therefore require a new
argument for the different effective weights. The result
here counts the primitive quotients once each.

The general small-denominator odd-moment bounds of Bloom
and Kuperberg price the numerical size of fraction
denominators. They do not supply a radical-sensitive
replacement when the output denominators and input powers
are numerically tall. Here it is the two linked polar
orders, together with the small-budget one-row price, that
allow recovery of complete values. Large families of
powerful representations are still needed to force a
large tuple count in the small-denominator output sector.
That incidence assertion has not been derived from
projected graph density or complete-target multiplicity.
[Bloom and Kuperberg, Theorem 3](https://arxiv.org/abs/2312.09021)
· [Fixed-target counts with arbitrary common divisor](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-fixed-target-estimate-with-a-repaired-target)
· [The minimally repaired primitive quotient equivalence](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)


## Six moving reciprocal sums on the actual least-repaired letters

The variable-target result also has a form on the actual repaired
values. Let \(Z(x)\) be the prescribed least repair of \(x\), for
fixed orders at each budget. More generally the following proof
allows any one prescribed multiplier per quotient, supported on
\(abc\), with \(Z(x)\le X^J\) for a fixed \(J\). For least repairs
this height bound is automatic when the product of \(L\) and the
largest prescribed order is bounded. In role \(W\) write
\(\widetilde W=Z(x)W\), and for each role retain the two sums
\[
        t_{W,s}=\sum_{i=1}^k w_i\widetilde W_i^{-s},
                             \qquad s=q,q+1 .
\]
Assume the union of the prime supports of the six reduced
denominators, outside the prescribed coefficient pool, has
product at most \(n^{\mu_0}\). For \(q\) allowed to depend on
\(J\), the bounds (VP1) and (VP2) hold with the same original
radical hypotheses, and with constants permitted to depend
on \(J\). Here \(T\), if fixed for (VP2), is fixed on the
primitive quotient. Each primitive quotient is counted once.

The pricing is on the original row supports throughout.
A singleton at an original prime in role \(W\) survives in the
corresponding output denominator, since it divides the repaired
value in that role and is absent from every repaired value on
other rows. Thus (VP7) and the discard of rows without (VP3)
are unchanged. For the remaining row \(i\), use the same
exclusive partner \(j\), eligible original compact bands,
and exposure order. If an eligible prime on \(i\) in role
\(W\) lies in the exclusive packet with \(j\), integrality
of the repaired role sum at that prime forces \(p\) to
divide \(\widetilde W_j\). The difference of consecutive
weighted reciprocal powers then gives (VP4) for
\(\widetilde W_i,\widetilde W_j\) in that role. Its order
on the repaired \(i\)-value is at least the original order.
On every occupied role this therefore forces
\(\widetilde W_i=\widetilde W_j\), with \(q\) chosen so that
\((q+1)r_0\gamma/12>2+J\). Two occupied roles force agreement
of the entire repaired triples and their primitive quotients.

There are at least two unoccupied roles in the global case
and at least one free unoccupied role in the fixed-letter
case. On each such role the eligible compact primes absent
from \(j\) are absent from all its repaired values as well.
Use (VP5) with the known prefix sums for that repaired role,
whose denominator again has polynomial height. The resulting
modulus exceeds the height of the repaired value. Two recovered
values give the whole repaired row and its quotient. On a fixed
original \(T\)-fibre one recovered free repaired value suffices
up to \(X^{o_X(1)}\) choices. Its multiplier must be a divisor
of this value. Guess that divisor and divide to recover the
original free value. It and the fixed \(T\) determine the row.
Thus the replacement of the charged cost of the intermediate
row by a subpower count survives. The two signed pooled sums
on the original primitive letters use fewer coordinates. The
six separate sums on repaired letters allow the multiplier to
vary without changing the original prime budget.

