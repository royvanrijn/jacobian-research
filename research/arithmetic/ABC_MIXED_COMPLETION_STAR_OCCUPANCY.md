# Mixed orientation in actual completion stars

A same-target star can fix complete summands from both free roles.
One orientation can be sampled many times to recover its available
numerical modulus. A single anchor in the other orientation can
then supply the missing part. The graph cost charges the two
dictionaries separately, and the modulus retains the actual
target's own high-valuation tail.

## Mixed stars with exact expected modulus

Use the parameters of
[PL5](ABC_POPULAR_COMPOSITE_PACKET_LCM.md#forced-high-valuation-tails-in-common-completion-stars).
In particular \(n=X^L\), all rows are positive primitive
\(A+B=C\) with \(X/2<C\le X\), each radical is at most \(n^{C_0}\),
and \(L\) is sufficiently small and fixed before the height limit.
Use \(\theta,\sigma\) from PL5 at a fixed \(r>0\). Fix \(R>0\) and
write \(Q(W)\) for the product of full powers \(p^e\parallel W\)
with \(Le\ge R\). For a target \(C\), let \(m_C\) be its number of
rows and put \(\gamma_C=\log_X Q(C)\). In its fibre set
\[
 f_{A,C}(q)=\frac{\#\{(A,B,C):q\mid Q(A)\}}{m_C},\qquad
 f_{B,C}(q)=\frac{\#\{(A,B,C):q\mid Q(B)\}}{m_C}.
\]
Here \(q=p^j\) runs over all depths, including the lower depths of
each prime power selected in the tail. Set \(w_q=\log_X p\).
Primitivity gives \(f_{A,C}(q)+f_{B,C}(q)\le1\) and makes the
target's tail coprime to every tail used on its fibre.

Fix nonnegative integers \(k,\ell\) with \(s=k+\ell\ge1\).
Sample independently \(k\) rows from which to fix the \(A\) values
and \(\ell\) rows from which to fix the \(B\) values. A fixed
\(A\) exposes the moving \(B\)-tail, and a fixed \(B\) exposes
the moving \(A\)-tail. The exact expected numerical LCM weight
of the exposed tails together with the target tail is
\[
 H_{k,\ell}(C)=\gamma_C+
 \sum_q w_q
 \bigl[1-(1-f_{B,C}(q))^k(1-f_{A,C}(q))^\ell\bigr].
 \tag{MSO1}
\]

**MSO2.** For fixed \(\zeta,u>0\), consider targets with
\(H_{k,\ell}(C)\ge1-\theta+\zeta\) and \(m_C\ge n^u\).
If the \(A,B\) dictionaries have sizes at most \(n^a,n^b\),
respectively, the number of these targets is
\[
 O\left(
 n^{k(a-u)+\ell(b-u)+r-\sigma+2/R}X^\varepsilon
 \right).
 \tag{MSO2}
\]
More generally their sum of \(m_C^s\) is
\(O\bigl(\left|\mathcal T_{k,\ell}\right|n^{r-\sigma+2/R}
X^\varepsilon\bigr)\), where \(\mathcal T_{k,\ell}\) is the set
of distinct mixed ordered anchor tuples occurring in their
stars. The constants are uniform in the dictionaries and
depend on the fixed star sizes, the fixed parameters of PL5,
\(u,\zeta,\varepsilon\). They do not require any assumption
on the prime supports of the individual high tails.

For (MSO1), sum the probability of capture of each depth over
all depths of each prime. The LCM weight is at most \(s+1\), and
the expectation hypothesis guarantees that a fixed positive
fraction of stars have weight at least \(1-\theta\).
Removing repeated rows loses \(O_s(1/m_C)\), which tends to zero
uniformly on the cohort. A fixed target together with a typed
anchor tuple specifies its chosen rows injectively. Now fix a
tuple with \(A\) entries \(x_1,\ldots,x_k\) and \(B\) entries
\(y_1,\ldots,y_\ell\). The moving values are
\[
 W,\quad W-x_1,\ldots,W-x_k,\quad W-y_1,\ldots,W-y_\ell,\qquad W=C.
\]
They are integral affine forms of slope one. The first of the
\(s\) anchors puts the candidate targets in one original fixed
complete-letter fibre. All three original caps still hold in
that fibre. Apply PL5 to these forms, counting only completions
of the mixed star. This gives the bound on \(\sum_C m_C^s\).
There are at most \(n^{ka+\ell b}\) possible tuples. Dividing
by \(n^{su}\) proves (MSO2).

One can use any fixed \(\eta>0\) in place of \(r-\sigma+2/R\)
on the right as an upper exponent, by choosing \(r>0\) small,
then \(R\) large and \(L\) sufficiently small. This choice
may be made before the fixed star sizes, with counting constants
allowed to depend on them. For a target cohort of exponent \(c\)
there is a strict saving if \(k(a-u)+\ell(b-u)+\eta<c\) by a
fixed margin. Use \(X\to\infty\) at each fixed \(L\), choosing
the inner allowance \(\varepsilon\) in terms of that budget.

## Union bounds with two depth moments

In this section suppress the target subscript on the frequencies.
Put
\[
 \omega_D=\sum_q w_q f_D(q),\qquad
 \xi_D=\sum_q w_q f_D(q)^2\quad(D=A,B),\qquad
 \chi=\sum_q w_q f_A(q)f_B(q).
\]
Thus \(\chi\) is the mean numerical common weight of an \(A\)-tail
and a \(B\)-tail from independent rows on this same target.

**MSO3.** For \(k=\ell\ge1\), write
\(\Omega=\omega_A+\omega_B\) and
\(\Xi=\xi_A+\xi_B+2\chi\). If \(\Omega>0\), then
\[
 H_{k,k}(C)\ \ge\
 \gamma_C+\frac{\Omega^2}{\Xi}
     \left[1-\left(1-\frac{\Xi}{2\Omega}\right)^{2k}\right].
 \tag{MSO3}
\]
Indeed \((1-f_A)(1-f_B)\le(1-(f_A+f_B)/2)^2\).
Apply the convexity argument of
[PL8](ABC_POPULAR_COMPOSITE_PACKET_LCM.md#a-sharp-occupancy-criterion-retaining-the-target-tail)
to frequency \((f_A+f_B)/2\). Equality holds in this lower bound
when \(f_A=f_B\) and their positive sum is constant across depths.
To see the limiting scope, take abstract blocks of equal weight
with total weight \(0<J\le1\). In an independent abstract row, assign
each block to \(A\) with probability \(t/2\), to \(B\) with
probability \(t/2\), and to neither with probability \(1-t\).
This gives \(\Omega=Jt,\ \Xi=Jt^2\) for \(0<t<1\).
The union of any fixed number of samples has weight at most
\(J=\Omega^2/\Xi\). The example treats depth blocks as abstract
sets and gives no arithmetic realization.

There is also a one-expensive-anchor bound. For \(k\ge1\) and
\(\ell=1\), set
\[
 \widetilde\omega_B=\omega_B-\chi,\qquad
 \widetilde\xi_B=\sum_q w_q(1-f_A(q)) f_B(q)^2.
\]
The exact gain over \(\gamma_C+\omega_A\) in (MSO1) is
\(\sum_q w_q(1-f_A(q))[1-(1-f_B(q))^k]\).
If \(\widetilde\omega_B>0\), it is at least
\[
 \frac{\widetilde\omega_B^2}{\widetilde\xi_B}
 \left[1-\left(1-\frac{\widetilde\xi_B}{\widetilde\omega_B}\right)^k\right].
 \tag{MSO4}
\]
This is the same Jensen bound with weights \(w_q(1-f_A(q))\).
In particular, let \(J=1-\theta-\gamma_C-\omega_A>0\).
Suppose \(J\ge J_0>0\), \(\widetilde\omega_B\ge w_0>0\)
and \(\widetilde\xi_B\le(1-d)\widetilde\omega_B^2/J\)
for a fixed \(0<d<1\). Since \(J\le1\), take any fixed integer
\(k\ge2/(d w_0)\). The elementary inequality
\(1-(1-t)^k\ge kt/(1+(k-1)t)\), followed by weighted
Cauchy–Schwarz, gives gain at least
\(J/(1-d/2)\). Thus (MSO2) applies with
\(\zeta=J_0d/(2-d)\). This uses \(k\) anchors from the cheap
dictionary and one from the other, even when the individual
opposite tail overlaps \(\xi_A,\xi_B\) are large.

## Asymmetric exponent comparison and cross shifts

Randomly orienting each sampled row gives the averaged law
\((f_A+f_B)/2\), with both orientations charged at the cap of
the union dictionary. At equal dictionary caps this gives the
two-moment lower bound in (MSO3) by the PL7 and PL8 argument.
The exact fixed-mix expectation (MSO1) may be larger. The
separate costs in (MSO2) matter when the two dictionary caps
differ. As an exponent comparison, take
\(u=3/20\), \(c=1/3-u=11/60\), \(a=1/6\), and \(b=1/4\).
Consider abstract depth frequencies with \(\gamma_C=0\),
one \(A\)-block of numerical weight \(63/100\) present on every
row, and three disjoint \(B\)-blocks, one present on each row
with equal frequencies. Their total available \(B\)-weight is
\(29/50\) and they are disjoint from the \(A\)-block. Neither
one-sided union has weight above \(63/100\). The mixed sampling
\((k,\ell)=(3,1)\) gives
\[
       H_{3,1}=\frac{63}{100}+\frac{29}{50}
                  \left[1-\left(\frac23\right)^3\right]
                =\frac{2803}{2700}>1 .
\]
The cost \(3(a-u)+(b-u)=3/20\) is below \(c\). Charging
the larger dictionary cap on the symmetrized union would already
cost \(2(b-u)=1/5>c\) for two anchors, and one anchor cannot
reach numerical weight one. These are strict comparisons for
the supplied parameters, with \(\theta\) chosen small and the
counting allowances below the margins. No realization of these
frequencies by a three-cap arithmetic graph is asserted.

**MSO5.** Large \(\chi\) puts a deep divisor on many normalized
cross shifts. For two independently sampled rows \((A_i,B_i,C)\)
and \((A_j,B_j,C)\), let
\[
       Q_{ij}=\gcd(Q(A_i),Q(B_j)),\qquad
       h_{ij}=|C-A_i-A_j|=|B_i-A_j|=|B_j-A_i|.
\]
If \(\chi\ge\eta>0\), a fixed positive proportion of pairs
has \(\log_X Q_{ij}\ge\eta/2\), since this weight is at most
one. The diagonal has \(Q_{ij}=1\) by primitivity. The pairs
with \(h_{ij}=0\) are swapped rows, with at most one such \(j\)
per \(i\), and can be deleted at \(O(1/m_C)\) proportional cost.
For a remaining pair \(Q_{ij}\mid h_{ij}\).
Original primitivity makes \(Q_{ij}\) coprime to \(B_i\), to
\(A_j\), and to \(C\). It therefore divides
\(h_{ij}/\gcd(B_i,h_{ij})\), the shift on the normalized
primitive difference of the opposite endpoints \(B_i,A_j\).
The completed full divisor on those primes has the same radical
at most \(n^{1/R}\), retains depth at least \(R/L\) at its primes,
and retains numerical weight at least \(\eta/2\).
This implication preserves an actual high-valuation part after
normalization. A signed radical cap on its remaining primes is
still needed for the
[three-letter shift count](ABC_COMPRESSED_SHIFT_TARGET_TAILS.md#endpoint-support-and-original-target-discounts).
