# Endpoint congruences for fully shared mirrored factorial moments

Let a fully shared ordinary mirrored profile contain groups of
four rows
\[
             (a_s,b_s),\ (b_s,a_s),\
             (Q_s-a_s,Q_s-b_s),\ (Q_s-b_s,Q_s-a_s),
 \tag{1}
\]
with \(Q_s>\max(a_s,b_s)\) and all \(a_s,b_s\) positive integers.
The remaining rows form a fixed mirrored core, possibly empty.
Choose an integer \(J\ge3\) at least the
[absolute eventual ordinary cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Take an odd prime \(p>J\) that exceeds \(J\) times the largest
core entry, when the core is nonempty. Suppose:

- For every group, \(a_s,b_s\) are \(p\)-units and
  \(v_p(Q_s)>\lfloor\log_p(J\max(a_s,b_s))\rfloor\).
- For each \(2\le j<J\) and \(1\le k\le\lfloor j/2\rfloor\),
  at least one group satisfies
  \[
                       j a_s+(b_s-a_s)k\equiv0\pmod p .
  \tag{2}
  \]

Then the actual profile satisfies SFC(2) in every degree.
This is a condition on the integer rows of the profile itself.
It does not require estimating a real or Gaussian comparison.
For every mirrored core such groups can be chosen fully shared
with arbitrarily small total added curvature.

## Exact paired complementary factors

In degree \(m\), put \(x=ma_s+(b_s-a_s)k\) and
\(y=mb_s-(b_s-a_s)k\). Both lie strictly between zero
and \(mQ_s\). A complementary row pair has the normalized identity
\[
             \frac{(-1)^{x-1}mQ_s}{\binom{mQ_s}{x}}
                  =\frac{x}{\prod_{h=1}^{x-1}(1-mQ_s/h)} .
\]
[Complementary row identity](MULTIVARIATE_COMPLEMENTARY_ROW_EULER_REDUCTION.md#a-complementary-pair-over-a-retained-core)
Applying it to \(x\) and \(y\) accounts for all four rows
in (1). Their two signs combine into \((-1)^{m(a_s+b_s)}\),
which is independent of \(k\). After one scalar normalization
the group therefore multiplies coefficient \(k\) by
\[
 \frac{xy}{\prod_{h=1}^{x-1}(1-mQ_s/h)
              \prod_{h=1}^{y-1}(1-mQ_s/h)}\equiv xy\pmod p
                  \qquad (m\le J).
 \tag{3}
\]
Each denominator factor is a \(p\)-unit congruent to one,
since \(h\le J\max(a_s,b_s)\) and \(Q_s\) satisfies the
specified valuation bound.

At \(k=0,m\) the value in (3) is \(m^2a_sb_s\) modulo \(p\),
which is a nonzero unit. At each internal index, one group
kills the coefficient by (2), either directly or by reflection
\(k\mapsto m-k\). The core has equal endpoints which are
\(p\)-units by the prime bound. Thus every normalized ordinary
polynomial of order \(m<J\) reduces to a unit multiple of
\(1+z^m\). These reductions have their original degrees and
are coprime at consecutive orders: a common root would
satisfy \(z^m=z^{m+1}=-1\), giving \(z=1\), impossible at
odd \(p\). The rational polynomials are therefore coprime
at every terminal below \(J\), and the eventual theorem
covers the remaining terminals. Their rational Bézout
identities give the result over every characteristic-zero field.

## Realizing the congruences

One group per reduced fraction \(k/j\) in (2) suffices.
Choose any large \(a_s\equiv k\pmod p\), and take
\(b_s-a_s=p-j>0\). Then \(b_s\equiv k-j\pmod p\) is a unit,
and (2) also holds for every occurrence of the same reduced
fraction. Choose \(Q_s\) as a sufficiently large power
of \(p\), in particular at least \(2b_s\). The four rows
are ordinary and shared. Their curvature is at most
\(4(b_s-a_s)^2/a_s\), which tends to zero as the
\(a_s\) increase in the chosen residue classes.
