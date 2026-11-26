# A dyadic analytic interpolation for quadratic moment tails

Let
\[
 f=\sum_i(\alpha_iX_i^2/2+d_iX_i)+\sum_{i<j}\beta_{ij}X_iX_j,
 \qquad       P_n=\Lambda_m(f^n)/n!,
\]
with zero constant and integral natural coefficients in a
finite extension of \(\mathbb Q_2\). Let \(v(2)=1\).
There are a positive integer \(L\), constants \(A,B\),
and locally analytic functions \(F_a\) on \(\mathbb Z_2\)
for \(a=0,\ldots,L-1\), with values in a finite extension,
such that
\[
       P_n\equiv F_a((n-a)/L)\pmod{2^k}
               \quad(n\equiv a\bmod L,\ n\ge A k+B).        \tag{1}
\]
The constants and finite-precision coefficients of the
functions come from the native quadratic coefficients.
This separates the normalized moments into the
interpolation of their unit-root contribution and an
error with valuation increasing at least linearly in \(n\).

For a proposed window length \(w\), write
\(a+j=a_j+L b_j\) with \(0\le a_j<L\), for \(0\le j<w\).
If the functions
\[
                       F_{a_j}(z+b_j),\qquad 0\le j<w      \tag{2}
\]
have no common zero on \(\mathbb Z_2\) for every \(a\),
then all sufficiently late windows of length \(w\)
contain a nonzero exact moment. Conversely, infinitely
many exact zero windows force a common zero in (2)
for at least one \(a\). A finite-precision
certificate for the exclusion is given below.

## Degree bounds for the rational approximations

Use the [fixed quadratic pole theorem](QUADRATIC_FIXED_DYADIC_POLE_SUPPORT.md).
Let \(U_i=1-\alpha_iT-d_i^2T^2\), and choose the
common denominator \(D_k=\prod_i U_i^{2k-1}\) in the
rational approximation modulo \(2^k\).
If its numerator is \(N_k\), then we can take
\[
 d_k=\deg D_k\le 2m(2k-1),\qquad
 \deg N_k\le m(4k-3)+q+2(k-1),                             \tag{3}
\]
where \(q\) is the number of nonzero cross monomials.

Here is the numerator bound. In the vertex construction
write \(\mathcal L^a F=P_a/U^{2a+1}\).
The differential operator
\(\mathcal L=(\alpha T^2/U)\partial_T\)
increases the numerator degree by at most three.
Thus \(\deg_T P_a\le1+3a\). Each increment in denominator
power in the weighted differentiation can add at most
two to the numerator degree, so a retained term
with denominator \(U_i^{2a_i+1+b_i}\) has numerator
degree at most \(1+3a_i+2b_i\).
Clearing that vertex against \(U_i^{2k-1}\)
gives degree at most \(4k-3\).
The cross shift is at most
\(q+2\sum_e v_2(j_e!)\le q+2(k-1)\)
in every retained array, because
\(j_e\le1+2v_2(j_e!)\) for its nonzero counts.
Summing proves (3). Division by scalar factors in the
construction does not change degrees, and integrality
of the resulting numerator is part of the pole theorem.

Let \(\mathscr S\) be the forward shift on sequences.
With \(D_k^\sharp(X)=X^{d_k}D_k(1/X)\), equation (3)
gives the recurrence
\[
                       D_k^\sharp(\mathscr S)P_n=0\bmod2^k
                                  \quad(n\ge n_k),
\]
where \(n_k=\max(0,\deg N_k-d_k+1)\).
Indeed this is the coefficient at \(n+d_k\) in the
product of \(D_k\) with the generating series.
Both \(n_k\) and \(d_k\) have a linear bound in \(k\).

## The two reciprocal-root parts and their compatibility

In a finite splitting extension, split \(D_k^\sharp\) into
its monic factors \(U_k^\sharp,V_k^\sharp\) for unit and
nonunit roots respectively, with multiplicities.
All roots belong to the fixed set of reciprocal roots of
the \(U_i\). Their residues separate the two groups, so the
resultant is a unit. Thus the monic factors have a
Bézout relation with integral coefficients. For \(n\ge n_k\),
apply its two summands, as polynomials in \(\mathscr S\),
to \(P_n\bmod2^k\). This decomposes the sequence as
\(u_k+z_k\), with
\[
 U_k^\sharp(\mathscr S)u_k=0,\qquad
                      V_k^\sharp(\mathscr S)z_k=0 .
\]
The unit sequence extends backwards to all nonnegative
indices by the unit constant coefficient of its recurrence.
If one of the monic factors is \(1\), the corresponding
sequence is zero.

We use a direct recurrence estimate. For a monic
\(F(X)=\prod_{i=1}^d(X-\mu_i)\) with
\(v(\mu_i)\ge\epsilon>0\), the coefficient of \(X^i\)
in the remainder of \(X^t\) modulo \(F\) has valuation
at least \((t-i)\epsilon\) for \(i<d\).
The coefficients of \(F\) have the corresponding bounds
\((d-i)\epsilon\), and repeated substitution of
\(X^d\) proves the assertion by induction on \(t\).
In particular all those coefficients are zero modulo
\(2^k\) when \(t\ge d+\lceil k/\epsilon\rceil\).

Choose \(\epsilon>0\) a lower bound for the valuations
of all nonunit reciprocal roots of the fixed vertex
denominators. If there are none, the nonunit sequence is
already zero. Otherwise the estimate shows
\(z_k(n)=0\bmod2^k\) for
\(n\ge n_k+\deg V_k^\sharp+\lceil k/\epsilon\rceil\).
By (3) this is a linear threshold \(Ak+B\).
In particular \(u_{k+1}\equiv u_k\bmod2^k\) on a
sufficiently late interval, since both agree with \(P_n\).
Their difference satisfies the product of their two
unit recurrences and so propagates backwards modulo
\(2^k\). Hence the congruence holds at every nonnegative
index. Define \(u(n)\) to be the compatible limit of
these unit sequences in the complete valuation ring.

## Finite differences and local analyticity

Choose \(L\) so that every fixed unit reciprocal root
\(\mu\) satisfies \(v(\mu^L-1)\ge1\).
A residue-field order first gives positive valuation;
successive squaring raises this valuation to at least
one by \(v((1+x)^2-1)\ge\min(1+v(x),2v(x))\).
This requires only a finite splitting extension.

For \(Y=\mathscr S^L-1\), the unit sequence \(u_k\)
is annihilated modulo \(2^k\) by
\[
               \prod_{\mu\ {\rm unit}}
                         (Y-(\mu^L-1))^{m_{\mu,k}},
\]
where \(m_{\mu,k}\) is its multiplicity in \(U_k^\sharp\).
Indeed after substitution of \(Y\), the polynomial
in \(\mathscr S\) is a multiple of \(U_k^\sharp\).
The preceding remainder estimate with \(\epsilon=1\)
shows that \(Y^r u_k=0\bmod2^k\) whenever
\(r\ge\deg U_k^\sharp+k\).
By (3), \(r\ge(4m+1)k\) suffices, so for every
nonnegative \(a\),
\[
                v\bigl((\mathscr S^L-1)^r u(a)\bigr)
                               \ge \left\lfloor
                                           \frac r{4m+1}\right\rfloor .
\]
Define
\[
 F_a(z)=\sum_{r\ge0}(\mathscr S^L-1)^r u(a)\binom zr .
\]
At nonnegative integer \(z\) this interpolates
\(u(a+Lz)\) by finite binomial inversion.
For analytic convergence, on a disk \(z=b+2^j y\),
with \(b\) an integer representative and \(|y|\le1\),
the \(r\) numerator factors of \(\binom zr\) have
total Gauss valuation at least
\(\sum_{l=1}^j\lfloor r/2^l\rfloor\), since among \(r\)
consecutive integers at least \(\lfloor r/2^l\rfloor\)
are divisible by \(2^l\).
By Legendre's formula, division by \(r!\) therefore gives
Gauss valuation at least
\(-\sum_{l>j}\lfloor r/2^l\rfloor\ge-r2^{-j}\).
Choosing \(2^{-j}<1/(4m+1)\) makes the displayed series
converge in the algebra of restricted power series in
\(y\). This proves local analyticity, and the comparison
between \(u(n)\) and \(P_n\) proves (1).

## The window test

For each \(a\), the absence of a common zero in (2)
gives a positive lower bound for the maximum of their
absolute values on the compact domain \(\mathbb Z_2\).
There are finitely many \(a\), so choose \(k\) with
\(|2^k|\) smaller than all these bounds and then use (1).
A putative zero window sufficiently late would make
every value in its vector divisible by \(2^k\).
Conversely, an infinite set of zero windows has a
subsequence in one class \(a\), whose start parameters
have a limit in \(\mathbb Z_2\). By (1) their vector
values tend to zero, giving a common zero by continuity.

The test can be certified at finite precision:
compute \(u_k\) by the split recurrence and truncate
the binomial series at \(r<(4m+1)k\).
For the resulting finite list of binomial polynomials,
their values modulo \(2^k\) are fixed by the input
residue modulo \(2^{\,k+\max_r v_2(r!)}\).
A simultaneous nonvanishing test on those finitely many
residues certifies the desired bound. The interpolation
may vanish identically, for example when there are no
unit reciprocal roots. The criterion separates a
concrete analytic obstruction from that degenerate case.

## Continuous interpolation for arbitrary support

The continuity and eventual window conclusions apply at
every prime to the arbitrary finite supports in the
[fixed singleton pole theorem](FACTORIAL_MOMENT_SERIES_PRIME_POWER_RATIONALITY.md#which-poles-can-occur-at-higher-precision).
In its notation, assume the natural coefficients are
integral and \(c=0\) or \(v(c)>1/(p-1)\).
For each \(k\), use the rational expression modulo \(p^k\)
with denominator \(D_k=\prod_i U_i^{D_i(k)}\).
Every nonconstant \(U_i\) has constant coefficient one
and unit leading coefficient. The latter comes from the
unique singleton at the smallest retained position.
Consequently the reciprocal roots are all units and
belong to a fixed finite set. Indeed their monic
polynomials are integral and have unit constant
coefficients.

For each precision the recurrence with characteristic
polynomial \(D_k^\sharp\) holds beyond the numerator
threshold, and its constant coefficient is a unit.
It extends the tail uniquely backwards to a sequence
\(u_k\) modulo \(p^k\). Comparing two precisions using
the product recurrence shows that the extensions are
compatible, since they agree on their common tail.
Thus they define a limit \(u\) that agrees with \(P_n\)
modulo \(p^k\) for sufficiently large \(n\) at each \(k\).

Choose a fixed \(L\) such that \(v(\mu^L-1)\ge1\)
for every reciprocal root \(\mu\). A common residue order
followed by enough powers of \(p\) has this property,
using \(v((1+x)^p-1)\ge\min(1+v(x),p v(x))\).
The same recurrence estimate gives
\[
       (\mathscr S^L-1)^r u(a)\equiv0\bmod p^k
                              \quad(r\ge \deg D_k+k).
\]
The binomial series therefore converges uniformly on
\(\mathbb Z_p\), as \(\binom zr\) is integral there.
It gives continuous functions \(F_a\), with
\(P_n\equiv F_a((n-a)/L)\pmod{p^k}\) for sufficiently
large \(n\) in every class \(a\bmod L\).

The same compactness argument gives the eventual window
criterion (2) over \(\mathbb Z_p\). At any precision \(k\)
a finite certificate uses \(r<\deg D_k+k\) and input
residues modulo \(p^{\,k+\max_r v_p(r!)}\).
This general case gives continuity and a threshold
at each precision. The quadratic degree bounds above
supply the linear threshold and local analyticity.
