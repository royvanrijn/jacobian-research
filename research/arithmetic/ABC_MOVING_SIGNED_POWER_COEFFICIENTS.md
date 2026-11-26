# Moving signed power coefficients in a full-modulus chart

Two letters can have moving coefficients without a count for their
coefficient dictionaries. The argument applies when their two bases
have one exact common order and their signed coefficient ratio has
sufficiently small numerical height. A rectangular determinant gives
an exact curve in each residue chart. Low-degree factors reduce to
unit root classes with controlled multiplicity, while high-degree
factors have few rational points at the relevant height.

## A count with unrecorded coefficients

For a reduced rational number \(z\), write \(H(z)\) for the larger
of the absolute numerator and the positive denominator. Fix
\(0<\kappa_-\le\kappa_+\), \(\sigma,w>0\), \(0<\mu_0\le\mu\),
a positive integer \(a\), and a counting allowance \(\eta>0\).
Write \(q_0=a\mu_0\).
Put \(n=X^L\) with the height limit taken at each fixed \(L>0\).
Use positive primitive \(abc\) rows in an original shell of height
at most \(X\). Orient two chosen complete letters \(U,V\) as
\[
        U=C_U x^k,\qquad V=C_V y^k,\qquad
        \kappa_-\le kL\le\kappa_+ .
\]
The integer roots \(x,y\) have prime support in their own complete
roles. The positive rational coefficients \(C_U,C_V\) also have
support in their own roles. This includes signed valuation
remainders. Define \(r=x/y\) and \(s=\pm C_V/C_U\), with the sign
chosen so that \(r^k\equiv s\pmod D\) when \(D\) divides the third
letter. The common order \(k\) can vary by row.

Assume that the third letter has a divisor \(D\) such that
\[
 D\ge X^w,\quad n^{\mu_0}\le \operatorname{rad}(D)\le n^\mu,\quad
 p^e\parallel D\ \Longrightarrow\ Le\ge\sigma,\ p>P_0(L).
\]
The cutoff \(P_0(L)\) is fixed at each budget and is specified
below. All coefficients and roots on the other roles are units
at \(D\) by primitivity. Retain any original radical and
eligibility conditions.

**MPC1.** There are \(c>0,L_0>0\) and a cutoff \(P_0(L)\)
depending on the indicated fixed parameters such that, for every
\(0<L<L_0\) fixed, the extra condition
\[
                            H(s)\le X^{c kL}
\]
gives the count
\[
                    O\!\left(n^{(a+1)\mu+\eta}X^\varepsilon\right)
                           \qquad(\varepsilon>0).
\]
No coefficient value or coefficient support is recorded. In
particular, if \((a+1)\mu<1/3\) by a fixed margin, this gives a
strict small-budget saving on the stated moving-coefficient
sector. The local counting proof also works directly with
unit rational pairs \(r,s\), given \(H(r)\le X^{1/k+cL}\)
and \(r^k\equiv s\bmod D\).

To deduce the bound on \(H(r)\) in the original application,
the prime supports of \(C_U,C_V\) are disjoint. Every signed
numerator and denominator of either coefficient divides one of
the numerator and denominator of \(s\). Thus
\(x^k,y^k\le X H(s)\). Also \(U/V=\pm r^k/s\) with the
prescribed signs, so the coprime complete pair \(U,V\), and
hence the full row, is determined by \(r,s\).

## Choosing a rectangle and forcing an equation

Choose a positive integer \(J\) with
\[
                    q_0\kappa_-(J+1)/2>4,
\]
and then choose \(0<t<1/2\) with
\(a t\kappa_+(J+1)<\sigma/2\).
For each \(k\) put \(I=\lfloor tk\rfloor\) and
\(N=(I+1)(J+1)\). Take
\[
                   P_0(L) = \big\lceil\kappa_+/L\big\rceil(J+2).
\]
All of these choices except \(I,N\) and \(P_0(L)\) are made
before the budget is chosen. For each \(k,D\), refine by the
residue of \(r\) modulo \(\operatorname{rad}(D)^a\). There are
at most \(n^{a\mu}\) such residues. There are at most
\(n^\mu X^{o_X(1)}\) choices for \(D\), since its squarefree
support has the given cap and its powers on fixed support cost
height-subpower uniformly in \(D\le X\).
[Fixed-support enumeration (EB3–EB4)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)

In one of the resulting charts, take an \(N\)-row determinant of
the monomials \(r^i s^j\), \(0\le i\le I\), \(0\le j\le J\).
Clear each row by the product of the denominators of \(r^I\)
and \(s^J\). The cleared integer determinant has absolute value
at most
\[
 N!\,X^{\,N\{I/k+cLI+ckLJ\}}.
\]
At each \(p^e\parallel D\), put \(r_0\) for a representative
of the fixed residue of \(r\bmod p^a\).
The row polynomials agree modulo \(p^e\) with
\(r^{i+kj}\). Expand around \(r_0\) through order \(N-2\).
Because \(r-r_0\in p^a\mathbb Z_p\) and
\(e\ge a(N-1)\) once \(L\) is small, the row is a sum of the
fixed Taylor directions, of depths \(0,a,\ldots,a(N-2)\),
and an error of depth \(a(N-1)\).
A nonzero determinant term uses distinct Taylor directions.
With \(b\) error rows the depth is at least
\[
        a\{(N-b)(N-b-1)/2+b(N-1)\}\ge aN(N-1)/2 .
\]
Clearing denominators does not change local valuations because
the denominators are units at \(D\).

The lower depth summed over \(p\mid D\) is at least
\(Nq_0L(N-1)/2\) in units of \(\log X\). For sufficiently
small fixed \(L\) the lower
depth per row exceeds \(3t\), whereas the real upper bound per
row has exponent at most \(2t\) if \(c\) is small enough that
\(c\kappa_+(t+J)<t\). Thus the determinant is zero at all late
heights, uniformly in the record. Every chart lies on a
nonzero integer equation \(F(r,s)=0\) of bidegrees at most
\((I,J)\). No bound on the coefficients of that equation
is required below.

## Uniform completion on each factor

Fix a small constant \(\delta_0>0\) such that
\[
 \delta_0<t/2,\qquad
               \frac{w}{3\delta_0(J+1)}>4+2\mu\kappa_+ .
\]
Work with each nonconstant irreducible factor \(G\in\mathbb
Z[r,s]\), primitive as a polynomial. There are finitely many
factors at each fixed budget.

If \(\deg_r G\ge\delta_0 k\), the projective closure of this
factor has degree at least \(\delta_0 k\). The primitive
projective height of \([1:r:s]\) is at most \(H(r)H(s)\):
use a common denominator formed as the product of the two
reduced denominators. The coefficient-independent bound
of Binyamini–Cluckers–Kato, Theorem 1 for \(K=\mathbb Q,n=2\),
therefore counts the factor at exponent in base \(n\) at most
\[
 \frac{2(1/k+cL+ckL)}{\delta_0 kL}+o_X(1)
                  =\frac{2c}{\delta_0}+O(L)+o_X(1).
\]
This remains uniform when the equation depends on \(D\) and
\(X\). Alternatively, a rational factor that is not
geometrically irreducible has at most \(O((I+J)^2)\)
rational points by intersecting two conjugate components.
[Binyamini–Cluckers–Kato, Theorem 1](https://arxiv.org/html/2401.03982v2#S1.SS1)

If \(\deg_rG<\delta_0 k\), the one-variable polynomial
\(P(Z)=G(Z,Z^k)\) is nonzero and primitive. Its monomial
exponents \(i+kj\) are distinct because \(I<k\).
The degree is less than \(P_0(L)\), and the number of
nonzero monomials is at most
\((\delta_0 k+1)(J+1)\).
For any \(p\mid D\), its reduction is not zero.
At a nonzero residue its root multiplicity is at most one
less than the number of nonzero monomials. In fact, applying
successive powers of \(Z\,d/dZ\) at such a residue gives a
Vandermonde matrix in the surviving exponents. These exponents
are distinct modulo \(p\), so the derivatives cannot
all vanish. We may use the uniform multiplicity bound
\(m_0=\lceil 2\delta_0 k(J+1)\rceil\) at sufficiently
small \(L\).

Here is a direct root-class estimate that does not depend on
coefficient sizes or on nonsingularity of the reduction.
Factor \(P\) over an algebraic closure of \(\mathbb Q_p\).
Normalize each linear factor to Gauss norm one. The
remaining scalar is a unit, because the primitive
polynomial has Gauss norm one and nonzero reductions of
the normalized factors multiply to a nonzero polynomial.
Factors whose roots have absolute value greater than one
reduce to constants and evaluate to units at \(p\)-adic
integers. At a \(p\)-adic unit \(z\), at most \(m_0\)
of the remaining roots \(\beta\), counted with multiplicity,
can satisfy \(|z-\beta|_p<1\). They are precisely the
factors over the given residue. If \(|P(z)|_p\le p^{-e}\),
at least one of these roots satisfies
\(|z-\beta|_p\le p^{-e/m_0}\).
Any two unit solutions assigned to the same root agree
modulo \(p^{\lfloor e/m_0\rfloor}\) by the ultrametric
inequality. Thus there are at most \(\deg P\) relevant
classes at that precision.

Apply this to \(z=r\) using \(P(r)\equiv0\bmod D\).
There are at most \((\deg P)^{\omega(D)}=X^{o_X(1)}\)
joint classes modulo
\[
                \prod_{p^e\parallel D}p^{\lfloor e/m_0\rfloor}
                         \ge D^{1/m_0}/\operatorname{rad}(D).
\]
The lower bound exceeds \(2 X^{2/k+2cL}\) at late heights:
after multiplying its logarithmic exponent by \(k\), the
lower exponent is at least
\(w/(3\delta_0(J+1))-\mu\kappa_+\), and the upper
exponent is \(2+2c\kappa_+\) before accounting for the
constant \(2\). Thus a joint class determines the reduced
rational \(r\) by its cross determinant. Once \(r\) is
known, \(s\) is determined by \(s\equiv r^k\bmod D\) and
\(D>2H(s)^2\), which follows by choosing
\(2c\kappa_+<w\). Each low-degree factor therefore
costs only height-subpower.

Choose \(c\) also so that \(2c/\delta_0<\eta/4\) and
take \(L\) sufficiently small for the displayed margins.
Summing the factor counts, order choices and residue
charts proves (MPC1). All orders,
Taylor degrees and prime cutoffs are fixed at each
budget before taking the height limit. In the diffuse
lock application the lower radical bound \(n^{\mu_0}\) is produced by
the valuation band, rather than inferred from the
numerical height of a divisor with arbitrary high prime
orders.

## Which valuation profiles the count covers

For the exact signed decomposition, at each prime
\(p^e\parallel U\) choose \(e=kq_e+r_e\) with
\(|r_e|\le k/2\), and do the same independently on
\(V\). Put \(x=\prod p^{q_e}\) on \(U\); its
coefficient has valuation \(r_e\). Both powers and
coefficients are supported on their own complete letter.
The ratio \(s\) is therefore reduced by disjoint support,
and its numerator and denominator logarithms are the
two actual signed sums of these prime valuation
remainders. A sufficient geometric condition for
(MPC1) is that all degrees except a small fixed raw mass
lie within a sufficiently small fraction of \(k\) of an
integer multiple of \(k\) on the two roles. More precisely,
if the two roles have total normalized raw mass at most
\(R\), all but mass \(\rho\) have \(|r_e|\le \tau k\), and
\(\tau R+\rho/2\le c\), then
\(H(s)\le X^{ckL}\), with no restriction on the exact
coefficient dictionary. The scaled degree measure may be diffuse inside bands
of width \(2\tau kL\) about the scaled grid. For the
stated range \(kL\in[\kappa_-,\kappa_+]\), the grid
has positive spacing at the small-budget limit. The
condition does not follow from diffuseness alone.
## Centred remainder benchmark on diffuse spectra

For a valuation degree \(e\), the centred remainder divided by
\(k\) is \(e/k-\lfloor e/k+1/2\rfloor\). Its positive and
negative magnitudes are \(g_+(\{e/k\})\) and
\(g_-(\{e/k\})\), where
\[
       g_+(u)=u{\bf 1}_{[0,1/2)}(u),\qquad
       g_-(u)=(1-u){\bf 1}_{[1/2,1)}(u).
\]
Their means are both \(1/8\). For any periodic bounded-variation
function \(g\), the mean over integer \(K\le k<2K\) of
\(g(e/k)\) differs from the circle mean by
\(O_g(K/e+e/K^2+1/K)\).
For the sum-to-integral error use the variation of \(g(e/x)\)
over the integration interval, which is \(O_g(1+e/K)\).
For the integral substitute \(v=e/x\) and integrate the
mean-zero part against \(v^{-2}\), using its bounded periodic
primitive.

For \(K\asymp L^{-\alpha}\), \(1/2<\alpha<1\), these bounds
apply uniformly on every compact positive scaled band
\(Le\in[\sigma_1,\sigma_2]\).
Suppose the two complete letters have wholly nonatomic limiting
scaled-degree laws with bounded first moments. In base \(n\)
write \(r_U,r_V\) for their limiting radical masses. For a proportion \(1-o(1)\) of these choices of \(k\), the common
signed coefficient ratio satisfies
\[
                     \log_X H(s)
                    \le \left(\frac{r_U+r_V}{8}+o(1)\right)kL.
                                                               \tag{MPC2}
\]
The two signed sums concentrate separately, as can be seen
from their second moments. On each compact band approximate
\(g_+\) and \(g_-\) by finite trigonometric polynomials in
\(L^2\). The preceding uniform mean estimate applies to the
squared approximation error, and the triangle inequality in
\(L^2\) controls its average against any bounded radical-mass
law. In the variance of a centred trigonometric sum only
frequencies \(\ell e+\ell'e'\) with \(\ell,\ell'\ne0\) arise.
If \(|L(\ell e+\ell'e')|\ge\xi\) for fixed \(\xi>0\),
the same integration and variation estimates make their
mean \(o(1)\) uniformly. In the complementary region the
product measures converge weakly and the limiting resonance
\(\ell u+\ell'v=0\) has measure zero by nonatomicity. Let
\(\xi\downarrow0\) and then remove the low and high bands
using nonatomicity at zero and the first-moment bound. This
proves the variance estimate and (MPC2) simultaneously for
both signs. At balanced masses the benchmark is \(1/12\).
Here \(kL\to0\), so (MPC2) is a separate finer-grid benchmark
and does not meet the bounded-positive-\(kL\) hypothesis of
(MPC1).

## A limitation of the rational rectangular estimate

The \(1/12\) benchmark does not by itself make the
moving-coefficient rectangle work. Consider the formal
Taylor estimate that grants independent directions at the
smallest weights in the full two-coordinate Taylor lattice.
At a prime \(p^e\parallel D\), with known root residue
modulo \(p^{a_p}\), use coordinates \(r-r_0\) and \(s-r^k\).
Write \(q_p=a_p\log_n p\) and \(w_p=e\log_X p\).
The assigned numerical weights are
\(q_p L j+w_p\ell\) for \(j,\ell\ge0\). Put
\(b=NL\), \(q=\sum q_p\), and \(w=\sum w_p\).
For any \(v>0\), the number of local weights at most \(v\)
is at least \(v^2/(2q_p w_p L)\), by adding
\(\lfloor(v-w_p\ell)/(q_p L)\rfloor+1\) over the
nonnegative \(\ell\) for which \(w_p\ell\le v\)
and then using the monotone integral. Their first \(N\)
weights thus have mean at most
\((2/3+O(1/N))\sqrt{2q_p w_p b}\).
Summing and using Cauchy–Schwarz gives the formal per-row
gain at most
\[
                     (2/3+O(1/N))\sqrt{2qwb}.
\]
Even allowing \(q\le1/3\) and \(w\le1\), its leading term
is at most \((2\sqrt6/9)\sqrt b\).
The rational rectangle with \(I,J\to\infty\),
\(N=(I+1)(J+1)\), and coefficient height \(X^{ckL}\)
has ordinary cleared archimedean price
\[
      I/k+ckLJ
      \ \ge\ (1-o(1))\,2\sqrt{cb}.
\]
At \(c=1/12\) its leading minimum is
\(\sqrt{b/3}\), larger than the Taylor-lattice gain.
For the benchmark \(k\asymp L^{-\alpha}\), a bounded
\(I\) and positive \(J\) also fail because
\(ckLJ/\sqrt{NL}\to\infty\).
For bounded \(J\) and positive \(I\), the first \(N\)
weights have mean at most \(qL(N-1)/2\) after summing
the one-coordinate choice \(\ell=0\). This is smaller
than \(I/k\) as \(kL\to0\). The same one-coordinate
bound handles an identically zero side degree.

This compares the complete Taylor lattice calculation
with the usual cleared-height estimate. It is not an
upper bound for the actual valuation of an arithmetic
determinant. Additional dependencies among the local
Taylor directions or a stronger archimedean
normalization would require a separate argument.

## Comparing sufficient width bounds on a faint lock

Consider a collection of fixed outer parameter choices with
\(\mu\downarrow0\), \(\mu_0\asymp\mu\), and
\[
                   w_-\mu\le\log_X D\le w_+\mu
\]
for fixed positive \(w_-,w_+\). Keep the original lower bound
\(\sigma\) on \(Le\) and the compact positive range for \(kL\)
fixed. At each choice of these parameters take the small-budget
and height limits in the order used in (MPC1). Hold a chart
price \(q_*<1/3\) and analytic allowance \(\eta>0\) fixed
with \(q_*+\eta<1/3\). Set \(a=\lfloor q_*/\mu\rfloor\).
Then \(q_0=a\mu_0\) has a fixed positive lower bound, so \(J\)
can be fixed for this comparison. The parameters
\(t,\delta_0,c\) in the proof of (MPC1) can all be chosen to
be sufficiently small constant multiples of \(\mu\).
Indeed the constraints \(at\kappa_+(J+1)<\sigma/2\) and
\(w/(3\delta_0(J+1))>4+2\mu\kappa_+\) allow \(t,\delta_0\)
of this order, taking \(w=w_-\mu\). The upper determinant
height and the high-degree count then allow \(c\) of the
same order. Thus this certificate accommodates signed
coefficient ratios of numerical height exponent comparable
to \(\mu\), with total count charge at most
\(q_*+\mu+\eta\).

The existing two-base estimate SP3 is another certificate
for such moving ratios. Apply its rational congruence to
\(r^k s^{-1}\). In its proof, if the auxiliary fullness
parameter is held at a fixed multiple of \(1/L\), then
the product-height condition for \(H(r)\le X^{O(L)}\),
\(H(s)\le X^{ckL}\) requires \(c\ll\mu^2\) in this regime.
That is a restriction of that choice of parameters, not a
limit of SP3. Indeed the same \(D\) can instead be treated
as \(K\)-full with \(K=\lfloor u\mu/L\rfloor\), for any
fixed \(u>0\) at sufficiently small \(\mu\).
Use amplification order one. Its height coefficients can
be taken as
\[
 C_r=\max\{1,\ K(1/k+cL)/(w_-\mu)\},\qquad
 C_s=\max\{1,\ KckL/(w_-\mu)\}.
\]
Then \(C_r=O_u(1)\) and \(C_rC_s/K=O_u(L/\mu+c/\mu)\).
For a sufficiently small constant multiple \(c\) of \(\mu\)
this meets the determinant condition in SP3, with \(L\)
small at the fixed \(\mu\). Moreover the count
\(D^{1/K+\varepsilon}\) costs at most
\(2w_+/u+o(1)\) in base \(n\). Choosing \(u\) large also
gives a cut below \(1/3\) after adding the raw support charge.
Thus both certificates permit linear numerical width in
this fixed-lock comparison. They use different completion
arguments.
[Signed two-base rational congruence (SP3)](ABC_SIGNED_VALUATION_PHASES.md#from-shifted-phase-concentration-to-small-height-product)

For a formal profile satisfying the linear width condition,
take \(kL\to1\) and spread radical mass \(1/3\) on each of the
two roles over the exact positions
\(0<|e-3k|\le\theta k\), with \(\theta\) a sufficiently small
fixed multiple of \(\mu\). Centered rounding gives a moving
coefficient ratio with logarithmic height exponent at most
\(2\theta kL/3\). The scaled profile can be wholly diffuse
and every occupied position is off the exact grid. Within
disjoint prime palettes all those positions contribute to
the canonical coefficient record; small numerical
coefficient height does not bound the cardinality of that
record. This is a profile and dictionary comparison. It
does not build primitive additive triples carrying the
prescribed third-letter lock and all three original radical
caps.
