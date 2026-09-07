# Valuation fan for a fourth suspension search

## 0. Scope

Consider the first root-changing ansatz which mixes the old root with the
reciprocal coefficient:

\[
 T=\frac{U_0(S)+P\,U_1(S)}
         {V_0(S)+P\,V_1(S)}.                           \tag{0.1}
\]

At the old reciprocal boundary,

\[
 v(t)=1,\qquad v(P)=1,\qquad v(S)=-1.                  \tag{0.2}
\]

This note classifies the possible leading contacts of (0.1) before any
coefficient elimination.  It is deliberately only an old-boundary theorem.
The finite valuation theorem still requires every prime in the denominator
of \(T\), every prime in \(\partial_ST\), and every exceptional divisor over
their intersections to be put in the complete polar ledger.

The conclusion is useful:

* all leading contacts are indexed by two integers and two three-valued
  support flags;
* away from three exceptional contacts, the old-boundary test is a single
  inequality;
* two exceptional contacts are impossible;
* the only exceptional contact which can cancel is one explicit lattice
  ray, followed by one initial-residue identity; in degree two that residue
  identity forces degeneration to \(T\in k(P)\).

Thus the old boundary produces a small finite input for a subsequent
prime-by-prime or Gröbner calculation.  It does not by itself prove that a
fourth suspension exists or that none exists.

## 1. Effective degrees and the nine support contacts

Use the convention \(\deg 0=-\infty\), and put

\[
 \lambda_U=\max\{\deg U_0,\deg U_1-1\},\qquad
 \lambda_V=\max\{\deg V_0,\deg V_1-1\}.                \tag{1.1}
\]

Then

\[
 v(U_0+PU_1)=-\lambda_U,\qquad
 v(V_0+PV_1)=-\lambda_V.                               \tag{1.2}
\]

There is no generic cancellation in (1.2).  In the tied case the two
residues have different powers of the boundary residue \(z=\operatorname
{in}_v(PS)\).

For each row record whether its two entries attain the maximum:

\[
 \epsilon_U=
 \bigl({\bf1}_{\deg U_0=\lambda_U},
       {\bf1}_{\deg U_1-1=\lambda_U}\bigr),
\quad
 \epsilon_V=
 \bigl({\bf1}_{\deg V_0=\lambda_V},
       {\bf1}_{\deg V_1-1=\lambda_V}\bigr).            \tag{1.3}
\]

Each flag is one of

\[
 (1,0),\qquad(0,1),\qquad(1,1).                        \tag{1.4}
\]

Consequently, for each pair \((\lambda_U,\lambda_V)\), there are only nine
leading support contacts.  If \(\deg U_i,\deg V_i\le2\), then

\[
 \lambda_U,\lambda_V\in\{-1,0,1,2\},                  \tag{1.5}
\]

where value \(-1\) means that only \(P\) times a nonzero constant occurs.
This is already a finite \(4\cdot4\cdot9\) box, before removing zero
denominators, common factors, and nonbirational strata.

## 2. The derivative contact operator

Put

\[
 r=S^{-1},\qquad z=PS,\qquad m=\lambda_U-\lambda_V.
\]

The ansatz has a unique expansion

\[
 T=r^{-m}F(r,z),\qquad F(0,z)=F_0(z)=\frac{A(z)}{B(z)}, \tag{2.1}
\]

where \(A\) and \(B\) are nonzero affine polynomials in \(z\).  Their
constant and linear supports are exactly the flags (1.3).  At fixed \(P\),

\[
 \partial_S=r(-r\partial_r+z\partial_z),
\]

and hence

\[
 \boxed{\partial_ST=
 r^{\,1-m}K(r,z),\qquad
 K=mF-rF_r+zF_z.}                                     \tag{2.2}
\]

Let

\[
 c=\operatorname{ord}_r K(r,z)\in\mathbb Z_{\ge0}
       \cup\{\infty\}.                                 \tag{2.3}
\]

The value \(c=\infty\) says that \(T\in k(P)\) and therefore is not a root
coordinate.  Otherwise

\[
 v(T)=-m,\qquad
 \boxed{v(\partial_ST)=\delta=1-m+c.}                  \tag{2.4}
\]

For a cotangent lift, the \(Q\)-linear term is

\[
 \widehat Q=\frac{L}{\partial_ST}+H(P,T),               \tag{2.5}
\]

where \(L\) is a boundary unit.  The direct determinant-one lift has
\(L=Q\); the controlled-divisor normalization used in the reciprocal chart
has \(L=Q-PS=y\).  Both have value zero at the generic point of \(t=0\).
Thus the first summand has order \(-\delta\).  A positive \(\delta\) is a
pole which must be killed by the initial residue of \(H\).

Formula (2.4), rather than the four raw polynomial degrees, is the contact
vector relevant to polynomiality:

\[
 \boxed{(\lambda_U,\lambda_V;\epsilon_U,\epsilon_V;c).}\tag{2.6}
\]

## 3. Complete classification of leading derivative cancellation

At \(r=0\),

\[
 K_0(z)=mF_0(z)+zF_0'(z).                              \tag{3.1}
\]

Thus \(c>0\) exactly when

\[
 F_0(z)=a z^{-m}.                                      \tag{3.2}
\]

Because \(F_0=A/B\) is a ratio of affine polynomials, after cancelling a
common factor (3.2) has only the following possibilities:

\[
\begin{array}{c|c|c|c}
m&F_0&(\epsilon_U,\epsilon_V)&\text{leading meaning}\\ \hline
0&a&\text{compatible proportional supports}
  &\text{constant leading root},\\
1&a/z&((1,0),(0,1))
  &\text{leading root }a/P,\\
-1&az&((0,1),(1,0))
  &\text{leading root }aP.
\end{array}                                             \tag{3.3}
\]

Every other contact has \(c=0\).  This proves that the support enumeration
does not hide an unbounded collection of derivative cancellations.

The generic rows are:

\[
\begin{array}{c|c|c|c}
m&v(T)&v(\partial_ST)&\text{old-boundary result}\\ \hline
m\le0&-m\ge0&1-m>0&\text{pole; impossible to cancel},\\
m=1&-1&0&\text{no pole},\\
m\ge2&-m\le-2&1-m<0&\text{cotangent term is regular}.
\end{array}                                             \tag{3.4}
\]

The first row is impossible because every monomial \(P^iT^j\) has
nonnegative value when \(m\le0\).

In the exceptional constant contact \(m=0\), one has \(c\ge1\) and hence
\(\delta=1+c\ge2\).  Again \(P\) and \(T\) both have nonnegative value, so
the pole cannot cancel.  The root-only Möbius obstruction is the first
member \(c=1,\delta=2\) of this row.

In the exceptional \(m=-1\) contact, \(\delta=2+c>0\) and \(T\) is regular,
so it is also impossible.

Only the exceptional \(m=1,\ F_0=a/z\) contact remains.  Here

\[
 \delta=c.
\]

If \(c>0\), the possible cancelling monomials are exactly

\[
 v(P^iT^j)=i-j=-c
 \quad\Longleftrightarrow\quad
 \boxed{j-i=c}.                                       \tag{3.5}
\]

This is the unique exceptional cancellation ray.

## 4. The residue equation on a surviving lattice row

Valuation equality is necessary but not sufficient.  Write

\[
 K(r,z)=r^cK_c(z)+O(r^{c+1}),\qquad
 \eta=\operatorname{in}_v(L).
\]

For a general pole \(\delta=1-m+c>0\), the complete first residue equation
is

\[
 \boxed{
 \frac{\eta}{K_c(z)}
 +
 \sum_{\substack{i,j\ge0\\i-mj=-\delta}}
 h_{ij}\,z^iF_0(z)^j=0
 }\quad\text{in }k(E).                                 \tag{4.1}
\]

With a fixed degree bound on \(H\), the sum is finite.  If (4.1) holds, the
finite valuation theorem requires the same calculation in the successive
negative graded pieces down to order \(-1\).  There are exactly \(\delta\)
such slots.

For the exceptional surviving row \(m=1\), (4.1) uses only the ray
\(j-i=c\).  Since \(F_0=a/z\), every term on that ray is a scalar multiple
of \(z^{-c}\).  Therefore the entire first residue test reduces further to

\[
 \boxed{\eta/K_c(z)\in k^\times z^{-c}}
 \quad\Longleftrightarrow\quad
 K_c(z)\in k^\times\eta z^c.                           \tag{4.2}
\]

It is a one-dimensional residue-shape test before any nonlinear coefficient
equations are formed.

For the standard reciprocal source jet,
\[
 q|_{t=0}=h_0y^2,\qquad z=\operatorname{in}_v(PS)=-h_0y,
\]
so \(\eta=y=-z/h_0\).  In that chart (4.1) is an identity in the single
residue parameter \(z\), and (4.2) says simply

\[
 \boxed{K_c(z)\in k^\times z^{c+1}.}                   \tag{4.3}
\]

This can be cleared and compared coefficient by coefficient.

## 5. The bounded degree-two lattice

Under \(\deg U_i,\deg V_i\le2\), the rows which pass the old-boundary
inequality are particularly small.

1. **Generic primitive row:** \(m=1,c=0\), with
   \[
   (\lambda_U,\lambda_V)=(0,-1),(1,0),(2,1).
   \]
   It has no old-boundary cotangent pole.  Known affine root changes live
   here, so reconstruction and orbit filters are essential.

2. **Higher-pole rows:** \(m=2,c=0\), with
   \[
   (\lambda_U,\lambda_V)=(1,-1),(2,0),
   \]
   and \(m=3,c=0\), with
   \[
   (\lambda_U,\lambda_V)=(2,-1).
   \]
   The cotangent term is regular at the old boundary, but \(T\) has pole
   order two or three.  These rows need a primitivity/birationality test.

3. **Exceptional cancellation candidate:** \(m=1,c>0\), necessarily
   \[
   \epsilon_U=(1,0),\qquad\epsilon_V=(0,1),
   \]
   with the same three degree pairs as item 1.  Its pole can only cancel on
   \(j-i=c\), subject to (4.1).

Every \(m\le0\) row, including both exceptional contacts in (3.3), is
eliminated at the old boundary.  There are no further leading cases.

### Proposition 5.1 -- the exceptional degree-two candidate is empty

Assume the standard reciprocal source residue

\[
 \eta\in k^\times z.
\]

Then none of the three exceptional degree pairs in item 3 satisfies (4.3)
unless \(T\in k(P)\).

#### Proof

For \((\lambda_U,\lambda_V)=(0,-1)\), the support conditions force

\[
 T=\frac{a+bP}{dP}\in k(P)
\]

immediately.

For \((1,0)\), write

\[
\begin{aligned}
N&=a_1S+a_0+P(b_1S+b_0),\\
D&=P(c_1S+c_0).
\end{aligned}
\]

The first derivative residue is

\[
 K_1=-\frac{a_0c_1-a_1c_0}{c_1^2z}.
\]

It cannot be proportional to \(z^2\), as (4.3) requires for \(c=1\).
Killing it gives \(a_0c_1=a_1c_0\).  The next residue is then

\[
 K_2=-\frac{b_0c_1-b_1c_0}{c_1^2},
\]

which cannot be proportional to \(z^3\).  Killing it makes both numerator
coefficient vectors proportional to \((c_1,c_0)\), and hence
\(T\in k(P)\).

For \((2,1)\), write

\[
\begin{aligned}
N={}&a_2S^2+a_1S+a_0+P(b_2S^2+b_1S+b_0),\\
D={}&d_0+P(c_2S^2+c_1S+c_0).
\end{aligned}
\]

The first residue is

\[
K_1=-\frac{(a_1c_2-a_2c_1)z-2a_2d_0}{c_2^2z^2}.
\]

It cannot have the required shape \(z^2\).  Its vanishing forces

\[
d_0=0,\qquad a_1c_2=a_2c_1.
\]

After these substitutions,

\[
K_2=
\frac{(-2a_0a_2c_2+2a_2^2c_0)
 +
 (a_1b_2c_2-a_2b_1c_2)z}
{a_2c_2^2z}.
\]

It cannot have shape \(z^3\).  Its vanishing forces

\[
a_2c_0=a_0c_2,\qquad a_1b_2=a_2b_1.
\]

The next residue is

\[
 K_3=\frac{2(a_0b_2-a_2b_0)}{a_2c_2},
\]

which cannot have shape \(z^4\).  If it too vanishes, \(U_0,U_1,V_1\)
are all proportional, \(d_0=0\), and again \(T\in k(P)\).  QED

Therefore the exceptional \(P\)-dependent contact does not repair the
root-only Möbius pole in the degree-two ansatz.  The complete old-boundary
survivor list is:

\[
\boxed{
\begin{array}{c|c}
m&(\lambda_U,\lambda_V)\\ \hline
1&(0,-1),(1,0),(2,1),\quad c=0,\\
2&(1,-1),(2,0),\quad c=0,\\
3&(2,-1),\quad c=0.
\end{array}}                                           \tag{5.1}
\]

## 6. Search order forced by the fan

The proposed experimental order can now be made precise.

### A. One denominator prime and one \(P\)-dependent coefficient

Start with the three \(m=1\) degree pairs.  Proposition 5.1 removes the
exceptional support
\[
\epsilon_U=(1,0),\quad\epsilon_V=(0,1).
\]
The remaining \(c=0\) strata have no old-boundary cotangent pole.  Factor
their denominators and Wronskians, remove common-factor and \(T\in k(P)\)
strata, and test whether any is outside the known affine-root orbit.

### B. Two denominator primes, one primitive root variable

Retain \(m=1\), but factor the complete denominator of (2.5):

\[
(V_0+PV_1)\cdot
\bigl((U_0'+PU_1')(V_0+PV_1)
      -(U_0+PU_1)(V_0'+PV_1')\bigr).                  \tag{6.1}
\]

Compile every irreducible factor and their intersection valuations.  The
old \(t\)-row is already classified; this stage tests the new Wronskian
prime rather than repeating it.

### C. One root plus one primitive reconstruction variable

Only after A and B should one enlarge the residue algebra.  If a primitive
variable \(R\) has value \(\rho\), replace (4.1)'s lattice equation by

\[
 i-mj+k\rho=-\delta,\qquad i,j,k\ge0.                  \tag{6.2}
\]

This is again a finite lattice slice under fixed degree bounds.  It records
exactly how an additional primitive can supply a residue unavailable in
\(k[P,T]\).

### D. Genuinely nonmonomial denominator

Finally retain irreducible factors of (6.1) as actual primes.  A
nonmonomial prime cannot be replaced by its old-boundary exponent vector:
its residue field and the zeros of \(K_0(z)\) carry information lost by that
replacement.  Resolve intersections, add the exceptional valuations, and
apply the successive residue test of the finite valuation theorem.

## 7. What an empty computation would prove

An empty calculation in stages A--D would be strong evidence for the MLE
trichotomy only if it carries the following certificate:

1. the degree/support box (1.5) and all nine flags were enumerated;
2. the three derivative contacts (3.3) were treated separately;
3. every factor of (6.1), not only \(t\), was placed in the polar ledger;
4. all negative residue slots, not only the first tie, were killed;
5. common-factor, denominator, inseparable, and nonbirational strata were
   removed by saturation after the residue calculation.

Without items 3--5, emptiness is only an old-boundary obstruction.  With
them, it is a finite-valuation/Hilbert-basis certificate for the declared
degree-two ansatz.
