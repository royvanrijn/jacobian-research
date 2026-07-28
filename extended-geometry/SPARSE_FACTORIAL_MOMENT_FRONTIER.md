# Sparse factorial-moment frontier

## 1. Statement and scope

Let

\[
 \mathcal L(U^\alpha)=\prod_i\alpha_i!
 \qquad
 \left(\mathcal L:\mathbb C[U_1,\ldots,U_n]\to\mathbb C\right). \tag{1.1}
\]

The exact searches in this note prove three finite results.

### Theorem 1.1 (three-term binary search)

Let

\[
 f=c_0x^{a_0}y^{b_0}+c_1x^{a_1}y^{b_1}
   +c_2x^{a_2}y^{b_2},                                       \tag{1.2}
\]

where the three monomials are distinct,
\(a_i+b_i\leq6\), and \(c_0c_1c_2\ne0\).  Then

\[
 \bigl(\mathcal L(f),\mathcal L(f^2),\mathcal L(f^3)\bigr)
 \ne(0,0,0).                                                   \tag{1.3}
\]

The search exhausts all

\[
 \binom{\binom{6+2}{2}}3=\binom{28}{3}=3276                 \tag{1.4}
\]

supports over \(\mathbb C\).

### Theorem 1.2 (paired four-term search)

Let

\[
 \sigma(U_0,U_1,U_2,U_3)=(U_2,U_3,U_0,U_1)
\]

and let \(M,N\) represent two distinct nontrivial monomial orbits under
\(\sigma\), with \(\deg M,\deg N\leq6\).  For \(c\ne0\), put

\[
 f=(M-\sigma M)+c(N-\sigma N).                                \tag{1.5}
\]

Then

\[
 \bigl(\mathcal L(f),\mathcal L(f^2),
        \mathcal L(f^3),\mathcal L(f^4)\bigr)\ne(0,0,0,0).     \tag{1.6}
\]

There are 100 nontrivial monomial orbits in the declared degree range, and
the search exhausts all

\[
 \binom{100}{2}=4950                                          \tag{1.7}
\]

two-orbit supports.  The involution makes every odd moment vanish
identically, so the substantive certificate is the absence of a common
nonzero root of the second and fourth moments.

### Theorem 1.3 (sharp binary homogeneous cutoffs through degree four)

Let \(1\leq d\leq4\), and let

\[
 f(x,y)=\sum_{j=0}^d c_jx^{d-j}y^j                            \tag{1.8}
\]

be a homogeneous binary form.  If

\[
 \mathcal L(f^m)=0\qquad(1\leq m\leq d+1),                    \tag{1.9}
\]

then \(f=0\).  The cutoff is sharp: for every \(d=1,2,3,4\),
there is a nonzero form for which the first \(d\) moments vanish.

Theorem 1.3 is an exact fixed-dimension, fixed-degree result.  Theorems
1.1--1.2 are finite support theorems.  None proves the Factorial Conjecture
or the Strong Factorial Conjecture outside the displayed ranges.

The canonical checker is
[`verify_sparse_factorial_moment_frontier.py`](../scripts/verify_sparse_factorial_moment_frontier.py).

## 2. Relation to the existing factorial literature

Van den Essen--Wright--Zhao introduced the Factorial Conjecture and proved
its one-variable case, the linear-form case, and the initial-moment result
for at most two monomials in
[*On the Image Conjecture*](https://arxiv.org/abs/1008.3962).
Edo--van den Essen formulated the Strong Factorial Conjecture in
[*The Strong Factorial Conjecture*](https://arxiv.org/abs/1304.3956).
Rocks subsequently proved the strong statement for powers of linear forms
and sums of prime powers, gave a construction principle, and obtained
partial two-monomial results in
[*On the Incompatibility of Diophantine Equations Arising from the Strong
Factorial Conjecture*](https://doi.org/10.1016/j.jalgebra.2017.01.039).

Theorem 1.1 is therefore aimed at the next sparse coefficient count rather
than repeating the established linear-form calculations.  It is only a
bounded exponent search; no literature-wide minimality claim is made.

## 3. Exact three-term elimination

Fix a support \(A_i=(a_i,b_i)\), normalize \(c_0=1\), and write

\[
 f=x^{A_0}+u x^{A_1}+v x^{A_2}.                               \tag{3.1}
\]

Here \(x^{(a,b)}=x^ay^b\).  Put \(L_i=a_i!b_i!\).  The first moment gives

\[
 v=-\frac{L_0+L_1u}{L_2}.                                    \tag{3.2}
\]

For every \(m\), direct multinomial expansion gives

\[
 \mathcal L(f^m)=
 \sum_{k_0+k_1+k_2=m}
 \binom{m}{k_0,k_1,k_2}
 u^{k_1}v^{k_2}
 (k_0A_0+k_1A_1+k_2A_2)!.                                   \tag{3.3}
\]

Substitution of (3.2) turns the second and third moments into exact
univariate polynomials

\[
 P_2(u),P_3(u)\in\mathbb Q[u].                                \tag{3.4}
\]

The checker computes their gcd, takes its square-free part, and removes the
factors \(u\) and \(L_0+L_1u\), which are precisely the forbidden loci
\(c_1=0\) and \(c_2=0\).  Every one of the 3276 resulting saturated gcds is
one.  This proves Theorem 1.1 over the algebraic closure of \(\mathbb Q\),
and hence over \(\mathbb C\).

If a support had survived, each irreducible gcd factor would have been
tested against moments through order twelve before being recorded.

## 4. Involution compression of the Dvorsky shadow

The first diagonal shadow in
[`FACTORIAL_MOMENT_WITNESSES.md`](FACTORIAL_MOMENT_WITNESSES.md) is

\[
 U_t(U_aU_d-U_bU_c).                                         \tag{4.1}
\]

After dropping the disjoint factor \(U_t\), its two monomials form an orbit
difference for the involution exchanging the two variable pairs.  Formula
(1.5) is the smallest extension by one further orbit difference.

Since \(\mathcal L\) is invariant under coordinate permutations and
\(\sigma(f)=-f\),

\[
 \mathcal L(f^m)=(-1)^m\mathcal L(f^m),                       \tag{4.2}
\]

so every odd moment is zero.  For a fixed orbit pair, multinomial expansion
produces

\[
 Q_2(c),Q_4(c)\in\mathbb Z[c].                                \tag{4.3}
\]

The checker removes the forbidden factor \(c\) from the square-free gcd of
\(Q_2,Q_4\).  All 4950 saturated gcds are one.  Thus no four-term polynomial
in this Dvorsky-aligned family has four initial zero moments.

Any survivor would have violated the Strong Factorial Conjecture because
the polynomial has exactly four nonzero monomials.  The negative result
therefore closes the first symmetry-compressed attack, but not arbitrary
four-term supports.

## 5. Binary homogeneous forms as interval moments

Let \(X,Y\) be independent mean-one exponential variables.  Put

\[
 S=X+Y,\qquad T=\frac{X}{X+Y}.
\]

Then \(S\) and \(T\) are independent, \(S\) has the Gamma\((2,1)\)
distribution, and \(T\) is uniform on \([0,1]\).  If \(f\) is homogeneous
of degree \(d\), then

\[
\boxed{
 \mathcal L(f^m)
 =(dm+1)!\int_0^1 f(t,1-t)^m\,dt.}                            \tag{5.1}
\]

Consequently Theorem 1.3 is equivalently a finite complex polynomial-moment
theorem on the unit interval for degrees at most four.

For each \(d\), the checker forms the moment polynomials

\[
 M_m(c_0,\ldots,c_d)=\mathcal L(f^m).                          \tag{5.2}
\]

It covers projective coefficient space by the \(d+1\) charts \(c_j=1\).
On every chart, the exact rational Gröbner basis of

\[
 (M_1,\ldots,M_{d+1})                                        \tag{5.3}
\]

is \([1]\).  Hence the only affine common zero of the first \(d+1\)
moments is the origin.

Sharpness is also exact.  On the chart \(c_0=1\), the ideal

\[
 (M_1,\ldots,M_d)                                             \tag{5.4}
\]

is proper and zero-dimensional.  Its lexicographic terminal eliminant has
degree

\[
\begin{array}{c|rrrr}
d&1&2&3&4\\ \hline
\deg E_d&1&2&5&24.
\end{array}                                                    \tag{5.5}
\]

Thus (5.4) has a complex point for every displayed \(d\), while adding
\(M_{d+1}\) removes every point.

Two sharp witnesses are especially small.  In degree two,

\[
\begin{aligned}
 f_2={}&x^2+
 \left(-\frac52+\frac{i\sqrt{15}}2\right)xy
 +\left(\frac14-\frac{i\sqrt{15}}4\right)y^2,\\
 \mathcal L(f_2)&=\mathcal L(f_2^2)=0,\\
 \mathcal L(f_2^3)&=180+108i\sqrt{15}\ne0.                    \tag{5.6}
\end{aligned}
\]

In degree three,

\[
\begin{aligned}
 f_3={}&x^3-(6+i\sqrt{21})x^2y
 +(6+i\sqrt{21})xy^2-y^3,\\
 \mathcal L(f_3^m)&=0\qquad(m=1,2,3),\\
 \mathcal L(f_3^4)
 &=252\,564\,480-42\,301\,440i\sqrt{21}\ne0.                  \tag{5.7}
\end{aligned}
\]

The degree-three example is anti-invariant under \(x\leftrightarrow y\),
which explains its odd moment vanishings.

## 6. Reproduction

Run

```bash
.venv/bin/python scripts/verify_sparse_factorial_moment_frontier.py
```

The command:

- exhausts the 3276 three-term supports using exact rational gcds;
- exhausts the 4950 paired four-term supports using exact integer gcds;
- directly checks odd-moment cancellation on every paired support;
- computes every projective Gröbner chart for binary homogeneous degrees
  one through four;
- verifies the two explicit sharp witnesses; and
- writes
  [`sparse_factorial_moment_frontier.json`](../artifacts/generated-results/sparse_factorial_moment_frontier.json).

The attempted degree-five projective Gröbner calculation is not part of the
certificate: its first chart did not finish inside the short exploratory
window and was terminated.

## 7. Next exact frontier

The computations point to a more focused program.

1. **Arbitrary four-term binary supports.**  The paired family is closed, but
   a general four-term polynomial is the next possible sparse Strong
   Factorial counterexample.
2. **Degree-five binary homogeneous forms.**  The exact pattern suggests the
   bounded conjecture that the first \(d+1\) moments suffice for every
   homogeneous binary degree \(d\), and that \(d+1\) is sharp.  Degree five
   needs modular Gröbner discovery or a recurrence/residue proof rather than
   the direct rational calculation used here.
3. **Interval-moment structure.**  Formula (5.1) replaces multivariate
   factorial bookkeeping by the period
   \[
   \int_0^1\frac{dt}{1-zf(t,1-t)}.
   \]
   Picard--Fuchs recurrences or inverse-branch residues are the natural
   route to an all-degree theorem.

Cyclotomic tensor factors should not be the next priority: they lengthen a
finite zero prefix only by increasing the monomial count at the same rate.
