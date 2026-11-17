# A minimal even Long jet and its fixed-lower obstruction

The least common spatial degree for an even graded full Long jet
satisfying all six commutators through \(h^5\) is **27**. An explicit
rational polynomial point attains it.

Fixing its retained second coefficient gives a one-parameter family
of fourth coefficients. None of that entire fixed-lower family has
a degree-27 sixth coefficient satisfying the \(h^7\) equations.
This last obstruction is for the specified second coefficient;
it is not a classification of all degree-27 first corrections.

## The sharp positive stage

Use source order \((X,Q,P,Z)\), brackets
\(\{P,X\}=\{Z,Q\}=1\), weights \((-1,1,4,2)\), and \(h\)-weight 3.
Let \(f=(R,T,D,S)\) be the full canonical Long polynomial symbols from
[the exact transcription](../dixmier-audit/verify_long_pbw.py).
The image weights are \(w=(-1,2,4,1)\). All products displayed below
are commuting Weyl-symbol products, followed by Weyl quantization.

Let \(g_0\) be the retained 59-coefficient first particular in
[the cap-23 witness](graded_long_cap23_witness.json),
under the key first_particular. It satisfies the entire first equation
\[
 \mathcal A_2g_0=-\Pi^3(f,f)/24,\qquad
 (\mathcal A_mU)_{ij}=\{U_i,f_j\}+\{f_i,U_j\},\quad i<j.       \tag{1}
\]
Pair notation on the right means component \(i,j\).
Every coefficient of \(g_0\) has degree at most 23 and weight
\(w_i-6\).

**Positive theorem and minimum.** There is an actual finite rational
polynomial tuple
\[
                         F=f+h^2g_0+h^4K_0                  \tag{2}
\]
with common spatial cap 27 and the stated grading, satisfying all six
canonical commutators modulo \(h^7\). No such even graded jet has
common cap at most 26.

For the lower bound, use the separate
[complete cap-26 obstruction](GRADED_LONG_CAP26_EVEN_JET_OBSTRUCTION.md).
It allows every first and fourth coefficient at that cap, including
all first-image deformations. Thus the lower bound is not restricted
to \(g_0\). The positive point (2) attains the bound over \(\mathbb Q\)
and hence after every characteristic-zero scalar extension.

The [witness](minimal_even_long_jet_witness.json) retains every
coefficient of \(K_0\). Its component term counts and degrees are
\[
 \begin{array}{c|rrrr}
                 &R&T&D&S\\ \hline
 \hbox{terms}     &0&12&59&3\\
 \hbox{degree}    &-&19&27&15 .
 \end{array}                                                   \tag{3}
\]
A dash means zero. In particular the first image is \(R_0\) through
this stage. Two compact formulas for the point are
\[
\begin{split}
 B&=3X^2P+2(1-3XQ)Z,\qquad
 A_0=Q+3XQ^2-XB/3,\\
 (K_0)_T&=\frac{49}{64}X^8R^2(1+XA_0),\\
 (K_0)_S&=-\frac{49}{16}X^9R^2,\qquad (K_0)_R=0.
\end{split}                                                    \tag{4}
\]
The entire 59-term polynomial \((K_0)_D\), of momentum order three,
is supplied in the witness. All four components are global
polynomials, with no position denominators.

The point is checked by substituting it into every equation
\[
\begin{split}
 \mathcal A_4K_0={}&-\mathcal P_1(g_0,g_0)
 -\frac{\mathcal P_3(f,g_0)+\mathcal P_3(g_0,f)}{24}
 -\frac{\mathcal P_5(f,f)}{1920},                              \tag{5}\\
 \mathcal P_k(U,V)_{ij}&=\Pi^k(U_i,V_j).
\end{split}
\]
Together with (1) and the classical frame equations, this proves
the six commutators of (2) through \(h^5\). This is an actual
polynomial point, rather than a rank-only or formal existence result.

## Entire fourth family with the lower coefficient fixed

At cap \(C\), a coefficient at \(h^m\) in component \(i\) has the
ENTIRE monomial basis
\[
 X^aQ^bP^cZ^d,\qquad
 a=b+4c+2d+3m-w_i\ge0,\quad a+b+c+d\le C.                    \tag{6}
\]
For \(m=4,C=27\), its counts are
\((40,62,80,54)\), totaling 236. The matrix of (5) has 1687 support
rows. The witness supplies a 235-square coefficient minor with
determinant
\[
                         326829\pmod {1000003}.                \tag{7}
\]
All denominators are invertible at this prime. There is the exact
nonzero admitted homogeneous direction
\[
       K_{R^9}=\{R^9,f\}=(0,0,-9R^8,0).                       \tag{8}
\]
Hence the rational rank is exactly 235 and the kernel has dimension
one. The ENTIRE cap-27 solution family for the fixed \(g_0\) is
\[
                              K_0+sK_{R^9}.                   \tag{9}
\]
No source or target gauge removes \(s\) in the proof below.

There is also a direct fixed-\(g_0\) check of the preceding degree
threshold. At cap 26, the fourth space has 206 columns and the
nonzero kernel (8), so its coefficient rank is at most 205.
A retained augmented 206-square minor, including the fixed right
side of (5), has determinant
\[
                         964926\pmod {1000003}.                \tag{10}
\]
Thus cap 26 is inconsistent for this fixed lower coefficient.
This separate finite check proves the fixed-\(g_0\) minimum directly;
the universal lower bound in the theorem uses the stronger cap-26
source cited above.

## The next equation fails for every member of the fixed family

Permit EVERY sixth coefficient \(L\) at cap 27, of weight \(w_i-18\).
Its component counts are
\[
                         (13,24,34,20),\qquad 91.              \tag{11}
\]
In particular all thirteen first-image corrections are allowed.

For \(K=K_0+sK_{R^9}\), the coefficient of \(h^7\) requires
\[
                         \mathcal A_6L=\mathcal B_7(K),
\]
where
\[
\begin{split}
 \mathcal B_7(K)={}&-\mathcal P_1(g_0,K)-\mathcal P_1(K,g_0)\\
 &-\frac{\mathcal P_3(f,K)+\mathcal P_3(K,f)
                           +\mathcal P_3(g_0,g_0)}{24}\\
 &-\frac{\mathcal P_5(f,g_0)+\mathcal P_5(g_0,f)}{1920}
  -\frac{\mathcal P_7(f,f)}{322560}.                           \tag{12}
\end{split}
\]
This expression includes all contributions at the indicated order.
It is affine in \(s\), since a product of two fourth coefficients
first contributes at \(h^9\).

Exact polynomial calculation gives
\[
                   \mathcal B_7(K_0+sK_{R^9})
                         =\mathcal B_7(K_0).                  \tag{13}
\]
Thus the entire free parameter is retained and its coefficient is
verified to vanish; no parameter sampling or normalization premise
is used. The constant right side has 114 nonzero terms.

The complete row set has 974 entries. A 91-square coefficient
minor has determinant \(395212\) modulo \(1000003\), so
\(\mathcal A_6\) is injective. More decisively, the witness retains
a rational functional \(\lambda\) supported on 72 specified rows:
\[
                   \lambda\mathcal A_6=0,\qquad
                   \lambda\mathcal B_7(K_0)=1.                \tag{14}
\]
The replay verifies the first equality on all 91 coefficient
columns and the second exactly over \(\mathbb Q\). Equations
(12)--(14) would imply \(0=1\), for every \(s\) and every admitted
\(L\), after any characteristic-zero scalar extension.

**Fixed-lower obstruction.** No cap-27 sixth coefficient completes
any member of (9) through commutator order \(h^7\). Higher
coefficients cannot change that equation. This conclusion keeps
\(g_0\) fixed and does not exclude other degree-27 first families.

## Replay and interpretation

[The checker](check_minimal_even_long_jet.py) uses Python integers
and Fraction, [the complete column construction](check_graded_long_cap23.py),
and [the independent polynomial arithmetic](check_long_graph_curvature.py).
The canonical producer supplies polynomial input only.
[The receipt](minimal_even_long_jet_controls.json) pins the
positive point, full fourth kernel, direct fixed-cap lower check,
all sixth-coefficient columns, and the 72-term obstruction
functional. It separately pins the cap-26 theorem needed for the
universal minimum; that theorem is imported rather than rerun by
this checker.

Evenness here refers to Weyl symbols. Normal PBW reordering preserves
the spatial cap but can mix powers of \(h\). The positive finite
jet is a compatible partial tuple, not an asserted full lift.
The fixed-lower \(h^7\) obstruction supplies no universal cap-27
higher-order conclusion, no normalization theorem and no DC1/DC2
resolution.
