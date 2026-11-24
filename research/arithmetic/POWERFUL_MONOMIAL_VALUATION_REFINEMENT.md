# Prime-adic Vandermonde bounds for short monomial congruences

Fix \(A\ge1\). For distinct integers \(e,f\in[r,r+A]\), an
\(r\)-full modulus \(n\), a unit \(K\bmod n\), and
\(n^{1/(3r)}\le Y\le n^{2/r}\), consider
\[
     U^e/V^f\equiv K\pmod n                                           \tag{1}
\]
with \(U,V\) positive reduced rationals of height at most \(Y\),
whose numerators and denominators are units modulo \(n\). For
sufficiently large \(r\) depending on \(A\), and every
\(\varepsilon>0\), the number of pairs is
\[
 O_{A,r,\varepsilon}\!\left(
 Y^{\,O_A((\log r/r)^{1/8})+\varepsilon}\right).                         \tag{2}
\]
For fixed \(u\ge v\ge0\), with \(p=r+u\), \(q=r+v\), the
consequent primitive-triple bound is
\[
 N_{p,q,r}(X)\ll_{u,v,r,\varepsilon}
 X^{\,1/p+1/q-2/r^2+
          O_{u,v}(r^{-17/8}(\log r)^{1/8})+\varepsilon}.                 \tag{3}
\]
Here primitive triples are ordered positive \(p\)-full, \(q\)-full,
and \(r\)-full summands \(a+b=c\le X\) with \(\gcd(a,b,c)=1\).
These refine the errors in the manuscript's short rational
monomial congruence and all-fixed-offset theorems. The new
step measures the Vandermonde coefficient prime by prime.
[Manuscript, Short rational congruences and The doubled leading saving](../../unpublished/powerful-number-representations/main.tex)

## The valuation estimate

We may interchange \(e,f\) and invert \(K\) so that \(e>f\).
Set \(g=\gcd(e,f)\le A\). Every solution of (1) gives a unit
\(W=U^{e/g}/V^{f/g}\) with \(W^g\equiv K\bmod n\). At each odd
prime power the unit power map has fibers of size at most \(g\)
by cyclicity of the unit group. At powers of two the corresponding
bound is \(2g\), using the two cyclic factors. The Chinese
remainder theorem and the divisor bound therefore give
\(O_{A,\delta}(n^\delta)\) possible residues \(W\) for every
\(\delta>0\). Since \(n\le Y^{3r}\), this partition costs
\(O_{A,r,\varepsilon}(Y^{\varepsilon/2})\), by taking
\(\delta=\varepsilon/(6r)\). In each part we have coprime
exponents \(E=e/g>F=f/g\), difference \(d=E-F\le A\), and
\(F\ge r/A\). We continue to use \(K\) for the new residue.

Let \(h\ge1\) with \(h<F\) and \(L_0=h+1\). Use the
\(M=L_0^2\) grid weights \(Fi+Ej\) for \(0\le i,j\le h\).
They are distinct, and subtracting them from their maximum
and sorting gives nonnegative distinct weights
\(w_1<\cdots<w_M\). For \(J=\binom M2\), put
\[
 \kappa=\frac{\prod_{i<j}(w_j-w_i)}{\prod_{i=0}^{M-1}i!}
       =\det\!\left(\binom{w_j}{i-1}\right)_{1\le i,j\le M}.
\]
The second identity is the ordinary Vandermonde identity in a
polynomial basis with leading coefficients \(1/i!\).
In particular \(\kappa\) is a positive integer.

Fix a prime \(\ell\) and \(Q=\ell^\nu\). One of \(E,F\) is a
unit modulo \(Q\). In each of \(L_0\) rows, use the corresponding
variable to order the grid weights. Each row is \(L_0\)
consecutive terms of a progression with unit step.
For residue \(a\bmod Q\) let \(c_a\) count the grid weights
in that class, and let \(b_a\) count \(0,\ldots,M-1\) in that
class. Subtracting the valuation of the Vandermonde on those
consecutive weights gives the exact identity
\[
 v_\ell(\kappa)
   =\sum_{\nu\ge1}\left(\sum_{a\bmod\ell^\nu}\binom{c_a}{2}
                    -\sum_{a\bmod\ell^\nu}\binom{b_a}{2}\right).
                                                                        \tag{4}
\]
For \(Q\le L_0\) we have \(|c_a-M/Q|\le L_0\), by adding
the discrepancy of at most one from each row. Since the totals
of each count vector are \(M\), the summand in (4) is
\(\tfrac12\sum_a[(c_a-M/Q)^2-(b_a-M/Q)^2]\). It is at most
\(Q L_0^2/2\). Summing these powers \(Q\le L_0\) gives
\(O(L_0^3)\). If \(Q>L_0\), each row contains at most one
weight in a fixed class. Then \(\sum_a\binom{c_a}{2}\le M L_0/2\).
There are at most \(O_A(\log(rL_0))\) contributing larger
powers, since for \(Q>\max_j w_j\) the grid counts have zero
pair coincidences, and \(M-1\le\max_j w_j\).
Thus, uniformly in the prime,
\[
                  v_\ell(\kappa)\ll_A L_0^3\log(rL_0).                \tag{5}
\]

## Amplification and the determinant argument

Choose
\[
 h=\left\lfloor c_A (r/\log r)^{1/8}\right\rfloor,\qquad
 \delta_0=6A/h,\qquad
 t=\lceil4/\delta_0\rceil,\qquad B_*=\lceil(4J+2)/\delta_0\rceil .
\]
Here \(c_A>0\) is a sufficiently small fixed constant and the
onset for \(r\) depends on \(A\). These choices have \(h<F\)
and \(0<\delta_0\le1\). Since \(B_*\ll_A h^5\), (5) gives
\[
 B_*\bigl(1+\max_\ell v_\ell(\kappa)\bigr)
       \ll_A h^8\log r\le r/2.
                                                                        \tag{6}
\]
We also have \(B_*\le r/20\) and \(4dt<r\) for large \(r\).

Use the manuscript's multiplicative amplification in each
fixed residue branch. For a solution family \(\mathcal S\),
its \(t\)-fold coordinatewise products have height at most
\(Z=Y^t\), satisfy (1) with \(K^t\), and number
\(\gg_{t,\eta}|\mathcal S|^tY^{-t\eta}\) for every \(\eta>0\).
Indeed independently color each prime in each coordinate with
two signs, keeping a pair when the two numerators use positive
signs and the two denominators use negative signs. The number
of distinct primes across the four integers is
\(O(\log Y/\log\log(3Y))\). Some coloring therefore keeps
\(|\mathcal S|Y^{-o(1)}\) pairs and their products have no
cross-cancellation. A reduced product has four numerator and
denominator integers at most \(Z\), whose factorizations into
\(t\) integers number \(Y^{o(1)}\) by the fixed-order divisor
bound. This gives the stated amplification. We work with the
amplified family and rename its residue \(K\).

Write \(n=\prod_\ell \ell^{a_\ell}\) with \(a_\ell\ge r\).
By (6) the cap \(\lfloor a_\ell/B_*\rfloor\) exceeds or
equals \(v_\ell(\kappa)+1\) for each \(\ell\mid n\). The
divisor using all capped exponents is at least
\(n^{1/(2B_*)}\), hence at least \(Z^{\delta_0}\), because
\(Z^{\delta_0}\le n^{2t\delta_0/r}\le n^{10/r}\).
Build a divisor \(L\) using an initial exponent
\(v_\ell(\kappa)+1\) when activating a prime, then unit
exponent steps up to its cap. Stop when \(L\ge Z^{\delta_0}\).
Each prime \(\ell\mid n\) is at most \(n^{1/r}\le
Z^{\delta_0}\). Every jump has factor at most
\(\kappa Z^{\delta_0}\), giving
\[
 Z^{\delta_0}\le L\le\kappa Z^{2\delta_0},\qquad
 L^{B_*}\mid n_0,\qquad n_0\ge Z^{4J+2},                               \tag{7}
\]
where \(n_0\) is the product of full powers in \(n\) of the
primes activated in \(L\).

For \(T=U/V\) the congruence gives
\[
        U^d\equiv K T^{-F},\qquad V^d\equiv K T^{-E}\pmod n.
\]
Any two amplified pairs with the same rational \(T\) would
have congruent \(d\)-th powers of \(U\) and \(V\). Their
cross-numerators are bounded by \(Z^{2d}<n\) since \(4dt<r\).
The congruences are equalities, and positivity yields the
same pair. Thus the values \(T_i\) for distinct pairs are
distinct.

Suppose \(M\) distinct pairs have \(T_i\) in one unit class
modulo \(L\). The alternant factors as
\[
 \det(T_i^{w_j})_{i,j=1}^M
   =\left(\prod_{i<j}(T_j-T_i)\right)\Phi(T_1,\ldots,T_M)
\]
with \(\Phi\in\mathbb Z[T_1,\ldots,T_M]\). Taylor expansion
on the diagonal gives
\(\Phi(z,\ldots,z)=\kappa z^{\sum_j w_j-J}\). Since every
activated \(\ell\) occurs in \(L\) with exponent at least
\(v_\ell(\kappa)+1\), and all \(T_i\) are units, the value
\(\Phi(T_1,\ldots,T_M)\) has valuation exactly
\(v_\ell(\kappa)\) at those primes. A numerator of
\(T_j-T_i\) is bounded by \(Z^4\), since each \(T_i\) has
height at most \(Z^2\). If \(n_0\) divided the alternant, it
would divide the nonzero product of those numerators times
\(\kappa\), whose absolute value is at most \(\kappa Z^{4J}\).
This is smaller than \(n_0\) for large \(Z\), by (7).

The matrix of monomials \(U^{di}V^{dj}\), with denominators
cleared row by row, gives an integer determinant congruent
modulo \(n_0\) to the alternant times a unit. The column
multipliers are \(K^{i+j}\), and the row multipliers include
the fixed power of \(T_i\) that turns the negative grid
weights into the \(w_j\). The integer determinant is therefore
nonzero. It is divisible by \(L^J\), since the alternant is,
and \(L^J\mid n_0\). Clearing denominators to exponent \(dh\)
in each rational variable bounds the entries by \(Z^{2dh}\).
The resulting upper bound \(M! Z^{2dhM}\) is smaller than
\(L^J\) for large \(Z\), as
\[
                 \delta_0 J/M=3A(h+2)>2dh .
\]
This contradiction bounds the number of pairs in each class
by \(M-1\). There are \(O_{A,r}(Z^{2\delta_0})\) amplified
pairs by (7). Undoing amplification gives at most
\(O_{A,r,\eta}(Y^{2\delta_0+\eta})\) original pairs in the
fixed branch. The unit-root partition and choice
\(\eta=\varepsilon/2\) give (2). Bounded \(Y\), including
\(n=1\), is handled by the finite number of rational choices
below the proof's height cutoff. The case \(e=f\) is outside
(2) and is covered separately by the equal-exponent graph
where needed.

## Transfer to fixed-offset primitive triples

For \(u=v\), the manuscript's aligned gcd case already gives
the stronger error \(O_u(r^{-5/2})\). Suppose \(u>v\).
For a canonical dyadic box write \(S_A,S_B,S_C\) for the
sums of the logarithmic coordinate sizes on the three sides,
to base \(X\). The manuscript reduces to competitive boxes
where the available counting exponent \(\Delta\) is at least
\(2/r-(u+v+2)/r^2\), and hence
\(S_T=1/r+O_{u,v}(r^{-2})\) on every side. Normalize the
coordinate masses to probability vectors \(A,B,C\). For
offsets \(u,v,0\), their costs are
\(m_A=u+\sum_j jA_j\), \(m_B=v+\sum_j jB_j\), and
\(m_C=\sum_j jC_j\). All three costs are bounded, and with
\(D=\max(m_A+m_B,m_A+m_C,m_B+m_C)\) the manuscript's pair
bounds give
\[
 \Delta\le 2/r-D/r^2+O_{u,v}(r^{-3}).                                  \tag{8}
\]

Fixing all coordinates except a maximum mass from each of
the three sides leaves a primitive ternary equation with
nonzero coefficients and powers at least \(r\), in a box of
side \(O(X^{1/r})\). By Browning–Verzobio, Theorem 1.2, its
count is \(O_{u,v,r,\varepsilon}(X^{O(r^{-3/2})+\varepsilon})\),
uniform in the coefficients. Their distinguished-variable
parameter satisfies \(W\ll_r X^{r^{-3/2}}\) here, and the
remaining factors contribute \(O(r^{-3/2})\) in the exponent.
Writing \(\alpha,\beta,\gamma\) for the maxima of \(A,B,C\)
gives the constraint
\[
                     \alpha+\beta+\gamma\le1+O_{u,v}(r^{-1/2}).
                                                                        \tag{9}
\]
[Browning–Verzobio, arXiv:2608.24512v1, Theorem 1.2](https://arxiv.org/pdf/2608.24512v1)

Fix the \(r\)-full target \(n=c\), with
\(n\ge X^{rS_C}\ge X^{1-O_{u,v}(1/r)}\).
Omit coordinate \(i=0\) or \(1\) from both other sides,
of exponents \(e=p+i\) and \(f=q+i\).
For two edges sharing a right vertex, the manuscript's
reconstruction uses reduced rational quotients of the
omitted coordinates on each side and gets (1) with a unit
\(K\) determined by the distinct left vertices. All factors
are units by pairwise coprimality of the primitive summands.
The rational heights are \(O(X^{1/r})\), so
\(Y=n^{2/r}\) satisfies the required height inequalities
and includes them for large \(X\). The exact subtraction
recovers the canceled scalar on the left and the right
vertex from the two fractions. The recovery coefficient
vanishes only if the complete left summands coincide,
which by the unique canonical factorization would identify
the left vertices. Applying (2) gives codegrees
\(X^{O_{u,v}(\theta_r/r)+\varepsilon}\) with
\(\theta_r=(\log r/r)^{1/8}\).
[Exact reconstruction in the manuscript](../../unpublished/powerful-number-representations/main.tex)

Edge multiplicity has bound \(O(X^{1/r^2+\varepsilon})\).
At a fixed pair of vertices the plane curve has equation
\(P a^e+Q c^f=n\). Viewed as a polynomial in \(c\) over
\(\mathbb C[a]\), it is Eisenstein at a simple zero of
\(n-P a^e\), since \(PQn\ne0\). Its degree is
\(\max(e,f)\ge r\), so the coefficient-uniform curve count
in the enclosing square of side \(O(X^{1/r})\) gives the bound.
[Bombieri–Pila, Theorem 5](https://people.maths.ox.ac.uk/pila/Ovals.pdf)

Write \(a=A_i\), \(b=B_i\) for the two omitted masses.
After summing the fixed-side choices and including edge
multiplicity, the two orientations of the degree count give
\[
 r\Delta\le 1+\max\left\{1-a,1-b,
       \tfrac32-\max(a+b/2,b+a/2)\right\}
       +O_{u,v}(\theta_r).                                            \tag{10}
\]
In detail the vertex counts are \(X^{S_A(1-a)}\),
\(X^{S_B(1-b)}\). For bounds \(m,n_1\) on these counts and
\(K_1\) on codegrees, Cauchy–Schwarz bounds the simple edges by
\(m+n_1+\sqrt{K_1}\min(m\sqrt{n_1},n_1\sqrt m)\).
The factor \(X^{S_C}\), the height and multiplicity estimates,
and \(rS_T=1+O_{u,v}(r^{-1})\) give (10).
Thus either \(a\) or \(b\) is \(O_{u,v}(\theta_r)\), or
both \(2a+b\) and \(a+2b\) are at most \(1+O_{u,v}(\theta_r)\).

Call a mass small at the threshold in this alternative. If
\(A_0\) is small, truncating the two costs gives
\(m_A+m_B\ge u+v+3-2A_0-A_1-B_0\); (9) gives the
desired \(u+v+2-O_{u,v}(\theta_r)\). The \(B_0\) case is
symmetric. Otherwise the two inequalities are available at
\(i=0\). If \(A_1\) is small, truncating at indices three
and two gives
\[
 m_A+m_B\ge u+v+5-3A_0-2A_1-A_2-2B_0-B_1.
\]
Subtracting \(2A_0+B_0\le1+O_{u,v}(\theta_r)\) from the
deficit leaves \(A_0+A_2+B_0+B_1+O_{u,v}(\theta_r)\le
2+O_{u,v}(\theta_r)\). A small \(B_1\) uses the other
orientation. In the remaining case, both indices satisfy
\(A_i+B_i\le2/3+O_{u,v}(\theta_r)\), by adding their two
inequalities. Truncating both sides at index two gives
\(m_A+m_B\ge u+v+4-2(A_0+B_0)-(A_1+B_1)\).
All cases therefore give \(D\ge u+v+2-O_{u,v}(\theta_r)\).
Insert this in (8), expand \(1/p+1/q\) to order \(r^{-2}\),
and absorb the finitely many logarithmic box factors into
\(X^\varepsilon\) at each fixed \(r\). This gives (3).

The improving parameter is the minimum exponent of a sufficiently
powerful modulus. The balanced square-full kernel-pair problem and
the low-power cube-modulus collision bounds have separate modulus
and exponent conditions.
[Square-full kernel targets and the conic bound](TRIFONOV_2026_TERNARY_CONIC_SOURCE_SCOPE.md)

## Limit of the current determinant degree budget

The remaining target in
\(\texttt{GAP-ARITH-POWERFUL-FIXED-OFFSET-SHARP-ERROR}\) needs a
monomial-counting loss of order \(r^{-1/2}\) in this route. The
current square grid and greedy-divisor conditions have a stricter
limit even if every prime valuation of \(\kappa\) could be taken to
be zero. The size comparison requires
\(\delta_0 J>2dhM\), hence \(\delta_0\gg 1/h\) for \(d\ge1\).
The nonvanishing check in (7) uses \(B_*\delta_0\ge4J+2\), while
the capped divisor construction for exponents as small as \(r\)
requires \(B_*\le r\). Thus
\(\delta_0\gg h^4/r\) as well. The resulting lower bound
\(\max(1/h,h^4/r)\gg r^{-1/5}\) concerns these sufficient
conditions and this grid. A move towards \(r^{-1/2}\) needs to
change the determinant dimension or the nonvanishing argument
that uses the product of all pair differences.
