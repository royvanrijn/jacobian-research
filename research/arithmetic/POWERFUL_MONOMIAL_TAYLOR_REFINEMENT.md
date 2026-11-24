# Taylor stability of the Vandermonde factor over powerful moduli

Fix \(A\ge1\), and let \(e\ne f\) lie in \([r,r+A]\).
For an \(r\)-full modulus \(n\), a unit \(K\bmod n\), and
\(n^{1/(3r)}\le Y\le n^{2/r}\), let \(\mathcal S\) consist of
positive reduced rational pairs \(U,V\) of height at most \(Y\)
with numerators and denominators coprime to \(n\) and
\(U^e/V^f\equiv K\bmod n\). For sufficiently large \(r\)
depending on \(A\), and every \(\varepsilon>0\), we prove
\[
 |\mathcal S|\ll_{A,r,\varepsilon}
     Y^{\,O_A((\log r/r)^{1/5})+\varepsilon}.                         \tag{1}
\]
If the full-power part of \(n\) supported on primes
\(\ell\ge r\) is at least \(n^{1/2}\), the exponent in (1)
improves to \(O_A(r^{-1/5})+\varepsilon\). This gives the
unrestricted primitive-triple bound for fixed
\(p=r+u,\ q=r+v,\ u\ge v\ge0\)
\[
 N_{p,q,r}(X)\ll_{u,v,r,\varepsilon}
 X^{\,1/p+1/q-2/r^2+O_{u,v}(r^{-11/5})+\varepsilon}.                  \tag{2}
\]
Here \(N_{p,q,r}(X)\) counts ordered positive primitive triples
\(a+b=c\le X\) with the indicated fullness orders.

The change is in choosing a divisor on which the Vandermonde
quotient has stable prime-adic valuation. The preceding
valuation refinement bounds the valuation of the quotient
on the diagonal. The Taylor expansion below controls the
change relative to that diagonal value.
[Baseline determinant method and fixed-offset reconstruction](../../unpublished/powerful-number-representations/main.tex)
· [Preceding valuation refinement](POWERFUL_MONOMIAL_VALUATION_REFINEMENT.md)

## The Taylor lemma

Let \(0\le w_1<\cdots<w_M=H\) with \(M\ge2\). Put
\(J=\binom M2\), \(S=\sum_jw_j\), and use the polynomial
\[
 \det(x_i^{w_j})_{i,j=1}^M=
      \prod_{i<j}(x_j-x_i)\,\Phi(x_1,\ldots,x_M),\qquad
 \kappa=\det\left(\binom{w_j}{i-1}\right)_{i,j=1}^M
       =\frac{\prod_{i<j}(w_j-w_i)}{\prod_{i=0}^{M-1}i!}>0 .
\]
The alternant is divisible by the Vandermonde over
\(\mathbb Z[x_1,\ldots,x_M]\), so \(\Phi\) is integral.
The value of the binomial determinant follows from the
leading coefficients of the binomial polynomials.

**Lemma.** Fix a prime \(\ell\), a unit \(z\in\mathbb Z_\ell\),
and \(s\ge 1+\lfloor\log_\ell H\rfloor\). If
\(x_i\equiv z\bmod \ell^s\) for every \(i\), then
\[
 \frac{\Phi(x_1,\ldots,x_M)}{\kappa z^{S-J}}
                    \ \in\ 1+\ell\mathbb Z_\ell.                     \tag{3}
\]

Write \(x_i=z(1+\ell^s t_i)\) with \(t_i\in\mathbb Z_\ell\).
Expand \((1+\ell^s t_i)^{w_j}\) and apply Cauchy–Binet.
After extracting \(z^{S-J}\) and dividing by the
Vandermonde, this gives a finite sum for \(\Phi\):
\[
 z^{-(S-J)}\Phi(x_1,\ldots,x_M)
 =\sum_{0\le k_1<\cdots<k_M\le H}
       \ell^{s(\sum k_i-J)}
       \frac{\det(t_i^{k_j})}{\prod_{i<j}(t_j-t_i)}
       \det\left(\binom{w_j}{k_i}\right).                             \tag{4}
\]
The quotient in this formula is an integral polynomial, so
the formula holds even when some \(t_i\) coincide.
The term for \(k_i=i-1\) is \(\kappa\).

For another index set put \(a=\sum k_i-J\ge1\).
Since \(\det((w_j)_{k_i})\), with falling factorial entries,
is alternating in the \(w_j\), it is an integer multiple of
their Vandermonde. Therefore
\[
 \frac{\det\!\left(\binom{w_j}{k_i}\right)}{\kappa}
 \in \frac{\prod_{i=1}^M(i-1)!}{\prod_{i=1}^M k_i!}\,\mathbb Z .
\]
Each \(k_i\ge i-1\). The denominator ratio on the right
contains exactly \(a\) consecutive-factor contributions,
all at most \(H\), and each has \(\ell\)-adic valuation at
most \(\lfloor\log_\ell H\rfloor\). After division by
\(\kappa\), the corresponding term of (4) thus has valuation
at least \(a(s-\lfloor\log_\ell H\rfloor)\ge1\).
Summing proves (3).

The threshold in (3) can be attained for a general weight
list. Take \(M=\ell^a\), \(a\ge1\), and weights
\(0,1,\ldots,M-2,M\). The Vandermonde quotient is symmetric
and homogeneous of degree one. Its diagonal value at
\((1,\ldots,1)\) is \(\kappa=M\) by the displayed product
formula, so \(\Phi(x)=\sum_i x_i\). Choose distinct integers
\(t_i\) with \(\sum_i t_i\equiv-1\bmod\ell\), and take
\(x_i=1+\ell^a t_i\). The quotient value has valuation at
least \(a+1\), while \(v_\ell(\kappa)=a\), even though
the \(x_i\) coincide modulo \(\ell^a\). This example
concerns general weight lists.

## The older diagonal valuation can be large

Take \(r=2^b+1\) with \(b\) large and use \(e=r+1,\ f=r\).
For the \((h+1)\)-by-\((h+1)\) grid of weights \(f i+e j\),
assume \(6\le h\le 2^{b/2}\). The exponents are consecutive
and coprime, and \(e\equiv2f\bmod 2^b\) with \(f\) odd.
For every \(Q=2^\nu\) with \(16h\le Q\le2^b\),
the weight residues coincide exactly when the values of
\(i+2j\) agree, since those values lie in \([0,3h]\).
Writing \(M=(h+1)^2\), the number of grid-weight pairs
equal modulo \(Q\) is at least
\(M^2/(2(3h+1))-M/2\ge M^2/(16h)\).
For the consecutive weights \(0,\ldots,M-1\), the
corresponding count is at most \(M^2/(2Q)\). Indeed each
residue class of size \(a\) has \(a-1\le M/Q\).
The difference is at least \(M^2/(32h)\) at each of
these levels. Differences at other levels are
nonnegative because consecutive weights distribute the
fixed total among residue classes as evenly as possible.
Summing the level differences in the pair-valuations
identity yields
\[
 v_2(\kappa)\ge \frac{(h+1)^4}{32h}
        \bigl(b-\lceil\log_2(16h)\rceil+1\bigr)
        \gg h^3\log r                                                \tag{5}
\]
for sufficiently large \(b\). Here reversing the weight
order and translating by a constant to get the \(w_j\)
preserves all pair differences. This matches the order of
the earlier general primewise upper bound along these
parameters. For \(n=2^r\), a divisor-selection rule with
\(\lfloor r/B_*\rfloor\ge v_2(\kappa)+1\) and
\(B_*\asymp_A h^5\) forces \(h^8\log r\ll_A r\).
Lemma (3) replaces that condition by one depending only
on the largest weight \(H\).
[Pair-valuations identity and prior divisor cap](POWERFUL_MONOMIAL_VALUATION_REFINEMENT.md)

## Selecting the divisor and counting rational pairs

Assume \(e>f\), interchanging coordinates and inverting
\(K\) if needed, and put \(d=e-f\le A\).
Choose \(h\) to be \(\lfloor c_A(r/\log r)^{1/5}\rfloor\)
with a sufficiently small fixed \(c_A>0\). In the case
\(n_{\ge r}\ge n^{1/2}\), where \(n_{\ge r}\) denotes the
full-power part over primes at least \(r\), instead take
\(h=\lfloor c_A r^{1/5}\rfloor\). In each case use
\(M=(h+1)^2\), \(J=\binom M2\), and
\[
 \delta=6A/h,\qquad t=\lceil4/\delta\rceil,\qquad
 B_*=\lceil(4J+2)/\delta\rceil .
\]
All statements below concern \(r\) sufficiently large in
terms of \(A\). In particular \(\delta\le1\),
\(4\le t\delta\le5\), and \(4dt<r\).

Use multiplicative amplification from the manuscript.
For each \(\eta>0\), the coordinatewise \(t\)-fold
products of a solution set have size
\(\gg_{t,\eta}|\mathcal S|^tY^{-t\eta}\).
Indeed prime-sign coloring in each coordinate retains
\(|\mathcal S|Y^{-o(1)}\) pairs whose products have no
numerator-denominator cancellation, and ordered
factorizations of the four product integers have
multiplicity \(Y^{o(1)}\) by the divisor bound for fixed
\(t\). The products have height at most \(Z=Y^t\) and
satisfy the same congruence with \(K^t\). Work with that
amplified set and rename its new residue \(K\).

Write \(n=\prod_\ell\ell^{a_\ell}\). In the general
case take all primes of \(n\) as eligible. In the case
\(n_{\ge r}\ge n^{1/2}\), use only \(\ell\ge r\).
For grid pairs \(0\le i,j\le h\), choose \(w_j\) to be
the sorted list obtained by subtracting each weight
\(fi+ej\) from their maximum \(H=(e+f)h\).
These weights are distinct since
\(\gcd(e,f)\le A\) and \(h<\min(e,f)/A\). Define
\(s_\ell=1+\lfloor\log_\ell H\rfloor\).
We have \(B_*\ll_A h^5\), and we can choose \(c_A\)
so that
\[
                       B_*s_\ell\le r/2
                       \quad\text{on every eligible prime}.
                                                                        \tag{6}
\]
In the first case this uses \(s_\ell\ll_A\log r\).
In the second \(H<r^2\) for large \(r\), so
\(s_\ell\le2\) on the eligible primes. Also take
\(B_*\le r/40\) in both cases.

Cap exponents of eligible primes at
\(\lfloor a_\ell/B_*\rfloor\), which is at least
\(s_\ell\) by (6), and at least \(a_\ell/(2B_*)\).
Their total capped divisor is at least
\(n^{1/(4B_*)}\), since the eligible full-power part
is at least \(n^{1/2}\). This exceeds or equals
\(Z^\delta\) by \(Y\le n^{2/r}\) and \(t\delta\le5\).
Activate each selected prime at exponent \(s_\ell\),
with subsequent increases one at a time up to the cap,
and stop at the first divisor \(L\ge Z^\delta\).
A prime \(\ell\mid n\) is at most \(n^{1/r}\le Z^\delta\)
by \(Y\ge n^{1/(3r)}\). The activation factor satisfies
\(\ell^{s_\ell}\le \ell H\). Hence
\[
 Z^\delta\le L\le H Z^{2\delta},\qquad
 L^{B_*}\mid n_0,\qquad n_0\ge Z^{4J+2},                              \tag{7}
\]
where \(n_0\) contains the full powers of all activated
primes.

For an amplified solution put \(T=U/V\). Then
\(U^d\equiv K T^{-f}\) and \(V^d\equiv K T^{-e}\pmod n\).
Distinct pairs give distinct rational \(T\). Indeed equal
\(T\) implies congruent \(d\)-th powers of the two
coordinates, which are actually equal by the
cross-numerator bound \(Z^{2d}<n\), using \(4dt<r\).
Positivity gives equality of the pairs.

If \(M\) distinct pairs lie in one residue class of \(T\)
modulo \(L\), apply (3) at each prime of \(L\) to their
\(T_i\), taking \(z=T_1\). The value of \(\Phi\) has exactly
valuation \(v_\ell(\kappa)\) at those primes. The
alternant \(\det(T_i^{w_j})\) therefore has the same selected
valuations as \(\kappa\) times the product of all
cross-numerators of \(T_j-T_i\). Each \(T_i\) has height at
most \(Z^2\), so this nonzero product has magnitude at most
\(\kappa Z^{4J}\). For \(Z^2>\kappa\) it is smaller than
\(n_0\). Thus \(n_0\) does not divide the alternant.

Evaluate the \(M\) monomials \(U^{di}V^{dj}\) and clear the
two denominators in each row to exponent \(dh\). The
resulting integer determinant is congruent modulo \(n_0\)
to the alternant times a unit, since the column multipliers
are \(K^{i+j}\) and the row multipliers remove the common
negative weight. It is nonzero and divisible by \(L^J\).
Its magnitude is at most \(M! Z^{2dhM}\). But
\(\delta J/M=3A(h+2)>2dh\), so the upper bound is smaller
than \(L^J\) for sufficiently large \(Z\). For both
determinant inequalities an explicit height cutoff
independent of \(n,K\) is
\[
                     Z>(H+1)^J M!+1,
\]
using \(\kappa\le H^J\) and
\(\delta J-2dhM\ge A(h+6)M\ge1\).
The complementary bounded-height range has at most the
fourth power of this cutoff possible rational pairs,
uniformly in \(n,K\).

Each class contains at most \(M-1\) amplified pairs.
By (7) there are \(O_{A,r}(Z^{2\delta})\) such pairs.
Undoing amplification gives
\(O_{A,r,\varepsilon}(Y^{2\delta+\varepsilon})\).
The two choices of \(h\) prove (1) and its log-free
restricted-modulus variant. The unit-only modulus \(n=1\)
has \(Y=1\) and a single possible pair.

## The log-free bound for primitive triples

Take \(n=c\), the \(r\)-full target in the primitive triple.
If \(c_{\ge r}\ge c^{1/2}\), use the strengthened modulus
bound in the manuscript's aligned graph argument, with
\(Y=c^{2/r}\). For \(u>v\) the two aligned omissions
\(i=0,1\) have distinct exponents \(r+u+i,r+v+i\), which
lie in an interval of fixed length about \(r\). All
fractions are units by primitivity, and the short-height
condition on the ratios follows in the competitive boxes
from \(c\ge X^{1-O_{u,v}(1/r)}\). This gives codegree
\(X^{O_{u,v}(r^{-6/5})+\varepsilon}\) by the unchanged
integer common-neighbor recovery. The coefficient-uniform
edge and surface bounds have the same scope as before.
The aligned graph inequality in normalized masses now
has error \(O_{u,v}(r^{-1/5})\), since multiplying the
codegree exponent by \(r\) yields that order. The two
aligned cost certificates yield
\(D\ge u+v+2-O_{u,v}(r^{-1/5})\). The count exponent is
at most \(2/r-D/r^2+O_{u,v}(r^{-3})\), giving (2) on
these targets. For \(u=v\) the manuscript's aligned-gcd
case already gives the smaller error \(O_u(r^{-5/2})\).
[Graph inequalities and cost certificates in the earlier note](POWERFUL_MONOMIAL_VALUATION_REFINEMENT.md#transfer-to-fixed-offset-primitive-triples)

If \(c_{\ge r}<c^{1/2}\), the part of \(c\) supported on
primes below \(r\) has \(O_r((1+\log X)^r)\) possibilities.
The other part is \(r\)-full and at most \(X^{1/2}\),
so has \(O_r(X^{1/(2r)})\) possibilities. The latter
bound follows from the canonical factorization
\(m=a_r^r\prod_{j=r+1}^{2r-1}a_j^j\), by summing the
\(a_r\) count over the product of the convergent series
\(\sum_{a_j\ge1}a_j^{-j/r}\). For each such \(c\) choose
the \(p\)-full summand \(a\le X\) in \(O_p(X^{1/p})\)
ways by the same counting argument. Then \(b\) is
determined. Thus this portion contributes at most
\[
 O_{u,v,r,\varepsilon}(X^{1/p+1/(2r)+\varepsilon}),
\]
whose exponent apart from \(\varepsilon\) is below
\(1/p+1/q-2/r^2\) for all sufficiently large \(r\)
at fixed offsets. Combining the portions proves the
unrestricted estimate (2) with no logarithmic factor
in its exponent error.

The sharp-error obligation
\(\texttt{GAP-ARITH-POWERFUL-FIXED-OFFSET-SHARP-ERROR}\)
asks for an error \(O_{u,v}(r^{-5/2})\), beyond (2).
The small fixed powers in the balanced square-full
kernel graph fall outside the large-power hypotheses.

Exact signed-integer determinant checks confirm (3)
on 426 weight/prime examples, including 66 grid-weight
cases. Direct residue counting checks (5) at
\((b,h)=(12,6),(16,32),(20,64),(24,90)\).
