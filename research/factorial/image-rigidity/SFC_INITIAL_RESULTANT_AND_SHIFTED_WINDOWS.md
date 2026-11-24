# Exact initial resultant and two uniform shifted windows

## Scope

For a characteristic-zero field \(K\), parameters \(r_1,\ldots,r_s\in K\), and coefficient variables \(c_1,\ldots,c_s\), set

\[
g_m(\mathbf c;\mathbf r)=[z^m]\prod_{i=1}^s(1-c_i z)^{-mr_i-1},\qquad g_0=1.
\tag{1}
\]

Formal binomial expansion makes this a homogeneous polynomial of degree \(m\) in the \(c_i\), with coefficients in \(\mathbb Q[\mathbf r]\). For nonnegative integer \(r_i\), these are the actual normalized factorial moments of

\[
F=\Bigl(\prod_i X_i^{r_i}\Bigr)\sum_i c_iX_i,
\qquad
g_m=\frac{\Lambda(F^m)}{m!\prod_i(mr_i)!}.
\tag{2}
\]

The main result here is an exact resultant for the **initial** window \(g_1,\ldots,g_s\), for every \(s\ge1\) and all parameters. Two further complete polynomial certificates prove the windows \(g_2,g_3,g_4\) for exponents \((a,a,b)\), and \(g_3,g_4,g_5\) for exponents \((1,1,b)\), for all integers \(a,b\ge0\). They do not settle all shifted windows, even for the latter ternary families.

## 1. A triangular relation to weighted elementary polynomials

Write

\[
P(z)=\prod_i(1-c_i z),\quad A(z)=\prod_i(1-c_i z)^{r_i},\quad a(z)=zA(z),
\]

where \(A\) is a formal series when the parameters are not nonnegative integers. Define

\[
E_{\mathbf r}=\sum_i r_i c_i\partial_{c_i},\qquad
\delta_k=(1+E_{\mathbf r})e_k
=\sum_{|I|=k}\left(1+\sum_{i\in I}r_i\right)\prod_{i\in I}c_i.
\tag{3}
\]

Then the polynomial

\[
D(z)=P(z)+zP(z)A'(z)/A(z)
=1+\sum_{k=1}^s(-1)^k\delta_kz^k
\tag{4}
\]

satisfies the exact formal identity

\[
D(z)\sum_{m\ge0}g_m a(z)^m=1.
\tag{5}
\]

For completeness, take \(H(0)=0,H'=A/P\) and \(B=H\circ a^{-1}\). Formal Lagrange inversion gives

\[
[t^{m+1}]B=\frac1{m+1}[z^m]H'(z)A(z)^{-m-1}
=\frac{g_m}{m+1}.
\]

Since \(a'=H'D\), differentiation gives \(B'(a)=1/D\), proving (5). All units used have constant term one. Equivalently, the Lagrange coefficient identity follows by substituting \(t=a(z)\) in the formal coefficient residue and integrating by parts; no analytic convergence or choice of branches is involved.

Taking the coefficient of \(z^k\) in (5) gives, for \(1\le k\le s\),

\[
g_k=(-1)^{k+1}\delta_k-
\sum_{j=1}^{k-1}g_j[z^k](D a^j).
\tag{6}
\]

The multiplier of \(g_j\) is homogeneous of coefficient-degree \(k-j\). Thus

\[
(g_1,\ldots,g_k)=(\delta_1,\ldots,\delta_k)
\quad\text{in }\mathbb Q[\mathbf r,\mathbf c].
\tag{7}
\]

This is a polynomial triangular change of generators; no subset sum or coefficient is inverted.

## 2. The all-support exact resultant

Use the homogeneous resultant in the fixed degrees \(1,2,\ldots,s\), normalized by

\[
\operatorname{Res}(c_1,c_2^2,\ldots,c_s^s)=1.
\]

**Theorem 1.** As a polynomial identity in \(\mathbb Q[\mathbf r]\),

\[
\boxed{\operatorname{Res}_{\mathbf c}(g_1,\ldots,g_s)
=\prod_{\varnothing\ne I\subseteq[s]}
\left(1+\sum_{i\in I}r_i\right)^{(|I|-1)!(s-|I|)!}.}
\tag{8}
\]

In particular, the sign and constant are exactly \(+1\).

**Proof.** Let \(\mathcal R_s\) denote the left side. Equation (6) permits replacing the final generator by

\[
(-1)^{s+1}\delta_s
=(-1)^{s+1}\left(1+\sum_i r_i\right)\prod_i c_i.
\]

Homogeneous addition of multiples of earlier generators does not change the resultant. Homogeneity in the final equation extracts its scalar to power \((s-1)!\). Multiplicativity in the final equation then splits the product of coordinate variables. Setting \(c_j=0\) in each earlier \(g_i\) gives exactly the same moment polynomial in the remaining \(s-1\) variables and parameters. The coordinate restriction has sign

\[
\operatorname{Res}(g_1,\ldots,g_{s-1},c_j)
=(-1)^{(s-j)(s-1)!}\mathcal R_{s-1}(\mathbf r_{\widehat j}).
\tag{9}
\]

Indeed, move the \(j\)-th coordinate to the last position, using \(s-j\) transpositions; coordinate covariance contributes that sign to the power \(1\cdot2\cdots(s-1)\). With the last equation equal to the last coordinate, restriction introduces no further factor under the stated normalization.

The total sign from these operations is

\[
(-1)^{(s-1)!\{s+1+s(s-1)/2\}}=+1.
\]

This is direct for \(s=1,2\); for \(s\ge3\), \((s-1)!\) is even. Hence

\[
\mathcal R_s(\mathbf r)
=\left(1+\sum_i r_i\right)^{(s-1)!}
\prod_{j=1}^s\mathcal R_{s-1}(\mathbf r_{\widehat j}).
\tag{10}
\]

Starting with \(\mathcal R_1=1+r_1\), a proper subset \(I\) of size \(k\) appears in \(s-k\) factors, each with exponent \((k-1)!(s-k-1)!\). Its exponent is therefore \((k-1)!(s-k)!\). The full subset has exponent \((s-1)!\). This proves (8), including multiplicities, without a generic simplicity or local-intersection assumption. ∎

The resultant properties used here are its normalization, covariance, homogeneity, multiplicativity, and invariance under homogeneous elementary transformations. Their exact statements are in G. Staglianò, *A package for computations with classical resultants*, JSAG 8 (2018), 21–30, Theorem 1.1 and Proposition 1.5, pp. 22–23, [journal PDF](https://msp.org/jsag/2018/8-1/jsag-v8-n1-p03-s.pdf), DOI [10.2140/jsag.2018.8.21](https://doi.org/10.2140/jsag.2018.8.21). The coordinate restriction follows from the same normalization and resultant definition. These are standard imported resultant identities; the subset-product application is derived above.

### Consequences and resonances

The initial system has no nonzero projective zero over an algebraic closure exactly when

\[
1+\sum_{i\in I}r_i\ne0
\quad\text{for every nonempty }I.
\tag{11}
\]

For nonnegative integer exponents these factors are positive integers, so the initial window is regular. Its nonvanishing was already a consequence of the retained contact identity; (8) supplies an exact resultant and exact multiplicities. It does not establish a new arbitrary-start SFC theorem.

For integer parameters and any prime \(p\), whenever the factors in (8) are nonzero,

\[
v_p\mathcal R_s
=\sum_{\varnothing\ne I\subseteq[s]} (|I|-1)!(s-|I|)!
v_p\left(1+\sum_{i\in I}r_i\right).
\tag{12}
\]

A minimal resonant example has \(k\) active coefficients equal to all \(k\)-th roots of unity, and \(r_i=-1/k\) on those active indices. Then \(P=1-z^k\), \(A=(1-z^k)^{-1/k}\), and

\[
g_m=[z^m](1-z^k)^{m/k-1}=0\qquad(m\ge1).
\tag{13}
\]

If \(k\nmid m\), there is no term in that congruence class. If \(m=k\ell\), the polynomial on the right has degree \(k(\ell-1)<m\). No proper nonempty active subset is resonant. Inactive coefficient variables may be set to zero. These examples have negative formal exponents and are not counterexamples to ordinary factorial SFC.

In fact, for any parameters, an initial projective zero makes every \(\delta_k\) vanish by (7), hence \(D=1\). Equation (5) and invertibility of \(a\) show that all positive \(g_m\) vanish at that coefficient vector. This conclusion is specific to a zero starting at \(m=1\).

## 3. Why the initial formula cannot simply be shifted

For two variables, put \(c_1=t,c_2=1\), with parameters \(r,s\). Direct polynomial elimination gives

\[
\operatorname{Res}_t(g_2,g_3)
=\frac{(r+1)(s+1)(r+s+1)^2}{4}\,Q(r,s),
\tag{14}
\]

where

\[
\begin{aligned}
Q={}&81r^3s^2+74r^3s+17r^3+81r^2s^3+257r^2s^2\\
&+181r^2s+37r^2+74rs^3+181rs^2+115rs+22r\\
&+17s^3+37s^2+22s+4.
\end{aligned}
\]

All coefficients of \(Q\) are positive, but it is a further nonlinear factor, absent from the initial subset product. At \(r=s=-1/3\), all initial subset factors are nonzero, whereas

\[
g_3\equiv0,\qquad g_2(t)=\frac{2t^2+t+2}{9}.
\tag{15}
\]

The latter has two nonzero complex roots. Thus the same parameter resonance criterion fails for a shifted window in the formal parameter extension. This is a negative control on the proposed extrapolation, not a counterexample with nonnegative integer exponents.

## 4. A complete uniform ternary window at orders 2–4

**Theorem 2.** For all integers \(a,b\ge0\), the three homogeneous polynomials \(g_2,g_3,g_4\) for exponents \((a,a,b)\) have no nonzero common zero over any characteristic-zero field extension.

On the chart \(c_3\ne0\), scale it to one and put \(v=c_1+c_2,w=c_1c_2\). Then

\[
g_m(v,w)=[t^m](1-vt+wt^2)^{-am-1}(1-t)^{-bm-1}
\tag{16}
\]

and, with \(A_2=2a+1\),

\[
g_2=A_2\{(a+1)v^2+(2b+1)v-w\}+(2b+1)(b+1).
\]

Consequently \(g_2=0\) forces

\[
w=W=(a+1)v^2+(2b+1)v+\frac{(2b+1)(b+1)}{2a+1}.
\tag{17}
\]

The denominator is nonzero throughout the stated parameter range. Define the integer polynomials

\[
F=6(2a+1)g_3(v,W),\qquad
G=24(2a+1)^2g_4(v,W).
\]

Their exact resultant is

\[
\operatorname{Res}_v(F,G)
=41472(b+1)(a+1)^2(2a+1)^7(2a+b+1)(a+b+1)^3Q_{70}(a,b).
\tag{18}
\]

The full \(F,G,Q_{70}\) coefficients are retained in [the exact witness](sfc_ternary_orders_2_4_witness.json). The polynomial \(Q_{70}\) has 70 strictly positive coefficients, total degree 12, and constant coefficient 432. Therefore (18) never vanishes for \(a,b\ge0\).

This is a complete identity certificate: the independently reconstructed \(F,G\) have degrees in \((a,b,v)\) at most \((4,3,3)\) and \((6,4,4)\). The resultant consequently has degrees at most (34) in \(a\) and (24) in \(b\). The checker verifies equality of the exact integer Sylvester determinants and the displayed product at the entire \(35\times25\) grid, and checks those same degree bounds on the product. Polynomial interpolation proves (18) identically, rather than only for the tested parameters.

The boundary \(c_3=0\) is also included. There

\[
g_2=(2a+1)((a+1)v^2-w),
\]

and after its vanishing,

\[
g_3=-\frac{(3a+1)(3a+2)(a+1)}2v^3.
\]

Thus \(v=w=0\), which forces \(c_1=c_2=0\). This completes the proof. ∎

## 5. A complete uniform ternary window at orders 3–5

**Theorem 3.** For all integers \(b\ge0\), the homogeneous polynomials \(g_3,g_4,g_5\) for exponents \((1,1,b)\) have no nonzero common zero over any characteristic-zero field extension.

Use the same chart and variables as in (16), with \(a=1\). Set

\[
\begin{aligned}
N(v,b)={}&40v^3+(60b+20)v^2+(36b^2+36b+8)v\\
&+9b^3+18b^2+11b+2,\\
J(v,b)={}&8(3b+5v+1).
\end{aligned}
\]

Then \(2g_3=N-Jw\). A zero of \(J\) cannot be a zero of \(N\), since

\[
N(-(3b+1)/5,b)=\frac{(b+2)(3b+1)(3b+11)}{25}\ne0.
\tag{19}
\]

We may therefore substitute \(w=N/J\) and set

\[
F=24J^2g_4(v,N/J),\qquad G=120J^2g_5(v,N/J).
\]

They are integer polynomials, with \((b,v)\)-degrees at most \((6,6)\) and \((7,7)\). Their resultant is

\[
\begin{aligned}
\operatorname{Res}_v(F,G)={}&-185900678964646453482356736000000000000\\
&\cdot(b+1)(b+3)^2(3b+1)^4(3b+11)^4(b+2)^9Q_{22}(b).
\end{aligned}
\tag{20}
\]

The full polynomials are in [the exact witness](sfc_ternary_orders_3_5_witness.json). All 23 coefficients of the degree-22 polynomial \(Q_{22}\) are strictly positive. The checker reconstructs the elimination polynomials from (16) and verifies (20) at \(b=0,\ldots,84\). Both sides have degree at most \(7\cdot6+6\cdot7=84\), so these exact determinant values certify the full identity.

On \(c_3=0\),

\[
g_3=20v(v^2-w),\qquad
g_4=70v^4-105v^2w+15w^2.
\]

Their only common zero is \(v=w=0\): either \(v=0\), or \(w=v^2\) and then \(g_4=-20v^4\). This proves the theorem. ∎

## 6. Verification and remaining gap

Run from the repository root:

```sh
python3 research/factorial/image-rigidity/verify_sfc_initial_resultant.py
```

The standard-library checker pins this note and both polynomial witnesses. It independently constructs their eliminated polynomials, checks all 960 exact interpolation determinants, verifies the explicit boundary and shifted countercontrols, and checks the triangular formal-series identity and initial resultant using finite exact Macaulay determinant controls. The latter controls support the hand proof of the all-\(s\) formula; they do not replace its induction. No finite scan is presented as a proof of an arbitrary-start statement.

The all-start target remains:

\[
V(g_n,g_{n+1},\ldots,g_{n+s-1})=\{0\}
\quad\text{for every }n\ge1\text{ and nonnegative integer }\mathbf r.
\]

Equation (8) proves only \(n=1\), while Theorems 2 and 3 supply the stated later windows for the stated uniform ternary families. A proof for arbitrary shifted windows still needs a new structural control on their resultants or scalar contact; initial subset resonances alone do not provide it.
