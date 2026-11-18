# Fixed periodic factors in divisor boxes

Three interval factors remain uniformly distributed beyond the square-root
range when each is multiplied by a bounded function of period 12. The proof
uses completion modulo \(12q\). Its additional local factors separate in
the three dual variables and disappear before the incomplete
Kloosterman correlation is estimated.

## Theorem FixedPeriodicBoxes

Let \(\omega_1,\omega_2,\omega_3:\mathbb Z\to\mathbb C\) be 12-periodic
functions with \(|\omega_i|\le1\), and let
\[
 I_i=[N_i,N'_i]\cap\mathbb Z,\qquad
 1\le N_i\le N'_i\le2N_i,\qquad X=N_1N_2N_3.
\]
Write \(*\) for Dirichlet convolution, and put
\[
 f=(\omega_1\mathbf1_{I_1})*
   (\omega_2\mathbf1_{I_2})*
   (\omega_3\mathbf1_{I_3}),\qquad
 \Delta(f;q,a)=
 \sum_{n\equiv a\ (q)}f(n)
 -\frac1{\varphi(q)}\sum_{(n,q)=1}f(n).
\]
For
\[
 \delta=\frac1{10000},
\]
there is an absolute constant \(C\) such that
\[
 |\Delta(f;q,a)|\le C\frac{X^{1-\delta}}q
 \quad\text{if}\quad
 (q,12)=(a,q)=1,\qquad q\le X^{1/2+\delta}.
 \tag{1}
\]
The constant is uniform in the intervals, the residue \(a\), and the
periodic functions.

In particular, (1) applies to arbitrary triples of Dirichlet characters
modulo 12, including principal and imprimitive characters, with their
literal values zero on nonunits. It also applies to one or two factors:
pad the convolution with singleton factors at 1.

For additional functions \(v_i\) of bounded variation on
\([N_i,2N_i]\), the same estimate for
\((\omega_1v_1\mathbf1_{I_1})*
 (\omega_2v_2\mathbf1_{I_2})*
 (\omega_3v_3\mathbf1_{I_3})\)
has the extra factor
\[
 \prod_{i=1}^3
 \bigl(\|v_i\|_\infty+\operatorname{Var}_{[N_i,2N_i]}v_i\bigr).
 \tag{2}
\]
This includes smooth interval cutoffs with uniformly bounded variation.

The modulus condition \((q,12)=1\) is part of the theorem. There is no
claim here about moduli with unbounded powers of 2 or 3, or about a binary
representation problem.

## Primary input

The imported analytic source is D. R. Heath-Brown,
*The divisor function \(d_3(n)\) in arithmetic progressions*,
Acta Arithmetica **47** (1986), 29–56:
[DOI](https://doi.org/10.4064/aa-47-1-29-56),
[primary archive](https://matwbn.icm.edu.pl/ksiazki/aa/aa47/aa4713.pdf),
[cached primary scan](sources/heath-brown-d3-progressions-1986.pdf).

The proof below transfers the three box estimates in that paper,
Lemmas 5, 6 and 7, printed pp.44 and 47. It does so through their proofs,
rather than applying an indicator-only theorem to weighted indicators.
The deep input is its all-modulus complete-correlation Lemma 4, printed
p.40. The proof of that lemma occupies Section 3, printed pp.35–40, and
uses the Birch–Bombieri finite-field estimate cited there.

Here is its precise form in the notation needed below. Set
\[
 e_m(z)=\exp(2\pi iz/m),\qquad
 K_2(A,B,C;q)=
 \sum_{\substack{x,y,z\pmod q\\xyz\equiv1\ (q)}}
 e_q(Ax+By+Cz),
\]
and define the complete correlation
\[
 \mathcal C(k,t_1,t_2,\rho,\sigma;q)=
 \sum_{j\pmod q}^{*} e_q(kj)
 K_2(\rho,\sigma,jt_1;q)
 \overline{K_2(\rho,\sigma,jt_2;q)}.
\]
An asterisk restricts the summed variable to units. For some absolute
constant \(A>0\), Heath-Brown's Lemma 4 gives
\[
 \begin{split}
 |\mathcal C(k,t_1,t_2,\rho,\sigma;q)|
 \ll{}& d_A(q)q^{5/2}
 \bigl\{
 (k^2,(t_1-t_2)^2,q)
 (t_1,q)(t_2,q)(\rho,q)(\sigma,q)
 \bigr\}^{1/2}\\
 &{}\times ((t_1-t_2)^3,q)^{1/6}.
 \end{split}
 \tag{3}
\]
Parentheses denote greatest common divisors; \(d_A(q)\ll_\varepsilon
q^\varepsilon\) is the fixed divisor-function factor in that lemma.
There is no restriction on the additive parameter \(k\).

We also use the source's pointwise bound
\[
 |K_2(A,B,C;q)|\le q(A,B,C,q)d_3(q)
 \tag{4}
\]
from equation (3.1), printed p.35, the ordinary Weil–Estermann
Kloosterman bound used on printed p.44, and the elementary reciprocal
energy estimate proved in equations (6.6)–(6.13), printed pp.46–47.
These inputs concern the original unweighted sums.

## Completion with the fixed period

Fix \(c=12\), put \(Q=cq\), and define
\[
 W_i(h)=\sum_{n\in I_i}e_Q(-hn),\qquad
 G_i(k)=\sum_{r\pmod c}\omega_i(r)e_c(kr).
\]
Thus \(|G_i(k)|\le c\). Choose representatives
\(-Q/2<h\le Q/2\). The geometric-series estimate gives
\[
 |W_i(h)|\ll
 \min\bigl(\#I_i,Q/|h|\bigr)
 \qquad(0<|h|\le Q/2).
 \tag{5}
\]
Finite Fourier inversion reconstructs the number of interval integers
in each residue class, so this completion remains exact if an interval
is longer than \(Q\).

Let \(S(a)=\sum_{n\equiv a\ (q)}f(n)\). The Chinese remainder theorem gives
\[
 \begin{split}
 S(a)=Q^{-3}\!
 \sum_{h_1,h_2,h_3\pmod Q}
 &\prod_{i=1}^3
 \bigl(W_i(h_i)G_i(q^{-1}h_i)\bigr)\\
 &{}\times
 K_2(c^{-1}h_1,c^{-1}h_2,ac^{-1}h_3;q).
 \end{split}
 \tag{6}
\]
The inverse in \(G_i\) is modulo \(c\); those in \(K_2\) are modulo \(q\).
Indeed, for an integer \(x\),
\[
 e_Q(hx)=e_c(q^{-1}hx)e_q(c^{-1}hx).
\]
The product congruence constrains only the components modulo \(q\).
The components modulo \(c\) therefore give the three separate \(G_i\).

Changing the unit variables in \(K_2\) rewrites its factor in (6) as
\[
 K_2(h_1,h_2,a'h_3;q),
 \qquad a'\equiv ac^{-3}\pmod q.
 \tag{7}
\]
The residue \(a'\) is still a unit. All dual variables are therefore
centered at zero; the price of the fixed period is a bounded separable
factor and the larger frequency range \(|h_i|\le cq/2\).

If at least one \(h_i\) is zero modulo \(q\), its \(K_2\) factor is
independent of \(a\). For example, with \(c_q\) denoting a Ramanujan sum,
\[
 K_2(0,s,at;q)=c_q(s)c_q(at)=c_q(s)c_q(t).
\]
There are \(O(c)\) multiples of \(q\) in each chosen frequency range.
All terms having such a frequency can be included in a main term
independent of \(a\).

Different estimates below may produce different main terms \(M_j\).
They need only be independent of \(a\). If
\[
 |S(a)-M_j|\le E_j
\]
for every unit \(a\), averaging over those residues gives
\[
 |\Delta(f;q,a)|\le2E_j.
 \tag{8}
\]
Thus the estimates can be used separately without identifying any of
their main terms.

## Three estimates for weighted boxes

Put
\[
 U=N_1,\quad V=N_2,\quad W=N_3,\qquad
 I=\#I_1,\quad J=\#I_2,\quad K=\#I_3.
\]
For every \(\varepsilon>0\), each of the following bounds holds with
an \(a\)-independent main term:
\[
 E_{\mathrm{pr}}\ll_\varepsilon X^\varepsilon
 \left\{q^{5/6}+q^{1/4}(IJ)^{1/2}+q^{1/2}K^{1/2}\right\},
 \tag{9}
\]
\[
 E_{\mathrm{small}}\ll_\varepsilon X^\varepsilon q^{1/2}I,
 \tag{10}
\]
\[
 E_{\mathrm{long}}\ll_\varepsilon X^\varepsilon
 I^{3/4}J^{1/2}(K+U)^{1/4}
 \left(1+\frac{V^2}q+\frac{V^{3/2}}{q^{1/2}}\right)^{1/4}.
 \tag{11}
\]
The three variables may be permuted in each estimate. We prove that the
fixed-period factors preserve these three estimates.

### The principal estimate

Remove the frequencies zero modulo \(q\) from (6), and split the others
into dyadic ranges
\[
 R<|r|\le2R,\qquad S<|s|\le2S,\qquad T<|t|\le2T.
\]
Each range lies within \(|h|\le cq/2\). Write
\[
 F_I(R)=\min(I,q/R),\qquad
 F_J(S)=\min(J,q/S),\qquad F_K(T)=\min(K,q/T).
\]
Since \(c\) is fixed, these are envelopes for (5), up to fixed constants.

The triangle inequality over \(r,s\) bounds the contribution of one
dyadic box by
\[
 C_cq^{-3}F_I(R)F_J(S)
 \sum_{r,s}
 \left|\sum_t b_tK_2(r,s,a't;q)\right|,
 \qquad
 b_t=G_3(q^{-1}t)W_3(t),
 \tag{12}
\]
where \(|b_t|\ll_c F_K(T)\). The first two local factors have been
bounded in absolute value **before** combining \(r\) and \(s\).

As on printed pp.41–42 of the source, write
\[
 r=\rho r',\qquad s=\sigma s',\qquad (r's',q)=1,
\]
where the signed integers \(\rho,\sigma\) have prime factors only among
those dividing \(q\). The number of choices with
\(|\rho|,|\sigma|\le cq\) is \(O_{c,\varepsilon}(q^\varepsilon)\);
this is the smooth-number estimate used in source equation (4.2).
The change of unit variables gives
\[
 K_2(\rho r',\sigma s',a't;q)
 =K_2(\rho,\sigma,a'r's't;q).
 \tag{13}
\]
In particular, the original frequencies need not be units.

Fix \(\rho,\sigma\), set
\[
 R_0=R/|\rho|,\quad S_0=S/|\sigma|,\quad H=R_0S_0,
 \qquad \mathcal K(z)=K_2(\rho,\sigma,a'z;q),
\]
and absorb signs into \(\rho,\sigma\), so that \(r',s'>0\).
Combining \(h=r's'\) and applying Cauchy's inequality gives
\[
 \begin{split}
 &\sum_{r',s'}\left|\sum_t b_t\mathcal K(r's't)\right|\\
 &\quad\le
 \left(\sum_{H<h\le4H}d(h)^2\right)^{1/2}
 \left(
 \sum_{t_1,t_2}b_{t_1}\overline{b_{t_2}}
 \sum_{\substack{H<h\le4H\\(h,q)=1}}
 \mathcal K(ht_1)\overline{\mathcal K(ht_2)}
 \right)^{1/2}.
 \end{split}
 \tag{14}
\]
The expression in the second square root is nonnegative before
expansion. After expansion, triangle inequality bounds it by
\[
 C_cF_K(T)^2\sum_{t_1,t_2}|\Sigma(t_1,t_2)|,
 \qquad
 \Sigma(t_1,t_2)=
 \sum_{\substack{H<h\le4H\\(h,q)=1}}
 \mathcal K(ht_1)\overline{\mathcal K(ht_2)}.
 \tag{15}
\]
There is no periodic factor in the \(h\)-sum. The third local factor
remains outside this correlation in the coefficients \(b_t\).

Complete the \(h\)-sum:
\[
 \Sigma(t_1,t_2)=\frac1q
 \sum_{k\pmod q}
 \mathcal C(a'^{-1}k,t_1,t_2,\rho,\sigma;q)
 \sum_{H<h\le4H}e_q(-kh).
 \tag{16}
\]
This uses precisely the complete correlation (3). The interval sum is
\(O(H)\) at \(k=0\) and \(O(q/|k|)\) for centered nonzero \(k\). The
divisor estimate
\[
 \sum_{1\le k\le q}k^{-1}(k^2,D)^{1/2}
 \ll d(D)\log(2q)
\]
therefore gives the source's equation (4.7):
\[
 \begin{split}
 |\Sigma(t_1,t_2)|\ll_\varepsilon{}&
 q^{5/2+\varepsilon}
 \bigl\{(t_1,q)(t_2,q)(\rho,q)(\sigma,q)\bigr\}^{1/2}
 ((t_1-t_2)^3,q)^{1/6}\\
 &{}\times
 \left\{\frac Hq((t_1-t_2)^2,q)^{1/2}+1\right\}.
 \end{split}
 \tag{17}
\]
Grouping the noncongruent pairs by
\((t_1,q),(t_2,q),(t_1-t_2,q)\), as in the divisor summation on printed
p.43, yields
\[
 \begin{split}
 \sum_{t_1,t_2}|\Sigma(t_1,t_2)|
 \ll_\varepsilon{}&
 q^{5/2+\varepsilon}
 (\rho,q)^{1/2}(\sigma,q)^{1/2}T^2
 \left(Hq^{-5/6}+1\right)\\
 &{}+q^{2+\varepsilon}(\rho,q)(\sigma,q)HT.
 \end{split}
 \tag{18}
\]
Here the last term is the diagonal estimate from (4).
The divisor summation remains the one in the source: for fixed divisor
conditions, a residue class in an interval of length \(O_c(T)\)
contains \(O_c(1+T/d)\) integers. Also
\[
 ((t_1-t_2)^3,q)^{1/6}\le(t_1-t_2,q)^{1/2},
\]
and
\[
 ((t_1-t_2)^3,q)^{1/6}
 ((t_1-t_2)^2,q)^{1/2}
 \le q^{1/6}(t_1-t_2,q),
\]
which give its two off-diagonal terms.

For clarity, the extended frequency range introduces no additional
correlation assumption. Treat every pair \(t_1\equiv t_2\pmod q\) as
diagonal. Each \(t_1\) has at most \(c+1\) such \(t_2\) in the whole
range; the kernels coincide, giving \(O_c(T)\) diagonal pairs.
For the other pairs, \(|t_1-t_2|\ll_c T\), and the preceding divisor
count retains the same powers. The smooth parts extend only to \(cq\).
Finally, (16) is valid for every \(H\), including \(H>q\); its
\(k=0\) term has the factor \(H/q\) already present in (17).

Combining (12), (14), (15) and (18), and absorbing the number of dyadic
ranges and the divisor-square sum into \(X^\varepsilon\), reduces the
estimate to elementary envelopes. With
\(F_r=F_I(R_0|\rho|)\), \(F_s=F_J(S_0|\sigma|)\),
\[
 R_0(\rho,q)^{1/2}F_r\ll q,\qquad
 R_0^{1/2}(\rho,q)^{1/4}F_r\ll(Iq)^{1/2},
 \tag{19}
\]
and similarly for \(s\). Also
\[
 TF_K(T)\le q,\qquad T^{1/2}F_K(T)\le(Kq)^{1/2}.
\]
The first off-diagonal term in (18) gives \(q^{5/6}\);
the second gives \(q^{1/4}(IJ)^{1/2}\); and its diagonal term gives
\(q^{1/2}K^{1/2}\). This proves (9), the transferred form of
Heath-Brown's Lemma 5.

### The small-variable estimate

Complete the second and third variables modulo \(Q\), retaining
\(u\in I_1\), \((u,q)=1\). The component modulo \(q\) is
\[
 K_1(c^{-1}s,au^{-1}c^{-1}t;q),\qquad
 K_1(A,B;q)=\sum_{x\pmod q}^{*}e_q(Ax+Bx^{-1}),
\]
multiplied by the two separate local factors \(G_2,G_3\).
If \(s\) or \(t\) is zero modulo \(q\), the Kloosterman sum is a
Ramanujan sum independent of \(a\). Remove these terms into the main
term. For the others use
\[
 |K_1(A,B;q)|\le d(q)q^{1/2}(A,B,q)^{1/2}.
\]
The local factors and \(\omega_1(u)\) have bounded absolute value.
The Fourier envelopes leave a sum bounded by
\[
 C_cq^2\sum_{d\mid q}d^{1/2}
 \left(\sum_{\substack{1\le s\le cq\\d\mid s}}\frac1s\right)^2
 \ll_{c,\varepsilon}q^{2+\varepsilon}.
\]
Together with the normalization \(Q^{-2}\), the factor \(q^{1/2}d(q)\),
and the \(I\) possible values of \(u\), this proves (10). This is
the argument of source Section 5, printed p.44, without any assumption
on weighted Kloosterman sums.

### The long-variable estimate

Complete only the third variable modulo \(Q\), and remove \(h\equiv0
\pmod q\) into an \(a\)-independent main term \(M\). Then
\[
 \begin{split}
 S(a)-M=\frac1Q
 \sum_{\substack{h\pmod Q\\h\not\equiv0\ (q)}}
 &G_3(q^{-1}h)W_3(h)\\
 {}\times&
 \sum_{\substack{u\in I_1,v\in I_2\\(uv,q)=1}}
 \omega_1(u)\omega_2(v)e_q(ac^{-1}hu^{-1}v^{-1}).
 \end{split}
 \tag{20}
\]
Put \(k=ac^{-1}hu^{-1}\pmod q\). After triangle inequality over
\(u,h\), let \(n(k)\) be the sum of the resulting positive
coefficients at \(k\), and put
\[
 B(k)=\sum_{\substack{v\in I_2\\(v,q)=1}}
 \omega_2(v)e_q(kv^{-1}).
\]
Hölder's inequality gives
\[
 |S(a)-M|\ll_c q^{-1}
 \left(\sum_kn(k)\right)^{1/2}
 \left(\sum_kn(k)^2\right)^{1/4}
 \left(\sum_k|B(k)|^4\right)^{1/4}.
 \tag{21}
\]
The Fourier envelope and the bound for \(G_3\) imply
\[
 \sum_kn(k)\ll_c qI\log(2q),\qquad
 \sum_kn(k)^2\ll_{c,\varepsilon}X^\varepsilon qI(K+U).
 \tag{22}
\]
To check the second estimate, split \(h\) into dyadic ranges of size
\(T\). Collisions of the residues \(k\) require
\(h_1u_2\equiv h_2u_1\pmod q\). For fixed \(h_2,u_1\), the nonzero
integer \(h_1u_2\) lies in a range of length \(O_c(TU)\), in one
residue class modulo \(q\), and has \(O_\varepsilon(X^\varepsilon)\)
factorizations. Thus the collision count is
\[
 O_{c,\varepsilon}\bigl(
 X^\varepsilon(TI+T^2UI/q)\bigr).
\]
Weighting by \(\min(K,q/T)^2\) and summing the ranges proves (22),
as in source equations (6.3)–(6.4), printed pp.45–46.

Orthogonality and \(|\omega_2|\le1\) give
\[
 \sum_{k\pmod q}|B(k)|^4\le q\,\mathcal H(I_2;q),
\]
where \(\mathcal H(I_2;q)\) counts unweighted quadruples of units in
\(I_2\) satisfying
\[
 v_1^{-1}+v_2^{-1}\equiv v_3^{-1}+v_4^{-1}\pmod q.
\]
Every coefficient in the expanded fourth moment has absolute value at
most one. The elementary estimate proved in source equations
(6.6)–(6.13), printed pp.46–47, is
\[
 \mathcal H(I_2;q)\ll_\varepsilon X^\varepsilon J^2
 \left(1+\frac{V^2}q+\frac{V^{3/2}}{q^{1/2}}\right).
\]
It concerns this unweighted energy, so no periodic-weight extension is
needed. Substituting into (21)–(22) proves (11), the transferred form
of Heath-Brown's Lemma 7.

## Choosing a positive exponent

Order the lower interval scales as \(S\le M\le L\), with \(SML=X\).
Take
\[
 \eta=\frac1{100},\qquad \delta=\frac1{10000}.
\]
All finitely many logarithmic and divisor losses above can be absorbed
into \(X^\varepsilon\), with \(\varepsilon=1/1000\).
Since \(I_i\le2N_i\), constants from the numbers of interval integers
are harmless.

If \(S\le X^{1/4-\eta}\), apply (10) to the shortest variable. Then
\[
 qE_{\mathrm{small}}\ll
 X^\varepsilon q^{3/2}S
 \le X^{1-\eta+3\delta/2+\varepsilon}
 \le X^{1-\delta}.
\]

If \(L\le X^{1/2-\eta}\), apply (9) with \(K\) attached to the largest
variable. After multiplication by \(q\), its three terms are bounded
by \(X\) to the respective exponents
\[
 \frac{11}6\left(\frac12+\delta\right),\qquad
 \frac54\left(\frac12+\delta\right)+\frac13,\qquad
 \frac32\left(\frac12+\delta\right)
       +\frac12\left(\frac12-\eta\right).
\]
The middle bound uses \(L\ge X^{1/3}\).
These exponents are \(0.91685\), less than \(0.95846\), and
\(0.99515\), respectively. Adding \(\varepsilon=0.001\) leaves each
below \(1-\delta=0.9999\).

In the remaining case,
\[
 S>X^{0.24},\qquad L>X^{0.49},\qquad
 SM<X^{0.51},\quad M<X^{0.27},\quad L<X^{0.52}.
\]
Use (11) with \(U=S,V=M,W=L\). Since \(K+U\ll L\), expanding its
last factor bounds the three terms by
\[
 X^{1/4}S^{1/2}M^{1/4},\qquad
 X^{1/4}S^{1/2}M^{3/4}q^{-1/4},\qquad
 X^{1/4}S^{1/2}M^{5/8}q^{-1/8}.
\]
The first is at most \(X^{0.44125}\): use \(S\le M\) and
\(SM\le X^{0.51}\). The other two use
\[
 S^{1/2}M^{3/4}\le X^{0.3225},\qquad
 S^{1/2}M^{5/8}\le X^{0.28875},
\]
which follow from \(S\ge X^{0.24}\) and \(SM\le X^{0.51}\).
After multiplying by \(q\), the three exponents are therefore at most
\[
 0.44125+(0.5+\delta)=0.94135,
\]
\[
 0.5725+\tfrac34(0.5+\delta)=0.947575,
 \qquad
 0.53875+\tfrac78(0.5+\delta)=0.9763375.
\]
Again all remain below \(0.9999\) after the loss \(X^{0.001}\).

Thus one of (9)–(11), followed by (8), proves (1). Bounded \(X\) is
covered by increasing the absolute constant, and \(q=1\) has zero
discrepancy. Singleton padding gives the assertions for fewer factors.
Tensor-product partial summation, using the estimate uniformly in
all interval endpoints, proves (2). \(\square\)

## Evidence and scope

The primary scan has SHA256
\[
 \texttt{7102515f18e1498eb69b5c8dcc0dc094c5c94d73b4a710ec84fefb94d91e9fd0}.
\]
It contains 15 PDF spread-pages. Printed pp.41–43 are PDF pages 7–8;
printed pp.44–47 are PDF pages 9–10. The complete-correlation statement,
the three box arguments, their larger fixed frequency ranges, and the
weighted fourth-moment step were checked against those primary pages.

The theorem and its transfer argument have received internal mathematical
review, including a fresh review of the standalone proof. This is not
external review or formal verification. Its mathematical inputs are the
specific primary estimates identified above; a finite numerical experiment
is not used as evidence for the analytic bound.

