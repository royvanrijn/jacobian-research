# A prime plus a sum of two squares for every target under ZF

Let \(b\) indicate positive sums of two squares, \(K\) be the
Landau–Ramanujan constant, and \(C(T)\) be the positive product
\[
 C(T)=K\prod_{\substack{q\equiv3\ (4)\\q\nmid T}}
                 \left(1-\frac{1}{q(q-1)}\right)
       \prod_{\substack{q\equiv3\ (4)\\q\mid T}}(1+1/q).
\]
The starting estimates are the target-grid major-arc bound (PN.5)
under the fixed Dirichlet zero-free input ZF, and the uniform
auxiliary-weight removal (PN.7).
[Major arcs, weight removal, and the model](PRIME_PLUS_GAUSSIAN_NORM.md)
· [Scope of ZF](PRIMITIVE_NORM_BINARY_IMPORTED_INPUTS.md)

**Theorem ReflectedPrimeNormEveryTarget.** Under ZF, as \(T\) tends to
infinity through the positive integers,
\[
 \sum_{p<T}b(T-p)=C(T)\frac{T}{(\log T)^{3/2}}(1+o(1)).
\]
The relative error is uniform on every sufficiently large dyadic block
of targets. The proof uses the reflected Type I calculation below and
the reflected minor-arc comparison (PN.8). It uses the arithmetic and
operator estimates in Jiang–Lau–Lü for the specified purposes below.
[Jiang–Lau–Lü, v2, Sections 2–6](https://arxiv.org/abs/2611.12467v2)

## 1. Reflected Type I terms

Fix \(X<T\le 2X\), \(L=\log X\), \(U_0=X^{1/5}\),
\(Z=\exp(L^{0.24})\), and \(B_0=\exp(\sqrt L)\).
Use a family of \(k_0\) banks of split primes from Jiang–Lau–Lü,
Section 2, deleting all factors of \(T\). Write \(V_\nu\) for
the sum of \(1/q\) in a trimmed bank, \(\omega_\nu(n)\) for the
number of its members dividing \(n\), and put
\(W=\prod_\nu\omega_\nu/V_\nu\), \(f=bW\). The bank masses are
bounded above and below uniformly for fixed \(k_0\), since the
deleted mass is \(O(L\exp(-L^{0.1}))\).
Use a fixed smooth \(\Phi\) supported in \((0,1)\) as a function of
the prime variable divided by \(T\).

Apply Vaughan's identity to \(\Lambda(m)\), initially restricted to
\(Z\)-rough \(m\). In the notation of Jiang–Lau–Lü (2.8), let \(I_1\)
and \(I_2\) denote the terms with outer variables \(d\le U_0\),
and \(d,\ell\le U_0\), respectively, with \(b(m-1)W(m-1)\)
replaced by \(f(T-m)\) and cutoff \(\Phi(m/T)\). Thus \(m=dj\)
in \(I_1\), \(m=d\ell j\) in \(I_2\), the coefficients are
\(\mu(d)\log j,\mu(d)\Lambda(\ell)\), and all factors of \(m\)
are \(Z\)-rough. The excluded prime powers have total contribution
\(O_{A,k_0,\Phi}(XL^{-A})\) for each fixed \(A\), by the divisor
bound. The Type II term is the same reflected replacement in
the term with \(d,j>U_0\).

**Theorem ReflectedPrimeNormTypeI.** With these banks and cutoffs,
for each fixed \(A>0\), uniformly in \(T\),
\[
 I_1-I_2=\frac{1}{V(Z)}
             \sum_{n\ge1} f(n)w_{Z,T}(n)\Phi(1-n/T)
          + O_{A,k_0,\Phi}(XL^{-A}),                         \tag{1}
\]
where
\[
 V(Z)=\prod_{q\le Z}(1-1/q),\quad
 w_{Z,T}(n)=
 \prod_{\substack{q\le Z\\q\nmid T}}
       \left(1-\frac{\mathbf1_{q\nmid n}}{q-1}\right)
 \prod_{\substack{q\le Z\\q\mid T}}\mathbf1_{q\nmid n}.
\]
The weighted mean on the right of (1) equals
\[
 \frac{C(T)T}{\sqrt{\log T}}\int_0^1\Phi(u)\,du
      +o_{k_0,\Phi}\!\left(\frac{C(T)T}{\sqrt{\log T}}\right).
                                                               \tag{2}
\]
The main term includes a bank factor \(1+O_{k_0}(\exp(-cL^{0.1}))\),
which is absorbed in the remainder. This Type I theorem uses
classical distribution results.

### Distribution after trimming

For a function \(F\) write \(\Delta_F(y;s,a)\) for its sum on the
reduced class \(a\bmod s\) up to \(y\), minus the mean over all
reduced classes. With all bank parameters fixed at \(X\), uniformly
for \(X^{0.98}\le y\le 2X\),
\[
 \sum_{\substack{s\le X^{0.41}\\4\nmid s}}
       \max_{(a,s)=1}|\Delta_f(y;s,a)|\ll_{A,k_0} yL^{-A}.
                                                               \tag{3}
\]
First take the untrimmed banks, or any subset of them. Insert
parameters \(t_\nu\) with \(|t_\nu-1|=1/2\) into
\(b(n)\prod_\nu t_\nu^{\omega_\nu(n)}\). The value at primes is
constant on every block of the untrimmed bank partition in
each odd residue class modulo four. The endpoint ratios of that
partition satisfy Fouvry–Tenenbaum's spacing condition, and the
multiplicative function is bounded by \(\tau(n)\). Theorem 1.8
of Fouvry–Tenenbaum therefore applies in \(\mathcal F(4,2)\),
as in Jiang–Lau–Lü, Lemma 3.1. For \(4\nmid s\) the reference
term in their discrepancy is the reduced-class mean. Cauchy's
formula for the first derivatives in the selected parameters
gives (3) for the untrimmed weight. This argument applies at any
length \(\ge X^{0.97}\) with the same level \(X^{0.41}\).
[Fouvry–Tenenbaum, Definition 1.1 and Theorem 1.8](https://arxiv.org/abs/2004.04766v4)

For trimming, expand the product of
\(\omega_\nu^{\rm untrim}(n)-\sum_{q\in\mathcal P_\nu,\ q\mid T}\mathbf1_{q\mid n}\).
A selected deleted prime in each of some banks gives one
divisibility \(a_0\mid n\), where \(a_0\) is their product. After
\(n=a_0 v\), the remaining summand is \(b(v)\) with the untrimmed
weights in the other banks. This uses the exact split-prime
invariance of \(b\) and the disjointness of the banks. If
\((a_0,s)>1\), both reduced-class sums vanish. Otherwise the
discrepancy is the one at \(y/a_0\) in the residue \(a/a_0\bmod s\).
Each \(a_0\le \exp(O_{k_0}(L^{0.2}))\) and there are
\(O_{k_0}(L^{k_0})\) choices. The preceding untrimmed bound,
with a larger log saving, proves (3). The trimmed prime sets
themselves need no interval partition in \(\mathcal F(4,2)\).

### Gcd reduction and the roughness sieve

In either Type I term denote the outer product by \(j_0\)
(\(d\) or \(d\ell\)). Before expanding the roughness of the
remaining variable, we can restrict to \((j_0,T)=1\).
Any discarded \(m\) is a multiple of a prime \(q\mid T\),
\(q>Z\). There are \(O(XL/Z+L)\) such integers. The total
absolute weight at a given \(m\) is \(\ll L\tau(m)^2\tau(T-m)\)
on \(m<T\). Cauchy–Schwarz and fixed divisor moments give
\(O_A(XL^{-A})\). The same argument can discard \((m,T)>1\)
from the rough Type II terms.

Now expand the roughness of the remaining variable using
\(\sum_{c\mid \prod_{q\le Z}q}\mu(c)\). The additional congruence is
\(n=T-m\equiv T\pmod{j_0c}\), where \((j_0,c)=1\).
Truncate at \(c\le B_0\). To bound the tail, an integer \(m\)
with a discarded squarefree divisor has one divisor
\(B_0<a\le ZB_0\) made of primes \(\le Z\), by multiplying
its prime factors one at a time. The count of such integers is
\[
 O\left(X\sum_{\substack{a>B_0\\a\mid\prod_{q\le Z}q}}\frac1a
            +ZB_0\right)
    \ll_A XL^{-A}.
\]
The last estimate has an arbitrarily large log saving by
Rankin tilting with \(1/\log Z\). There are at most \(\tau(m)\)
additional signed divisor terms even before cancellation.
Another divisor-moment Cauchy–Schwarz bound proves the
truncation error, uniformly for both Type I sums.

Write \(g=(c,T)\), \(c_0=c/g\), \(s=j_0c_0\), and
\[
 g_-=\prod_{\substack{q\mid g\\q\equiv3\ (4)}}q,\qquad h=g g_- .
\]
In a nonzero summand \(b(n)=1\), the condition \(g\mid n\)
already implies \(h\mid n\). The change \(n=hv\) is exact:
\[
 b(hv)W(hv)=b(v)W(v),\qquad
 v\equiv (T/g)g_-^{-1}\pmod{s}.                              \tag{4}
\]
The norm equality follows from the even inert valuations
of \(h\), and all bank primes avoid \(T\). The residue in
(4) is reduced, \(4\nmid s\), and \(h\le B_0^2\).
The only small even divisor of \(s\) comes from the
squarefree \(c_0\). All statements hold when \(g=1\).

For fixed \(g\), apply (3) after partial summation at lengths
\(\asymp_\Phi X/h\) to replace (4) by its reduced-class mean.
The derivatives of the weights in the scaled variable are
\(O_\Phi(Lh/X)\), since \(m/T\) and \(n/T\) are bounded away
from zero on the cutoff. The modulus \(s\le U_0^2 B_0<X^{0.41}\)
determines \(j_0,c_0\) by its rough and smooth factors. At
\(j_0=d\ell\) use \(\sum_{d\ell=j_0}|\mu(d)|\Lambda(\ell)
\le\log j_0\). Thus the summed error for each \(g\) is
\(O_{A,k_0,\Phi}((X/h)L^{-A})\) after adjusting \(A\).
The sum of \(1/h\) over the squarefree target divisors is
at most a fixed power of \(L\).

In the reduced-class main term, changing back to \(n\)
gives the average with density \(1/(\varphi(j_0)\varphi(c_0))\)
on \(g\mid n\), \((n,j_0 c_0)=1\). Replace
\(1/\varphi(j_0)\) by \(1/j_0\) and drop \((n,j_0)=1\).
The relative totient change is \(O(L/Z)\). For the
coprimality error use, for \(q\mid j_0\),
\[
 \sum_{\substack{n\le 2X\\h q\mid n}} f(n)
       \ll_{k_0}\frac{\tau(h)XL}{h q}.
\]
Indeed \(q>Z\) is coprime to \(h\) and \(\tau(n)\) bounds \(f(n)\)
up to a constant. The sums of \(\tau(h)/h\),
\(1/\varphi(c_0)\) and \(1/j_0\) cost powers of \(L\),
whereas \(\sum_{q\mid j_0}1/q\ll L/Z\). Extending the
\(c=c_0g\) sum in the main term to all squarefree small-prime
products also costs \(O_A(XL^{-A})\). Apply the same Rankin
tilt to \(\tau(h)/(h\varphi(c_0))\), whose local nonconstant
terms after tilting are \(O(1/q)\). These bounds apply in
both Type I terms.

Restore the outer factors divisible by target primes in the
separated main coefficients. Each such prime exceeds \(Z\) and
\(\sum_{j_0\le U_0^2,\ q\mid j_0}1/j_0\ll L/q\); the
logarithmic coefficient and divisor sums already retained cost
fixed powers of \(L\). Thus the restoration costs \(O_A(XL^{-A})\).

The resulting signed \(c\)-sum is exactly \(w_{Z,T}(n)\).
The other factors separate as in (3.6)–(3.9) of Jiang–Lau–Lü.
In detail, the first term multiplies the weight sum by
\(M_0\log(T-n)-M_1\), where
\[
 M_j=\sum_{\substack{d\le U_0\\d\ Z\text{-rough}}}
               \frac{\mu(d)(\log d)^j}{d}\quad(j=0,1).
\]
The second term multiplies it by
\(M_0\sum_{\ell\le U_0,\ \ell\ Z\text{-rough}}\Lambda(\ell)/\ell\).
The usual rough Möbius calculation gives
\(M_0=O_A(L^{-A})\) and \(M_1=-1/V(Z)+O_A(L^{-A})\).
For this calculation convolve \(\mu\) with the indicator of
\(Z\)-smooth integers, use the classical logarithmic savings
for the ordinary Möbius sums at lengths \(\ge\sqrt{U_0}\),
and apply Rankin to the remaining smooth divisors. This
proves (1).

### Recovering the constant

Here are the details for evaluating the signed weighted mean
uniformly. For \(q\ge3,\ q\le Z,\ q\nmid T\), extract the
constant \((q-2)/(q-1)\) from \(w_{Z,T}\), leaving
\(1+\mathbf1_{q\mid n}/(q-2)\). Put \(A_Z(T)\) for the
product of the extracted constants. At a prime \(q\mid T\)
the factor instead removes all \(q\)-divisible norms.
The factor at two imposes even \(n\) for odd \(T\), and odd
\(n\) for even \(T\).

After factoring out \(A_Z(T)\), the Dirichlet series of \(f w_{Z,T}\)
is \(B(s)H_T(s)\), where \(B(s)=\sum b(n)n^{-s}\) and \(H_T\)
is a finite Dirichlet polynomial. To see this, for \(q\nmid T\),
\(q\ge3\), the quotient of Euler factors is
\(1+q^{-s}/(q-2)\) at split primes and \(1+q^{-2s}/(q-2)\)
at inert primes. At \(q\mid T,\ q\ge3\), it is \(1-q^{-s}\)
or \(1-q^{-2s}\), respectively. At two the quotient is
\(2^{-s}\) or \(1-2^{-s}\). In each bank, expand the
\(\omega_\nu\) sum by selecting one prime. Its local
split-prime factor \(1+q^{-s}/(q-2)\) is replaced by
\((q-1)q^{-s}/(q-2)\), and the bank sum is divided by \(V_\nu\).
All selected primes are distinct and avoid \(T\).

Write \(H_T(s)=\sum_d h_T(d)d^{-s}\). With
\(\sigma=1-1/\log Z\), these formulas give
\(\sum_d |h_T(d)|/d^\sigma\ll_{k_0}(\log L)^{O(1)}\).
For the target factors use
\(\sum_{q\mid T}1/q\le\sum_{q\le L^2}1/q+O(L^{-1})\),
and \(q^{1/\log Z}\le e\) for all retained primes. The
unselected non-target factors contribute a bounded product,
and the selected sums have bounded reciprocal mass.
The same formulas give \(H_T(1)\gg_{k_0}(\log L)^{-O(1)}\).
The moment with a factor \(\log d\) at \(s=1\) costs at
most an extra \(O(\log Z)\).

Convolve the Landau estimate
\(\sum_{n\le u}b(n)=Ku/\sqrt{\log u}+O(u/(\log u)^{3/2})\)
with \(h_T\). For \(d\le\sqrt X\), replacement of the
logarithm of \(u/d\) by that of \(u\) costs
\(O((1+\log d)/L)\) relatively on \(u\asymp_\Phi X\).
The remaining \(d\) have an exponentially small reciprocal
tail by the \(\sigma\)-bound. Partial summation gives the
leading term of (1) as \(K T A_Z(T)H_T(1)/(\sqrt{\log T}V(Z))\)
times \(\int_0^1\Phi\), with an error on the scale \(C(T)T/\sqrt{\log T}\) bounded by \(O_{k_0,\Phi}((\log L)^{O(1)}L^{-0.76})\).
This is the same Landau estimate used in Jiang–Lau–Lü (3.15)
and the following paragraph.

At \(s=1\), divide the mean factor at each prime by \(1-1/q\).
The norm's local probability of being \(q\)-divisible is
\(1/q\) at a split prime and \(1/q^2\) at an inert prime.
The normalized inert factors are therefore
\(1-1/[q(q-1)]\) for \(q\nmid T\), and \(1+1/q\) for \(q\mid T\).
The split and dyadic factors are one. At a selected bank
prime \(q\) the prime-facing mean of \(\mathbf1_{q\mid n}\)
is \(1/(q-1)\), so the bank contributes
\(\sum_{q\in\mathcal P_\nu}1/(q-1)\,/V_\nu
=1+O(\exp(-cL^{0.1}))\).
Finally \(\sum_{q>Z,\ q\mid T}1/q\ll L/(Z\log Z)\);
the non-target tails converge absolutely. This recovers (2).

## 2. Averaging determinant-\(T\) frames in the physical box

If \(P=(u,v)^{\mathsf t}\) is primitive and \((uv,T)=1\),
let \(g_T(P)\) be the integral matrix
\[
 g_T(P)=
 \begin{pmatrix}u&c\\ v&d\end{pmatrix},\qquad
 vc\equiv-T\pmod u,\quad 0\le c<u,\quad d=(T+vc)/u .
\]
Its columns are a basis of the index-\(T\) lattice
\(\Lambda_\lambda=\{(u',v'):v'\equiv\lambda u'\pmod T\}\)
where \(\lambda=vu^{-1}\bmod T\). The first column is \(P\),
and \(r=c/u\) is the normalized complement parameter.
Define the corresponding real frame by
\[
 g_T(x,y,r)=
 \begin{pmatrix} Ux&Ux r\\ Vy&Vy r+T/(Ux)\end{pmatrix},
 \qquad x=u/U,\quad y=v/V .
\]
Let \(S\) be squarefree, \((S,T)=1\), and define
\(G_{S,T}=\{g\bmod S:\det g=T\bmod S\}\).

**Lemma DeterminantTargetFrameAveraging.** Suppose \(U,V\ge X^\delta\),
\(U,V\le X^C\), \(UV\ge T\), \(X<T\le2X\), and \(S=X^{o(1)}\).
Let \(\Omega\) be a region in \([1,16]^2\times[0,1)\) whose boundary
is contained in \(N=X^{o(1)}\) real algebraic hypersurfaces of
bounded degree. For any function \(0\le F\le1\) on \(G_{S,T}\),
\[
 \frac1{UV}\sum_{\substack{(u,v)=1\\(uv,T)=1}}
     \mathbf1_\Omega(u/U,v/V,c/u)F(g_T(P)\bmod S)
 = \rho(T)|\Omega|\frac1{|G_{S,T}|}\sum_{g\in G_{S,T}}F(g)
        + O(X^{-c_\delta+o(1)}),                             \tag{5}
\]
with some \(c_\delta>0\), uniformly in \(\Omega,F,T\) at these
complexity bounds. Here
\[
 \rho(T)=\frac1{\zeta(2)}
                 \prod_{p\mid T}\frac{1-1/p}{1+1/p}.
\]
In the main term the kernels of the first row modulo the prime
divisors of \(S\) are independent uniform projective lines.

First take \(U\le V\), a rectangular grid cell \(I_1\times I_2
\times I_3\), and one residue matrix \(\bigl(\begin{smallmatrix}
 u_0&c_0\\v_0&d_0
 \end{smallmatrix}\bigr)\) in \(G_{S,T}\). For fixed
\(u\in UI_1\), \(u\equiv u_0\bmod S\), \((u,T)=1\), impose
\[
 c\in uI_3,\quad v\in VI_2,\quad
 c\equiv c_0,\ v\equiv v_0\pmod S,\quad
 cv\equiv ud_0-T\pmod{uS},\quad (v,T)=1.                  \tag{6}
\]
Use \(c<u\) at the end of \(I_3\). These are exactly the
canonical-frame conditions, since the congruence mod \(u\)
forces \((u,v)=1\). Write \(u=d_S q\) with \(d_S\) the full
\(S\)-primary part of \(u\), and \(B=d_SS\). As in
Jiang–Lau–Lü, Lemma 6.3, there are exactly \(d_S\) admissible
pairs of classes \(c_*,v_*\bmod B\). At \(p\mid S,\ p\mid u\),
the classes \(c_0,v_0\bmod p\) are units, since
\(c_0 v_0\equiv -T\bmod p\). Each lift of the first
fixes a unique lift of the second. At \(p\mid S,\ p\nmid u\)
the pair is already prescribed.

Use Möbius inversion for \((v,T)=1\), setting \(v=a v'\)
for \(a\mid T\), with sign \(\mu(a)\). Both \(a\) and \(B\)
are units modulo \(q\). For each of the \(d_S\) pairs, the
remaining congruence reduces, by parametrizing \(c,v'\) in their
classes modulo \(B\) and translating, to a modular hyperbola with
unit right-hand side modulo \(q\). The interval lengths are
\(u|I_3|/B\) and \(V|I_2|/(aB)\). The rectangle estimate for
\(xy\equiv k\bmod q\), \((k,q)=1\), is its area times
\(\varphi(q)/q^2\) plus \(O_\epsilon(q^{1/2+\epsilon}
(1+|I|/q)(1+|J|/q))\), including integer endpoints. Use
Boca–Zaharescu, Proposition 4.4 with constant test function,
and elementary counting for \(q=1\).
[Boca–Zaharescu, Proposition 4.4](https://arxiv.org/abs/math/0404112)

The main contribution at \(u\) is
\((\varphi(T)/T)V|I_2||I_3|\varphi(q)/(S^2q)\).
The error at \(u\) is
\[
 O_\epsilon\!\left(\tau(T)d_S q^{1/2+\epsilon}
                      (1+V/(uS))\right).
\]
With \(\sum_{u\ll U}d_S(u)^{1/2}\ll US^2\), the error summed
over \(u\) is \(O_\epsilon(\tau(T)VU^{1/2+\epsilon}S^2)\)
when \(U\le V\). Here \(\tau(T)=X^{o(1)}\) uniformly.
For the remaining main sum expand
\(\varphi(q)/q=\sum_{r'\mid u,(r',S)=1}\mu(r')/r'\),
retaining only \((r',T)=1\). Möbius inversion on \(u\)
and CRT in the class \(u_0\bmod S\) give
\[
 \sum_{\substack{u\in UI_1,\ u\equiv u_0\bmod S\\(u,T)=1}}
      \frac{\varphi(q)}q
 =\frac{U|I_1|}{S}\frac{\varphi(T)}T
      \prod_{p\nmid ST}(1-p^{-2})+O(\tau(T)\log(2U)+1).
\]
For the error, sum the unit endpoint errors against \(1/r'\)
and bound the convergent \(r'>16U\) main tail. Since
\(|G_{S,T}|=S^3\prod_{p\mid S}(1-p^{-2})\), the cube main term
is \(\rho(T)UV|I_1||I_2||I_3|/|G_{S,T}|\).
Dividing the error by \(UV\) gives a fixed power saving, for
example \(X^{-\delta/4+o(1)}\).

If \(U>V\), exchange the physical coordinates first. The canonical
frame at \((v,u)\) is
\[
 \begin{pmatrix}v&jv-d\\ u&ju-c\end{pmatrix},
 \qquad r'=j-r-T/(uv),\qquad j=\lceil r+T/(uv)\rceil .
\]
There are at most two values of \(j\), since \(T/(uv)\le1\).
On each piece the transformation on matrix residues is an
invertible fixed linear change, and \((x,y,r)\mapsto (y,x,r')\)
has absolute Jacobian one. The pieces and transformed regions
are bounded by \(O(N+1)\) hypersurfaces of bounded degree.
Indeed one substitutes \(r=j-r'-T/(UVxy)\) and clears the
positive denominator \(xy\). This handles the frame change
exactly, including its correction \(T/(uv)\).

In either orientation, approximate the region by inner and
outer unions of cubes of side \(\varepsilon\), as in
Jiang–Lau–Lü, Lemma 6.2. A bounded-degree hypersurface
meets \(O(\varepsilon^{-2})\) cubes in this fixed box by
the tube estimate used there, so the volume error is
\(O((N+1)\varepsilon)\). The summed rectangle error is
\(O(\varepsilon^{-3}|G_{S,T}|X^{-\delta/4+o(1)})\).
Taking \(\varepsilon\asymp X^{-\delta/16}\) proves (5).
Since \(T\) is a unit modulo \(S\), the determinant-\(T\) residue
matrices are obtained from \(\mathrm{SL}_2(\mathbb Z/S\mathbb Z)\)
by scaling the second row by \(T\). The first-row kernels are
unchanged, and transitivity and CRT give the asserted product law.

## 3. The reflected high moment and (PN.8)

Use the \(\alpha_m,\beta_n\) of (PN.4), with dyadic scales
\(M,R\ge X^\delta\), \(MR\asymp X\). In the minor-arc
comparison both sequences can be arbitrary bounded rough
sequences with the same unit support and fixed logarithmic
envelopes. Take \(\eta\) and \(\psi\) as in Jiang–Lau–Lü
Section 5, in particular \(\operatorname{supp}\psi\subset[-5,5]\).
Use the arcs of (PN.4), with \(B=A_0\). Lemma 6.8 of
Jiang–Lau–Lü also holds with this radius. In its Fourier
inversion truncate at \(|s|\le L^{A_0}/2\) and take the
Dirichlet ceiling \(\lceil4Y/L^{A_0}\rceil\). A denominator
at most \(L^{A_0}\) would put \(\theta\) on an arc of radius
\(L^{A_0}/Y\); the upper denominator bound changes only by a
constant in the large-sieve estimate.

**Proposition ReflectedPrimeNormMinorArcs.** For each fixed \(B>0\)
and sufficiently many fixed banks, with \(A_0\) chosen as in
the moment argument,
\[
       |Q_Y(T)-Q_Y^{\mathfrak M}(T)|\ll_B XY L^{-B}.
                                                               \tag{PN.8}
\]
The estimate is uniform in \(T\asymp X\) and in the dyadic
parameters of (PN.4). It also holds for the signed-quotient
version with \(k\in\mathbb Z\) and the smooth cutoff removed.

### Integral reconstruction and the choice of state space

For the signed-quotient version, put \(t=(cmn-ars)/T\) when
the numerator is divisible by \(T\). For \((a,c)=1\), equality of
the signed integral quotients is equivalent to \(t=c-a\).
For the major expression, averaging over \(j\bmod T\) in
(PN.4) imposes \(t\in\mathbb Z\), and integration in \(\xi\)
gives \(\psi(t/Y)\). The remaining factor is the integral of
\(e(\theta(t-c+a))\) over \(\mathfrak M\). Thus the difference
kernel is the minor-arc kernel of Jiang–Lau–Lü (5.10) with
the integer \(t\) just defined. If \(0<mn,rs<T\) on the
cutoff, the signed quotient is automatically positive.

Adjoin the ordered shared labels with weights \(1/(V_\nu p)\)
as in Jiang–Lau–Lü (5.5)–(5.6). Denote their product by \(D\).
For a dyad \(d_0\le D<2d_0\), put \(H=d_0Y\), \(U=HM\),
\(V=R\). The physical endpoints are \(P=(Dcm,s)\),
\(Q=(Dar,n)\), both contained in \([U,16U]\times[V,2V]\).
They satisfy
\[
 \det(P,Q)=TDt,\qquad |t|\le 5Y.
\]
For squarefree \(Dac\) the nonzero terms have \(t\ne0\);
otherwise \(a\mid mn\) contradicts roughness. A possible
ambient common divisor at either endpoint is \(Z\)-rough,
coprime to \(T\), and divides \(t\). It is therefore one.
The endpoints are units in both coordinates modulo \(T\) and
lie in the same \(\Lambda_\lambda\).

Use the disjoint union over \(\lambda\in(\mathbb Z/T\mathbb Z)^\times\)
of the labelled state spaces in these lattices. In each lattice
use points primitive **in that lattice**, with the same physical
box and the same divisibility by bank labels as in the source
argument. Intermediate points in an operator power are allowed
to have nonunit coordinates modulo \(T\). Along edges keep the
lattice fixed and use the intrinsic determinant \(\det(P,Q)/T\).
The intrinsic complement parameter of \(P\) is the ratio of
first coordinates of a lattice complement to \(P\), modulo one.
Use this parameter in the separation tests of Jiang–Lau–Lü (5.11).

This state space includes every arithmetic pair above and introduces
no extra term with nonzero input vectors. Indeed the input is
zero unless both coordinates are units, and then the lattice is
unique. With complex coefficients the fibre function is
\(\overline{\alpha_{P_1/\Pi(L)}}\beta_{P_2}\);
its conjugate at \(P\) times the value at \(Q\) gives the desired
crossed product. Roughness makes this function constant on each
full tuple fibre, since a nonzero value requires that the tuple
account for all auxiliary prime divisors of \(P_1\), each to the
first power. Thus the ordering factor and state-mass calculation
of Proposition 5.3 are unchanged. The two exceptions are
overlapping labels and failure of separation at an endpoint.

For clarity, the signed-quotient comparison also holds in polarized
form with distinct bounded rough coefficients at the two sides.
The absolute exceptional estimates and the fibre pairing apply to
the two crossed vectors. Equivalently, the minor expression is a
Hermitian form in the coefficient on the pair \((m,n)\), and its
polarization between two products follows by the same crossed
pairing. One can therefore insert
\(\Phi(mn/T)\overline{\Phi(rs/T)}\) by Mellin separation with a
finite \(L^1\) norm, without any loss in the minor-arc bound.
The signed and positive quotient versions coincide on this support.
The Mellin separation in this step needs no bound on the heights of
the twists, since the minor estimate uses only absolute envelopes.
The major-arc estimate (PN.5) already allows the smooth cutoff, and
its proof also applies without it to the signed-quotient version.

### Rooted moments

Use an even power \(r_0\asymp\sqrt L\). Roots with a nonzero
fibre function are primitive unit points, so only these roots are
needed on the diagonal in the fibre Hölder bound. We may sum
over all primitive unit roots, since the even diagonal moments
are nonnegative. For a fixed lattice and distinct intrinsic
primitive points in the positive box, either their slopes differ
by at least \(T/(16U)^2\), or the points coincide: primitive points
on the same positive ray in a lattice are unique. A path of
\(r_0\) steps changes its slope by \(O(r_0 T H/U^2)\).
Thus the support parameter in Jiang–Lau–Lü, Lemma 5.4, is
\(O(1+r_0 H)\), also for the matrix restricted to the unit roots.
The fibre norms are \(O(UV L^{O(1)})\) in squared norm. To see
this, for a fixed ordered tuple count first coordinates divisible
by its product and use the probability sums
\(\sum_{p\in\mathcal P_\nu}(V_\nu p)^{-1}=1\).
There is only one copy of any unit root. This is precisely the
norm and row-count input of Proposition 5.5.

At a unit root \(P_0\), expand the walk in the basis
\(g_T(P_0)\), writing its points \(g_T(P_0)z_j\) with
\(z_0=z_{r_0}=e_1\). Lattice primitivity is now ordinary
primitivity. All intrinsic determinants are \(\det(z_j,z_{j+1})\).
The slope bound gives \(|z_{j,2}|\ll r_0 H\); the first
coordinate gives \(|z_{j,1}|\ll1+r_0 H\).
For every bank prime \(p\), divisibility of \((g_T z_j)_1\)
is exactly the condition \([z_j]_p=\ker(\operatorname{row}_1
g_T\bmod p)\). This is valid at the intermediate points as
well, since \(p\nmid T\) and the \(z_j\) are primitive.

Now the coefficient of a fixed labelled path is independent of
the frame. The number of paths and the largest absolute
coefficient are \(\exp(O(L^{0.73}))\), using the same count as
(5.19)–(5.22). The product \(S\) of its distinct labels has
\(\log S=O(L^{0.70})\). Its box and separation conditions give
a region of complexity \(\exp(O(L^{0.21}))\) as in Lemma 6.2.
In fact the normalized coordinates in the frame are
\[
 x(z_1+r z_2),\qquad
 y(z_1+r z_2)+\frac{T}{UV}\frac{z_2}{x}.
\]
The box boundaries are polynomials of degree at most three
after clearing \(x\). The separation ratios depend only on the
first row, and their derivatives with respect to \(r\) on the
box are bounded by the reciprocal square of \(z_1+r z_2\).
Thus the separation-boundary count of Lemma 6.2 is unchanged.
The averaging (5), applied path by path, replaces the initial
point sum divided by \(UV\) by \(\rho(T)\) times the independent-line
continuous average with a power-saving error after summing over
all paths.

We give the scale check for this continuous model. Replace the
second physical coordinate by its value divided by \(T\), and
the second physical scale \(V\) by \(V/T\). The real frame now
has determinant one. In Jiang–Lau–Lü, §§6.2–6.4, the estimates
for a **fixed** real frame allow this smaller second scale.
The construction and adjoint identity of the recorded-label
operator (6.6)–(6.8) use only intrinsic determinants, projective
lines and the fixed bank masses. In Lemma 6.9, the slope intervals
in the rescaled box have width \(H/U^2\). The determinant bound
forces neighboring intervals, and the first-coordinate box
gives \(|h+Dlr_\circ|\le16\) in the index-\(D\) lattice.
The other constraint is \(|l|\ll H/D\asymp Y\), obtained from
the determinant. Consequently the same separated basis of
Lemma 6.7 and the same large-sieve rectangle of Lemma 6.8
apply at both endpoints. The second-coordinate box is an
endpoint projection. No lower bound for \(V/T\) enters.

The absolute row and column estimate of Lemma 6.10 likewise
counts \(O(1)\) values of \(h\) for fixed intrinsic determinant.
In Lemma 6.11 the fresh-label term uses the preceding common
rectangle and Parseval in the label products. For a returning
label the equality term fixes an intrinsic congruence on \(h\);
the separation test and \(|h+d_*r_\circ|\le16\) give its small
row or column sum. The corresponding major-arc term counts
primitive points in the strip \(|h+Dt r_\circ|\le16,\ |t|\le5Y\)
in an index-\(p\) lattice, again with the second-coordinate
box just a restriction. These are all the geometry inputs to
that norm proof. The trimmed banks retain their lower masses
and progression gaps; the measures \(\pi_\nu(p)=1/(V_\nu(p+1))\)
still obey \(\sum\pi_\nu\le1\) and the comparison with
\(1/(V_\nu p)\) in the separation test. Thus Lemma 6.11 and
the Fejér restriction at the end of Proposition 5.7 give
\(O_E(L^{-E r_0})\) for the model rooted moment, with \(k_0\)
large enough for \(E\). This is a bound for each frame.
The polynomial physical side lengths are used in (5) for
the **arithmetic** average.

The power-saving transfer error from (5) is smaller than
\(L^{-E r_0}\), for each fixed \(E\). Also \(\rho(T)\ll1\).
The rooted diagonal bound divided by \(UV\) is therefore
\(O_E(L^{-E r_0})\). Lemma 5.4, the slope support bound,
and \((1+r_0 H)^{1/r_0}\ll1\) give the necessary pairing
bound \(UV L^{-E+O(1)}\). The dyadic prefactor is \(1/d_0\),
so \(UV/d_0\asymp XY\), and the number of dyads is \(O(L)\).
This handles the separated squarefree-label contribution.

### Exceptions

For a squarefree edge with failed separation at its source, use
(5) for the one-edge family with the source separation indicator
replaced by its complement. The coordinate size and coefficient
counts are \(\exp(O(L^{0.21}))\); the transfer again has power
saving. In the model use the absolute kernel bound
\(\lambda_{a,c}(t)=\mathbf1_{t=c-a}
      +O(L^{3A_0}/Y)\mathbf1_{|t|\le5Y}\),
with its positive upper interpretation. Write \(z'=(h,Dt)\).
For fixed \(t,r\) the first-coordinate box allows \(O(1)\) choices
of \(h\). Each label has its independent line probability, so
the joint probability is at most \(1/(Dac)\). This also holds
for shared labels because they divide the intrinsic
determinant. The ordered-label sum with these weights is at
most one after the reciprocal-bank normalizers. The set of
\(r\) failing separation has measure \(O(\exp(-cL^{0.1}))\),
by the circle union bound and Markov argument of Lemma 6.13,
which applies to the trimmed sets. The contribution after
the dyadic prefactor is \(O_B(XYL^{-B})\). Reverse the edge
for a failure at its target.

Overlapping labels are bounded before the graph. For fixed
\(a,c\asymp Y\) their arithmetic term has absolute size
\(\ll XL^{O(1)}/Y\). One way to see this is to count
\(\tau(mn)^2\le\tau_4(mn)\) in the reduced progression
\(mn\equiv T\bmod a\), and likewise for \(rs\bmod c\).
In each four-factorization the product of the three factors
other than a largest one is \(O(X^{3/4})\). Counting the
largest factor in one residue class gives
\(O((X/a+X^{3/4})L^{O(1)})\); use Cauchy–Schwarz over the
common quotients. The major term also costs
\(XL^{O(1)}/Y\) per pair: the target congruence leaves \(O(1)\)
values of \(rs\) for each \(mn\) and conversely because \(a,c\)
are target units, and the arc measure is \(O(L^{3A_0}/Y)\).
The ordinary divisor second moment pays for multiplicities.

There are \(O(Y^2\sum_{p\ {\rm bank}}p^{-2})\) pairs \(a,c\)
sharing a prime. For fixed distinct \(a,c\), a repeat or a
collision with their labels among the \(O_{k_0}(L^{0.01})\)
inserted labels has probability
\(O_{k_0}(L^{0.02}/p_{\min})\) under the independent reciprocal
bank measures. Both bounds give an arbitrarily large log
saving times \(XY\). This supplies the other exception and
proves (PN.8).

## 4. Completing the pointwise asymptotic

Combine (PN.5) and (PN.8) under ZF. In the Type II term from
Vaughan's identity, split \(d,j>U_0\) into \(O(L^2)\)
boxes with \(M,R\ge X^{1/5}\) and \(MR\asymp_\Phi X\).
The coprimality deletion of Section 1 applies. Mellin
inversion of \(\Phi(dj/T)\) gives \(\mu(d)d^{iv}\) on a
rough interval and \(\beta_j=c_{U_0}(j)j^{iv}\) with
\(|\beta_j|\ll L\), as in Jiang–Lau–Lü, Corollary 4.3.
Heights \(|v|>L^C\) cost an arbitrary logarithmic saving
by smooth decay and the divisor-moment bound
\(\sum_m\tau(m)\tau(T-m)\ll XL^{O(1)}\) on \(1\le m<T\).
For the retained heights use the signed-quotient mean square.
Expanding \(f(T-dj)\) by the bank products uses
\(b(ak)=b(k)\) for \(a\) split and \(k>0\), with \(b=0\) for
\(k\le0\). In a single dyad \(a\asymp Y\), Cauchy–Schwarz
therefore costs at most \((X Q_Y(T)/Y)^{1/2}\).
The mean-square bounds, with sufficient saving and bank count,
pay for all dyads and Mellin costs. Thus the Type II term is
\(O_{A,k_0,\Phi}(XL^{-A})\) for each fixed \(A\).

By the Type I theorem the von Mangoldt sum with \(bW\) has
the smoothed main term \(C(T)T/\sqrt{\log T}\) times
\(\int_0^1\Phi\). All conclusions are uniform over the
\(J\asymp_{k_0}\log L\) disjoint families from Section 2 of
Jiang–Lau–Lü. Prime powers contribute \(X^{1/2+o(1)}\).
On the fixed bulk cutoff \(\log p=\log T+O_\Phi(1)\).
This gives the smoothed unweighted prime asymptotic.
Changing \(V_\nu\) to \(\sum_{p\in\mathcal P_\nu}1/(p-1)\)
changes each family by a scalar \(1+O(\exp(-cL^{0.1}))\).
Average the families and apply (PN.7). Its weighted-variance
proof controls \(\sum_{p<T}b(T-p)|\overline W(T-p)-1|\),
so the same removal estimate applies with the smooth cutoff.

Finally the positive sieve proof of (PN.2) works on any
subinterval \(I\subset(0,T)\), with the interval length in
the sieve main term and the same CRT error. Thus its
upper bound is \(O(C(T)|I|/L^{3/2}+X^{1/10})\), including
the omitted small primes. The two endpoint intervals of length
\(\varepsilon T\) can therefore be removed uniformly at
relative cost \(O(\varepsilon)+o(1)\). Use bulk cutoffs
equal to one off these intervals, take \(X\) to infinity,
and then let \(\varepsilon\) tend to zero. This proves the
every-target theorem with the stated \(1+o(1)\) and explicit
constant.
