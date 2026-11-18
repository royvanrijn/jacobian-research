# Power saving correlations for simultaneous norms

The raw Gaussian--Eisenstein coordinate weight has a pointwise main term outside
a power-sized exceptional set. The same argument treats the fully sifted weight
and a growing family of individual bad-prime tests. This gives the actual raw
lower count and bad-prime marginal interface almost everywhere, including
targets with arbitrarily high prime-power divisors.

The proof uses integer powers of four fixed Dirichlet \(L\)-functions. It
adapts the Fourier method in
[D12Power](PRESCRIBED_PRIME_POWER_SAVING_EXCEPTIONS.md), with a separate choice
of the major denominator cutoff for the unrestricted correlation. It does not
remove the exceptional set or prove the every-target binary conjecture.

## Statement

Let
\[
 h(n)=\frac{r_G(n)r_E(n)}{24},\qquad
 h_6(n)=h(n)\mathbf1_{(n,6)=1},\qquad
 C_{12}=\{c\ge1:p\mid c\Longrightarrow p\equiv1\pmod{12}\}.
\]
The representation counts use all ordered integer coordinates, and all value
summands below are positive. Put
\[
 h_C(n)=h_6(n)\mathbf1_{C_{12}}(n),\qquad
 v_p(n)=h_6(n)\mathbf1_{p\mid n}
\]
for forbidden primes \(p\ge5\), meaning \(p\not\equiv1\pmod{12}\).
The support of \(h_6\) makes \(v_p=h_6\mathbf1_{p^2\mid n}\).
For arithmetic functions \(u,v\), write
\[
 C_{u,v}(N)=\sum_{a+b=N}u(a)v(b).
\]
Fix the absolute constants
\[
 \rho=10^{-4},\qquad \kappa=\rho/10000,
 \qquad \varepsilon_*=\rho/100=10^{-6},
\]
and define the finite family
\[
 \mathcal U_X=\{h_6,h_C\}\cup
       \{v_p:5\le p\le X^\kappa,\ p\text{ forbidden}\}.
\]

**Theorem SimNormPower.** The local-model coefficients
\(\mathfrak S_{u,v}(N)\) defined below satisfy, uniformly for
\(u,v\in\mathcal U_X\),
\[
 \sum_{X<N\le2X}
 \left|C_{u,v}(N)-N\mathfrak S_{u,v}(N)\right|^2
 \ll X^{3-10\varepsilon_*}.
 \tag{1}
\]
Consequently, for every sufficiently large \(X\), one exceptional set
\(\mathcal E_X\subset(X,2X]\), of cardinality
\(O(X^{1-7\varepsilon_*})\), has the following property. For every integer
\(N\notin\mathcal E_X\), simultaneously for every pair in the family,
\[
 C_{u,v}(N)=N\mathfrak S_{u,v}(N)+O(X^{1-\varepsilon_*}).
 \tag{2}
\]
The constants are independent of the individual target, marked primes and
all target valuations.

For \(N\equiv2\pmod{12}\), let
\[
 \mathcal T(N)=
 \prod_{\substack{p\mid N\\p\equiv1\ (12)}}(1+9/p)
 \prod_{\substack{p\mid N\\p\ge5,\ p\not\equiv1\ (12)}}(1+1/p).
\]
The main coefficients obey
\[
 \mathfrak S_{h_6,h_6}(N)\asymp\mathcal T(N),
 \tag{3}
\]
\[
 2\mathfrak S_{v_p,h_6}(N)-\mathfrak S_{v_p,v_p}(N)
 =\beta_p(N)\mathfrak S_{h_6,h_6}(N),
 \tag{4}
\]
and
\[
 \frac{\mathfrak S_{h_C,h_C}(N)}{\mathfrak S_{h_6,h_6}(N)}
 =\Gamma(N):=\prod_{\substack{p\ge5\\p\text{ forbidden}}}
                      (1-\beta_p(N))
 \ge10-\pi^2>0.
 \tag{5}
\]
Here \(\beta_p(N)\) is the exact lifted local bad probability from the
[local-density theorem](SIMULTANEOUS_NORM_LOCAL_DENSITIES.md#the-local-theorem).

In particular, for admissible targets outside the same exceptional set, the
unrestricted coordinate counts satisfy
\[
 \begin{split}
 R(N)&=576C_{h_6,h_6}(N)\asymp N\mathcal T(N),\\
 \frac{B_p(N)}{R(N)}&=\beta_p(N)+O(X^{-\varepsilon_*})
                 \qquad(p\le X^\kappa),\\
 \frac{R_{\mathrm{good}}(N)}{R(N)}
       &=\Gamma(N)+O(X^{-\varepsilon_*}),
 \end{split}
 \tag{6}
\]
where \(R_{\mathrm{good}}=576C_{h_C,h_C}\). These are full unrestricted
counts. They are not ratios for a selected subpopulation in one short window.
No numerical onset or priority claim is attached to this theorem.

## Four character factorizations

Use the four Dirichlet characters \(\psi\) literally modulo \(12\),
including the principal character, and let \(a_4\) be their Dirichlet
convolution. The norm Euler coefficients are
\[
 h(p^j)=\begin{cases}
 (j+1)^2,&p\equiv1\pmod{12},\\
 \mathbf1_{2\mid j}(j+1),&p\equiv5,7\pmod{12},\\
 \mathbf1_{2\mid j},&p\equiv11\pmod{12}.
 \end{cases}
 \tag{7}
\]
For every \(u\in\mathcal U_X\), write
\[
 \sum_{n\ge1}\frac{u(n)}{n^s}
 =E_u(s)\prod_{\psi\ (12)}L(s,\psi),\qquad
 E_u(s)=\sum_{k\ge1}\frac{e_u(k)}{k^s}.
 \tag{8}
\]
With \(z=p^{-s}\), the local factors of \(E_u\) outside \(12\) are
\[
 \begin{array}{c|c}
 u&\text{local correction}\\\hline
 h_6&1-\chi_{12}(p)z^2\\
 h_C&1-z^2\ (p\text{ allowed}),\quad(1-z^2)^2\ (p\text{ forbidden})\\
 v_p&3z^2-z^4\ (p\equiv5,7\ (12)),\quad z^2-z^4\ (p\equiv11\ (12)),
 \end{array}
 \tag{9}
\]
where for \(v_p\) only its marked prime changes from the \(h_6\) row.
All factors at \(2,3\) are one. The marked-prime formula multiplies the
four-character local factor \((1-z^2)^{-2}\) by the norm factor with its
constant coefficient removed. It never divides by \(h(p)=0\).

These square-supported corrections give, for a fixed integer \(A\),
\[
 \sum_{k\ge1}\frac{|e_u(k)|}{k^{3/5}}\ll1,
 \qquad |e_u(k)|\le d_A(k),
 \tag{10}
\]
uniformly throughout the family. The convergent Euler product has factors
\(1+O(p^{-6/5})\); the finitely altered coefficients are bounded by three.
All family members are nonnegative and bounded by \(h_6\le d_4\).

The positive residue of the \(h_6\) series is
\[
 \rho_6=\frac{\log(2+\sqrt3)}{3\sqrt3}.
\]
The Euler factors and this normalization are also recorded in the
[weighted AP note](SIMULTANEOUS_NORM_WEIGHTED_AP_MEANS.md).

## Shifted means and nonreduced residue classes

The major arcs require polynomially uniform prefix means, not a short interval
inference from a logarithmic prefix error. We give the shifted raw-weight
argument, including bases for which \(h_6(B)=0\).

For \((B,6)=1\), define
\[
 j(B)=\prod_{\substack{p\text{ forbidden}\\v_p(B)\text{ odd}}}p.
\]
Then
\[
 h_6(Bn)=0\quad\text{unless }j(B)\mid n,
 \qquad
 h_6(Bn)=h_6\!\left(Bj(B)\frac n{j(B)}\right)
       \quad(j(B)\mid n).
 \tag{11}
\]
The stable base \(B'=Bj(B)\) has even bad-prime valuations and
\(h_6(B')>0\). For a good \(p^v\Vert B'\), put \(a=v+1\). The
local correction to the four-character product for
\(h_6(B'n)/h_6(B')\) is
\[
 (1-z)P_v(z),\qquad
 P_v(z)=1+\left(-2+\frac2a+\frac1{a^2}\right)z
            +\left(1-\frac2a+\frac1{a^2}\right)z^2.
 \tag{12}
\]
For \(p\equiv5,7\pmod{12}\), with even \(v\), it is
\(1+(1-v)z^2/(v+1)\). For \(p\equiv11\pmod{12}\) it is
\(1-z^2\). At unchanged primes use the \(h_6\) row of (9).
Their coefficients are uniformly bounded, their absolute \(3/5\)-norm
is \(B^{o(1)}\), and \(h_6(B')=B^{o(1)}\). If \((B,6)>1\), the
whole shifted sequence is zero instead.

For a character \(\chi\pmod D\) whose primitive character does not
factor through \(12\), all four \(\chi\psi\) are nonprincipal, of
period at most \(12D\). In a four-index dyadic box with product at most
\(Y\), sum the largest index last. The other three lower endpoints have
product at most \(Y^{3/4}\), and the final periodic sum is \(O(D)\).
The \(O((\log(2Y))^4)\) boxes, convolution with the correction kernel,
and (11) give
\[
 \sum_{n\le Y}h_6(Bn)\chi(n)
 \ll D(BY)^{o(1)}Y^{3/4}.
 \tag{13}
\]
If \(\chi(j(B))=0\), the left side is identically zero, so that case is
included.

Untwisted, peel the fixed principal character's \(\zeta\)-factor out of
(8). The remaining three-character sequence has cumulative bound
\(B^{o(1)}Y^{2/3}(\log(2Y))^3\) by the same largest-index argument.
Its absolute mass to \(Y\) is
\(O(B^{o(1)}Y(\log(2Y))^3)\). Hyperbola splitting at \(Y^{3/4}\)
gives
\[
 \sum_{n\le Y}h_6(Bn)=C_BY+O((BY)^{o(1)}Y^{3/4}),
 \tag{14}
\]
\[
 C_B=\rho_6\prod_{p^v\Vert B}
       \frac{\sum_{i\ge0}h(p^{v+i})p^{-i}}{D_p},
 \qquad D_p=\sum_{i\ge0}h(p^i)p^{-i}.
 \tag{15}
\]
Define \(C_B=0\) when \((B,6)>1\). For completeness, if \(A(t)\)
is the cumulative bound for the remaining sequence, the split at
\(R=Y^{3/4}\) has floor error \(O(B^{o(1)}R\log^3(2Y))\).
Partial summation bounds the main-series tail by
\(O(B^{o(1)}YR^{-1/3}\log^3(2Y))\). The other terms are
\(\sum_{\ell\le Y/R}(A(Y/\ell)-A(R))\), with the same bound.
This proves (14); (15) follows from residues and then (11).

An exceptional character has underlying primitive character \(\psi\)
factoring through \(12\). On supported terms its value is
\(\psi(B)^{-1}\), subject to coprimality with its modulus, because
\(Bn\equiv1\pmod{12}\). This phase is necessary when \(B\) is not
\(1\pmod{12}\). Inclusion-exclusion at divisors of the modulus removes
coprimality using (14). In particular a zero value \(h_6(B)\) does not
make \(h_6(Bn)\) identically zero when \((B,6)=1\).

Now split an arbitrary class \(c\pmod Q\) by \(g=(c,Q)\), writing
\(n=gn'\) with \(n'\) in a unit class modulo \(Q/g\). Character
orthogonality, (13), and exceptional inclusion-exclusion use bases at most
\(Q^2\). A deliberately loose budget of \(Q^4X^{3/4+o(1)}\) bounds
the error in each residue class, uniformly for \(0\le Y\le2X\).
The main coefficient is its exact local mass, described in the next section.

For \(h_C\), a base with a forbidden factor makes the shifted sequence
zero; allowed bases have correction (12), and the unchanged forbidden
factors are \((1-z^2)^2\). Thus the same polynomial budget applies.
For \(v_p\) at rational denominator \(D\le X^\rho\), take
\(Q=\operatorname{lcm}(D,p)\le X^{\rho+\kappa}<X^{2\rho}\), retain
only classes divisible by \(p\), and sum their additive phases. Even the
larger budget
\[
 Q^6X^{3/4+o(1)}\le X^{3/4+12\rho+o(1)}<X^{0.80}
 \tag{16}
\]
therefore covers every class sum and rational prefix for the whole family.
This polynomial budget is retained before summing classes; one does not sum
\(Q\) copies of a coarse \(X^{0.80}\) bound and call it unchanged.

## Local measures and rational coefficients

For \(p\ge5\), put \(M_p=(1-p^{-1})D_p\) and
\[
 W_p(t)=\frac{h(p^{v_p(t)})}{M_p}.
 \tag{17}
\]
Its value at zero is irrelevant. At \(p^j\Vert12\), use instead
\(W_p(t)=p^j\mathbf1_{t\equiv1\ (p^j)}\). Each has Haar mass one.
The carrier constant for \(h_6\) is \(\rho_6\).
For \(v_p\), retain that same carrier constant, and replace only its
marked density by \(W_p\mathbf1_{v_p(t)\ge2}\), without renormalizing.
The marked measure has mass less than one; its mass contributes even when
that prime does not divide the frequency denominator.

For \(h_C\), use the carrier constant
\[
 \rho_C=\rho_6\prod_{\substack{p\ge5\\p\text{ forbidden}}}D_p^{-1}>0,
 \tag{18}
\]
replace bad local densities by
\((1-p^{-1})^{-1}\mathbf1_{\mathbb Z_p^\times}\), and retain the
good and fixed-12 densities. The product in (18) converges absolutely because
\(D_p=1+O(p^{-2})\) at forbidden primes. Denote a family's carrier constant
by \(c_u\), and its densities by \(\Omega_{u,p}\).

These are the arithmetic local laws of the preceding progression calculation.
Indeed the raw tail mass at \(p^v\) is
\[
 p^{-v}\frac{\sum_{i\ge0}h(p^{v+i})p^{-i}}{D_p}
       =\frac{C_{p^v}}{\rho_6p^v}.
 \tag{19}
\]
Exact valuations are differences of successive tails. Unit residue classes
are uniform by the exceptional character sum. The fixed-12 support gives the
stated classes at two and three. These facts reconstruct all nonreduced
classes by CRT; the local product is not an additional independence hypothesis.

Choose the local additive characters so their product on an integer \(t\)
is \(e(at/D)\), where \(e(y)=\exp(2\pi iy)\). For a reduced rational,
define
\[
 B_u(a/D)=c_u\prod_p
       \int_{\mathbb Z_p}\Omega_{u,p}(t)e_{p,a/D}(t)\,dt.
 \tag{20}
\]
Formula (16) proves that these are the actual rational prefix coefficients.
They are uniformly \(O(1)\). For a frequency of conductor \(p^j\)
outside \(12\), the raw radial coefficient is
\[
 -\frac{f(j-1)p^{-j}}{M_p}
 +\frac{1-p^{-1}}{M_p}\sum_{i\ge j}f(i)p^{-i},
 \qquad f(i)=h(p^i).
 \tag{21}
\]
Use the masked \(f\) at a marked prime. For the mask, \(j=1,2\)
give \(O(p^{-2})\), hence \(O(p^{-j})\); for \(j\ge3\) the bound
is \(O((j+1)^2p^{-j})\). The other local laws have that same bound.
Multiplication gives
\[
 |B_u(a/D)|\ll D^{-1+o(1)},
 \tag{22}
\]
uniformly even in the marked prime, without a factor \(p^2\).
Only conductors dividing \(12\) survive at its two primes. Outside \(12\),
the coefficients are radial and independent of the primitive numerator.

## Long interval Fourier estimates

Let
\[
 F_u(\beta)=\sum_{n\le2X}u(n)e(n\beta).
\]
Choose one common measurable reduced Dirichlet approximant \(a/D\) for each
\(\beta\), for all family members, with
\[
 D\le T=\lfloor X^{0.98}\rfloor,\qquad
 |\beta-a/D|\le1/(DT),\qquad R=X^\rho.
 \tag{23}
\]
Differences on the circle use signed representatives. The major cutoff
\(R\) is separate from the maximum denominator \(T\).
For \(D<R\), partial summation of (16) gives
\[
 F_u(\beta)=B_u(a/D)J(\beta-a/D)+O(X^{0.83}),
 \qquad J(v)=\int_0^{2X}e(vt)\,dt.
 \tag{24}
\]
The prefix error is \(O(X^{0.80})\), and
\(1+|\beta-a/D|X\le1+X^{0.02}/D\). This proves (24) directly,
without substituting a short interval into a prefix logarithmic error.

For \(R\le D\le T\), we claim
\[
 F_u(\beta)\ll X^{1-\rho/8+o(1)}.
 \tag{25}
\]
Discard \(n\le X^{0.99}\), costing \(X^{0.99+o(1)}\) by the divisor
bound. Expand \(u=e_u*a_4\). For \(k>X^{0.01}\), (10) and the count
\(2X/k+1\) give total cost
\[
 X^{o(1)}\left\{X\,X^{-0.004}+X^{0.6}\right\},
\]
which is acceptable. For smaller \(k\), partition the output product into
dyadic strips \((P,2P]\), \(P\ge X^{0.99}\), then split all five
variables dyadically. The largest free character index exceeds \(X^{0.24}\).
All combined coefficients are bounded by fixed divisor functions uniformly
throughout the family.

If that largest index exceeds \(X^{0.7}\), its opposite combined index
\(U\) is at most \(X^{0.30}\). Split the fixed character into residue
classes modulo \(12\), giving geometric frequencies \(12\ell\beta\).
Reduction of \(12a/D\) changes the denominator by a bounded factor. If
\(U<X^{\rho/2}\), rational spacing gives geometric bound \(O(D)\);
the relative contribution is at most
\[
 X^{o(1)}UD/P\le X^{-0.01+\rho/2+o(1)}.
\]
For \(U\ge X^{\rho/2}\), Cauchy--Schwarz and spacing in blocks of
\(cD\), for a sufficiently small fixed \(c>0\), give the relative
squared bound
\[
 X^{o(1)}\{1/D+1/U+U/P+D/P\}.
 \tag{26}
\]
For example the needed squared geometric-sum estimate is
\((U/D+1)(V^2+DV)\), with \(V\asymp P/U\). The fixed rational
approximation constant after reduction of \(12\beta\) is handled by
choosing \(c\) small. Since \(D/P\le X^{-0.01}\), (26) saves
\(X^{-\rho/4+o(1)}\) before weakening the exponent.

Otherwise isolate the largest free index and combine the others. Both sides
are between \(X^{0.23}\) and \(X^{0.77}\); let \(M\) be the smaller.
Cauchy--Schwarz on the larger variable leaves differences \(O(M)\) and
geometric intervals of length \(O(P/M)\). For each difference the pair
coefficient mass is \(O(MX^{o(1)})\). Any remaining periodic factor from a
fixed character or its unit indicator is split into the twelve residue
classes. The geometric frequency is then \(12j\beta\); its reduced
denominator is comparable to \(D\), so the same spacing estimate applies.
One can alternatively enlarge a nonnegative squared sum to all residue classes.
The relative squared bound is
\[
 X^{o(1)}\{1/M+1/D+M\log X/P+D\log X/P\}.
 \tag{27}
\]
The first term includes the diagonal. The strip intersection and separate
dyadic restrictions leave genuine intervals for the geometric sum. Since
\(M\ll\sqrt P\), \(D/P\le X^{-0.01}\), and \(D\ge X^\rho\),
this gives a stronger saving. Summing the boxes and strips, whose output
scales sum to \(O(X)\), proves (25). All logarithmic and divisor losses
are \(X^{o(1)}\), with fixed positive exponent margins.

Define \(P_u(\beta)=B_u(a/D)J(\beta-a/D)\) when the common assigned
\(D<R\), and zero otherwise. Then
\[
 \|F_u-P_u\|_\infty\ll X^{1-\rho/8+o(1)},\qquad
 \int|F_u|^2\ll X^{1+o(1)},\qquad
 \int|P_u|^2\ll X^{1+o(1)}.
\]
The last estimate uses (22), disjoint assignment regions,
\(\int_{-1/2}^{1/2}|J(v)|^2dv\ll X\), and
\(\sum_{D\le R}D^{-1+o(1)}=X^{o(1)}\). Therefore
\[
 \int_0^1|F_uF_v-P_uP_v|^2\,d\beta
       \ll X^{3-\rho/4+o(1)}.
 \tag{28}
\]

## Complete model coefficients and their tail

Every major rational \(a/D\) owns the entire core
\(|\beta-a/D|<1/(2DT)\). A second valid approximant would imply
\[
 1\le D'/(2T)+D/T<1,
\]
which is impossible. Its assigned region lies inside \(1/(DT)\).
Extend \(J(v)^2\) to the real line at cost \(O(DT)\). Even summing
crudely over at most \(R^2\) modes, using \(|B_uB_v|\ll1\), costs
\[
 O(R^3T)=O(X^{0.98+3\rho}).
 \tag{29}
\]
This separation of \(R,T\) supplies a power error at length \(X\).
The same-cutoff parameters of the short-window argument cannot simply be
substituted with window length \(X\) in its crude model-tail estimate.

For \(X<N\le2X\), real Fourier inversion gives
\[
 \int_{\mathbb R}J(v)^2e(-Nv)\,dv=N.
\]
Define the completed coefficient by grouping full primitive numerator sums at
each denominator:
\[
 \mathfrak S_{u,v}(N)=\sum_{D\ge1}\sum_{a\ (D)}^*
          B_u(a/D)B_v(a/D)e(-aN/D).
 \tag{30}
\]
The convention at \(D=1\) is the single zero frequency. Write a nonzero
denominator as \(D=D_0D_1\), with \(D_0\mid12\),
\((D_1,12)=1\). Radial numerator dependence outside \(12\) gives a
Ramanujan sum after CRT, while the fixed-12 part has bounded size. A complete
denominator group is at most
\[
 D^{-2+o(1)}(D_1,N).
 \tag{31}
\]
For fixed \(0<a<1/2\),
\[
 \sum_D D^{-3/2+a}\sum_{s\mid(D,N)}s
 \le\zeta(3/2-a)\sum_{s\mid N}s^{-1/2+a}\ll\tau(N).
\]
Extracting \(R^{-1/2}\) from (31) thus bounds the tail of (30) by
\(O(R^{-1/2}X^{o(1)})\), uniformly for \(N\asymp X\). This is absolute
convergence over denominator groups, not a claim of absolute summability over
individual numerator modes.

Parseval in (28), then (29) and the completed tail at the \(O(X)\) target
integers, prove (1). The savings are respectively
\[
 \rho/4,\qquad 0.04-6\rho,\qquad\rho,
\]
each strictly exceeding \(10\varepsilon_*\). All family uniformities
and \(X^{o(1)}\) losses have fixed positive margins. Chebyshev with error
\(X^{1-\varepsilon_*}\), followed by a union over
\(O(X^{2\kappa+o(1)})\) pairs, leaves at most
\(O(X^{1-7\varepsilon_*})\) targets. This proves (2).

## Exact local main terms and survival ratios

Each local density above belongs to \(L^2(\mathbb Z_p)\). Its convolution
with another such density is continuous, and the additive Fourier series is
absolutely convergent by Cauchy--Schwarz. Finite-prime expansions and CRT
therefore identify (30) with the product of the local integrals. Passage to
all primes is justified by the absolute denominator-group tail, including
increasing exponents within a finite prime set. No Hardy--Littlewood
correlation hypothesis is used here.

For the raw pair put \(f_p(j)=h(p^j)\), \(x=p^{-1}\), and
\[
 I_{p,r}=\int_{\mathbb Z_p}f_p(v_p(t))f_p(v_p(N-t))\,dt,
 \qquad r=v_p(N).
\]
The exact valuation partition is
\[
 \begin{split}
 I_{p,r}={}&(1-x)\sum_{j<r}f_p(j)^2x^j
       +(1-2x)f_p(r)^2x^r\\
       &+2(1-x)f_p(r)\sum_{j>r}f_p(j)x^j.
 \end{split}
 \tag{32}
\]
The fixed-12 local pair factor is \(12\) for an admissible target. Thus
\[
 \mathfrak S_{h_6,h_6}(N)
       =12\rho_6^2\prod_{p\ge5}\frac{I_{p,v_p(N)}}{M_p^2}.
 \tag{33}
\]
At unit targets the factors in this normalized product are \(1+O(p^{-2})\).
For every positive target valuation, (32) gives uniformly
\[
 \frac{I_{p,r}}{I_{p,0}}
       =1+\frac{(f_p(1)-1)^2}{p}+O(p^{-2}).
 \tag{34}
\]
The tails use only polynomial valuation growth and geometric decay. Unit
strata give positive lower bounds \(I_{p,r}\ge1-p^{-1}\) for \(r\ge1\)
and \(I_{p,0}\ge1-2p^{-1}\). Local Cauchy--Schwarz gives a finite upper
bound independent of \(r\). These bounds also handle the finitely many small
primes uniformly. Formulae (33)--(34) prove (3), in agreement with the
[natural target normalization](SIMULTANEOUS_NORM_POINTWISE_SIEVE.md#the-normalization-is-the-natural-local-one).

At a marked prime, the integrand of
\(2\mathfrak S_{v_p,h_6}-\mathfrak S_{v_p,v_p}\) is the raw product
weight times the indicator that at least one component is divisible by \(p\).
Its local integral is exactly \(\beta_p(N)\) times the raw local integral,
for every target valuation. This proves (4). The corresponding identity
before any asymptotic is
\[
 B_p(N)=576\{2C_{v_p,h_6}(N)-C_{v_p,v_p}(N)\}.
 \tag{35}
\]
No finite-valuation or selected-subpopulation restriction is present.

For the fully sifted pair, let
\[
 g_r=\begin{cases}1-2/p,&r=0,\\1-1/p,&r\ge1.\end{cases}
\]
At a forbidden prime its normalized local numerator is
\(g_r/(1-p^{-1})^2\). The raw numerator is
\(I_{p,r}/[(1-p^{-1})D_p]^2\). The residue ratio (18) cancels the
extra \(D_p^2\), so
\[
 \frac{\mathfrak S_{h_C,h_C}(N)}{\mathfrak S_{h_6,h_6}(N)}
       =\prod_{p\text{ forbidden}}\frac{g_{v_p(N)}}{I_{p,v_p(N)}}
       =\prod_{p\text{ forbidden}}(1-\beta_p(N)).
 \tag{36}
\]
The [lifted local bound](SIMULTANEOUS_NORM_LOCAL_DENSITIES.md#the-local-theorem)
is \(\beta_p(N)\le9/p^2\), including \(p\mid N\). Hence the product
converges uniformly and
\[
 \Gamma(N)\ge1-\sum_{p\text{ forbidden}}\beta_p(N)
 \ge1-9\sum_{\substack{m\ge5\\\gcd(m,6)=1}}m^{-2}
       =10-\pi^2>0.
\]
This proves (5). Equations (2)--(5), together with the lower bound in (3),
prove (6). For admissible targets the raw positive candidate pairs are already
coprime to \(6\): an even common norm is divisible by four, whereas
\(N\equiv2\pmod4\), and the Eisenstein norms are zero or one modulo
three, whereas \(N\equiv2\pmod3\). Thus using \(h_6\) loses no
canonical coordinate count.

The survival proportion in (36) depends on the target. It is distinct from
the previously established averaged weighted survival proportion of about
\(63.3\%\).

## Dependencies and the remaining pointwise problem

The proof dependencies are
\[
 \begin{gathered}
 \text{norm Euler factors and four-character periodicity}
       \Longrightarrow\text{uniform shifted and rational means},\\
 \text{those means, Type I and Type II sums, and the repaired model}
       \Longrightarrow\text{power mean square (1)},\\
 \text{(1), simultaneous family control and lifted local densities}
       \Longrightarrow\text{the raw and fully sifted ratios (6)}.
 \end{gathered}
\]
The integer-four-character Fourier argument is developed in
[D12Power](PRESCRIBED_PRIME_POWER_SAVING_EXCEPTIONS.md). The Euler factors and
residue are established in the
[weighted AP note](SIMULTANEOUS_NORM_WEIGHTED_AP_MEANS.md); the uniform local
normalization is proved in the
[pointwise sieve note](SIMULTANEOUS_NORM_POINTWISE_SIEVE.md). Its prime harmonic
input is [Loughran--Matthiesen, Lemma 2.4](https://arxiv.org/html/1904.12845v2#S2).
The lifted bad-prime theorem and its exact prime-power scope are in the
[local-density note](SIMULTANEOUS_NORM_LOCAL_DENSITIES.md). The present shifted
means and strip estimates are proved here using periodic character sums,
convolution, geometric sums and Cauchy--Schwarz.

The global error is mean square. The same estimates give only a crude
pointwise product-coefficient error of order \(X^{3/2-c}\), using the
supremum error and an \(L^2\) norm, rather than \(o(X)\). They do not remove
the exceptional set or establish the every-target H1 and fixed-prime marginal
hypotheses. A further pointwise signed-correlation argument is needed there.
The unweighted semigroup-scale lower bound with power exceptions is already a
separate result of D12Power; this theorem supplies exact raw and fully sifted
coordinate-weight main terms and their actual marginal interface.

## Fixed macroscopic balance intervals

**Corollary SimNormBalancedPower.** Fix \(I=[\alpha,\beta]\subset(0,1)\)
with \(\alpha<\beta\), and define
\[
 C^I_{u,v}(N)=\sum_{\substack{a+b=N\\a/N\in I}}u(a)v(b).
\]
Outside one set of \(O_I(X^{1-7\varepsilon_*})\) targets, simultaneously
for every ordered pair \(u,v\in\mathcal U_X\),
\[
 C^I_{u,v}(N)=|I|N\mathfrak S_{u,v}(N)
                +O_I(X^{1-\varepsilon_*/20}).
 \tag{37}
\]
For admissible \(N\equiv2\pmod{12}\), the balanced raw count is
\(\asymp_I N\mathcal T(N)\), and its bad-prime and fully sifted ratios are
\(\beta_p(N)+O_I(X^{-\varepsilon_*/20})\) and
\(\Gamma(N)+O_I(X^{-\varepsilon_*/20})\), uniformly for
\(p\le X^\kappa\). In particular the actual balanced H1 and marginal
interface hold outside a power exceptional set. They are not asserted for
every target.

Here are the interval and uniformity details. All preceding Fourier estimates
hold uniformly for \(u\mathbf1_{J_1},v\mathbf1_{J_2}\), where
\(J_1,J_2\subset[0,2X]\) are arbitrary intervals. The major estimate is a
difference of two prefix estimates, followed by the same partial summation.
In the minor estimate, intersecting the product cutoff with an interval still
leaves an interval for every geometric sum. The same \(L^2\) bounds hold:
sequence restrictions decrease their squared sums, and each interval kernel
has real \(L^2\) mass at most \(2X\).
Replace \(J(v)^2\) by the product of the two interval kernels; its absolute
tail is still bounded by \(\min(2X,1/|v|)^2\). Its real Fourier coefficient
is
\[
 K_{J_1,J_2}(N)=|J_1\cap(N-J_2)|.
\]
The same Parseval and grouped-tail proof therefore gives the uniform mean
square (1), with \(N\mathfrak S\) replaced by
\(K_{J_1,J_2}(N)\mathfrak S\). No additional short-interval theorem is
being used.

Put \(\mu=\varepsilon_*/10\), \(\Delta=X^{1-\mu}\), and partition
\((X,2X]\) into \(O(X^\mu)\) blocks with left endpoints \(N_0\).
For each block choose
\[
 J_1=[\alpha N_0,\beta N_0],\qquad
 J_2=[(1-\beta)N_0,(1-\alpha)N_0].
\]
Both lie in \([0,2X]\). Chebyshev at error \(X^{1-\varepsilon_*}\),
union over these blocks and all ordered family pairs, excludes at most
\[
 O\!\left(X^{1-8\varepsilon_*+2\kappa+\mu+o(1)}\right)
       =O(X^{1-7\varepsilon_*})
\]
targets. For \(N\) in the block the fixed-interval correlation differs from
\(C^I_{u,v}(N)\) at only \(O_I(\Delta+1)\) possible values of \(a\).
Since \(u,v\le d_4\), their total weight is
\(O_I(\Delta X^{o(1)})\). Also
\[
 K_{J_1,J_2}(N)=|I|N+O_I(\Delta),\qquad
 \mathfrak S_{u,v}(N)\ll X^{o(1)}.
\]
The last bound follows uniformly from (31) and the same divisor sum as in
the completed tail. Both changes are absorbed by
\(O_I(X^{1-\mu/2})\), proving (37).

For a general interval \(I\), the exact bad-count identity is
\[
 B_{p,I}(N)=576\{C^I_{v_p,h_6}(N)+C^I_{h_6,v_p}(N)
                                      -C^I_{v_p,v_p}(N)\}.
 \tag{38}
\]
It reduces to twice the first mixed term only when the interval is invariant
under reflection \(t\mapsto1-t\). The two ordered local main coefficients
are equal, so (4) still gives the stated marginal ratio for every fixed
\(I\). The ratio errors can also be summed over \(p\le X^\kappa\),
since \(\varepsilon_*/20>\kappa\). This includes the balanced interval
\([1/3,2/3]\) with its full moving value range, not a selected central
subpopulation.

## A gap compatible with strong Fourier norm control

The need for further pointwise information is not merely a weakness of a
particular error estimate. There are nonnegative congruence-supported weights
with a complete gap at one compatible target despite much stronger Fourier
and mean-square closeness to a positive model.

For every integer \(m\ge0\), set \(N=24m+14\), put
\(M=(N-2)/12=2m+1\), and consider
\[
 a_j=12j+1,\qquad 0\le j\le M,\qquad N-a_j=a_{M-j}.
\]
On \(1\le a<N\) let \(v(a)=\mathbf1_{a\equiv1\ (12)}\), with zero
outside this interval. Choose independent signs \(\epsilon_j\) for \(0\le j\le m\), set
\(\epsilon_{M-j}=-\epsilon_j\), and define
\(u(a_j)=1+\epsilon_j\), with zero elsewhere. Thus \(u\) takes values
zero and two and has the same total mass as \(v\). For every choice of signs,
\[
 C_{u,u}(N)=0,\qquad C_{v,v}(N)=M+1\asymp N.
 \tag{39}
\]

Let \(F=\sum_a u(a)e(a\beta)\), \(D=\sum_a v(a)e(a\beta)\), and
\(E=F-D\). For a fixed \(\beta\),
\[
 E(\beta)=\sum_{j=0}^m\epsilon_j
              \{e(a_j\beta)-e((N-a_j)\beta)\}.
\]
Each real or imaginary coefficient has modulus at most two. The bound
\(\cosh t\le\exp(t^2/2)\), exponential Markov, and a union over real
and imaginary parts give
\[
 \mathbb P(|E(\beta)|>T)
       \le4\exp\{-T^2/[16(m+1)]\}.
\]
Apply this on the \(N^3\)-point grid with
\(T=16\sqrt{N\log N}\). Since \(N/(m+1)>12\), the probability of
any failure is at most \(4N^{3-192}<1\). There is therefore a choice of
signs with that bound at every grid point. Also
\[
 |E'(\beta)|\le2\pi\sum_{a:v(a)=1}a\le2\pi N^2,
\]
so interpolation from the grid costs at most \(2\pi/N\). This choice
satisfies \(\|F-D\|_\infty\ll\sqrt{N\log N}\).
Parseval gives
\[
 \int|F|^2=4(m+1),\qquad\int|D|^2=2(m+1),
\]
and hence
\[
 \sum_t|C_{u,u}(t)-C_{v,v}(t)|^2
   =\int|F^2-D^2|^2
   \le\|F-D\|_\infty^2\int|F+D|^2\ll N^2\log N.
 \tag{40}
\]
**Explicit-family corollary.** For every integer \(k\ge0\), the compatible
target \(N=24\cdot2^k-10\) admits such weights with the sharper bounds
\[
 \|F-D\|_\infty\le2\sqrt{2\cdot2^k},\qquad
 \sum_{t\in\mathbb Z}|C_{u,u}(t)-C_{v,v}(t)|^2\le80\cdot2^{2k}.
 \tag{41}
\]
The weight \(u\) has an explicit recursion, takes values zero and two, and has
the same mass as the periodic baseline on \(1\le a<N\). Their target counts are
\(C_{u,u}(N)=0\) and \(C_{v,v}(N)=2^{k+1}\asymp N\).

To prove this, define the Rudin--Shapiro polynomials by
\[
 P_0(z)=Q_0(z)=1,\qquad
 \begin{cases}
 P_{k+1}(z)=P_k(z)+z^{2^k}Q_k(z),\\
 Q_{k+1}(z)=P_k(z)-z^{2^k}Q_k(z).
 \end{cases}
\]
Inductively, both polynomials have degree less than \(L=2^k\) and all \(L\)
coefficients are signs: the two summands in each recursion have disjoint
coefficient ranges. On \(|z|=1\), cancellation of the cross terms gives
\[
 |P_{k+1}(z)|^2+|Q_{k+1}(z)|^2
       =2\{|P_k(z)|^2+|Q_k(z)|^2\}.
\]
The initial sum is two, so induction gives
\(|P_k(z)|^2+|Q_k(z)|^2=2L\), and in particular
\(|P_k(z)|\le\sqrt{2L}\).

Write \(P_k(z)=\sum_{j=0}^{L-1}\epsilon_jz^j\), put \(M=2L-1\), and
extend the signs by \(\epsilon_{M-j}=-\epsilon_j\) for \(0\le j<L\).
Use the weights above with \(m=L-1\) and these signs. Each reflected pair
has one weight two and one weight zero, proving equality of masses and the
target gap. With \(z=e(12\beta)\), the Fourier difference is exactly
\[
 E(\beta)=e(\beta)\{P_k(z)-z^MP_k(z^{-1})\}.
\]
Both polynomial values have modulus at most \(\sqrt{2L}\), proving the
supremum bound in (41). Exactly \(L\) coefficients of \(F+D\) are three
and \(L\) are one, so Parseval gives \(\int_0^1|F+D|^2=10L\).
Applying the product estimate in (40) yields the full correlation bound
\((8L)(10L)=80L^2\), proving (41). Thus this infinite family has Fourier
error \(O(\sqrt N)\) and full correlation mean square \(O(N^2)\). Since
\(u(1)=2\), the constructed weight is not multiplicative.

Thus positivity, congruence support, square-root Fourier closeness and this
strong mean-square estimate alone do not force a representation at every
compatible target. These weights are not multiplicative. The construction
is not a counterexample to \(C_{12}\), its local Euler model, or a method
using additional arithmetic information.

## Scope and review

For comparison, [Loughran--Matthiesen, Theorem 1.9](https://arxiv.org/html/1904.12845v2#S1.SS4)
gives an asymptotic with an o(1) normalized error for fixed nonnegative
Frobenian multiplicative functions on linear forms whose nonconstant parts
are pairwise independent, averaged over a fixed positive-volume convex body.
Its statement does not give the uniform power error for a growing family of
bad-prime masks used here. The proof above derives that error from the special
integer-four-character structure. This comparison is about theorem scope;
priority for the specialized quantitative result remains unaudited.

The shifted-base parity and character phases, polynomial AP budgets, strip
geometry, periodic indicators, common major assignment, grouped singular-series
tails, local normalizations and power budgets received separate internal
mathematical reviews. These do not constitute external specialist review,
formal verification or an established priority determination. No new finite
coverage search is used, and no practical cutoff is claimed.
