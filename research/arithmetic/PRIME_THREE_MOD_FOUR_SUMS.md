# Binary sums with prime factors congruent to three modulo four

The subgroup restriction can be removed for this particular unfiltered
prime prescription. Two changes are essential: retain the negative
character branch in the rough density, and reproduce characters on a
fixed prime coset with their slot-dependent phases.

The argument uses exactly inputs ZF and G of the
[existing input audit](PRIMITIVE_NORM_BINARY_IMPORTED_INPUTS.md):
the fixed all-conductor Dirichlet zero-free half-plane and the physical
weighted-dilation root-moment theorem. The statement below makes that
dependence explicit. The source record supplies those inputs in the
current research framework; this application has internal written
review, not external human review or formal verification.

## 1. Pointwise count and local criterion

Put
\[
 C^-_4=\{n\ge1:p\mid n\Longrightarrow p\equiv3\pmod4\},\qquad
 R_-(N)=\#\{(a,b)\in(C^-_4)^2:a+b=N\}.
\]
We include 1 and count ordered numerical pairs. Define
\[
 A_-=\frac12\prod_{p\equiv1(4)}
                  \left(1-\frac2{p(p-1)}\right)>0,\qquad
 J_-(N)=\prod_{\substack{p\mid N\\p\equiv1(4)}}\frac{p-1}{p-2}.
                                                               \tag{M.1}
\]

**Theorem PrimeThreeModuloFourCorrelation.** Under the same ZF and G
inputs as the subgroup correlation theorem, uniformly over sufficiently
large even integers \(N\),
\[
 R_-(N)=A_-J_-(N)\frac N{\log N}
       \left(1+O((\log N)^{-1/16384})\right).                    \tag{M.2}
\]
Odd targets have zero count identically. Evenness is also exactly the
condition for representability in every finite quotient. Thus every
sufficiently large locally admissible target has a representation.

The error is uniform in all target prime valuations. There is no
numerical cutoff or exact small-exception classification. No assertion
for prescribed residue states of the individual summands, or for all
non-subgroup prime prescriptions, is part of this theorem.

Here the generated group is
\(H=\langle3\bmod4\rangle=(\mathbb Z/4)^\times\), but the prime density
is \(\rho=1/2\). In particular \(\rho\ne1/[(\mathbb Z/4)^\times:H]\).
The subgroup density-index relation must not be used.

## 2. The marginal and the negative secondary branch

Write \(f=\mathbf1_{C^-_4}\), \(\chi=\chi_{-4}\), and
\[
 F(s)=\prod_{p\equiv3(4)}(1-p^{-s})^{-1},\qquad
 P_3(s)=F(2s),\quad P_3=P_3(1).
\]
Euler factors give
\[
 F(s)^2=(1-2^{-s})\frac{\zeta(s)}{L(s,\chi)}P_3(s),\qquad
 F_\chi(s):=\sum_n\chi(n)f(n)n^{-s}=\frac{F(2s)}{F(s)}.
                                                               \tag{M.3}
\]
Near one, with positive determinations on the real line to its right,
\[
 F(s)=(s-1)^{-1/2}h_0(s),\qquad
 F_\chi(s)=(s-1)^{1/2}h_\chi(s),
\]
where
\[
 h_0(1)^2=\frac{2P_3}{\pi},\qquad
 \frac{h_\chi(1)}{h_0(1)}=\frac\pi2,\qquad
 c_-=\frac{h_0(1)}{\Gamma(1/2)}=\frac{\sqrt{2P_3}}\pi.
\]
The fixed-function Selberg--Delange theorem therefore gives
\[
 \sum_{n\le X}f(n)\sim c_-\frac X{\sqrt{\log X}},\qquad
 \sum_{n\le X}\chi(n)f(n)
 =-\frac{\pi c_-}{4}\frac X{(\log X)^{3/2}}
       +O\left(\frac X{(\log X)^{5/2}}\right).                  \tag{M.4}
\]
The corrected author edition of de la Bretèche--Tenenbaum,
Theorem 1.2, printed p. 3, permits a complex parameter, including
\(-1/2\); its normalization (1.10) supplies
\(\Gamma(-1/2)^{-1}=-1/(2\sqrt\pi)\).
Alternatively the fixed-cutoff Hankel calculation below gives these
two leading coefficients directly.

Both odd states have leading density \(c_-/2\). Their equality at
leading order does not justify omitting the secondary branch when
arbitrary logarithmic centering is required.

## 3. Two-branch rough density

Set
\[
 L=\log X,\quad y=\exp(L^{2/5}),\quad S=\{p:p\le y\},\quad
 V=\prod_{p\le y}(1-1/p).
\]
Write \(n=A v\), with \(A=n_S\) and \((v,\prod_{p\le y}p)=1\).
Removing the allowed small Euler factors gives
\[
 F_R(s)=(s-1)^{-1/2}h_{0,R}(s),\qquad
 F_{\chi,R}(s)=(s-1)^{1/2}h_{\chi,R}(s).
\]
Exactly,
\[
 F_{\chi,R}(s)=\frac{F_R(2s)}{F_R(s)},\qquad
 \frac{h_{\chi,R}(1)}{h_{0,R}(1)}
   =\frac\pi2\prod_{\substack{p\le y\\p\equiv3(4)}}\frac{1+1/p}{1-1/p}
   \sim e^\gamma\log y.                                      \tag{M.5}
\]
In particular \(h_{0,R}(1)\asymp(\log y)^{-1/2}\) and
\(h_{\chi,R}(1)\asymp(\log y)^{1/2}\).
A disk of radius \(c/\log y\) gives, for each fixed \(j\),
\[
 |h_{\sigma,R}^{(j)}(1)|/j!
       \ll_j h_{\sigma,R}(1)(C\log y)^j,\qquad \sigma=0,\chi.
                                                               \tag{M.6}
\]
Indeed the logarithm of each finite Euler product relative to its
value at one is bounded by
\(O(|s-1|\sum_{p\le y}\log p/p)=O(c)\); apply Cauchy's estimate.

For a fixed sufficiently deep \(K\), put
\[
 \begin{split}
 M_0(t)&=\sum_{j=0}^K\frac{h_{0,R}^{(j)}(1)}{j!}
             \frac{(\log t)^{-1/2-j}}{\Gamma(1/2-j)},\\
 M_\chi(t)&=\sum_{j=0}^K\frac{h_{\chi,R}^{(j)}(1)}{j!}
             \frac{(\log t)^{-3/2-j}}{\Gamma(-1/2-j)}.
 \end{split}                                                  \tag{M.7}
\]
Uniformly on every fixed power bulk \(X^a\le t\le X^C\),
\[
 M_0(t)>0,\qquad |M_\chi(t)|\ll L^{-3/5}M_0(t).                 \tag{M.8}
\]
The first negative-branch coefficient is negative, but the combined
density is positive. Extend \(M_0\) positively and \(M_\chi\) smoothly
outside the bulk, preserving (M.8) and the global bound
\(M_0(t)\ll h_{0,R}(1)L^{-1/2}\) on all relevant arguments.
For example first make \(M_\chi\) zero below a smaller bulk and
extend \(M_0\) by a bounded positive smooth function. Define
\[
 D(n)=\frac{f(A)}V\{M_0(v)+\chi(v)M_\chi(v)\},\qquad
 D_0(n)=\frac{f(A)}V M_0(v),\qquad \xi(n)=f(n)-D(n).             \tag{M.9}
\]

**Lemma MinusFourInvariantModel.** The model is nonnegative,
\(D(n)\ll L^{-3/10}\), and for every allowed \(p\le y\),
\[
 f(pn)=f(n),\quad D(pn)=D(n),\quad \xi(pn)=\xi(n)               \tag{M.10}
\]
exactly, including repeated factors. It has arbitrary logarithmic
rough centering against bounded periodic tests of period at most
a fixed power of \(L\), with fixed polynomial logarithmic smooth
costs, after \(K\) is chosen sufficiently large.

For invariance, multiplication changes \(A\) to \(pA\) and leaves
\(v\) unchanged. The character in (M.9) is \(\chi(v)\), not
\(\chi(n)\). In \(\chi(v)=\chi(n)\chi(A)\), both signs change.
The size bound follows from
\(h_{0,R}(1)\asymp L^{-1/5}\), \(V\asymp L^{-2/5}\).
The negative branch is smaller by (M.8).

Here is the centering argument with the required uniformity.
Inflate characters \(\psi\) to a modulus divisible by four. Literal
Euler factors satisfy
\[
 F_\psi(s)^2=\frac{L(s,\psi)}{L(s,\psi\chi)}
       \prod_{p\equiv3(4)}
                (1-\psi(p)^2p^{-2s})^{-1}.                    \tag{M.11}
\]
Only \(\psi=1,\chi\) have a zeta branch at one; their exponents are
\(1/2,-1/2\). All conductor primes at polylogarithmic moduli are
already below \(y\).

Negative powers require reciprocal \(L\)-bounds. Replace a principal
\(L(s,\psi)\) by
\(L^\#(s,\psi)=(s-1)L(s,\psi)/(s+1)\), and otherwise put \(L^\#=L\).
Under ZF, uniformly on \(23/24\le\Re s\le2\),
\[
 |L^\#(s,\psi)|+|L^\#(s,\psi)|^{-1}
       \ll(q(2+|t|))^C.
\]
For precision, take a fixed disk about \(2+it\) contained in
\(\Re s>11/12\). Convexity bounds the real part of \(\log L\)
by \(O(\log(q(2+|t|)))\), the center logarithm has the same bound, and
Borel--Carathéodory on a smaller disk gives the same bound for
\(|\log L|\). Exponentiating controls both signs.

Use smooth Mellin inversion to height \(\exp(L^{1/10})\), and move
to \(1-c/\log y\). Direct and reciprocal factors cost
\(\exp(O(L^{1/10}))\), finite Euler products cost fixed powers of \(L\),
and the shifted bulk integral saves \(\exp(-cL^{3/5})\).
Smooth tails are negligible after a fixed sufficiently high derivative
order. Hankel integration about one, with (M.6), gives (M.7) for
both exponents. Their relative Taylor remainders are
\(O_K(L^{-3(K+1)/5})\); the negative branch is additionally
\(O(L^{-3/5})\) relative to the positive model.

On the comparison side the unrestricted rough population has only
its principal density \(V\). Thus its \(M_0\) and \(\chi(v)M_\chi\)
pieces match precisely the two original branches. Character
orthogonality proves centering against every permitted periodic
test; no secondary branch is left over. This proves the lemma.
A fixed-function Selberg--Delange theorem alone would not establish
this moving-\(y\) assertion.

Retained smooth records satisfy \(A\le E=\exp(\sqrt L)\). The
discarded-record bound is the same tilted Euler-product estimate
as Module D, (W8), with positive prime density \(1/2\). Recordwise
centering therefore allows the bounded coefficients used by the
comparison construction. For reflected endpoints use \(D(|n|)\)
and \(\chi(|v|)\); a character evaluated on a signed rough argument
adds only its fixed value at \(-1\).

## 4. Character frames on prime banks in the coset 3 modulo 4

Choose the small and big banks from primes \(p\equiv3\pmod4\)
in the location ranges of Module D, deleting \(p\mid N\).
The imported graph source, Definition 3.1, pp. 10--11, and
Theorem 3.5, p. 12, permits arbitrary finite prime groups with the
specified size and reciprocal masses. Its own interval-prime
ideal-kernel theorem is not being applied literally.

**Lemma MinusFourCosetFrame.** The corrected quotient-character
construction of Module I works for these banks, with exactly the
same operator, reproduction and comparison bounds.

For proof, quotient primitive characters of conductor at most \(Q\)
by the relation
\(\operatorname{cond}(\psi\overline\eta)\mid4\), and choose one
representative per class. For a fixed selected slot set \(\Lambda\),
the input product \(D_\Lambda\) and comparison product \(X_\Lambda\)
have the same slot type. Hence
\[
 \chi(D_\Lambda)=\chi(X_\Lambda)=(-1)^{|\Lambda|}.
\]
If \(\psi=\eta\chi^e\) on the relevant units, both products obey
\[
 \psi(D_\Lambda)=c_{\psi,\Lambda}\eta(D_\Lambda),\quad
 \psi(X_\Lambda)=c_{\psi,\Lambda}\eta(X_\Lambda),\quad
 c_{\psi,\Lambda}=(-1)^{e|\Lambda|}.                            \tag{M.12}
\]
For inequivalent representatives, the character sum over the coset
\(3\bmod4\) is a unit phase times its sum over the kernel
\(1\bmod4\), and is zero. The prime theorem at polylogarithmic
moduli, followed by partial summation on arbitrary subintervals,
therefore proves exactly the cell estimate (AP5).
Deleting target divisors costs
\(O(\log N/(P_{\min}\log P_{\min}))\), smaller than every
negative logarithmic power, uniformly in target valuations.

Use logarithmic product cells \(I_j\), with
\(\nu_j=\mathbb P(\log D_\Lambda\in I_j)\), and keep the conjugated frame
\[
 K_\Lambda(D,X)=
 \sum_{\nu_j\ge\zeta}
 \frac{\mathbf1_{\log D\in I_j}\mathbf1_{\log X\in I_j}}{\nu_j}
       \sum_{\eta\in{\cal R}_Q}\overline{\eta(D)}\eta(X).
                                                               \tag{M.13}
\]
Its Gram diagonal is one and its off-diagonal entries are
\(O_B(\zeta^{-1}L^{-B})\). Thus it is positive, has operator norm
at most two, and obeys
\[
 |K|\le Q^2\zeta^{-1},\qquad \sup_X\mathbb E_D|K(D,X)|\le Q^2.
\]
Its matching term reproduces \(\psi\):
\[
 \mathbb E_D K(D,X)\psi(D)
 =c_{\psi,\Lambda}\eta(X)+O_B(Q^2\zeta^{-1}L^{-B})
 =\psi(X)+O_B(Q^2\zeta^{-1}L^{-B}).                             \tag{M.14}
\]
In rational-phase expansions, grouping equivalent characters only
multiplies coefficients by unit phases. Their total absolute
coefficient mass does not increase. The proof of (AP13), including
the fixed ideal/probe construction and symmetric-probe cancellation,
therefore retains every bound and parameter inequality of Module I.
This proves the lemma. A growing number of slots causes no loss:
the phases in (M.12) always have absolute value one.

## 5. Target-grid comparison and reconstruction

We specify all uses of the old subgroup hypotheses in transferring
the remaining proof modules.

| Interface | Verification for \(C^-_4\) |
|---|---|
| Endpoint dilation and extraction | (M.10) gives exact invariance at every selected allowed prime. The rough bracket is the full two-branch bracket in (M.9). |
| Large-shift and reflected pairing | Module G uses the physical root moment, bounded endpoint/list moments and the corrected ideal frame. All are available; reflection uses absolute arguments as above. |
| Nonexceptional target characters | A bank sum over \(p\equiv3(4)\) is one half the difference of the prime sums for \(\eta\) and \(\eta\chi\). If \(\eta\) is not induced modulo four, neither is principal. The same ZF estimate saves \(P^{-1/24}\) up to fixed logarithmic factors. |
| Exceptional target characters | Both characters induced modulo four are centered by (M.7)--(M.9), after the exact smooth record and gcd split. This step would fail at arbitrary accuracy if \(M_\chi\) were dropped. |
| Target gcd profiles | The residual periodic test is \(\eta(v)\), and the record contributes the unit phase \(\eta(A/d)\). Restoring the large-prime avoidance condition costs \(O_B((Y/d)L^{-B})\), exactly as in Module E. The record sum is \(1/(dV)\), retaining its \(1/d\). |
| Positive bank envelope | The forbidden good primes are \(1\bmod4\), of density \(1/2\). The two-endpoint sieve dimension is one. The bank labels are allowed and do not divide \(N\). Formula (W13), the CRT counts and smooth-record tails use these facts, not that the residue is the identity. |

In particular the cofactor/modulus ratio in the energy argument is
still \(M/Q_d=Y/(QP)\), independently of \(d\). At the actual target
dyad its lower bound has the form
\[
 \frac{Y'}{qP}\gg
 \exp(cL^{1/3}-CL^{1/4})L^{-B_{\rm out}-\delta}.                 \tag{M.15}
\]
It beats every fixed logarithmic cost. The nonexceptional energy
and the centered exceptional energy are therefore (G13) and (G16).
Summing \(1/d\) and using discrete additive orthogonality gives
(G17) at every requested fixed precision, without a leftover
factor of \(N\).

For the positive envelope, \(D\) is bounded by a fixed multiple
of \(D_0\), so the same bound applies to all four actual/model
pair types. The exact reconstruction (W15)--(W17) uses
\(f(pn)=f(n)\) and \(D(pn)=D(n)\), not unchanged residue states
of the whole endpoint. The two lemmas above consequently supply
the comparison and reconstruction required by the main proof.

## 6. Leading model and complete error budget

On the retained bulk and records,
\[
 |D-D_0|\ll L^{-3/5}D_0,\qquad
 |D*D-D_0*D_0|\ll L^{-3/5}(D_0*D_0).                          \tag{M.16}
\]
This pointwise domination controls all secondary mixed terms.
Replacing the positive rough logarithm by \(L\) costs
\(O(L^{-1/2})\) relatively. With
\(E_1(y)=\prod_{p\le y,\ p\equiv1(4)}(1-1/p)\),
the remaining leading model is
\[
 \lambda_X f(n_S),\qquad
 \lambda_X=\frac{c_-L^{-1/2}}{(1/2)E_1(y)}.
\]
The presieve restricts \(a,b\) to odd integers, contributing
\(1/2\) on every even target. Dividing by the squared marginal
\( (1/2)^2\) gives finite factor two. At \(p\equiv1(4)\), its
normalized pair factor is
\[
 \frac{1-(2-\mathbf1_{p\mid N})/p}{(1-1/p)^2};
\]
allowed good primes have factor one. Thus the model coefficient is
\(2c_-^2 C_1J_-(N)\), where
\(C_1=\prod_{p\equiv1(4)}(1-(p-1)^{-2})\).
The identity
\(\prod_{p\ {\rm odd}}(1-p^{-2})^{-1}=\pi^2/8\) gives
\[
 2c_-^2C_1
 =\frac{4P_3 C_1}{\pi^2}
 =\frac12\prod_{p\equiv1(4)}\left(1-\frac2{p(p-1)}\right)=A_-.
\]
The finite local criterion follows either directly from odd unit
pairs modulo four and positive unit-pair counts at \(p\equiv1(4)\),
or from the finite-image theorem in the
[Hardy--Littlewood note](MULTIPLICATIVE_HARDY_LITTLEWOOD.md#5-prime-factor-semigroups-beyond-subgroups).
CRT combines the conditions. Oddness is the only target exclusion.

Here is an explicit precision order, using the estimates of Module J
whose hypotheses were verified in Sections 3--5. Put
\[
 \phi=\varphi(4)=2,\quad \delta=1/(8192\phi)=1/16384,\quad
 {\cal S}_N=NJ_-(N)/L.
\]
Bank density is \(1/\phi\); it is unrelated to the index of \(H\).
The small and big location ranges have reciprocal masses
\((1/(8\phi)+o(1))\log L\) and
\((1/(15\phi)+o(1))\log L\). Take
\(b=\lfloor\log L/(32\phi)\rfloor\) groups of each type, each of
mass \(1+o(1)\). Formula (W13) gives success probability tending to
\(e^{-2}>1/8\). Chernoff bounds the rare-bank probability by
\(2L^{-1/(4096\phi)}=2L^{-2\delta}\), and the exact reconstruction
coefficient mass is at most \(L\).

Use the two-sided cutoff \(\epsilon=L^{-\delta}\), external phase
range \(L^3\), and total localization reserve \(L^{8\delta}\).
Request graph error \(L^{-22}\) per source block before the at most
\(L^2\) blocks, then fix the ideal/probe data and eventual mark length.
Let \(C_{\rm cmp},C_w,C_{\rm low}\) be the resulting fixed coefficient,
endpoint and centering exponents. They do not depend on the final
Hankel depth. If \(m_*\) bounds the free draws per group, take
\[
 K_*=\frac{\log(m_* !)+m_*\log2}{16\phi},\qquad
 C_{\rm prod}=2K_*+2,
\]
and choose
\[
 \begin{split}
 A_{\min}&>C_{\rm cmp}+2C_w+40+8\delta,\\
 B_{\rm out}&>2A_{\min}+C_{\rm prod}+4,\\
 T_E&>2(C_{\rm cmp}+C_w+40+8\delta).
 \end{split}                                                  \tag{M.17}
\]
The elementary two-free-product multiplier bound of (J.8) pays the
minor comparison. Equations (M.15), (G15)--(G17), and Parseval pay
the major comparison at precision \(T_E\).
Choose the fixed density depth last, with
\(3(K+1)/5>T_E+C_{\rm low}+40\), and choose the Mellin derivative
order with \((3-\delta)(j-1)>C_{\rm tail}+40\).
Extra branches only change fixed constants. For each fixed depth,
\(D\le1\) eventually; its smooth seminorms introduce no
depth-dependent exponent in the earlier comparison family.

The resulting complete error list is:

| Error | Bound |
|---|---|
| Rare banks | \(O({\cal S}_N L^{-2\delta})\) |
| Graph comparison and reconstruction | \(O(NL^{-19+8\delta})=o({\cal S}_N L^{-\delta})\) |
| Nonexceptional prime sums, discarded smooth records, deleted target primes and bank-square records | Exponential savings times fixed logarithmic powers |
| CRT endpoints | \(N^{0.45+o(1)}\) |
| Positive rough density variation | \(O({\cal S}_N L^{-1/2})\) |
| Negative rough branch in model evaluation | \(O({\cal S}_N L^{-3/5})\) |
| Fundamental-lemma and completed Euler tails | \(O_B({\cal S}_N L^{-B})\) for every fixed \(B>0\) |
| Endpoint strips | \(O(\epsilon{\cal S}_N+N^{0.45+o(1)})\) |

Every entry is \(O({\cal S}_N L^{-\delta})\). This proves (M.2).
No new graph estimate, compound-prime bank, zero-free hypothesis or
restriction on target valuations has been inserted.

## 7. Historical scope and limits

Dietmann--Elsholtz, *Sums of two numbers having only prime factors
congruent to one modulo four*, Banach Center Publications 118 (2019),
177--179, proves the plus-class existence theorem. Its final remark,
printed p. 179, explicitly poses the prime-\(3\bmod4\) sum problem
on \(N\equiv2\bmod4\). The
[published primary PDF](sources/dietmann-elsholtz-prime-factors-2019-published.pdf)
is retained. The problem is an established historical question.
That 2019 statement does not determine its current literature frontier.

The Gaussian/Greaves argument cannot simply exchange prime classes:
a prime \(1\bmod4\) can divide a primitive sum of two squares,
as \(5=1^2+2^2\) shows. The extension here is a new application proof
within the already recorded G/ZF framework. A bounded literature
review and internal agent agreement do not establish historical
priority, external mathematical review or formal verification.
The source graph theorem and zero-free theorem are not contributions
of this application. General residue sets and fixed summand-state
filters require their own verified comparison.
