# Binary sums supported on primes 3 modulo 10

This case is also the order-four specialization of
[the fixed prime-coset self-sum theorem](PRIME_COSET_SELF_SUMS.md), using
the [shared arithmetic-model comparison](SPARSE_ARITHMETIC_MODEL_COMPARISON.md).
The standalone flagship gives the common proof once, followed by its coset
checks and this explicit example in Section 9.5. The calculations below retain
the quartic-character normalization and its local factor at five.

The order-four prime class gives a further unfiltered extension of
[the prime-\(3\bmod4\) proof](PRIME_THREE_MOD_FOUR_SUMS.md).
All global conclusions retain exactly the same ZF and G inputs,
specified in the [input audit](PRIMITIVE_NORM_BINARY_IMPORTED_INPUTS.md).
The new points are complex conjugate Hankel branches, their Fourier
orientation, control of the chosen analytic logarithms, and fourth-root
slot phases. This is an internally reviewed written application.
No external human or formal review, effective cutoff or historical
priority certification is asserted. General prime cosets, unions of
cosets and prescribed individual summand states remain outside this
specific theorem.

## 1. Statement and exact arithmetic factor

Let
\[
 C=C_{10,3}=\{n\ge1:p\mid n\Longrightarrow p\equiv3\pmod{10}\},
 \qquad f=\mathbf1_C,\qquad
 R(N)=\#\{(a,b)\in C^2:a+b=N\}.
\]
We include 1 and count ordered numerical pairs. Choose the order-four
character \(\chi\bmod10\) with \(\chi(3)=i\).
The generated group \(\langle3\rangle\) is all of
\(G=(\mathbb Z/10)^\times\), but the prime density is
\(\rho=1/4\). The density is not the reciprocal of this subgroup's index.

Define
\[
 c=\frac1{\Gamma(1/4)}
       \lim_{s\downarrow1}(s-1)^{1/4}
              \prod_{p\equiv3(10)}(1-p^{-s})^{-1}>0.
\]
Equivalently,
\[
 c=\frac1{\Gamma(1/4)}
    \prod_p(1-p^{-1})^{\,1/4-\mathbf1_{p\equiv3(10)}}.
\]
The latter positive product is in increasing-prime order. Its convergence
follows by the fixed-AP prime theorem; an absolute Euler correction is
given below. Put
\[
 {\cal E}=\{p:p\nmid10,\ p\not\equiv3\pmod{10}\},\qquad
 C_{\cal E}=\prod_{p\in{\cal E}}\left(1-\frac1{(p-1)^2}\right),
\]
\[
 A=\frac{15}{8}c^2 C_{\cal E},\qquad
 J(N)=\left(\frac43\right)^{\mathbf1_{5\mid N}}
         \prod_{\substack{p\mid N\\p\in{\cal E}}}\frac{p-1}{p-2}.
\]

**Theorem PrimeThreeModuloTenCorrelation.** Under the existing ZF and G,
uniformly for sufficiently large even integers \(N\),
\[
 R(N)=A J(N)\frac{N}{(\log N)^{3/2}}
          \left(1+O((\log N)^{-1/32768})\right).
\]
Odd targets have zero count identically. The error is uniform in all
prime valuations of N. There is no effective numerical onset.

The Euler factors can also be written by including5 among the forbidden
odd primes:
\[
 {\cal F}=\{p\text{ odd}:p\not\equiv3\pmod{10}\},\quad
 A=2c^2\prod_{p\in{\cal F}}\left(1-\frac1{(p-1)^2}\right),\quad
 J(N)=\prod_{\substack{p\mid N\\p\in{\cal F}}}\frac{p-1}{p-2}.
\]
These two normalizations agree exactly. It is essential not to omit
the target-dependent factor at5.

Every value in C is a unit modulo10. For even N the number of unit pairs
modulo10 with a+b=N is3 when5 does not divide N and4 when5 divides N.
Thus the finite normalized factor is
\[
 \gamma_{10}(N)=\frac{10}{16}\#\{(a,b)\in G^2:a+b=N\bmod10\}
 =\begin{cases}15/8,&5\nmid N,\\5/2,&5\mid N.\end{cases}
\]
At p in E the normalized pair factor is
\[
 \frac{1-(2-\mathbf1_{p\mid N})/p}{(1-1/p)^2}
 =\left(1-\frac1{(p-1)^2}\right)
   \left(\frac{p-1}{p-2}\right)^{\mathbf1_{p\mid N}}.
\]
Allowed good primes have factor1. All forbidden odd primes are at least5,
so every such local pair factor is positive. The finite-image theorem
for arbitrary prime residue sets, applied to S={3}, gives the exact
all-finite-quotient criterion: N is even. Dirichlet/CRT words in the
generator3 realize every unit state modulo10, after extracting allowed
valuation prefixes. Small admissible exceptions are not excluded:
for example8 has no actual pair.

## 2. All four Dirichlet-series branches

Set
\[
 F_j(s)=\sum_n f(n)\chi(n)^j n^{-s}
       =\prod_{p\equiv3(10)}(1-i^j p^{-s})^{-1},
 \qquad j=0,1,2,3,
\]
and
\[
 \alpha_j=i^j/4\in\{1/4,i/4,-1/4,-i/4\}.
\]
For literal characters modulo10,
\[
 F_j(s)=E_j(s)
          \prod_{r=0}^3 L(s,\chi^r)^{\,i^{j-r}/4},             \tag{T.1}
\]
where the logarithm of E_j is absolutely convergent in Re s>1/2.
The identity follows by cancelling prime coefficients:
\[
 \mathbf1_{p\equiv3(10)}
       =\frac14\sum_{r=0}^3 i^{-r}\chi(p)^r.
\]
At2 and5 every literal character vanishes; these omitted Euler factors
are carried by L(s,chi^0), with no invented primitive Euler factor.
The powers and logs in(T.1) are defined by Euler logarithms on Re s>1
and continued coherently to the zero-free domain.

Near one,
\[
 F_j(s)=(s-1)^{-\alpha_j}h_j(s),
\]
with h_j analytic and nonzero. The fixed complex-parameter
Selberg--Delange theorem gives, to any fixed depth,
\[
 \sum_{n\le X}\chi(n)^j f(n)
 =\frac{h_j(1)}{\Gamma(\alpha_j)}
        X(\log X)^{\alpha_j-1}
   +O_j\!\left(X(\log X)^{\Re\alpha_j-2}\right).
\]
Thus the nonreal secondary terms have absolute order X/log X,
oscillating through exp(±i loglog X/4); they are not negligible to
arbitrary logarithmic precision.

Conjugacy is exact:
h_3(s)=overline{h_1(bar s)}, while h_0,h_2 are real on the real axis.
Two useful exact identities are
\[
 F_0(s)F_2(s)=F_0(2s),\qquad F_1(s)F_3(s)=F_2(2s).             \tag{T.2}
\]
In particular h_2(1)=F_0(2)/h_0(1)>0 and
|h_1(1)|^2=F_2(2)>0. The complex-parameter theorem is the same retained
de la Breteche--Tenenbaum Theorem1.2 used for the minus4 case;
it expressly permits alpha in C. Its fixed-function conclusion is
not used as a moving-cutoff uniformity theorem.

For a character psi at an induced modulus divisible by10, the general
twist has
\[
 F_\psi(s)=E_\psi(s)
       \prod_{r=0}^3 L(s,\psi\chi^r)^{\,i^{-r}/4}.              \tag{T.3}
\]
Only psi induced from a power chi^j yields a principal L-factor.
Its exponent is alpha_j. This lists the four and only four rough
exceptional branches.

## 3. Rough factors, complex conjugacy and positive invariant model

Put L=log X, y=exp(L^(2/5)), S={p<=y}, and
V=prod_{p<=y}(1-1/p). Delete the allowed small factors:
\[
 F_{j,R}(s)=F_j(s)
       \prod_{\substack{p\le y\\p\equiv3(10)}}(1-i^j p^{-s})
       =(s-1)^{-\alpha_j}h_{j,R}(s).
\]
For fixed k,
\[
 |h_{j,R}^{(k)}(1)|/k!
      \ll_k |h_{j,R}(1)|(C\log y)^k.                          \tag{T.4}
\]
As in the minus4 proof, a disk of radius c/log y bounds the ratio
h_{j,R}(s)/h_{j,R}(1) and Cauchy's formula gives(T.4).

AP Mertens with coherent Euler logarithms gives the sharper universal
asymptotic
\[
 h_{j,R}(1)\sim e^{-\gamma\alpha_j}(\log y)^{-\alpha_j}.
                                                                    \tag{T.5}
\]
Indeed write the logarithm as its first prime-power term plus the
absolutely convergent higher-power terms. The higher-power terms
cancel between log h_j(1) and the removed product; the difference
between the prime-harmonic and Abelian constants is gamma/4.
At minimum, only the following bounds are needed:
\[
 h_{0,R}(1)\asymp(\log y)^{-1/4},\quad
 h_{2,R}(1)\asymp(\log y)^{1/4},\quad
 |h_{1,R}(1)|=|h_{3,R}(1)|\asymp1.
\]
The latter two modulus claims also follow exactly from(T.2) for rough
series:
h_0,R(1)h_2,R(1)=F_0,R(2)->1 and
|h_1,R(1)|²=F_2,R(2)->1.

For a fixed sufficiently deep K define
\[
 M_j(t)=\sum_{k=0}^K\frac{h_{j,R}^{(k)}(1)}{k!}
             \frac{(\log t)^{\alpha_j-k-1}}
                  {\Gamma(\alpha_j-k)}.                      \tag{T.6}
\]
Powers of log t use its real positive logarithm on the power bulk.
Then M_0>0, M_2 is real, and M_3=bar M_1 there. More precisely,
\[
 |M_j(t)|/M_0(t)
 \ll \left(\frac{\log y}{\log t}\right)^{1/4-\Re\alpha_j}.
\]
Uniformly for X^a<=t<=X^C this gives
\[
 |M_1|+|M_3|\ll L^{-3/20}M_0,\qquad
 |M_2|\ll L^{-3/10}M_0.                                     \tag{T.7}
\]

Write n=A v with A=n_S. The correct model is
\[
 D(n)=\frac{f(A)}V
          \sum_{j=0}^3\overline{\chi(v)}^{\,j}M_j(v),\qquad
 D_0(n)=\frac{f(A)}V M_0(v).                                  \tag{T.8}
\]
**Lemma FourBranchInvariantModel.** For a sufficiently deep fixed
truncation, (T.8) is real, nonnegative, bounded by \(O(L^{-9/20})\),
and exactly invariant under all allowed primes at most \(y\).
It has arbitrary logarithmic rough centering against the periodic
and smooth tests used by the existing graph comparison.

The conjugation in(T.8) is essential: testing against chi^k must
select M_k, not M_{-k}. The order-two minus4 example cannot detect
this orientation error.

The combined numerator equals
M0+chi(v)^2 M2+2Re(bar chi(v) M1), so is real.
By(T.7) it is positive for large X and
\[
 |D-D_0|\ll L^{-3/20}D_0,\qquad D(n)\ll L^{-9/20}.              \tag{T.9}
\]
Choose smooth extensions outside the bulk with
M0≪h0,R(1)L^-3/4 globally, M3=bar M1, M2 real, and(T.7)
preserved. Such extensions exist by first cutting the secondary
terms to zero below a smaller bulk. This justifies the global bound.

For every allowed p<=y, n_S changes to pA and v does not change.
Consequently
\[
 f(pn)=f(n),\quad D(pn)=D(n),\quad (f-D)(pn)=(f-D)(n)           \tag{T.10}
\]
exactly, including repeated prime powers. This remains true although
multiplication rotates the residue of the whole n by3 modulo10.
For reflected endpoints use D(|n|); characters of signed rough
arguments add only fixed phases chi(-1)^j=(-1)^j.

## 4. Uniform centering and complex powers

At q<=L^B, every conductor prime lies below y. Character expansion
reduces all rough periodic tests to(T.3). The four characters
chi^j have precisely the Hankel densities(T.6); every other
character is regular at one and has a negligible smooth mean.

Complex powers require the FULL chosen logarithm of an L-function,
not just bounds for its absolute value and reciprocal.
Under the same fixed ZF, define L#=L off the principal character,
and L#=(s-1)L(s)/(s+1) at a principal character. Then on
23/24<=Re s<=2,
\[
 |\log L^\#(s,\psi)|\ll \log(q(2+|t|)),                        \tag{T.11}
\]
where the analytic logarithm is Euler-normalized on Re s>1.
To prove it, use a fixed disk centered at2+it contained in
Re s>11/12; convexity bounds Re log L# on its boundary, the
center Euler log is bounded by the same logarithm, and
Borel--Caratheodory gives(T.11) on a smaller disk.
This controls both the real and imaginary parts of the chosen log,
including possible winding; hence every fixed complex power in(T.3)
has polynomial conductor-height size.

Smooth Mellin inversion can now be moved to1-c/log y at height
exp(L^.1). The L-powers cost exp(O(L^.1)); finite Euler factors cost
a fixed log power; the left bulk integral saves exp(-cL^.6).
High derivatives make the Mellin tails negligible.
Hankel integration of(s-1)^(-alpha_j), with(T.4), gives(T.6)
and relative Taylor error O_K(L^(-3(K+1)/5)).
The secondary errors are additionally smaller than the positive
branch by(T.7).

On the comparison side the unrestricted rough population has only
its principal mean V. Thus testing(T.8) by chi^k leaves exactly
its j=k term, with coefficient M_k. The other terms have negligible
rough character means. Every four-branch mean therefore matches.
After choosing K after the fixed character/smoothness costs,
the centered rough sum is O_B(T L^-B) for every fixed B.
This is the moving-y proof; no unproved uniform use of a fixed-function
Selberg--Delange theorem is made.

## 5. Order-four prime-coset frame

**Lemma FourthRootCosetFrame.** The quotient-character ideal frame
on banks \(p\equiv3\pmod{10}\) has the same positivity, norm,
reproduction and comparison bounds as the corrected frame in
[Module I](PRIMITIVE_NORM_BINARY_CORRELATION_SUPPORT.md#i-the-fixed-ap-ideal-kernels).
The fourth-root slot phases introduce no coefficient or slot-count loss.

Choose banks p=3 mod10 in precisely the old small and big size ranges,
and delete primes dividing the target. Their density is1/phi(10)=1/4.
The physical graph input G allows arbitrary finite prime groups;
its separate interval-prime ideal theorem is not applied literally.

Quotient primitive characters of conductor<=Q by the relation
cond(psi bar eta)|10. The quotient is one of1,chi,chi²,chi³.
All bank primes exceed10Q² for large X, so all relevant character
values are units.

For matching input/comparison slot sets Lambda,
\[
 \chi(D_\Lambda)=\chi(X_\Lambda)=i^{|\Lambda|}.
\]
If psi=eta chi^e, the same deterministic phase acts on both:
\[
 \psi(D_\Lambda)=i^{e|\Lambda|}\eta(D_\Lambda),\qquad
 \psi(X_\Lambda)=i^{e|\Lambda|}\eta(X_\Lambda).                 \tag{T.12}
\]
For inequivalent representatives, their character sum on the
3mod10 fiber is a fixed unit phase times the zero sum on the
1mod10 kernel fiber. Siegel--Walfisz and partial summation prove
the same log-cell orthogonality as old(AP5).

Use the conjugated kernel on LOG-product cells
\[
 K_\Lambda(D,X)=\sum_{\nu_t\ge\zeta}
  \frac{1_{\log D\in I_t}1_{\log X\in I_t}}{\nu_t}
      \sum_{\eta\in{\cal R}_Q}\bar\eta(D)\eta(X),\qquad
 \nu_t=P(\log D_\Lambda\in I_t).
\]
It is positive, has norm<=2, and obeys
|K|<=Q²/zeta and sup_X E_D|K|<=Q².
The matching representative reproduces psi because its phase
i^(e|Lambda|) is identical on D and X. Grouping rational-phase
coefficients changes them by unit phases and cannot enlarge their
absolute coefficient mass. The old AP12/AP13 errors, symmetric-probe
cancellation and accuracy order are unchanged. No slot-count
factor is introduced even when the number of slots grows.

## 6. Target energy and envelope

For a nonexceptional character eta, a selected bank prime sum is
\[
 \frac14\sum_{r=0}^3 i^{-r}
       \sum_p \eta(p)\chi(p)^r p^{-1-it}.
\]
None of the four characters is principal. Their absolute coefficient
mass is one, so the same ZF prime estimate saves P^-1/24 with only
the old log costs. Cofactor/modulus ratios and all d-gcd layers
therefore remain exactly those of Module E.

For the four exceptional characters induced modulo10, write n=A v
after the exact gcd split. The record contributes eta(A/d), a unit
phase, and the residual test is a fixed-period character eta(v).
Section4 centers all four such tests, with arbitrary precision.
Restoring large-target-prime avoidance has the same Y/d error
and sum_(A:d|A)1/A=1/(dV). Summing d|Q with1/d gives the same
target-grid energy(G17), without a leftover target factor.

The leading positive model contains no subgroup projection:
D0≈lambda f(n_S), with
lambda=c L^-3/4/[alpha0 E_forbidden(y)], alpha0=phi(10)/10=2/5.
Here the forbidden-good set has density3/4.
The positive bank envelope has two-endpoint sieve dimension
2(1-rho)=3/2; the generated-group index is irrelevant.
All allowed bank primes avoid the target, so their joint record law
(W13), CRT counts and reconstruction are unchanged. Since
D=(1+O(L^-3/20))D0, the envelope also bounds all actual/model pair types.

For the model main term, finite bad-prime roots give exactly gamma10
from Section1. The completed good-prime factors give C_E J_E.
Moreover
\[
 |D*D-D0*D0|\ll L^{-3/20}(D0*D0).
\]
Thus none of the oscillatory mixed branches changes the leading constant.

## 7. Explicit positive error budget

Take phi=phi(10)=4,
delta=1/(8192phi)=1/32768 and
S_N=N J(N)/L^(3/2).
The same AP-Mertens masses are
(1/(8phi)+o(1))logL and(1/(15phi)+o(1))logL.
Choose b=floor(logL/(32phi))=floor(logL/128) banks of each type.
Each has mass1+o(1) and success probability tending e^-2>1/8.
Chernoff yields rare-bank probability
2L^(-1/(4096phi))=2L^-2delta; reconstruction mass is at most L.

Use epsilon=L^-delta, external phase rangeL³ and localization
reserveL^(8delta). Request L^-22 graph error per source dyad,
retaining L^-20 after at most L² dyads. Fix the ideal/probe data
and eventual mark length before choosing density depth.

For the resulting actual fixed family exponents Ccmp,Cw,Clow, use
K_*=[log(m_* !)+m_*log2]/(16phi), Cprod=2K_*+2, and choose
\[
 A_{\min}>C_{\rm cmp}+2C_w+40+8\delta,\qquad
 B_{\rm out}>2A_{\min}+C_{\rm prod}+4,
\]
\[
 T_E>2(C_{\rm cmp}+C_w+40+8\delta).
\]
The same target-dyad cofactor ratio
Y'/(qP)≫exp(cL^(1/3)-CL^(1/4))L^(-B_out-delta)
pays every fixed log cost. Choose K LAST with
3(K+1)/5>T_E+C_low+40, then the external Mellin derivative order
j with(3-delta)(j-1)>C_tail+40.

The complete error list relative to S_N is:

- Rare banks: O(L^-2delta).
- Graph/comparison/reconstruction:
  O(N L^(-19+8delta))/S_N
  =O(L^(-35/2+8delta)/J(N)), smaller than L^-delta.
- Nonexceptional primes, large smooth records, deleted target primes,
  bank-square records: exponential savings times fixed log powers.
- CRT endpoints: N^(.45+o(1)), negligible compared with S_N L^-delta.
- Positive rough logarithm variation: O(L^-1/2).
- Sum of complex and negative secondary branches in model evaluation:
  O(L^-3/20).
- Fundamental-lemma and Euler tails: O_B(L^-B) for every fixed B.
- Endpoint strips: O(epsilon)+N^(.45+o(1))/S_N.

Every error is O(S_N L^-delta) absolutely.
All four branches are bounded by a fixed multiple of M0, and D<=1
eventually for each fixed depth; changing the final depth changes
constants/onset, not an earlier comparison exponent. This completes
the stated pointwise asymptotic under the same ZF/G.

The actual new checks beyond the minus4 proof are the conjugated
four-branch orientation, realness via conjugacy, the smaller
L^-3/20 secondary gap, normalized analytic-log bounds for imaginary
L-powers, order-four slot phases, the factor at5 and sieve dimension3/2.

## 8. Optional absolute-Euler calibration of c

The c-product in Section1 is sufficient to specify the exact constant.
A finite-L version follows from(T.1). Put x=p^-s and
\[
 E_0(s)^4=
 \prod_{p\equiv3(10)}
    \frac{e^{-2\arctan x}}{(1-x)^3(1+x)}
 \prod_{p\equiv7(10)}
    \frac{1-x}{1+x}e^{2\arctan x}.
\]
Each factor is1+O(p^-2s), positive on the real line, so this product
converges absolutely at one. Classes1 and9 give factor1.
With Euler-normalized log and theta=Im log L(1,chi),
\[
 h_0(1)^4=\frac25\,
       \frac{e^{2\theta}}{L(1,\chi^2)}E_0(1)^4,\qquad
 c=h_0(1)/\Gamma(1/4).
\]
The literal characters omit2 and5. One may evaluate
L(1,chi)=(pi/10)[cot(pi/10)+i cot(3pi/10)] by the paired
periodic-series cotangent identity, so its usual argument is
theta=arctan(sqrt5-2). The real and imaginary periodic-series
parts are both positive for real s>1, confirming that this argument
matches the Euler-normalized continuation without winding.

No historical novelty or current unsolved-status claim is made here.
