# Cofactor drift in a coupled Gaussian–Eisenstein sieve

The unweighted mixed sumset offers a half-dimensional sieve toward
simultaneous Gaussian–Eisenstein prime restrictions. Its fixed-prime
Haar density does not, however, give a valid main term for coefficients
of power size. The omitted cofactor-size correction accumulates to a
positive proportion of the candidate mass. This note proves that
obstruction, together with an unconditional bound for the two-large-prime
pairs left after sieving. The full binary C12+C12 problem remains open.

## 1. Exact mixed identities and the local model

Let
\[
\begin{split}
 C_4&=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod4\},\\
 C_3&=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod3\},\\
 C&=C_{12}=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod{12}\}.
\end{split}
\]
Each set includes 1. Throughout, target integers satisfy
\(N\equiv2\pmod{12}\), and
\[
 R_0(N)=\#\{a+b=N:a\in C_4,\ b\in C_3\}.
 \tag{1}
\]
Representations are ordered. Let \(\mathcal P\) consist of the primes
congruent to 5 or 7 modulo 12, and put
\[
 P(z)=\prod_{\substack{p<z\\p\in\mathcal P}}p.
\]
For squarefree \(d\) supported on \(\mathcal P\), write uniquely
\(d=d_5d_7\), separating its prime factors by their residue class, and set
\[
 A_d(N)=\#\{a+b=N:a\in C_4,\ b\in C_3,\ d\mid ab\}.
 \tag{2}
\]

**Lemma CoupledCoefficientIdentity.** One has the exact identity
\[
 A_d(N)=\#\{d_5x+d_7y=N:x\in C_4,\ y\in C_3\}.
 \tag{3}
\]
If \((d,N)>1\), both sides are zero. For the corresponding local
Haar model the event density is
\[
 g_N(p)=
 \begin{cases}(p-1)^{-1},&p\in\mathcal P,\ p\nmid N,\\
 0,&\text{otherwise},\end{cases}
 \qquad
 g_N(d)=\frac{\mathbf1_{(d,N)=1}}{\varphi(d)}.
 \tag{4}
\]
This local calculation does not assert a global unweighted marginal
asymptotic.

A prime congruent to 5 cannot divide the second component, and a prime
congruent to 7 cannot divide the first. Division by the indicated
allowed prime factors preserves the two membership indicators exactly,
which proves (3). If a prime of \(d\) also divides \(N\), the equation
forces it into the forbidden component, proving the zero assertion.
For a prime congruent to 5 and a unit target, the first local component
has arbitrary valuation, while the second must be a unit. The total
Haar mass is \(1-1/p\), and the first-component divisibility event has
mass \(1/p\), giving (4). At a nonunit target that event is impossible.
The other residue class is symmetric; the joint local model follows
by the Chinese remainder theorem.

The allowed residues of \(C_4\) are 1 and 5 modulo 12, and those of
\(C_3\) are 1 and 7. Their sum is 2 only for the pair (1,1). Consequently
\[
 a\equiv b\equiv1\pmod{12},\qquad
 \Omega_5(a)\equiv\Omega_7(b)\equiv0\pmod2,
 \tag{5}
\]
where \(\Omega_j\) counts prime factors in residue class \(j\), with
multiplicity. These are separate exact parities, not independence
assumptions.

## 2. Cofactor drift and failure at every power level

Fix \(0<\alpha<\beta<1\), and write
\[
\begin{split}
 \mathcal P_X&=\{p:X^\alpha\le p\le X^\beta,\ p\equiv5\pmod{12}\},\\
 \mathcal M(X)&=\sum_{\substack{X<N\le2X\\N\equiv2\ (12)}}R_0(N),\\
 r_p(N)&=A_p(N)-g_N(p)R_0(N).
\end{split}
\]

**Theorem CoupledCofactorDrift.** There is a constant \(C>0\) such that
\[
 \mathcal M(X)\sim C\frac{X^2}{\log X}.
 \tag{6}
\]
Uniformly for \(p\in\mathcal P_X\),
\[
 \sum_{\substack{X<N\le2X\\N\equiv2\ (12)}}A_p(N)
 =\frac{\mathcal M(X)}p
       \left(1-\frac{\log p}{\log X}\right)^{-1/2}(1+o(1)).
 \tag{7}
\]
In particular,
\[
 \sum_{\substack{X<N\le2X\\N\equiv2\ (12)}}
       \sum_{p\in\mathcal P_X}r_p(N)
 =\bigl(K_{\alpha,\beta}+o(1)\bigr)\mathcal M(X),
 \tag{8}
\]
where
\[
 K_{\alpha,\beta}
 =\frac14\int_\alpha^\beta
       \frac{(1-t)^{-1/2}-1}{t}\,dt
 =\frac12\log\frac{1+\sqrt{1-\alpha}}{1+\sqrt{1-\beta}}>0.
 \tag{9}
\]
All limits have \(\alpha,\beta\) fixed. The analogous assertions hold
for primes congruent to 7, acting on the second component.

**Corollary CoupledBarePowerLevelFailure.** Fix \(\rho,\zeta>0\) and
\(B\ge0\), and put
\[
 D_N=N^\rho(\log N)^{-B},\qquad z_N=N^\zeta.
\]
For every \(0<\alpha<\beta<\min(\rho,\zeta,1)\),
\[
\begin{split}
 &\sum_{\substack{X<N\le2X\\N\equiv2\ (12)}}
 \sum_{\substack{d<D_N\\d\mid P(z_N)}}
      |A_d(N)-g_N(d)R_0(N)|\\
 &\hspace{20mm}\ge
 \bigl(K_{\alpha,\beta}+o(1)\bigr)\mathcal M(X).
\end{split}
 \tag{10}
\]
Thus the inner absolute remainder cannot be \(o(R_0(N))\) uniformly
for all sufficiently large admissible targets, at any fixed positive
power level. This refutes a bare critical-level hypothesis with an
arbitrary negative power of \(\log N\), or with error
\(o(R_0(N)\mathcal V_N(z_N))\), where
\(\mathcal V_N(z)=\prod_{p<z}(1-g_N(p))\le1\).
The obstruction already uses prime divisors, without composite ones.
It does not refute the desired C12+C12 theorem or a distribution
statement with a corrected main term.

### 2.1. Fixed residue means

Put \(f_i=\mathbf1_{C_i}\). There are constants \(c_4,c_3>0\) such that
\[
 \sum_{n\le x}f_i(n)
 =c_i\frac{x}{\sqrt{\log x}}
       +O\!\left(\frac{x}{(\log x)^{3/2}}\right),\qquad i=3,4.
 \tag{11}
\]
Moreover, for every fixed \(A>0\),
\[
 \sum_{n\le x}f_4(n)\chi_{-3}(n),\quad
 \sum_{n\le x}f_3(n)\chi_{-4}(n)
       \ll_A x(\log x)^{-A}.
 \tag{12}
\]
To verify the input precisely, apply Granville–Koukoulopoulos,
Theorem 1, with divisor bound \(\tau_1=1\). The fixed-modulus prime
number theorem supplies its prime-sum hypothesis with every fixed
logarithmic error exponent. The untwisted prime means are \(1/2\),
and the twisted means are zero. In the latter case every reciprocal
\(\Gamma(-j)\) vanishes. In the former case the positive Euler factor,
divided by \(\Gamma(1/2)\), is \(c_i\). This proves (11)–(12).

Since the two possible residues of each set have opposite values of
the twisting character, these estimates give
\[
\begin{split}
 U_1(x)&:=\#\{n\le x:n\in C_4,\ n\equiv1\ (12)\}
      \sim\frac{c_4}2\frac{x}{\sqrt{\log x}},\\
 U_5(x)&:=\#\{n\le x:n\in C_4,\ n\equiv5\ (12)\}
      \sim\frac{c_4}2\frac{x}{\sqrt{\log x}},\\
 V_1(x)&:=\#\{n\le x:n\in C_3,\ n\equiv1\ (12)\}
      \sim\frac{c_3}2\frac{x}{\sqrt{\log x}}.
\end{split}
 \tag{13}
\]
Only modulus 12 is used. No progression estimate with modulus \(p\)
is needed.

### 2.2. Uniform product measures on the additive region

Let \(L=\log X\). For \(p\in\mathcal P_X\), define finite nonnegative
measures on \([0,2]\) by their distribution functions
\[
 \mu_{p,X}([0,u])=
   \frac{p\sqrt{\log(X/p)}}X U_5(Xu/p),\qquad
 \nu_X([0,v])=\frac{\sqrt L}X V_1(Xv).
 \tag{14}
\]
For every fixed \(u>0\), (13) gives
\[
 \mu_{p,X}([0,u])=\frac{c_4}2u+o(1)
 \tag{15}
\]
uniformly in \(p\): indeed \(X/p\ge X^{1-\beta}\) tends to infinity,
and \(\log(Xu/p)/\log(X/p)\to1\) uniformly. The corresponding limit
for \(\nu_X\) is \(c_3v/2\). The convergence of distribution functions
is uniform on all of \([0,2]\). To prove this, first restrict to
\([\delta,2]\) for fixed \(\delta>0\), where (13) is uniform, and then
bound the mass of \([0,\delta]\) by its value at \(\delta\), using
monotonicity. Let \(X\to\infty\) before \(\delta\to0\).

It follows that
\[
 (\mu_{p,X}\otimes\nu_X)
       \{(u,v)\in[0,2]^2:1<u+v\le2\}
 \longrightarrow \frac{c_4c_3}4\cdot\frac32
 \tag{16}
\]
uniformly in \(p\). Here is an explicit domain justification. Partition
\([0,2]^2\) into a fixed grid of squares. Product masses converge
uniformly on each square by (15). The squares meeting either boundary
line have total limiting area tending to zero with the grid size.
Their masses are bounded by those of their union, which has that same
limiting bound. First take \(X\to\infty\) with the grid fixed, then
refine the grid. The region has area \((2^2-1^2)/2=3/2\).

For \(p\equiv5\pmod{12}\), the substitution \(a=pm\) forces
\(m\equiv5\pmod{12}\). Multiplication by \(p\) preserves membership
in \(C_4\) exactly. The unnormalized left side of (16) is therefore
exactly the left side of (7), giving
\[
 \sum_{\substack{X<N\le2X\\N\equiv2\ (12)}}A_p(N)
 =\left(\frac{3c_4c_3}8+o(1)\right)
       \frac{X^2}{p\sqrt{L\log(X/p)}}
 \tag{17}
\]
uniformly on the whole power band.

Apply the same argument with \(U_1(Xu)\) and \(V_1(Xv)\), using the
normalization \(\sqrt L/X\) for both. By (5), precisely the residue
pair (1,1) contributes. Thus
\[
 \mathcal M(X)=\left(\frac{3c_4c_3}8+o(1)\right)\frac{X^2}L.
 \tag{18}
\]
This proves (6)–(7), with \(C=3c_4c_3/8\), without a pointwise binary
asymptotic. The prime class 7 case uses the equal mean of residue 7
in \(C_3\), obtained from the same argument as (13).

### 2.3. The accumulated drift and target-prime correction

Set \(\Phi(t)=(1-t)^{-1/2}-1\). Fixed-class prime Mertens and partial
summation give
\[
 \sum_{p\in\mathcal P_X}\frac{\Phi(\log p/L)}p
       \longrightarrow\frac14\int_\alpha^\beta\frac{\Phi(t)}t\,dt.
 \tag{19}
\]
The reciprocal prime sum on the band is bounded. Hence the uniform
\(o(1)\) in (7) stays \(o(1)\) after summation. Also
\(\sum_{p\in\mathcal P_X}(1/(p-1)-1/p)=o(1)\). Consequently
\[
 \sum_{p\in\mathcal P_X}\sum_N A_p(N)
       -\mathcal M(X)\sum_{p\in\mathcal P_X}\frac1{p-1}
       =\bigl(K_{\alpha,\beta}+o(1)\bigr)\mathcal M(X),
 \tag{20}
\]
where each target sum retains the range and congruence of (6).

Replacing \(1/(p-1)\) by \(g_N(p)\) adds a nonnegative term. Using only
\(R_0(N)\le N\), this term is negligible:
\[
\begin{split}
 &\sum_{p\in\mathcal P_X}\frac1{p-1}
        \sum_{\substack{X<N\le2X\\N\equiv2\ (12)\\p\mid N}}R_0(N)\\
 &\qquad\ll\sum_{p\in\mathcal P_X}
                 \left(\frac{X^2}{p^2}+\frac Xp\right)
 \ll X^{2-\alpha}+X=o(\mathcal M(X)).
\end{split}
 \tag{21}
\]
This proves (8). Substituting \(u=\sqrt{1-t}\) in the integral gives
the closed form (9).

For the corollary, the chosen prime band lies below \(D_N,z_N\)
uniformly for \(N\in(X,2X]\), once \(X\) is large. The absolute sum
in (10) dominates the signed prime sum (8). In particular each
sufficiently large dyadic interval contains an admissible target with
\[
 \sum_{\substack{d<D_N\\d\mid P(z_N)}}
       |A_d(N)-g_N(d)R_0(N)|
 \ge\frac12K_{\alpha,\beta}R_0(N).
 \tag{22}
\]
Otherwise summing would contradict (10). Targets with zero candidate
mass contribute nothing. This proves failure of the proposed uniform
bare remainder at every positive power level. No conclusion about
logarithmic or other subpower levels is made.

## 3. Candidate mass and the natural sieve scale

For primes at least 5 write \(j_p=(p-1)/(p-2)\), and define
\[
 J_{11}(N)=\prod_{\substack{p\mid N\\p\equiv11\ (12)}}j_p,
 \qquad
 J_{12}(N)=\prod_{\substack{p\mid N\\p\equiv5,7,11\ (12)}}j_p.
\]
The mixed specialization (R1) of
[SubgroupSemigroupScale](ASYMMETRIC_SUBGROUP_SPLITTING.md) gives
\[
 R_0(N)\asymp\frac{N}{\log N}J_{11}(N)
 \tag{23}
\]
for every sufficiently large admissible target. This is an order
comparison, not an identified unweighted leading constant or a marginal
distribution theorem. The averaged drift proof above did not use (23).

Fix \(0<\epsilon<1/100\), and put
\[
 z=N^{1/2-\epsilon},\qquad w=N^{1/10},\qquad
 Y=N^{2\epsilon}<w,\qquad
 \mathcal V_N(z)=\prod_{p<z}(1-g_N(p)).
 \tag{24}
\]
The baseline product has dimension \(1/2\): the density is
\(1/p+O(p^{-2})\) on a set of prime density \(1/2\). Precisely,
\[
 \mathcal V_N(z)=\mathcal V_0(z)
      \prod_{\substack{p\mid N\\p\in\mathcal P,\ p<z}}j_p,
 \qquad \mathcal V_0(z)\asymp(\log z)^{-1/2}.
 \tag{25}
\]
There are at most two distinct target primes in \(\mathcal P\) at
least \(z\), since \(z^3>N\); their omitted correction factors have a
bounded product. Combining (23)–(25) gives
\[
 R_0(N)\mathcal V_N(z)
       \asymp\frac{N J_{12}(N)}{(\log N)^{3/2}}.
 \tag{26}
\]
The comparison constants are independent of \(\epsilon\) in its stated
range; the onset may depend on \(\epsilon\). Thus (26) is a useful
scale for a prospective sieve. Equations (8)–(10) show why it cannot be
obtained from a bare power-level remainder using only (4).

## 4. An unconditional two-large-prime bound

Let \(T_z(N)\) consist of the candidate pairs in (1) with
\((ab,P(z))=1\), and let \(B_z(N)\) count those pairs for which at least
one component is outside \(C\).

**Theorem CoupledTwoLargePrimes.** There is a constant \(C_0\),
independent of \(\epsilon\), such that for every fixed
\(0<\epsilon<1/100\) and all sufficiently large admissible \(N\),
\[
 B_z(N)\le C_0\epsilon^{5/4}
        \frac{N J_{12}(N)}{(\log N)^{3/2}}
 \ll\epsilon^{5/4}R_0(N)\mathcal V_N(z).
 \tag{27}
\]
Every target valuation is allowed. This is a bad-mass upper bound;
it does not supply a lower bound for \(T_z(N)\).

### 4.1. Exact factorization and retention of the other sieve

If the first component is bad, (5) and \(z^3>N\) imply the unique
factorization
\[
 a=cpq,\qquad c\in C,\quad c\le Y,\quad
 z\le p\le q,\quad p,q\equiv5\pmod{12},
 \tag{28}
\]
where \(p,q\) are primes and \(p=q\) is permitted. The two prime
factors are counted with multiplicity. They cannot divide \(N\),
since the other component belongs to \(C_3\). A bad second component
has the symmetric factorization with prime class 7.

It is essential to retain the sieve condition on the other component.
Define completely multiplicative indicators on positive integers by
\[
\begin{split}
 f_{3,z}(n)&=\mathbf1_{n\in C_3,\ (n,P_7(z))=1},\\
 f_{4,z}(n)&=\mathbf1_{n\in C_4,\ (n,P_5(z))=1},
 \qquad
 P_j(z)=\prod_{\substack{\ell<z\\\ell\equiv j\ (12)}}\ell.
\end{split}
 \tag{29}
\]
Below \(z\) each permits only primes congruent to 1 modulo 12. At
larger primes each remains an indicator, with growth constants uniform
in \(z\).

For fixed \(c,p\) from (28), we claim
\[
 \sum_{\substack{z\le q<N/(cp)\\q\text{ prime},\ q\equiv5\ (12)}}
       f_{3,z}(N-cpq)
 \ll\frac{N}{cp(\log N)^{7/4}}J_w(N),
 \qquad
 J_w(N)=\prod_{\substack{\ell\mid N,\ \ell\le w\\
                         \ell\equiv5,7,11\ (12)}}j_\ell.
 \tag{30}
\]
The endpoint may equivalently be rounded to include only positive
values of the second component. The symmetric claim uses prime class
7 and \(f_{4,z}\). The constants in both estimates are independent of
\(\epsilon,c,p,N\).

### 4.2. Corrected Henriot input and exact local controls

Put \(h=(c,N)\in C\). Multiplication by this content leaves the other
indicator unchanged:
\[
 f_{3,z}(hn)=f_{3,z}(n)\qquad(n\ge1).
 \tag{31}
\]
Write \(q=12t+5\), and use the primitive linear polynomials
\[
 Q_1(t)=12t+5,\qquad
 Q_2(t)=N/h-(cp/h)(12t+5).
 \tag{32}
\]
They are primitive because \((N/h,cp/h)=1\) and
\(Q_2(t)\equiv1\pmod{12}\). Their product has degree 2 and coefficient
norm \(O(N)\).

Majorize the prime condition on \(q\) by the multiplicative rough
indicator
\[
 f_w(n)=\mathbf1_{(n,\prod_{\ell\le w}\ell)=1},\qquad
 F(u,v)=f_w(u)f_{3,z}(v).
 \tag{33}
\]
A prime \(q\ge z>w\) has \(f_w(q)=1\). No prime indicator is treated
as a multiplicative function. Both factors in (33) are completely
multiplicative indicators; the constants in Henriot's class
\(\mathcal M_2(A,B,\nu)\) can be taken to be \(A=B=1\), uniformly in
the moving parameters.

Cover the interval from \(z\) to the upper endpoint by geometric
intervals \([2T/3,T]\), allowing the last one to extend below \(z\)
for an upper bound. Take the associated half-open integer intervals in
\(t\) with upper endpoint \(\lfloor(T-5)/12\rfloor\). In the source
notation they have
\[
 x\asymp T,\quad y\asymp T/3\le x,\quad
 T\gg z\ge N^{0.48},
 \tag{34}
\]
for sufficiently large \(N\). Their total lengths are \(O(N/(cp))\).
All values of \(Q_2\) used are positive; equality \(cpq=N\) cannot
occur since \(p\nmid N\).

Apply New Theorem 5 of Henriot's 2014 erratum with the fixed parameters
\[
 \alpha_H=1/2,\qquad \Delta_H=1/4,\qquad \nu=1/10000.
\]
Equation (34) implies \(x\ge C\|Q_1Q_2\|^{\Delta_H}\) and
\(x^{\alpha_H}<y\le x\). The required growth restriction is
\(\nu<\alpha_H/[50\cdot2(2+\Delta_H^{-1})]\), which holds.
The source constants therefore do not depend on \(\epsilon,c,p,N\).

The corrected exact-valuation local factor admits a useful elementary
upper bound. For a prime \(\ell\ge5\), let
\(P_\ell(\nu_1,\nu_2)\) be the Haar probability of the exact
valuation event, imposing unit conditions on the zero coordinates
when the valuation vector is nonzero. In the source coefficient at
\((0,0)\) there is no condition: its value is 1. The other events
partition the complement of the unit-unit event. If \(\rho_Q(\ell)\)
counts roots of the polynomial product modulo \(\ell\), then
\[
\begin{split}
 &(1-\rho_Q(\ell)/\ell)
       \left(1+\sum_{\boldsymbol\nu\ne0}
          F(\ell^{\nu_1},\ell^{\nu_2})
                    P_\ell(\boldsymbol\nu)\right)\\
 &\hspace{20mm}\le
       \int_{\mathbb Z_\ell}F(Q_1(t),Q_2(t))\,dt.
\end{split}
 \tag{35}
\]
In this integral the local valuation profiles of \(F\) are used.
Indeed the constant part on the left equals the true unit-unit mass,
and its remaining nonnegative terms are decreased by the extra factor
\(1-\rho_Q(\ell)/\ell\le1\).

When enlarging the nonnegative source coefficient sum, remove the
product cutoff \(n_1n_2\le x\) but **retain the restriction that every
prime factor is at most \(x\)**. It then factors into a finite Euler
product over \(\ell\le x\), so that the source prefactors combine
with (35). Removing this prime support as well would be unjustified.
For \(w<\ell\le x\), the integral in (35) is at most 1 because both
weights are indicators. All such factors can be dropped, including
any exceptional leading or target primes in that range.

For \(5\le\ell\le w\), the remaining integrals are exactly:

- If \(\ell\equiv1\pmod{12}\), the other component permits every
  valuation, so only \(Q_1\) must be a unit. The mass is \(1-1/\ell\),
  independently of \(v_\ell(c)\) and \(v_\ell(N)\), after (31).
- If \(\ell\equiv5,7,11\pmod{12}\), both components must be units.
  The slope \(cp/h\) is a unit: its primes come from good primes below
  \(w\), or the outer prime \(p\ge z>w\). The roots are distinct for
  \(\ell\nmid N\), giving \(1-2/\ell\). They coincide for
  \(\ell\mid N\), giving \(1-1/\ell\), at every positive target
  valuation. The target correction is exactly \(j_\ell\).

At 2 and 3 both forms are units; the finitely many source factors are
bounded by an absolute constant. The generic small-prime product is
\[
 \prod_{\substack{\ell\le w\\\ell\equiv1\ (12)}}(1-1/\ell)
 \prod_{\substack{\ell\le w\\\ell\equiv5,7,11\ (12)}}(1-2/\ell)
 \ll(\log w)^{-7/4}\asymp(\log N)^{-7/4}.
 \tag{36}
\]
Its first-order prime average is \(-1/4-2(3/4)=-7/4\).
Fixed-modulus prime Mertens supplies the constants. The exceptional
small target primes contribute \(J_w(N)\). Summing the source
interval bounds proves (30).

For the symmetric case, use \(Q_1(t)=12t+7\), the corresponding
\(Q_2\), and \(F=f_w f_{4,z}\). The content identity (31) is unchanged,
\(Q_2\equiv1\pmod{12}\), and all local factors are identical.

### 4.3. Summing outer primes and good cofactors

Let \(B_5\) count pairs with a bad first component. Equations (28)
and (30) give
\[
 B_5\ll\frac{NJ_w(N)}{(\log N)^{7/4}}
       \sum_{\substack{c\le Y\\c\in C}}\frac1c
       \sum_{\substack{z\le p\le\sqrt N\\p\equiv5\ (12)}}\frac1p.
 \tag{37}
\]
For fixed \(\epsilon\) in the stated range, the final prime sum is
\[
 \frac14\log\frac{1/2}{1/2-\epsilon}+o(1)\ll\epsilon,
 \tag{38}
\]
with an absolute constant once \(N\) is sufficiently large for that
\(\epsilon\). Positivity and the good-prime Euler product give
\[
 \sum_{\substack{c\le Y\\c\in C}}\frac1c
 \le\prod_{\substack{\ell\le Y\\\ell\equiv1\ (12)}}(1-1/\ell)^{-1}
 \ll(\log(2Y))^{1/4}
 \ll\epsilon^{1/4}(\log N)^{1/4}.
 \tag{39}
\]
Thus \(B_5\ll\epsilon^{5/4}NJ_w(N)/(\log N)^{3/2}\). The same
bound holds for \(B_7\). Pairs bad on both sides are covered because
the other component's sieve was retained in (29); no intersection
estimate is needed. Since \(B_z\le B_5+B_7\) and \(J_w\le J_{12}\),
this proves (27). Henriot's parameters, the Mertens constants, and
(23) are fixed independently of \(\epsilon\).

## 5. Scope and the remaining binary problem

The exact coefficient identity (3), the local Haar model (4), and the
unconditional two-large-prime estimate (27) all remain valid. What
fails is substituting \(g_N(d)R_0(N)\) as a uniform main term for a
critical or other positive-power absolute distribution estimate.
Extracting a divisor changes the remaining logarithmic density;
(7) quantifies that change, and (8) proves that the omitted drift is
not an error tending to zero relative to \(R_0\).

A corrected approach must retain coefficient-size dependence or
reorganize the sieve to account for it. A sufficient lower-sieve
framework with that correction is not supplied here. The former bare
input is disproved, rather than an open analytical hypothesis. No
failure of C12+C12, of a differently normalized distribution theorem,
or of a separate signed cancellation mechanism follows from this
obstruction. The all-target binary theorem, effective cutoff and
balance refinement remain open in this route.

The unweighted density here is \(1/p+O(p^{-2})\) on the union of two
prime classes, giving dimension \(1/2\). The tau-weighted mixed route
has density \(2/p+O(p^{-2})\), giving dimension 1. Changing weights
changes the distribution problem. Neither tau-weighted marginal
estimates nor the order comparison (23) establish unweighted
coefficient asymptotics.

## 6. Dependencies, primary sources and verification scope

The proof dependencies are:

- CoupledCoefficientIdentity: exact membership and local Haar counting.
- CoupledCofactorDrift and CoupledBarePowerLevelFailure: the one-variable
  means from Granville–Koukoulopoulos Theorem 1, fixed-class prime
  Mertens, and the product-measure argument above. No all-target mixed
  theorem is needed for these two results.
- CoupledTwoLargePrimes: corrected Henriot New Theorem 5, the exact
  parity (5), and fixed-class Mertens. Its expression in terms of
  \(R_0\mathcal V_N\) additionally uses the mixed
  [SubgroupSemigroupScale specialization (R1)](ASYMMETRIC_SUBGROUP_SPLITTING.md).

Primary sources:

1. A. Granville and D. Koukoulopoulos, *Beyond the LSD method for the
   partial sums of multiplicative functions*, Ramanujan J. **49**
   (2019), 287–319, [DOI](https://doi.org/10.1007/s11139-018-0119-3).
   Theorem 1, equations (1.2), (1.4), (1.5), is checked in the exact
   [arXiv v3](https://arxiv.org/abs/1710.01389v3), PDF pp.2–3, and the
   [published primary text](https://pmc.ncbi.nlm.nih.gov/articles/PMC6555439/).
   The [retained v3 PDF](sources/granville-koukoulopoulos-1710.01389v3.pdf)
   has SHA-256 `8dd052fd0a3cc0167058a678450d0c70f999a7ed42aa8542025f840570d58187`.
2. K. Henriot, *Nair–Tenenbaum uniform with respect to the discriminant —
   Erratum*, Math. Proc. Cambridge Philos. Soc. **157** (2014), 375–377,
   [DOI](https://doi.org/10.1017/S0305004114000280).
   [New Theorem 5, printed p.377](sources/henriot-erratum-2014.pdf),
   uses the corrected exact-valuation function defined on p.375.
   Its SHA-256 is `bd56f8f411457549c5cc53f9b24bbdec68544dca09cfe0836e35cc1b25677d6f`.
3. D. Loughran and L. Matthiesen, *Frobenian multiplicative functions
   and rational points in fibrations*,
   [arXiv:1904.12845v2](https://arxiv.org/abs/1904.12845v2), Lemma 2.4,
   supplies fixed-class prime sums and Mertens estimates.
   The [retained PDF](sources/loughran-matthiesen-1904.12845v2.pdf)
   has SHA-256 `ffbb0c2e6a6a39e090ebd6102f2a3b8936208df33c11bdf7c892744756c10203`.

Separate internal mathematical reviews reconstructed the drift argument
and the corrected-Henriot upper estimate, including the uniform
prime-band limits, zero-mean Gamma factors, and finite Euler-product
restriction. These are internal reviews, not external human review,
formal verification or a priority determination.

The retained
[stdlib diagnostic checker](eisenstein-restricted-prime-sums/verify_coupled_gaussian_eisenstein_sieve.py)
and [receipt](eisenstein-restricted-prime-sums/coupled_gaussian_eisenstein_sieve_checks.json)
replay finite coefficient, parity, content, local-density and scaled-domain
identities with negative controls. They do not certify an analytic
asymptotic, its onset, a distribution theorem, or eventual binary coverage.
No manuscript, submission or public release is made by this note.
