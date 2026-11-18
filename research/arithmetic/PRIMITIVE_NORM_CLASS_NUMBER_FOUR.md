# Primitive x²+14y² sums: genus conditions and a C4 class correction

Let
\[
 P=\{x^2+14y^2>0:\gcd(x,y)=1\},\qquad
 R_P(N)=\#\{(a,b)\in P^2:a+b=N\}.
\]
These are ordered pairs of values, without coordinate multiplicity. Put
\[
 \chi=\chi_{-56},\quad \psi=\chi_8,\qquad
 J_{-56}(N)=
 \prod_{\substack{\ell\mid N\\\ell\nmid14,\ \chi(\ell)=-1}}
                                \frac{\ell-1}{\ell-2}.
\]

**Theorem CN4.** There are constants \(c_1,c_2>0,N_0\) such that every
\(N\ge N_0\) with \(49\nmid N\) satisfies
\[
 c_1\frac{N J_{-56}(N)}{\log N}
 \le R_P(N)\le c_2\frac{N J_{-56}(N)}{\log N}.              \tag{CN4.1}
\]
Every positive multiple of 49 has count zero. The condition \(49\nmid N\)
is the complete primitive local sum condition; there is no dyadic
exclusion, including at arbitrarily large 2-adic target valuation.
No numerical cutoff, complete unrestricted exception list or leading
constant is supplied.

Here the class group is \(C_4\). Genus conditions alone do not characterize
individual primitive principal norms. On a sufficient additive subset,
two fixed primes correct the remaining subgroup \(2C_4\) while a
periodic character controls the genus quotient. The analytic proof
combines the explicit finite-genus transfer in
[the class-number-two case](PRIMITIVE_NORM_CLASS_NUMBER_TWO.md) with
the exact-prime profiles in
[the class-number-three case](PRIMITIVE_NORM_CLASS_NUMBER_THREE.md).
It requires no new additive class-field distribution estimate. This is
an internally audited deduction; no external review, formal verification
or historical priority is claimed.

## The class group and its genus character

Write \(\vartheta=\sqrt{-14}\), so
\(\mathcal O_K=\mathbb Z[\vartheta]\), with discriminant \(-56\).
Minkowski's bound \((2/\pi)\sqrt{56}<5\) gives ideal-class representatives
of norms 1,2,3,4. Let
\[
 \mathfrak p_2=(2,\vartheta),\qquad
 \mathfrak p_3=(3,\vartheta-1).
\]
There is no element of norm 2, so the ramified class
\([\mathfrak p_2]\) is nontrivial of order two. The element
\(2+\vartheta\) has norm 18 and
\[
                 (2+\vartheta)=\mathfrak p_2\mathfrak p_3^2.
                                                               \tag{CN4.2}
\]
Its reduction vanishes at \(\vartheta=1\pmod3\) and not at the
conjugate root; its norm has 2-adic valuation one. Thus
\(g=[\mathfrak p_3]\) has order four and \(2g=[\mathfrak p_2]\).
The representatives of norm at most four leave only
\(0,2g,g,-g\), since the norm-four ideal is
\(\mathfrak p_2^2=(2)\). The full class group is \(C_4\).
Also
\[
       [\mathfrak p_7]=2g,\qquad
       (\vartheta)=\mathfrak p_2\mathfrak p_7.              \tag{CN4.3}
\]
The primary Minkowski and ideal/form correspondence statements are
Milne, *Algebraic Number Theory*, v3.08, Theorems 4.3 and 4.29;
the exact retained edition and locators are in the source ledger.

The forms for these four classes are
\[
 Q_0=x^2+14y^2,\quad Q_2=2x^2+7y^2,\quad
 Q_\pm=3x^2\pm2xy+5y^2.
\]
For example, the bases of \(\mathfrak p_2\) and
\(\mathfrak p_3\) give \(Q_2\) and one of \(Q_\pm\) after dividing
the element norm by the ideal norm. The inverse forms have the same
unoriented values. For odd represented primes the first two forms have
\(\psi=+1\), and \(Q_\pm\) have \(\psi=-1\), as direct reduction
modulo eight verifies. Thus \(\psi\) distinguishes class parity.

Equivalently, this character is well defined on odd ideal classes.
If two odd ideals differ by \((\alpha)\), then \(\alpha\) is a unit
in \(\mathcal O_K\otimes\mathbb Z_2\). Its norm is \(1\) or \(7\)
modulo eight, so the norm ratio has \(\psi=+1\).
The character is multiplicative, and \(\psi(3)=-1\) on the generator
\(g\). Its kernel is \(2G=\{0,2g\}\).
Every class has an odd representative, since \(g\) does.

## Exact individual principal norm criterion

For a positive integer \(n\), first require no inert prime factor and
write
\[
                    n=2^e7^f c,\qquad e,f\in\{0,1\},
\]
where \(c\) is supported on unramified split primes.
Let \(s\) be the number of distinct split primes dividing \(c\).
For an unoriented split prime class, call it principal, class two,
or odd class according to whether its classes are \(0\), \(2g\),
or \(\{g,-g\}\).

Let \(r\) count the odd-class prime powers in \(c\) whose exponents
are odd. Define the parity
\[
 b=e+f+
   \sum_{\substack{q\mid c\\q\text{ of class }2g}}v_q(c)
   +\#\{q\mid c:[q]\text{ odd},\ v_q(c)\equiv2\pmod4\}
                      \pmod2.                             \tag{CN4.4}
\]

**Lemma CN4Individual.** The integer \(n\) is in \(P\) if and only if
the prime-factor and ramified-valuation conditions just stated hold and
\[
                 r\text{ is even},\qquad
                 r\ge2\text{ or }b=0.                     \tag{CN4.5}
\]
Its exact signed coordinate representation count is
\[
 R^*(n)=
 \begin{cases}
 2^s,&r\ge2\text{ even},\\
 2^{s+1},&r=0,\ b=0,\\
 0,&\text{otherwise}.
 \end{cases}                                               \tag{CN4.6}
\]
The count is zero as well if an inert factor or ramified square occurs.

**Proof.** Primitive ideal factorization orients each entire split prime
power to one conjugate ideal. An odd-class power of odd exponent
contributes either \(g\) or \(-g\); an exponent two modulo four
contributes \(2g\); exponent zero modulo four contributes zero.
Class-two powers contribute \(2g\) times their exponent, and the
ramified prefix contributes \(2g(e+f)\). This gives the fixed class
\(2gb\) plus \(r\) signs of \(g\).

Odd \(r\) cannot give the principal class. For even \(r\ge2\), exactly
\(2^{r-1}\) sign choices give either specified element of \(2G\):
already two signs give classes zero and \(2g\), and the count follows
by induction, or the character filter on \(C_4\).
For \(r=0\), the fixed class must be zero. The other \(s-r\)
split powers each have two distinct ideal orientations of the same
class. Every principal ideal has two generators. This proves
(CN4.5)--(CN4.6), including \(n=1\). \(\square\)

The local/genus condition is just \(r\) even, equivalently
\(\psi(c)=+1\), in addition to the prime-factor conditions.
It is insufficient globally: for every \(a\ge0\),
\[
                         3^{4a+2}                         \tag{CN4.7}
\]
satisfies all primitive local norm conditions but has \(r=0,b=1\).
Thus even after genus constraints, genuine individual class obstructions
persist. The additive theorem does not assert their disappearance for
individual norms.

## Complete primitive local conditions

Let \(C\) contain 1 and the integers supported on odd split primes,
excluding 7, and put \(S_\alpha=\{c\in C:\psi(c)=\alpha\}\).
The all-prime individual local candidates are exactly the disjoint union
\[
        \mathcal L=S_+\ \dot\cup\ 2S_+\
                         \dot\cup\ 7S_+\ \dot\cup\ 14S_+.
                                                               \tag{CN4.8}
\]
This is not a classification of \(P\).

At 2, put \(U=\{u\in\mathbb Z_2^*:u\equiv1,7\pmod8\}\).
The primitive norm image is precisely \(U\cup2U\).
For odd values \(x\) is odd, and the residues are \(1,7\pmod8\).
Conversely, for a target \(n\equiv1\pmod8\) take \(y=0\) and
\(x=\sqrt n\); for \(n\equiv7\pmod8\) take \(y=1\) and
\(x=\sqrt{n-14}\). These square roots are 2-adic units.
For an even primitive value, \(x=2z,y\) is odd, and
\(n=2u\), \(u=2z^2+7y^2\), so \(u\equiv1,7\pmod8\).
If \(u\equiv7\), use \(z=0,y=\sqrt{u/7}\); if \(u\equiv1\),
use \(z=1,y=\sqrt{(u-2)/7}\). This also proves sufficiency at
every higher dyadic digit.

At 7, unit values are nonzero squares; a primitive value divisible
by 7 has valuation one and unit quotient twice a square. Since
\((2/7)=1\), that is again a nonzero square. These conditions lift
by unit-square lifting. The factorization
\(\chi_{-56}=\chi_8\chi_{-7}\) shows that, on \(c\in C\),
\(\chi_{-7}(c)=\psi(c)\). Since \(\psi(7)=1\) and 2 is a
square modulo 7, the single genus condition in (CN4.8) gives both
ramified local tests. At other primes, inert factors are excluded
and split valuations are unrestricted.

The primitive local sum condition is exactly \(49\nmid N\).
At 7, \(-1\) is nonsquare: a divisible sum must use two
valuation-one values, and their unit quotients cannot sum to zero.
For each nonzero residue the ordered nonzero-square pair count is
\((6-2(N/7))/4\), either one or two, giving all other 7-adic
targets by lifting and simultaneous multiplication by 7 locally.

At 2 there is no exclusion. Odd targets lie in \(U+2U\);
the four residues \(\pm1\pm2\) cover all odd classes modulo eight.
Even targets with residues \(0,2,6\pmod8\) lie in \(U+U\).
For \(N\equiv4\pmod8\), \(N/2\equiv2,6\pmod8\) lies in
\(U+U\), so \(N\in2U+2U\).
Each is an exact union of residue classes with all higher digits
free. Other odd primes give all targets by unit norm surjectivity,
or the split local algebra. This proves the full local statement.

## Two primes correct the subgroup invisible to the genus

Use the fixed pool
\[
                        F=3\cdot5\cdot13\cdot19=3705.
\]
Every prime is split, has \(\psi=-1\), and hence class \(g\) or
\(-g\). The form \(Q_+\) represents them at
\((1,0),(0,1),(-2,1),(-1,2)\), respectively.

Two such distinct primes at exponent exactly one have possible
class sums \(0,2g\). Therefore any value \(2^e7^f c\), with
\(c\in S_+\) containing two reserved primes at exact exponent one,
is a primitive principal norm. Orient all other whole split powers
arbitrarily; their combined class is even. The ramified prefix
adds an even class by (CN4.3). The two reserved orientations can
cancel the result. No individual prime is required to be principal.

For a reduced target \(T\), use the 16 profiles from CN3Profile:
choose two dividing pool primes common if possible; if exactly one
divides, use it common and one other prime exclusively on each side;
if none divides, use two exclusive primes on each side. Leave unused
pool primes absent from both values. Each side gets two exact-one
primes, and its assigned product \(a_0,b_0\) has \(\psi=+1\).
Residual cofactors must lie in \(S_+\) and remain coprime to all of \(F\).
The exact zero/one local masses are CN3.7 at these odd split primes;
all selected cases are positive uniformly in target valuations.

After canceling a common assigned factor, its genus sign need not
be positive. The coefficient genera and residual signs must be retained
through the analytic transfer below. No division preserving the
principal genus without its sign correction is assumed.

## The genus-filter estimate at periods 56 and 112

On units modulo \(L=56\), let \(H=\ker\chi\), of size 12. Its
genus parts are
\[
 A_+=\{1,9,15,23,25,39\},\qquad
 A_-=\{3,5,13,19,27,45\}.
\]
Let \(V=\tau\mathbf1_C\), with normalized shifted families
\(F_B(m)=V(Bm)/V(B)\), \(G_D(l)=V(Dl)/V(D)\).
For \(t=1,2\), both signs and good shifts/coefficients in \(C\)
within \(T^\theta\), the existing CN2Filter proof gives
\[
 \sum_{\substack{um+tvl=T\\\psi(m)=\alpha,\ \psi(l)=\beta}}
 F_B(m)G_D(l)
 =\frac{T}{utv}\Gamma_{t;u,v;\alpha,\beta}(T)\rho_{F_B}\rho_{G_D}
 \prod_{\ell\nmid14}\frac{D_\ell}{M_{A,\ell}M_{B,\ell}}
                  +O(T^{1-\eta}),                         \tag{CN4.9}
\]
uniformly in gcds, all target valuations and the small divisibility
predicates used in restoration. Good radial factors and their zeros
are exactly (H3)--(H4) of the subgroup theorem, on this linear equation.
The finite factor is the actual count
\[
 \Gamma=\frac{tL}{|H_{tL}|^2}
 \#\{(a,b)\in H_{tL}^2:ua+tvb=T\bmod tL,\
                             \psi(a)=\alpha,\psi(b)=\beta\}.
                                                               \tag{CN4.10}
\]

For precision, this transfers the old proof rather than applying a
literal fixed-discriminant-minus-20 statement outside its hypotheses.
The first-variable filter has period \(tL\):
\[
 \mathbf1_{\psi(m)=\alpha}\mathbf1_{t\mid T-um}
 \mathbf1_{(z,56)=1,\ \chi(z)=1,\ \psi(z)=\psi(v)\beta},
                 \qquad z=(T-um)/t.                       \tag{CN4.11}
\]
It removes the bad 2-part of the paired-divisor congruence, leaving
only the good varying modulus \(vka\), coprime to 112. Expansion on
the first unit variable leaves exactly the character factors
\(\xi,\xi\chi\). Unless \(\xi\in\{1,\chi\}\), both are nonprincipal.
Retain all smooth scalar phases; kernel indices need only be units
with \(\chi=1\), not numbers with split support. Principal terms
reassemble through \(F_B\xi=F_B\). The fixed-period box input and
the signed bilinear application are unchanged with fixed period 112.
All original (H17)--(H18) losses fit the same hierarchy after taking
smaller positive \(\eta,\theta\). No distribution input at a bad
varying modulus is used.

If \(\mu\) is the first-filter average on \(H_{tL}\), each surviving
first residue has \(t\) second lifts. Hence
\(\Gamma=2tL\mu/\varphi(L)\), as in CN2.8; this gives the
real factor \(T/(utv)\), including the required \(1/t\).
Every coefficient/gcd/cofactor substitution carries both genus signs
as in CN2Filter. Scaling the finite variables by the actual removed
units bijects the old and transformed pairs and preserves \(\Gamma\).
In particular a negative-genus common factor causes a sign change,
not loss of that entire family.

For \(u=v=1,\alpha=\beta=+\), put
\[
 q_t(T)=\#\{b\in A_+:T-tb\in A_+\pmod{56}\}.
\]
Then the exact finite densities are
\[
 \Gamma_t(T)=\frac7{18}q_t(T),\qquad
 q_1(T)=d_2(T)\frac{6-2(T/7)}4,\quad
 q_2(T)=\frac{6-2(T/7)}4\ \ (T\text{ odd}),                \tag{CN4.12}
\]
for \(7\nmid T\), where \(d_2=2\) on \(T\equiv0\pmod8\),
\(d_2=1\) on \(T\equiv2,6\pmod8\), and zero otherwise.
For \(t=2\) even targets have zero density; both densities vanish
when \(7\mid T\). CRT separates the dyadic residues and nonzero
squares modulo 7. At period 112 each pair modulo 56 has two lifts,
which gives the displayed normalization. Every chosen positive
stratum has \(\Gamma\ge7/18\).

## Positive unweighting with the exact-prime block

The family (CN4.9) is precisely what the positive section-6 argument
needs, not just an estimate on a single base line.
For each fixed profile write \(m=a_0x,l=b_0y\), \((xy,F)=1\).
Its weights factor exactly:
\(V(a_0x)=V(a_0)V(x)\), and similarly on the second side.
Keep \(F\) as a fixed local block and exclude its primes from
every actual roughness/directive/common-cofactor restoration factor.
Analytic kernel indices may still contain these primes, whose finite
Euler factors must be retained.

The genus signs travel with every restored factor. Factors with
\(\psi=-1\) are permitted, with the residual signs changed by
the CN2 bijections; discarding them would reduce the split-prime
restoration density incorrectly. All fixed genus densities remain
the original positive \(\Gamma\). Exact-prime masses have a uniform
positive lower bound at all target depths. Outside \(F\), the rough
dimension remains two and the restored exponent remains one.
Removing finitely many split primes changes constants only.
The existing normalized-shift error budget, positive cofactor size
bound and recoverable patterns are unchanged.

Thus the existing positive unweighting proof gives
\[
 \#\{\text{pairs in any selected positive core profile}\}
                  \gg T J_{-56}(T)/\log T.                 \tag{CN4.13}
\]
A weighted asymptotic alone is not being used to infer this claim.
The finite exact-prime and genus modifications are retained at every
restoration stage. Take a common cutoff and lower constant over the
16 profiles and finitely many selected strata.

## Target selection and nonprincipal ramification

For \(49\nmid N\), write \(N=7^f M\), \(f=0,1\), \(7\nmid M\).
Select one of these sufficient families:
\[
\begin{array}{c|c|c|c}
 M\text{ condition}&T&t&\text{actual value pair}\\\hline
 M\text{ odd}&M&2&(7^f m,\ 2\cdot7^f l)\\
 M\equiv0,2,6\pmod8&M&1&(7^f m,\ 7^f l)\\
 M\equiv4\pmod8&M/2&1&(2\cdot7^f m,\ 2\cdot7^f l).
\end{array}                                                \tag{CN4.14}
\]
Here \(m+t l=T\), \(m,l\in S_+\), with the two exact-one pool
primes per side. Every selected base density is positive by CN4.12.
Each actual prefix has ramified exponents at most one.
The two-prime subgroup correction proves that both actual values
are primitive principal norms, including the ramified prefix classes.

Neither 2 nor 7 is a principal norm: \(1\in P\), but \(2,7\notin P\).
Consequently multiplication by 2 or 7 is not a closure property of
all \(P\), and no generator of norm 2 or 7 is invoked.
The prefix operation in (CN4.14) is justified by a fresh choice of
ideal orientations using the reserved primes. On numerical value
pairs it is injective, and all scale factors are at most 14.
This suffices to transfer (CN4.13) to the lower bound in CN4.1.
Deleting factors 2 and 7 leaves \(J_{-56}\) unchanged and changes
\(N/\log N\) only by bounded constants.

For the upper bound drop the genus, full class, ramified and dyadic
restrictions. Each odd inert prime excludes residues zero and \(N\)
on the line of two candidate values, merging only when it divides
the target. The existing upper block sieve has dimension one and
gives \(O(N J_{-56}(N)/\log N)\). Together with the local zero
cases and the lower construction, this proves Theorem CN4.

## A pointwise comparison with all local norm candidates

Let \(\mathcal L\) be the individual local candidate set (CN4.8), and
let \(R_{\mathcal L}(N)\) count ordered pairs of its positive values.

**Proposition CN4ClassError.** Uniformly for integers \(N\ge3\),
\[
 0\le R_{\mathcal L}(N)-R_P(N)
       \ll \frac{N J_{-56}(N)}{(\log N)^{5/4}}.             \tag{CN4.17}
\]
Consequently, through every sufficiently large target with \(49\nmid N\),
\[
 R_P(N)=R_{\mathcal L}(N)
              \left(1+O((\log N)^{-1/4})\right).            \tag{CN4.18}
\]
The constants are fixed-field constants, uniform in every target-prime
valuation. This compares two unweighted value convolutions. It does not
evaluate a leading constant or singular series for either convolution.

Partition the primes \(q\nmid14\) into
\[
 \mathscr H=\{\chi(q)=1,\psi(q)=1\},\quad
 \mathscr O=\{\chi(q)=1,\psi(q)=-1\},\quad
 \mathscr I=\{\chi(q)=-1\}.
\]
They comprise 6,6,12 of the 24 reduced residue classes modulo 56.
Write \(C_{\mathscr H},C_{\mathscr O},C_\chi\) for the corresponding
prime-factor semigroups; \(C_\chi\) allows both split parts.

**The square decomposition.** By CN4Individual, a hidden failure
\(a\in\mathcal L\setminus P\) has no odd-class prime to odd exponent.
It therefore has a unique expression
\[
 a=r d^2 c,\qquad r\in\{1,2,7,14\},\quad
       d\in C_{\mathscr O},\quad c\in C_{\mathscr H}.        \tag{CN4.19}
\]
The additional fixed class-two parity selecting failures inside this
family may be dropped for an upper bound. The prime supports of \(d,c\)
are disjoint; in particular \(d\) has no inert prime factor.
Dropping the other summand's genus condition and charging each pair
with a bad summand to one of its two orientations gives
\[
 R_{\mathcal L}(N)-R_P(N)
 \le2\sum_{r,s\in\{1,2,7,14\}}\sum_{d\in C_{\mathscr O}}
 \#\{rd^2c+sv=N:c\in C_{\mathscr H},v\in C_\chi,\ c,v\ge1\}.
                                                               \tag{CN4.20}
\]
All counts here are counts of values.

**The uniform residue conditions.** Fix \(r,s,d\). Put \(g=(r,s)\);
since \((d,14)=1\), also \((rd^2,s)=g\). If \(g\nmid N\), the
count is zero. Otherwise set
\[
 A=(r/g)d^2,\quad B=s/g,\quad n=N/g.
\]
Then \((A,B)=1\), and \(B\) has only prime factors 2 and 7.
Solving \(Ac_0\equiv n\pmod B\) gives every integral solution as
\[
             c=c_0+Bk,\qquad v=v_0-Ak.                    \tag{CN4.21}
\]
Positivity restricts \(k\) to an interval of real length \(X\ll N/d^2\).
An empty interval contributes zero; endpoints introduce only the
constant errors below.

For \(q\nmid14\), forbid both \(c=0,v=0\) at inert primes, and
only \(c=0\) at \(\mathscr O\)-primes. The number of forbidden
residues is exactly
\[
 \nu_N(q)=
 \begin{cases}
 2-\mathbf1_{q\mid N},&q\in\mathscr I,\\
 1,&q\in\mathscr O,\\
 0,&q\in\mathscr H.
 \end{cases}                                               \tag{CN4.22}
\]
At inert primes both slopes in (CN4.21) are units, because \(d\)
has no inert factors. The roots coincide precisely when \(q\mid n\),
equivalently \(q\mid N\). At an \(\mathscr O\)-prime only the
\(c\)-root is forbidden, and its slope \(B\) is a unit even if
\(q\mid d\). Thus no coefficient-radical loss occurs, including at
every prime dividing \(d\). The conditions at 2 and 7 may be dropped.

For a squarefree product \(h\) of sieved primes, CRT gives
\(\nu_N(h)=\prod_{q\mid h}\nu_N(q)\) residue classes. In the
interval of length \(X\) their count is
\[
                    X\nu_N(h)/h+O(\nu_N(h)),              \tag{CN4.23}
\]
uniformly in all root positions, coefficients and target valuations.
This is a direct residue count. The coefficient \(d^2\) need not
belong to \(C_{\mathscr H}\), so a literal coefficient specialization
of (H5) would be outside its prime-support hypotheses.

**A self-contained Selberg upper bound.** Use only primes with
\(\nu_N(q)>0\). On squarefree products of them put
\[
 f(h)=\nu_N(h)/h,\qquad
 w(h)=\prod_{q\mid h}\frac{\nu_N(q)}{q-\nu_N(q)},\qquad
 G_N(z)=\sum_{h\le z}\mu^2(h)w(h).
\]
Every denominator is positive: these primes are odd and \(\nu_N(q)\le2<q\).
Real weights \(\lambda_h\) supported on squarefree \(h\le z\),
with \(\lambda_1=1\), majorize the sifted set after their divisor
sum is squared. Using (CN4.23) gives
\[
 S(X,z)\le XQ(\lambda)
   +O\!\left(\sum_{h,k\le z}|\lambda_h\lambda_k|\nu_N([h,k])\right),
 \quad
 Q(\lambda)=\sum_{h,k\le z}\lambda_h\lambda_k f([h,k]).
                                                               \tag{CN4.24}
\]
Here \(S(X,z)\) counts the integers avoiding every forbidden residue
at primes at most \(z\); every pair in (CN4.20) is a survivor.

Set \(y_u=\sum_{u\mid h,\ h\le z}\lambda_h f(h)\).
The elementary identity
\(1/f((h,k))=\sum_{u\mid(h,k)}1/w(u)\) gives
\[
 Q(\lambda)=\sum_{u\le z}y_u^2/w(u),\qquad
                 1=\sum_{u\le z}\mu(u)y_u.
\]
Cauchy's inequality is attained by \(y_u=\mu(u)w(u)/G_N(z)\).
Finite Möbius inversion yields the actual admissible weights
\[
 \lambda_h=\mu(h)
   \prod_{q\mid h}(1-\nu_N(q)/q)^{-1}
   \frac{G_{N,h}(z/h)}{G_N(z)},\quad
 G_{N,h}(u)=\sum_{\substack{k\le u\\\gcd(k,h)=1}}\mu^2(k)w(k).
\]
Hence \(Q=1/G_N(z)\), \(\lambda_1=1\), and
\(\lvert\lambda_h\rvert\ll_\epsilon h^\epsilon\), uniformly in \(N\).
Indeed the displayed quotient is at most one, and the finite product
is bounded by \(\prod_{q\mid h}q/(q-2)\ll_\epsilon h^\epsilon\).
Likewise \(\nu_N(h)\ll_\epsilon h^\epsilon\), while
\(\nu_N([h,k])\le\nu_N(h)\nu_N(k)\). Relabelling \(\epsilon\)
in these bounds proves
\[
               S(X,z)\le \frac{X}{G_N(z)}
                            +O_\epsilon(z^{2+\epsilon}). \tag{CN4.25}
\]
There is no distribution estimate or coefficient-size hypothesis in
this CRT endpoint budget.

**A uniform lower bound for \(G_N(z)\).** Define
\[
 E_N(y)=\prod_{q\le y}(1-\nu_N(q)/q)^{-1},
\]
omitting zero-density primes, and let \(E_0\) use two roots at
every inert prime. The target correction is exact:
\[
 E_N(y)=
 \frac{E_0(y)}
 {\prod_{\substack{q\le y\\q\in\mathscr I,\ q\mid N}}(q-1)/(q-2)}
                           \ \ge\ \frac{E_0(y)}{J_{-56}(N)}.
                                                               \tag{CN4.26}
\]
Williams, *Mertens' Theorem for Arithmetic Progressions*, Theorem 1,
printed p.356, gives each fixed reduced residue-class Euler product
a positive constant times \((\log y)^{-1/24}\), with a smaller error.
The exact retained primary PDF is
[sources/williams-mertens-ap-1974.pdf](sources/williams-mertens-ap-1974.pdf).
Since there are twelve inert and six odd-class split residue classes,
\[
 E_0(y)\asymp
 \prod_{\substack{q\le y\\q\in\mathscr I}}(1-1/q)^{-2}
 \prod_{\substack{q\le y\\q\in\mathscr O}}(1-1/q)^{-1}
                         \asymp(\log y)^{5/4}.             \tag{CN4.27}
\]
The comparison at inert primes has convergent positive correction
\((1-1/q)^2/(1-2/q)=1+1/(q(q-2))\).
This fixed-modulus Mertens theorem is the analytic input to the
absolute comparison; no class-field prime distribution is used.

Truncation is uniform in \(N\). Normalize the finite positive Euler
weights on all squarefree products of primes at most \(y\).
Prime inclusion probabilities are \(w(q)/(1+w(q))=\nu_N(q)/q\).
For an integer \(m\ge2\),
\[
 \log(m!)\ge\sum_{q\le m}\lfloor m/q\rfloor\log q
       \ge\frac m2\sum_{q\le m}\frac{\log q}{q},
 \qquad \log(m!)\le m\log m.
\]
Thus \(\sum_{q\le y}\log q/q\le2\log y\). The Euler measure has
\(\mathbb E(\log h)\le4\log y\).
Take \(y=z^{1/16}\). Markov's inequality retains at least three
quarters of its mass on \(h\le z\). These terms are included in
\(G_N(z)\), so for all sufficiently large \(z\),
\[
 G_N(z)\ge\tfrac34 E_N(z^{1/16})
           \gg \frac{(\log z)^{5/4}}{J_{-56}(N)}.           \tag{CN4.28}
\]
Combining (CN4.25) and (CN4.28) gives the uniform bound
\[
 S(X,z)\ll \frac{XJ_{-56}(N)}{(\log z)^{5/4}}
                                  +O_\epsilon(z^{2+\epsilon}). \tag{CN4.29}
\]

**Summing small and large squares.** Take \(D=z=N^{1/10}\),
with integer cutoffs, and \(\epsilon=1/2\).
For \(d\le D\), sum (CN4.29) using
\(\sum d^{-2}\le\zeta(2)\) and \(\#\{d\le D\}\le D\).
For each fixed prefix pair this gives
\[
 \ll \frac{N J_{-56}(N)}{(\log N)^{5/4}}+N^{7/20}.
\]
For \(d>D\), ignore all support and congruence conditions. Each positive
\(c\le N/(rd^2)\) determines \(v\) at most once, so the total is
\[
                     \le N\sum_{d>D}d^{-2}\ll N^{9/10}.
\]
The sixteen prefix pairs and two orientations in (CN4.20) are fixed
costs. Since \(J_{-56}(N)\ge1\), the power errors are absorbed by
the main bound for large \(N\); enlarging the constant covers the
remaining integers \(N\ge3\). This proves (CN4.17).

Finally \(P\subset\mathcal L\) and the independent absorber proof
of Theorem CN4 give \(R_{\mathcal L}(N)\gg N J_{-56}(N)/\log N\)
on sufficiently large \(49\)-admissible targets. Division of
(CN4.17) by this lower bound proves (CN4.18). A leading asymptotic
for the local-candidate convolution, if established at this scale,
would transfer to the principal norm convolution. No such leading
asymptotic is established here.

## Evidence and limits

The [individual norm checker](verify_primitive_norm_minus56.py) compares
the exact coordinate counts with CN4.6 for every integer through \(10^6\),
using independent form labels for prime classes and whole-power
factorization. All counts agree: 87187 primitive values, 510420 signed
coordinate pairs and 26066 individual class failures among 113253
individually locally admissible values. The [finite profile checker](verify_primitive_norm_class4_profiles.py)
checks the explicit primes, exact profile masses, all target templates,
period-56/112 densities, genus transport and dyadic strata.
Neither finite checker proves the imported infinite analytic estimates.

The [finite Selberg checker](verify_primitive_norm_class4_sieve.py) and
[receipt](primitive_norm_class4_sieve_checks.json) compare direct coordinates
with independent local flags through one million, checking every hidden
failure's unique square envelope. Separate finite suites check combined
CRT roots, positive progressions, exact rational Selberg optimization and
interval bounds, and Euler marginals and truncation. These finite
identities do not prove the imported Mertens theorem or extrapolate the
pointwise comparison to all integers; the unrestricted proof is above.

The [source ledger](primitive_norm_minus56_source_ledger.json) pins the
retained primary sources and canonical proof dependencies. Internal
reviews separately reconstructed the C4 arithmetic, all local
valuations, fixed-period normalization and finite-block restoration.
There is no external review or formal verification.

A numerical cutoff, completeness of the bounded exception list below,
and an unweighted leading constant remain open. This concrete C4 case does
not assert a theorem for every even class group, arbitrary quadratic
orders, mixed fields or simultaneous norm constraints, and it does not
address the \(C_{12}+C_{12}\) convolution barrier.


Reproduce the finite receipts with:

    python3 research/arithmetic/verify_primitive_norm_minus56.py --output research/arithmetic/primitive_norm_minus56_checks.json
    python3 research/arithmetic/verify_primitive_norm_class4_profiles.py --output research/arithmetic/primitive_norm_class4_profile_checks.json


## Exact bounded binary exceptions

**Proposition CN4Finite.** For every \(2\le N\le2^{24}=16777216\),
\[
 N\in P+P\quad\Longleftrightarrow\quad
                     49\nmid N\text{ and }N\notin E_{56},
\]
where \(E_{56}\) is the following 118-element set:

    3,4,5,6,7,8,9,10,11,12,13,14,17,18,20,
    21,22,23,25,26,27,34,35,39,42,43,47,50,52,55,
    56,59,61,63,67,70,74,76,84,85,90,91,94,97,103,
    105,112,116,121,124,127,133,139,140,154,161,163,175,182,203,
    212,217,231,252,259,266,292,301,308,322,329,349,357,364,371,
    399,406,434,455,476,525,553,581,623,693,700,714,742,749,791,
    847,889,917,1022,1036,1057,1134,1141,1169,1204,1260,1281,1316,1379,1428,
    1477,1484,1645,1869,1967,2086,2380,2604,2807,2996,5012,5257,8876

The largest is 8876. Every listed number is an actual global
nonrepresentation. Absence of further exceptions is proved only through
the stated bound, and the unrestricted exact-exception gap remains open.

The [bounded sumset checker](verify_primitive_norm_minus56_sumset.py) and
[receipt](primitive_norm_minus56_sumset_checks.json) enumerate every
primitive \(Q_0\) value through the bound, and independently reconstruct
its membership by CN4.5. Prime classes are labelled using all four
ideal classes, via \(Q_0,Q_2,Q_\pm\). The coordinate and ideal flags
agree at every integer, giving 1344933 primitive values.
A second complete comparison verifies the genus parity against the
dyadic local residue set, yielding 1721437 individual local candidates
and 376504 genuine individual class obstructions through this bound.
Trial factorization through 10000 and a signed square-box enumeration
through 1000 challenge both reconstructions.

The 2507 primitive anchors at most 25000 are shifted against the complete
packed value flag. This gives only valid sums, and covers every local
target except the 118 displayed integers. Each residual receives a
complete search over \(1\le a\le\lfloor N/2\rfloor\), finding no
representation. The final mask equals the full local mask with exactly
\(E_{56}\) removed: 16434705 represented targets, 16434823 local
targets and 342392 separate local exclusions divisible by 49.
Because all summands are positive, values above \(2^{24}\) cannot
contribute to this finite range. Corruption controls reject erased
representations, injected global/local exceptions, wrong genus, and
the hidden class/valuation failures 9 and 450. This proves the
finite proposition; it gives no analytic cutoff.

## Ramification switches the principal and class-two norm forms

Let \(P_0=P\) and let
\(P_2=\{2u^2+7v^2>0:\gcd(u,v)=1\}\).
Put \(U_i=P_i\cap\{n:7\nmid n\}\), and let \(O_i\) denote the
odd values in \(P_i\). There are exact unrestricted identities
\[
\begin{split}
 P_0\cap7\mathbb Z&=7U_2,&P_2\cap7\mathbb Z&=7U_0,\\
 P_0\cap2\mathbb Z&=2O_2,&P_2\cap2\mathbb Z&=2O_0.
\end{split}                                                \tag{CN4.15}
\]
Indeed
\[
\begin{split}
 Q_0(7v,u)&=7Q_2(u,v),&
 Q_2(7y,x)&=7Q_0(x,y),\\
 Q_0(2u,v)&=2Q_2(u,v),&
 Q_2(x,2y)&=2Q_0(x,y).
\end{split}
\]
For the first map a 7-unit quotient forces \(7\nmid u\), so
\(\gcd(7v,u)=1\) exactly when \(\gcd(u,v)=1\). Conversely a
primitive \(Q_0\) value divisible by 7 forces \(7\mid x,7\nmid y\),
giving that map in reverse. The other 7-map has the analogous
\(7\nmid x\) guard. For the two dyadic maps the odd quotient forces
the unchanged coordinate odd, preserving gcd one; divisibility of
the scaled coordinate gives the converses. These guards also exclude
imprimitive iterations with a ramified square.

Since \(-1\) is nonsquare modulo 7, both summands of a 7-divisible
target are themselves divisible by 7. Therefore
\[
           7M\in P_0+P_0\quad\Longleftrightarrow\quad
                                      M\in U_2+U_2.        \tag{CN4.16}
\]
Unlike CN3.8 for discriminant \(-23\), the quotient uses the other
ideal-class norm form. It is not a closure statement \(7P_0\subset P_0\).
The bounded list partitions into 44 exceptions coprime to 7, largest
349, and 74 multiples of 7 whose quotients have largest 1268.
By CN4.16 the latter classify the missing 7-unit targets of \(U_2+U_2\)
only up to \(\lfloor2^{24}/7\rfloor=2396745\).
The structural identities are unrestricted; completeness of these
exception lists remains bounded.

Reproduce with:

    python3 research/arithmetic/verify_primitive_norm_minus56_sumset.py --output research/arithmetic/primitive_norm_minus56_sumset_checks.json


Reproduce the finite class-error identities with:

    python3 research/arithmetic/verify_primitive_norm_class4_sieve.py --output research/arithmetic/primitive_norm_class4_sieve_checks.json
