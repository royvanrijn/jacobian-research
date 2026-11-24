# Exact sums of two primitive discriminant −19 norms

Let
\[
 Q(x,y)=x^2+xy+5y^2,\qquad
 P=\{Q(x,y)>0:x,y\in\mathbb Z,\ \gcd(x,y)=1\}.
\]
Signs and zero coordinates are allowed. Let \(C\) contain 1 and the
positive integers whose prime factors are quadratic residues modulo 19.

**Theorem M19 (exact primitive classification).**
\[
 C+C=\{n\ge2:2\mid n,\ 19\nmid n\}\setminus\{4,20\},
                                                        \tag{M19.1}
\]
\[
 P+P=\{N\ge2:2\mid N,\ 361\nmid N\}
                           \setminus\{4,76,380\}.         \tag{M19.2}
\]
These are identities of numerical value sets, with positive summands.

The proof combines a classical quaternionic Hecke estimate with exact
finite coverage. The analytic part covers every inert core at least
\(2^{34}\), and two thin families arising from exceptional quotients.
The finite certificate covers **every** eligible target through
\(2^{34}=17179869184\). Multiplicative descent joins the two parts;
no finite scan is extrapolated to an unrestricted assertion.

The closest proof in the repository is the
[discriminant −11 classification](PRIMITIVE_NORM_MINUS11_EFFECTIVE.md).
Here the principal quaternion class has smaller mass, and the two core
exceptions prevent the exception-free descent used there. We supply
the changed normalization, the exceptional-quotient argument, and the
complete new numerical certificate. The classical inputs are identified
below and pinned in the [source ledger](primitive_norm_minus19_source_ledger.json).
No historical priority, formal verification, or external review is claimed.

## 1. Primitive values and the ramified prime

Put \(R=\mathbb Z[\omega]\), where
\(\omega=(1+\sqrt{-19})/2\). Its ideal class number is one:
the imaginary quadratic Minkowski bound (Milne, Theorem 4.3) is
\((2/\pi)\sqrt{19}<3\), while 2 is inert, so the only possible
integral representative ideal norm is 1. Quadratic reciprocity gives
\(\chi_{-19}(p)=(p/19)\), including at 2, which is inert.
Primitive ideal factorization therefore gives
\[
                         P=C\ \dot\cup\ 19C.             \tag{M19.3}
\]
Indeed an inert prime dividing a norm divides both coordinates, and a
primitive element has ramified norm valuation at most one. Conversely,
at every split prime choose an entire prime power on one oriented prime
ideal. Since the ideals are principal, the product of their generators,
optionally multiplied by \(\sqrt{-19}\), has the desired norm and is
primitive. This is an existence statement for a primitive representation
of the value; later quaternion witnesses may have split coordinate
common factors.

Every member of \(C\), and hence of \(P\), is odd. Members of \(C\)
are nonzero squares modulo 19. As \(-1\) is nonsquare modulo 19,
\(C+C\) contains no multiple of 19. The same observation and (M19.3)
give, for every positive integer \(m\),
\[
 19m\in P+P\quad\Longleftrightarrow\quad m\in C+C.         \tag{M19.4}
\]
In particular \(P+P\) has no multiple of \(361\).

For reference, multiplication by \(\sqrt{-19}\) is
\[
 (x,y)\longmapsto(-x-10y,2x+y),\qquad
 M^2=-19I,\quad\det M=19,\quad Q(Mv)=19Q(v).
\]
It carries a 19-unit primitive generator to a primitive generator of
norm 19 times as large. Primitivity at 19 follows from
\(4Q=(2x+y)^2+19y^2\); at the other primes the matrix is invertible.

## 2. The complete ideal-class quotient and its measure

Use
\[
 B=(-19,-1)_{\mathbb Q},\qquad O=R\oplus Rj,\qquad
 j^2=-1,\quad ja=\bar a j.
\]
Its norm is \(F=Q\oplus Q\). The algebra has discriminant 19 and
the order is maximal by Voight, Example 15.5.7 and Theorem 15.5.5,
printed p.241 of the pinned edition. Its units are \(\{\pm1,\pm j\}\):
the only field elements of norm 1 are \(\pm1\), and the two
nonnegative field norms must sum to 1.

Retain **all** locally principal right-ideal classes. If \(\Gamma_i\)
is the projective unit group of the corresponding maximal left order,
Voight's Theorem 25.1.1, printed p.403, gives
\[
 \mathfrak m=\sum_i\frac1{|\Gamma_i|}=\frac{19-1}{12}=\frac32,
 \qquad
 w=\frac{1/|\Gamma_0|}{\mathfrak m}=\frac13.               \tag{M19.5}
\]
Here \(|\Gamma_0|=2\). Thus \(w\) is the Haar probability of the
principal class; no quaternion class-number-one assertion is made.

At an odd inert prime \(p\), impose the projective reduction kernel
with finite quotient \(G_p=\mathrm{PGL}_2(\mathbb F_p)\). At 19 use
\[
 G_{19}=\mathbb F_{19^2}^{\times}/\mathbb F_{19}^{\times},
                         \qquad |G_{19}|=20.
\]
Write \(X_p\) for the resulting full finite adelic quotient. Its fibers
are \(\Gamma_i\backslash G_p\). Every projective unit group injects
into the indicated finite group. A noncentral norm-one unit has integral
trace \(-1,0,1\), hence discriminant \(-3\) or \(-4\). At a split
odd prime, scalar reduction forces this discriminant to be divisible
by \(p^2\), impossible even for \(p=3\). At 19, scalar reduction
modulo the radical forces divisibility by 19, also impossible.
Consequently the fine stabilizers are trivial, and
\[
 |X_p|=|G_p|\mathfrak m,
 \qquad |X_{19}|=30.                                    \tag{M19.6}
\]
Uniform fine-point measure pushes forward to inverse-unit mass.

For an index \(m\) coprime to 19 and to any added unramified level
prime, let \(T_m\) be the all-integral-ideal Hecke operator, including
imprimitive ideals. Voight, Propositions 41.3.1 and 41.3.6, give degree
\(\sigma(m)\) and recurrence
\(T_{\ell^{a+2}}=T_\ell T_{\ell^{a+1}}-\ell T_{\ell^a}\).
It is self-adjoint for the preceding measure. Explicitly, reversing
the marked ideal flag \(J\subset I\) to \(mI\subset J\) is an
involution up to central scaling. At auxiliary primes the local ideals
are equal, since \(m\) is a unit there; the remaining marking change
is the projectively trivial scalar \(m^{-1}\). This proves symmetry
on fine points and the stated weighted symmetry on coarse classes.

For \(h\) supported on the principal fiber and its base point \(x_0\),
\[
 (T_mh)(x_0)=\frac14
       \sum_{\substack{\alpha\in O\\\operatorname{nrd}\alpha=m}}
                        h(\rho_p(\alpha)).               \tag{M19.7}
\]
Only principal returns contribute. Each ideal has four generators;
principalization initially gives inverse generators, but projectively
\(\alpha^{-1}=\bar\alpha\), and conjugation permutes the norm shell.

At a fine prime \(p\), the only global reduced-norm characters are
\(1\) and the quadratic residue character \(\psi_p\). Their
\(T_m\)-eigenvalues are \(\sigma(m)\) and \((m/p)\sigma(m)\).
To check completeness of this list, maximal units away from \(p\)
have norm group \(\mathbb Z_\ell^\times\); the kernel at \(p\)
has square residue norms. Central scalars contribute squares, and
rational quaternion norms are all positive rationals by
Hasse–Schilling, Voight Main Theorem 14.7.4. The remaining norm quotient
has order two. With no auxiliary prime its order is one.

On the complement of all these characters,
\[
                         \|T_m\|\le\tau(m)\sqrt m.        \tag{M19.8}
\]
We use the complete
[classical Hecke transfer lemma](../../results/eisenstein-restricted-prime-sums/sections/hecke-transfer.tex).
Its proof removes the norm characters by strong approximation
(Voight, Main Theorem 28.5.3), applies global Jacquet–Langlands
(Theorem 14.4 and Section 5, Proposition 5.20), and uses the weight-two
Ramanujan bound from Deligne (the exact editions and theorem numbers
are in the ledger). Directly, Deligne's *La conjecture de Weil. I*,
Theorem 8.2, printed p.302, states the requisite root bound for a
normalized primitive cusp form of arbitrary level, weight at least two,
and every prime outside that level. Apply it to the primitive form
underlying each cuspidal constituent; its unramified roots are unchanged
on old spaces at other primes. The infinite type here is trivial, hence transfers
to weight two; every prime dividing \(m\) is split in \(B\) and has
maximal local level. These are the required hypotheses. No bound is
being used at an auxiliary level prime dividing its own Hecke index.
No explicit level-19 theta/newform identity is needed in this proof.

## 3. Odd components and the filtered spectral bound

From now through Section 6, \(n=2^t m\) is even, \(m\) is odd,
and \(19\nmid n\). Set
\[
 S(n)=\sigma(n)-\sigma(n/4)=3\cdot2^{t-1}\sigma(m)
                                        \ge\tfrac32n,
\]
\[
 \Delta(n)=\tau(m)\sqrt m\,2^{t/2}\frac{t+3}{2},\qquad
 E(n)=\frac{2\Delta(n)}{S(n)}.                            \tag{M19.9}
\]
Nonintegral indices contribute zero. Since \(R/2R=\mathbb F_4\),
an even field norm has both coordinates even. The number \(G(n)\)
of quaternion shell points having both field norms odd is therefore
\(r_F(n)-r_F(n/4)\). Both components of every such point are nonzero.
At every quotient above it is counted by
\(D_n=T_n-T_{n/4}\); multiplication by 2 is projectively trivial.

The constant and quadratic-character eigenvalues of \(D_n\) are
\(S(n)\) and \((n/p)S(n)\). Its complementary operator norm is at
most \(\Delta(n)\). Indeed, for normalized Hecke roots \(u,v\) at 2,
\(|u|=|v|=1\) and \(uv=1\), the local difference is
\[
 2^{t/2}\left(u^t+v^t+
                    \frac12\sum_{a=1}^{t-1}u^av^{t-a}\right).
\]
For \(t=1\) the interior sum is empty. Its absolute value is at most
\(2^{t/2}(t+3)/2\); the odd factor contributes \(\tau(m)\sqrt m\).
This holds at arbitrarily large two-adic valuations and at auxiliary
levels, where one must retain all their cusp constituents.

The principal-class indicator on the coarse quotient has mean \(w\)
and nonconstant squared norm \(w(1-w)\). Its return correlation gives
\[
 \frac{G(n)}4\ge wS(n)-(1-w)\Delta(n)
                  =\frac13 S(n)-\frac23\Delta(n).        \tag{M19.10}
\]
The opposite inequality, with a plus sign before the error, holds too.
We will subtract absolute bad-point bounds from (M19.10), so no
division by \(G(n)\) or assumption \(E(n)<1\) is required.

## 4. Excluding inert and ramified factors

Let \(\mathcal I\) be the odd inert primes, and put
\(\mu_p=2/(p(p-1))\). For \(p\in\mathcal I\), let \(b_p(n)\)
count points in the odd-component shell for which at least one field
norm is divisible by \(p\).

If \(p\nmid n\), the bad invertible reductions are exactly the
normalizer \(H_p\) of the nonsplit torus from \(R/pR\). It has order
\(2(p+1)\), contains \(\Gamma_0\), and consists of the two
field-component axes. Put \(h=1_{\Gamma_0\backslash H_p}\), supported
in the principal fiber. Its mean is \(\nu=w\mu_p\); its quadratic
norm-character mean is zero, since the torus norm is onto
\(\mathbb F_p^\times\). Right-\(H_p\)-invariance gives
\((D_nh)(x_0)=\langle h,D_nh\rangle/\nu\). The complementary
squared norm is \(\nu(1-\nu)\). Equations (M19.7)–(M19.9) imply
\[
             \frac{b_p(n)}4\le w\mu_pS(n)+\Delta(n).
                                                        \tag{M19.11}
\]
If \(p\mid n\), anisotropy forces both field components divisible by
\(p\) whenever one is. Thus the bad count is zero if \(p^2\nmid n\),
and otherwise equals \(G(n/p^2)\). Now
\[
 S(n/p^2)\le S(n)/p^2,\qquad
 \Delta(n/p^2)\le\Delta(n)/p.
\]
The upper version of (M19.10), together with
\(p^{-2}\le\mu_p\) and \((1-w)/p\le1\), proves (M19.11) on this
branch as well. This explicitly treats target-dividing inert primes
without an operator at its own level prime.

At 19 the radical quotient is \(\mathbb F_{19}[j]\). The two
field-component axes project to \(\{1,j\}=\Gamma_0\), so the bad
subset of the principal fiber is the single point \(x_0\). Its
indicator \(\delta\) has mean \(1/30\); its coefficients on
\(1,\psi_{19}\) are both \(1/30\). Dividing its return correlation
by the point probability gives
\[
 \frac{b_{19}(n)}4
 \le\frac{1+(n/19)}{30}S(n)+\frac{14}{15}\Delta(n)
 \le\frac1{15}S(n)+\frac{14}{15}\Delta(n).                \tag{M19.12}
\]

For large inert primes we use a lattice estimate:
\[
 \#\{(x,y)\in\mathbb Z^2:0<Q(x,y)\le X\}\le2X
                       \qquad(X\ge0).                   \tag{M19.13}
\]
The row count and monotonicity of
\(\sqrt{\max(X-19y^2/4,0)}\) give
\[
 S_Q(X)\le\frac{2\pi}{\sqrt{19}}X+
                             (2+4/\sqrt{19})\sqrt X.
\]
For integer \(X\ge36\), use \(\pi<22/7\), \(\sqrt{19}>43/10\),
and \(\sqrt X\ge6\). The resulting coefficient is
\(440/301+21/43=587/301<2\). Exact lattice enumeration checks
\(1\le X\le35\). Passing to \(\lfloor X\rfloor\) proves the real
version of (M19.13).

Define
\[
 M(X)=\max_{1\le a\le X}
       \prod_{\substack{q\mid a\\(q/19)=1}}(v_q(a)+1).
\]
Class number one and ideal factorization give \(r_Q(a)\le2M(X)\)
for \(1\le a\le X\). If an inert prime \(p\) divides a component
norm, that component is \(p\) times a nonzero vector. By (M19.13)
there are at most \(2n/p^2\) such vectors and at most \(2M(n)\)
choices for the other component. Counting both orientations yields
\[
                       b_p(n)\le8M(n)n/p^2.              \tag{M19.14}
\]

Write
\[
 B_0=\sum_{p\in\mathcal I}\mu_p,\quad
 \Pi(P)=\#\{p\in\mathcal I:p\le P\},\quad
 T(P)=\sum_{\substack{p\in\mathcal I\\p>P}}p^{-2}.
\]
Subtract (M19.11) for \(p\le P\), (M19.12), and (M19.14) for
\(p>P\) from (M19.10). Since \(S(n)\ge3n/2\), the number of
remaining points is at least
\[
 \frac43S(n)\left\{
 \frac45-B_0-4M(n)T(P)
              -E(n)\left(\frac{12}{5}+\frac32\Pi(P)\right)
                  \right\}.                             \tag{M19.15}
\]
Each surviving point has two positive odd norms, neither divisible by
19 nor by any inert prime. Its two values are in \(C\), so positivity
of the brace proves \(n\in C+C\).

## 5. Exact constants and the infinite ranges

The odd inert primes lie in the nine classes
\[
             3,13,15,21,27,29,31,33,37\pmod {38}.
\]
Their sum is 209. Counting these progressions and bounding the tail of
each decreasing square reciprocal series by its first term plus its
integral gives, for every real \(P>0\),
\[
 \Pi(P)\le\frac9{38}P+\frac72,\qquad
 T(P)\le\frac9{38P}+\frac9{P^2},\qquad B_0<\frac9{25}.    \tag{M19.16}
\]
The last inequality is checked rationally: sum the actual inert primes
through 1024 and majorize the remainder by \(3T(1024)\), using
\(\mu_p\le3/p^2\).

The finite Euler certificates give
\[
 M(X)\le7X^{1/6}\quad(X\ge1),\qquad
 E(n)<\frac{81}{10}n^{-1/3}
               \quad\text{if every prime dividing }n\text{ is inert}.
                                                        \tag{M19.17}
\]
Here and below an inert core is even. For the first bound, maximize
\((a+1)^6/p^a\) for each split prime \(p<64\) and \(0\le a\le10\),
then multiply. The product is
\[
 \frac{2337302235907620864}{35724352769875}<7^6.
\]
Exponents beyond 10 cannot improve it, since
\((12/11)^6<2\) and subsequent ratios decrease. Primes at least 64
contribute at most 1, since \(a+1\le2^a\).

For the second bound, (M19.9) gives exactly
\[
 E(n)n^{1/3}
 =\frac{2(t+3)}3\,2^{-t/6}
   \prod_{p^a\parallel m}
           \frac{(a+1)p^{5a/6}}{\sigma(p^a)}.
\]
The dyadic factor has maximum 3 at \(t=6\): its consecutive ratios
decrease, with \((9/8)^6>2\) and \((10/9)^6<2\). For odd inert
\(p<64\), only exponents \(0\le a\le4\) need be tested. Indeed
the ratio of consecutive sixth powers for \(a\ge4\) is less than
\(((a+2)/(a+1))^6/p\le(6/5)^6/3<1\).
For \(p\ge64\), the sixth power is at most \((a+1)^6/p^a\le1\).
The product including the dyadic factor is
\[
 \frac{732525035467020016584703898578215969335561207}
      {2777251580551473425917614639931925200896}
                                      <(81/10)^6.
\]
These arguments prove the Euler bounds on the full infinite domains.

For inert cores choose \(P=n^{1/4}\) in (M19.15). By
(M19.16)–(M19.17), the sum of the subtracted terms is at most
\[
 \frac9{25}+\frac{7227}{760}n^{-1/12}
                         +\frac{62793}{200}n^{-1/3}.
\]
For every \(n\ge2^{60}\), this is at most
\[
 \frac9{25}+\frac{7227}{760\cdot32}
                     +\frac{62793}{200\cdot2^{20}}<\frac45.
                                                        \tag{M19.18}
\]
Hence every such inert core is in \(C+C\).

There are two additional families needed for descent. If \(p\) is a
split prime and \(n=4p\), then directly from (M19.9),
\(E(n)\le(20/3)n^{-1/2}\). If \(n=20p\), \(p\ne5\), then
\(E(n)\le(100/9)n^{-1/2}\); here the odd part is \(5p\).
Thus the latter bound covers both families once \(n\ge2^{30}\).
With \(P=n^{1/3}\), the subtracted terms in (M19.15) are at most
\[
 \frac9{25}+\frac{201}{19}n^{-1/6}+337n^{-1/2}
 \le\frac9{25}+\frac{201}{19\cdot32}+\frac{337}{2^{15}}
                                                   <\frac45.
                                                        \tag{M19.19}
\]
Every member of these two thin families at least \(2^{30}\) therefore
belongs to \(C+C\). The isolated case \(p=5\) in the second family
is \(n=100\), inside the finite coverage range.

## 6. The complete interval bridge

For each \(k=34,\ldots,59\), let
\[
 L_k=2^k,\quad U_k=2^{k+1}-1,\quad d_k=M(U_k),\quad
 e_k=\frac{81}{10\lfloor\sqrt[3]{L_k}\rfloor},\quad
 P_k=\left\lfloor\sqrt{\frac{8d_k}{3e_k}}\right\rfloor.
\]
Then \(M(n)\le d_k\) and \(E(n)<e_k\) for each inert core in
the interval. The divisor maximum \(M(U_k)\) is computed exactly:
remove nonsplit factors, and assign decreasing positive exponents to
the smallest split primes. Exchanging exponents \(a<b\) at primes
\(p<q\) reduces the integer by \((p/q)^{b-a}<1\) without changing
its divisor product. Thus all maxima are attained on initial prime
prefixes with nonincreasing exponents. The recursive enumeration
stops before exceeding \(2^{60}\); the retained prime prefix has
product greater than that bound. This accounts for every possible
maximizer, including a single large exponent.

Set \(R_0=65536\), \(D_0=2^{64}\), sieve the primes through \(R_0\)
exactly, and define
\[
 A_k=\#\{p\in\mathcal I:p\le P_k\},\qquad
 V_k=\frac1{D_0}
       \sum_{\substack{p\in\mathcal I\\P_k<p\le R_0}}
                    \left\lceil\frac{D_0}{p^2}\right\rceil
                    +\frac9{38R_0}+\frac9{R_0^2}.
\]
Every \(P_k<R_0\), so \(A_k=\Pi(P_k)\) and \(V_k\ge T(P_k)\).
The exact rational certificate checks all 26 inequalities
\[
 \frac9{25}+4d_kV_k
       +e_k\left(\frac{12}{5}+\frac32A_k\right)
                                <\frac45-\frac1{100}.    \tag{M19.20}
\]
The worst interval is \(k=34\), where \(d_k=288\), \(P_k=494\),
and \(A_k=47\). Each exact numerator and denominator is retained.
Together (M19.18) and (M19.20) prove that **every even inert core at
least \(2^{34}\)** belongs to \(C+C\).

## 7. Exhaustive finite coverage and descent

The checker covers all even 19-unit targets from 2 through \(2^{34}\),
apart from 4 and 20. It uses consecutive intervals of target span
\(2^{25}\), so no giant global bitset is required. The permitted
short summands are the 657 members of \(C\cap[1,4999]\).

Here is the membership lemma used to construct each segment. For an
odd positive integer \(b\le H\),
\[
 b\in C\quad\Longleftrightarrow\quad
 (b/19)=1\ \text{and no inert prime }p\le\sqrt H
                                      \text{ divides }b. \tag{M19.21}
\]
The forward implication is immediate. For the reverse implication, a
remaining inert prime factor must exceed \(\sqrt H\), so there is at
most one and its exponent is one. All other factors are split, while
19 is excluded by the character. Such an inert factor would force
\((b/19)=-1\), a contradiction. This proves exact membership on
every segment, without a primality assumption on its elements.

For a target interval \([L,H]\) with even endpoints, construct the
odd membership flags from \(b_0=\max(1,L-a_{\max})\) to \(H-1\),
where \(a_{\max}\) is the largest short summand. The 19-periodic
character mask followed by deletion of multiples of the small inert
primes gives precisely (M19.21). Pack these flags into an integer
\(B\), whose bit \(j\) records \(b_0+2j\in C\). For a short
summand \(a\), shift \(B\) right by \((L-a-b_0)/2\), using a left
shift for a negative value. Its bit \(i\) then records exactly
\(L+2i-a\in C\). The OR of these shifts therefore provides genuine
positive pairs \(a+b=L+2i\). Flags outside the interval are masked.
The checker verifies every required target bit and checks consecutive
endpoint coverage through \(2^{34}\). The short-summand cutoff is
a successful certificate, not an assumed theorem about short sums.

Membership is independently checked against complete primitive
coordinate enumeration on a smaller range; the checker also tests
deliberately false membership and missing target bits. The exact
per-segment endpoints, required-bit counts, number of used summands,
and membership/coverage hashes are retained. The algorithm recomputes
and compares the actual masks; hashes alone are not the evidence.

To conclude, suppose a least even 19-unit target \(n\notin C+C\)
outside \(\{4,20\}\) exists. Finite coverage forces \(n>2^{34}\).
If a split prime \(p\) divides \(n\), then \(n/p\) is a smaller
even 19-unit target, at least 2. Unless \(n/p\in\{4,20\}\),
minimality gives a \(C+C\) representation that can be multiplied by
\(p\), because \(C\) is multiplicatively closed. The exceptional
possibilities give \(n=4p\) or \(n=20p\), both excluded by
(M19.19). If no split prime divides \(n\), every prime dividing it
is inert, so (M19.18)–(M19.20) give the same contradiction.

Finally \(C\cap[1,19]=\{1,5,7,11,17\}\), proving directly that
4 and 20 are absent from \(C+C\). This establishes (M19.1).
For a 19-unit target the same representations lie in \(P+P\);
the extra core exception 20 is supplied by \(20=1+19\), while 4
remains impossible. For targets divisible by 19, (M19.4) transports
the two core exceptions to 76 and 380. The parity and 361 exclusions
were already necessary. This proves (M19.2). \(\square\)

## Evidence and scope

The [checker](verify_primitive_norm_minus19.py) and
[receipt](primitive_norm_minus19_checks.json) reconstruct the finite
Euler maxima, lattice inequality, residue-class tails, all 26 rational
interval bounds, the two infinite-range numerical margins, and every
finite target through \(2^{34}\). Run from the repository root:

```sh
python3 research/arithmetic/verify_primitive_norm_minus19.py \
  --output research/arithmetic/primitive_norm_minus19_checks.json
```

It uses exact integer/rational arithmetic and the Python standard
library; optimized Python execution is rejected. The finite scan is
segmented and prints progress. The written proof supplies the infinite
spectral input, its normalization, and the finite-to-infinite descent;
the checker does not establish the classical automorphic theorems.
The [source ledger](primitive_norm_minus19_source_ledger.json) records
the inspected primary statements and precise reused proof sources.
