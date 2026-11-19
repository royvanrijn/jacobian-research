# Primitive norm sums in all imaginary quadratic orders

Let \(K\) be an imaginary quadratic field of fundamental discriminant
\(D_0<0\), and fix the order
\(\mathcal O=\mathbb Z+f\mathcal O_K=\mathbb Z+\mathbb Z\omega\), of
discriminant \(\Delta=f^2D_0\). Write \(G=\operatorname{Pic}(\mathcal O)\)
and \(h=|G|\). Primitivity throughout means primitiveness in this order:
\[
 P_{\mathcal O}=\{Q(x,y)>0:(x,y)\in\mathbb Z^2,\ \gcd(x,y)=1\},
 \qquad Q(x,y)=\operatorname N(x+y\omega).
\]
Let \(\mathscr V_p=Q(\mathbb Z_p^2\setminus p\mathbb Z_p^2)\), and let
\(\Lambda_{\mathcal O}\) be the positive integers belonging to every
\(\mathscr V_p\). A target \(N\ge2\) is locally admissible if
\(N\in\mathscr V_p+\mathscr V_p\) for every prime \(p\).
Counts \(R_P(N)\) and \(R_\Lambda(N)\) count ordered pairs of positive
numerical values with sum \(N\), without coordinate or ideal multiplicity.

Put \(S=\{p:p\mid2\Delta\}\), \(\chi=\chi_{D_0}\), and
\[
 J_{\mathcal O}(N)=
 \prod_{\substack{q\mid N\\q\notin S,\ \chi(q)=-1}}\frac{q-1}{q-2}.
 \tag{O.1}
\]

**Theorem OrderValueCounts.** For every fixed imaginary quadratic order,
there are positive constants \(c_{\mathcal O},C_{\mathcal O},N_{\mathcal O}\)
such that every locally admissible \(N\ge N_{\mathcal O}\) satisfies
\[
 c_{\mathcal O}\frac{N J_{\mathcal O}(N)}{\log N}
 \le R_P(N)\le R_\Lambda(N)
 \le C_{\mathcal O}\frac{N J_{\mathcal O}(N)}{\log N}.
 \tag{O.2}
\]
Every nonadmissible target has count zero. There is no restriction on the
ring class group, its even part, or the conductor. The estimates are
uniform in the permitted valuations of \(N\); the order is fixed.

**Theorem OrderClassError.** If \(2G=0\), then
\(P_{\mathcal O}=\Lambda_{\mathcal O}\). Otherwise, for each fixed
\(0<\delta<\delta_*(2G)\), uniformly for all \(N\ge3\),
\[
 0\le R_\Lambda(N)-R_P(N)
 \ll_{\mathcal O,\delta}
        \frac{N J_{\mathcal O}(N)}{(\log N)^{1+\delta}}.
 \tag{O.3}
\]
Here, for a nonzero finite abelian group \(B\),
\[
 \delta_*(B)=\min_{\substack{\ell\mid |B|\\\ell\ {\rm prime}}}
                    \frac{\ell-1}{2\ell}\ge\frac14.
\]
In particular the exponent \(\delta=1/5\) works for every fixed order;
the implied constants still depend on the order. No endpoint estimate
at \(\delta_*(2G)\) is asserted.
Thus on sufficiently large admissible targets the ratio \(R_P/R_\Lambda\)
is \(1+O_{\mathcal O,\delta}((\log N)^{-\delta})\).

These theorems extend the
[maximal-order result](PRIMITIVE_NORM_MAXIMAL_ORDERS.md).
The analytic inputs are its finite-character positive unweighting and
exact-one prime-square sieve. The conductor extension consists of the
local-ideal/genus interface and a finite cover of additive local targets.
It does not supply an unweighted pointwise leading constant, a numerical
onset, uniformity in varying orders, or an exact global exception list.
The [leading-constant reduction](PRIMITIVE_NORM_LEADING_CONSTANT.md)
separates that further problem from (O.2).

## 1. Proper invertible ideals at conductor primes

For an invertible fractional \(\mathcal O\)-ideal \(I\), write
\(I_p=i_p\mathcal O_p\), with \(i_p\in K_p^\times\).
Local principality, including at conductor primes, is exactly the
idele description in Khayutin, Definition 2.8. A nonzero principal
ideal is invertible even when its norm is divisible by the conductor.

For integral \(I\), let \(n(I)=[\mathcal O:I]\), and extend to fractional
ideals by covolume. Multiplication by \(i_p\) on the rank-two lattice has
determinant \(\operatorname N(i_p)\); consequently
\[
 v_p(n(I))=v_p(\operatorname N(i_p)),\qquad
 n(IJ)=n(I)n(J),\qquad I\overline I=n(I)\mathcal O.
 \tag{O.4}
\]
The last equality holds locally because
\(\operatorname N(i_p)/n(I)\) is a rational \(p\)-adic unit. Conjugation
therefore inverts \(G\). These statements do not require prime
factorization of conductor ideals.

Finitely many locally principal lattices \(i_p\mathcal O_p\), with
\(\mathcal O_p\) elsewhere, glue to an invertible fractional ideal.
For an elementary construction, place all of them between
\(m\mathcal O_p\) and \(m^{-1}\mathcal O_p\), glue their residue submodules
by CRT, and take the inverse image in \(m^{-1}\mathcal O\).
Gluing the inverse lattices verifies invertibility.
Integrality and the condition \(I\not\subset p\mathcal O\) are local.
We call such an integral ideal primitive.

Define the order-unit norm group
\[
 U_p=\operatorname N(\mathcal O_p^\times)\subset\mathbb Z_p^\times.
\]
It contains all rational unit squares. The quotient by \(U_p\) is
finite of exponent at most two, and is trivial outside \(S\).

**Lemma OrderGenusKernel.** The homomorphism
\[
 \Theta:G\longrightarrow
       \prod_{p\in S}\mathbb Z_p^\times/U_p,\qquad
 \Theta_p([I])=\frac{n(I)}{\operatorname N(i_p)}\pmod{U_p}
 \tag{O.5}
\]
is well defined and has kernel \(2G\).

**Proof.** The ratio is a unit by (O.4). A change of local generator
multiplies it by a unit norm. Global multiplication by \((\alpha)\)
changes numerator and denominator by the same positive element norm.
Thus the homomorphism is well defined; its target kills \(2G\).

Conversely, suppose \(\Theta([I])=1\). At every finite prime choose a
unit \(u_p\in\mathcal O_p^\times\) with
\(\operatorname N(i_pu_p)=n(I)\). Such a choice also exists outside \(S\).
The positive rational \(n(I)\) is therefore a field norm everywhere.
The cyclic Hasse norm theorem gives \(\alpha\in K^\times\) of this norm.
The elements \(i_pu_p/\alpha\) have norm one. Quadratic Hilbert 90 writes
them as \(z_p/\overline z_p\): for \(c\ne-1\), take \(z=1+c\);
for \(c=-1\), take \(z=\sqrt{D_0}\).
At almost every prime a solution can be chosen to be a unit. In an
unramified field normalize its valuation by a rational uniformizer;
in the split algebra use \(z=(c_1,1)\).
The local \(z_p\mathcal O_p\) thus glue to an invertible fractional
ideal \(J\), with \(I/(\alpha)=J/\overline J\).
By (O.4), \([I]=2[J]\). This is the local-ideal form of Khayutin's
general-order principal genus theorem, Proposition A.1. \(\square\)

**Lemma OrderPrimitiveShadow.** For a positive integer \(n\),
\[
 n\in\Lambda_{\mathcal O}
 \quad\Longleftrightarrow\quad
 \text{there exists a primitive integral invertible ideal \(I\)
       with \(n(I)=n\) and \([I]\in2G\)}.
 \tag{O.6}
\]
Replacing \(2G\) by the principal class characterizes \(P_{\mathcal O}\).

**Proof.** Primitive local elements of exact norm \(n\) generate local
ideals that glue to an integral primitive ideal of numerical norm \(n\).
Its (O.5) invariant is one. Conversely, for an ideal on the right,
multiply its local generator by a unit of norm
\(n/\operatorname N(i_p)\), using (O.5). This gives a primitive local
element of exact norm \(n\). A principal primitive integral ideal has
an integral primitive generator, whose norm is its ideal norm.
\(\square\)

The existence quantifier matters. For \(\Delta=-100\), the two reduced
form classes are \((1,0,25)\) and \((2,2,13)\); they form a group \(C_2\)
and hence distinct genera. Both represent \(25\) primitively:
\[
 25=Q_{(1,0,25)}(0,1)=Q_{(2,2,13)}(2,1).
\]
The second class represents \(2\), while the principal class does not.
Thus a conductor-divisible numerical norm need not determine an ideal
genus. The local generator in (O.5) retains the missing information.

## 2. Saturation by good ring-class primes

Put \(B=2G\). For each inverse pair \(\{g,-g\}\subset G\) with
\(2g\ne0\), let \(\Omega_g\) be the rational split primes outside
\(S\) whose proper prime-ideal classes are \(g,-g\), and set
\[
 m_g=\operatorname{ord}(2g)-1,\qquad
 k_g(n)=\#\{q\in\Omega_g:v_q(n)=1\},\qquad
 B_n=\langle2g:k_g(n)\ge m_g\rangle\le B.
\]
Prime ideals here are prime to the conductor and invertible. Each
inverse pair is counted once. If \(B\ne0\), write \(\mathcal H_B\)
for its maximal proper subgroups and define
\[
 \begin{split}
 \ell_M&=[B:M],\qquad
 \Omega_M=\bigcup_{\{g,-g\}:\,2g\notin M}\Omega_g,\\
 K_M&=\sum_{\{g,-g\}:\,2g\notin M}(m_g-1),\qquad
 k_M(n)=\#\{q\in\Omega_M:v_q(n)=1\}.
 \end{split}
\]

**Lemma OrderClassSaturation.** If \(n\in\Lambda_{\mathcal O}\) and
\(B_n=B\), then \(n\in P_{\mathcal O}\). In particular, for
\(B\ne0\),
\[
 \Lambda_{\mathcal O}\setminus P_{\mathcal O}
 \subset\bigcup_{M\in\mathcal H_B}\{n:k_M(n)\le K_M\}.
 \tag{O.7}
\]
If \(B=0\), the two individual sets coincide. The sufficient condition
\(k_g(n)\ge m_g\) for every inverse pair also remains valid.

**Proof.** Start with a primitive integral invertible ideal from (O.6),
whose class belongs to \(B\). At every good split prime, primitivity
places the entire prime power in exactly one conjugate prime ideal.
From each saturated family reserve \(m_g\) exponent-one primes and
initially orient all of them as \(g\). This changes the class only
by an element of \(B\). Reversing \(j\) of these orientations gives
the increments \(-2jg\), for \(0\le j\le m_g\), which fill
\(\langle2g\rangle\). The families are disjoint, so their independent
choices fill \(B_n\). If \(B_n=B\), they cancel the residual class.
Every conductor-local ideal, numerical norm and primitivity condition
is preserved. The resulting principal ideal has a primitive integral
generator.

For a hidden failure, \(B_n\) is therefore proper and is contained in
some \(M\in\mathcal H_B\). Every family with \(2g\notin M\) must
be unsaturated. Their disjoint union gives \(k_M(n)\le K_M\), proving
(O.7). For \(B=0\), (O.6) gives equality directly. \(\square\)

The fixed-order class-prime Mertens input remains
\[
 \sum_{\substack{q\le y\\q\in\Omega_g}}\frac1q
       =\frac1h\log\log y+O_{\mathcal O}(1)
       \qquad(2g\ne0).
 \tag{O.8}
\]
Bourgain--Fuchs, Section 2, part (5), Lemma 2, equation (5.2),
printed p.15, applies to a fixed primitive form of negative discriminant.
Their introduction expressly includes nonfundamental discriminants
\(-4a^2\). The identities \(P_C=P_{C^{-1}}\) and
\(\epsilon(C)=1\) for a nonambiguous class give the inverse-pair
coefficient \(1/h\), where \(h=|\operatorname{Pic}(\mathcal O)|\).
For each fixed discriminant the prime-counting error is
\(O_{\mathcal O}(x e^{-c\sqrt{\log x}})\), for all sufficiently
large \(x\). Partial summation proves (O.8), since
\[
 \int_{x_0}^{\infty}e^{-c\sqrt{\log t}}\frac{dt}{t}<\infty.
\]
The finite initial range and deletion of \(S\) change only the bounded
constant. A qualitative Chebotarev density alone would not supply
this bounded remainder. No distribution result at growing nonabelian
moduli is invoked.

For \(M\in\mathcal H_B\), the prime index \(\ell_M\) divides
\(|B|\). The doubling map \(G\to B\) has equal fibres, so exactly
\(h(\ell_M-1)/\ell_M\) classes have \(2g\notin M\). None is
ambiguous. They form \(h(\ell_M-1)/(2\ell_M)\) inverse pairs, and
summing the finitely many instances of (O.8) gives
\[
 \sum_{\substack{q\le y\\q\in\Omega_M}}\frac1q
 =d_M\log\log y+O_{\mathcal O}(1),\qquad
 d_M=\frac{\ell_M-1}{2\ell_M}\ge\delta_*(B)\ge\frac14.
 \tag{O.8a}
\]
Every prime divisor of \(|B|\) occurs as a maximal-subgroup index,
by reducing a cyclic primary factor modulo that prime. Thus
\(\min_Md_M=\delta_*(B)\). The proper-subgroup mechanism in the
exceptional-value problem is classical; see Bourgain--Fuchs, Theorem 1
and Section 2, part (5). The elementary fixed-group proof above is
used here on a binary line; their squarefree exceptional-value theorem
is not imported as this class-error estimate.

To prove OrderClassError, use the direct-integer random-mark proof of
MixedPrimeSquareSieve in the maximal-order note, Section 3, with
\(\Omega_g\) replaced by \(\Omega_M\) and the fixed primes \(S\)
deleted. For \(a+b=N\) in the individual shadow, every remaining
inert prime excludes \(a=0,N\), with one root when it divides \(N\).
A marked prime \(q\in\Omega_M\) excludes exactly \(v_q(a)=1\),
namely \(q-1\) classes modulo \(q^2\). Its density is
\(f_*(q)=(q-1)/q^2=q^{-1}+O(q^{-2})\). Higher powers remain
permitted, even when \(q\mid N\). Neither rule requires maximality
at a deleted prime.

The full sieve budget is unchanged. Its finite Selberg support has
prime-index product \(d\le z\); the mixed CRT moduli are at most
\(d^2\), with total endpoint \(O_\epsilon(z^{4+\epsilon})\).
The positive Markov Euler truncation uses \(y=z^{1/16}\).
Independently mark the primes in \(\Omega_M\cap[1,y]\) with probability
\(t<1\). An integer with \(k_M(a)\le K_M\) survives with probability
at least \((1-t)^{K_M}\), a fixed positive constant. The finite
expectation, using (O.8a), bounds each deficiency by
\[
 \ll_{\mathcal O,t}
 \frac{N J_{\mathcal O}(N)}{(\log y)^{1+td_M}}
       +O_\epsilon(z^{4+\epsilon}).
\]
Choose \(z=N^{1/10}\), \(\epsilon=1/2\), and
\(t\in(0,1)\) so that \(\delta<t\delta_*(B)\).
The endpoint is \(O(N^{9/20})\), absorbed in (O.3) because
\(J_{\mathcal O}(N)\ge1\). Sum the finitely many maximal subgroups
and both orientations of the additive pair, then enlarge the constant
to cover all \(N\ge3\). All target valuations remain in the exact
inert-root product (O.1). This proves OrderClassError. In particular,
\(t=9/10\) gives \(td_M\ge9/40>1/5\), so the exponent \(1/5\)
is available for every fixed order, with order-dependent constants.
No new zero-free or graph input is used, and no estimate at the
endpoint \(\delta_*(B)\) is asserted.

## 3. A finite additive cover despite unbounded conductor valuations

Write \(Q(x,y)=x^2+bxy+cy^2\); its Hessian determinant is \(-\Delta\).
For \(p\in S\), put \(t_p=v_p(\Delta)\) and \(h_p=2t_p+1\).

**Lemma OrderLocalCosets.** The primitive local value image
\(\mathscr V_p\) is a union of full cosets modulo \(p^{h_p}\).
Every \(N\in\mathscr V_p+\mathscr V_p\) has a representation by
nonzero values \(a,b\in\mathscr V_p\) with
\[
                  v_p(a),v_p(b)\le 2t_p+2.
 \tag{O.9}
\]

**Proof.** At a primitive vector \(v\), the identity
\(\operatorname{adj}(\operatorname{Hess}Q)\nabla Q(v)=-\Delta v\)
implies some partial derivative has valuation \(s\le t_p\).
If \(n'\equiv Q(v)\pmod{p^{2t_p+1}}\), fix the other coordinate and
apply strong one-variable Hensel: the residual valuation exceeds
\(2s\). The lifted coordinate agrees with the old one modulo
\(p^{s+1}\), so primitivity persists. This proves the coset assertion.

Starting from \(a+b=N\), replace the values by
\(a+p^{h_p}u,b-p^{h_p}u\). They remain in their full value cosets.
A summand of valuation below \(h_p\) keeps its valuation. Otherwise
exclude its one cancelling residue for \(u\) modulo \(p^2\);
the new valuation is at most \(h_p+1\). At most two residues are
excluded, and \(p^2>2\). This gives (O.9). \(\square\)

Thus every local additive target has bounded summand valuations.
Individual valuations need not be bounded: at a split conductor prime
zero can belong to the primitive local image. For example,
\(x^2+18y^2\) at \(3\) allows primitive norms \(27\) and \(18\),
but not \(9\). Reusing the maximal-order ramified exponent rule
or allowing every unit part at the conductor would be incorrect.
The full, proved local tables and their dyadic cases are in
[the conductor local-profile note](PRIMITIVE_NORM_CONDUCTOR_LOCAL_PROFILES.md).

Let
\[
 L=8\prod_{\substack{p\in S\\p\text{ odd}}}p,\qquad
 Z=\prod_{p\in S}p^{\,2v_p(\Delta)+2},\qquad
 H=\ker(\chi:(\mathbb Z/L)^\times\to\{\pm1\}).
 \tag{O.10}
\]
For any positive \(r\) supported on \(S\), define the actual value filter
\[
 A_r=\{a\in(\mathbb Z/L)^\times:ra\in\mathscr V_p
                                      \text{ for every }p\in S\}.
 \tag{O.11}
\]
Any integer representative of \(a\) can be used. Multiplication of an
element by a rational local unit multiplies its norm by a unit square.
Equal unit residues modulo odd \(p\), or modulo \(8\) at two, differ by
a unit square. Thus (O.11) is well defined and is a union of states of
the full quadratic vector
\[
 \Xi(a)=\big((a/p)_{p\in S,\ p\text{ odd}},
                     \chi_{-4}(a),\chi_8(a)\big).
 \tag{O.12}
\]
These are unit-square states; they need not be independent Picard-genus
characters.

Every \(A_r\) lies in \(H\). Indeed put
\(\lambda_p(x)=(x,D_0)_p\). If \(a\in A_r\), then
\(\lambda_p(a)=\lambda_p(r)\) at every \(p\in S\), since \(ra\) is a
field norm there. Reciprocity gives
\(\prod_{p\in S}\lambda_p(r)=1\): all other symbols and the infinite
symbol are one. For a positive unit representative \(a\),
\(\prod_{p\in S}\lambda_p(a)=\chi(a)\). Hence \(\chi(a)=1\).

Let \(\mathcal C\) be the semigroup generated by the primes outside
\(S\) split in \(K\). Then
\[
        r\{c\in\mathcal C:c\bmod L\in A_r\}
                                  \subset\Lambda_{\mathcal O}.
 \tag{O.13}
\]
This uses the actual conductor value images; no generator of norm \(r\)
is assumed. The union over all \(S\)-supported \(r\) is an exhaustive,
disjoint decomposition of the individual shadow. Only \(r\mid Z\)
will be used for the additive lower bound.

For \(r,s\mid Z\), let \(g=(r,s)\), \(A=r/g\), \(B=s/g\). If \(g\mid N\),
put \(M=N/g\) and
\[
 \Gamma_{r,s}(N)=\frac L{|H|^2}
     \#\{(a,b)\in A_r\times A_s:Aa+Bb=M\pmod L\};
 \tag{O.14}
\]
otherwise set \(\Gamma_{r,s}(N)=0\).
The original sets \(A_r,A_s\) remain after cancelling \(g\).

**Lemma OrderLocalProfiles.** A target is locally admissible exactly
when some \(r,s\mid Z\) have \(\Gamma_{r,s}(N)>0\).
This is an explicit finite criterion modulo \(ZL\).

**Proof.** Use (O.9) at each bad prime. The exact valuations of the
two local summands select \(r,s\mid Z\). Divide each summand by its
full prefix; its local quotients are units. Their residues assemble
by CRT to \(a,b\) in (O.14), and unit-square stability puts them in
\(A_r,A_s\).

Conversely, a root of (O.14) lifts to an exact linear equation over
each bad \(\mathbb Z_p\), preserving the unit squareclasses, since
\(A,B\) are coprime and at least one is a \(p\)-unit. The two values
remain in \(\mathscr V_p\).
Outside \(S\) the order is maximal and unramified at an odd prime.
At a split prime every value is a primitive norm, using the product
coordinates \((1,n)\). At an inert prime every unit is a primitive
norm, and every residue of \(N\) is a sum of two units since \(p>2\);
lift and use surjectivity of the unramified unit norm.
Finally \(g\mid Z\), so \(N\bmod ZL\) determines \(g\mid N\) and
\(N/g\bmod L\). This proves the finite criterion. \(\square\)

The criterion can be computed from \(Q\) alone. For each bad prime,
enumerate primitive pairs modulo \(p^{h_p}\); the first part of (O.9)
proves that this residue image gives the full \(p\)-adic value image.
Then determine (O.11), and test the finitely many roots in (O.14).
There is no empirical stabilization assumption.

## 4. Applying the existing finite-character counting theorem

The proof of FiniteGenusFilter in Section 5 of the maximal-order note
works for any fixed finite vector of quadratic characters on the
unit group modulo \(L\). In particular it applies to \(\Xi\) in (O.12).
For fixed coprime \(S\)-supported coefficients \(A,B\), and a pair of
states \(\alpha,\beta\), it gives
\[
 \#\{Ax+By=M:x,y\in\mathcal C,\ \Xi(x)=\alpha,\Xi(y)=\beta\}
             \asymp\frac{M J_{\mathcal O}(M)}{\log M}
 \tag{O.15}
\]
whenever the corresponding finite pair density is positive.

Here are the interface checks relevant to the extra conductor characters.
Expand the first-variable filter, including the second state, modulo
\(BL\). On the splitting core the two characters \(1,\chi\) have the
same principal mean; all other twists have both character factors
nonprincipal. The weighted proof therefore gives the actual root-count
density (O.14), not a product of independent marginal means.
The varying kernel indices and every cofactor are units modulo \(BL\).
Multiplication by a removed cofactor transports the entire vector
\(\Xi\), bijects the finite solution pairs, and preserves their density.
Thus the full normalized-shift weighted family required for positive
restoration survives. Good split primes with nontrivial \(\Xi\)-state
are retained in the restoration.

Both fixed coefficients may share primes with \(L\). The source proof
already checks this through BGLM v5, Proposition 5.1: it does not require
\((a_2,c)=1\), and the varying short variables/moduli are \(BL\)-units.
Every root modulo \(L\) has exactly \(B\) pair lifts modulo \(BL\), and
\(|H_{BL}|=B|H|\); this is precisely the normalization (O.14).
Positive density forces \((M,AB)=1\).
No growing conductor, new class-group dispersion, or new distribution
estimate is being imported in (O.15).

A union such as \(A_r\) is a disjoint finite union of the \(\Xi\)-states.
If (O.14) is positive, at least one contained state pair has positive
density, so (O.15) supplies a lower bound. Multiply the cofactor values
numerically by \(r,s\) and use (O.13). The map is injective on numerical
value pairs, and
\[
 J_{\mathcal O}(N/g)=J_{\mathcal O}(N)
\]
because \(g\) is \(S\)-supported. There are only finitely many
\(r,s\mid Z\), with fixed coefficients and positive finite densities.
Their minimum positive constant and maximum threshold give
\[
 R_\Lambda(N)\gg_{\mathcal O}
                 N J_{\mathcal O}(N)/\log N
\]
on every sufficiently large locally admissible target.

For the upper bound, forget all bad local conditions and sieve the
line \(a+(N-a)=N\) only at the good inert primes. The same direct
Selberg estimate gives \(R_\Lambda(N)\ll_{\mathcal O}NJ_{\mathcal O}/\log N\).
No infinite conductor-prefix tail is truncated in this upper bound.
Combining this lower bound with (O.3), its error is a smaller logarithmic
order, so it can be absorbed in half the lower bound for sufficiently
large \(N\). This proves OrderValueCounts.

### 4.1. Primitive forms and a translated genus coset

The same arithmetic applies to every primitive positive definite integral
binary quadratic form of a fixed discriminant \(\Delta\), not just the
principal norm form. Write
\[
 P_Q=\{Q(x,y)>0:(x,y)\in\mathbb Z^2,\ \gcd(x,y)=1\},\qquad
 \mathscr V_{p,Q}=Q(\mathbb Z_p^2\setminus p\mathbb Z_p^2),
 \quad \Lambda_Q=\mathbb Z_{>0}\cap\bigcap_p\mathscr V_{p,Q}.
\]
Numerical values are counted once. For two sets, let
\(R_{T_1,T_2}(N)=\#\{(a,b)\in T_1\times T_2:a+b=N\}\), and write
\(R_{Q_1,Q_2}=R_{P_{Q_1},P_{Q_2}}\). The term genus means local
integral equivalence at every prime, with the positive real condition.
It preserves primitive local value sets. For the comparisons in this
section, set \(\delta_*(0)=+\infty\); the corresponding class error is
exactly zero, rather than a new quantitative density estimate.

For \(Q=ax^2+bxy+cy^2\), put
\[
 \beta=(b+\sqrt\Delta)/2,\qquad
 A=\mathbb Za+\mathbb Z\beta,\qquad c_Q=-[A]\in G.
 \tag{F.1}
\]
Here \(\mathcal O=\mathbb Z[\beta]\), \(n(A)=a\), and
\(Q(x,y)=N(ax+\beta y)/a\). The lattice \(A\) is an ideal and
\[
 A\overline A=a(a,\beta,\overline\beta,c)\mathcal O=a\mathcal O,
\]
since \(\beta+\overline\beta=b\) and \(\gcd(a,b,c)=1\). Thus
\(A^{-1}=\overline A/a\), so it is proper invertible even when
\(a\) shares a prime with the conductor. This also fixes the class sign.
Khayutin's Lemma 10.2, printed p.240, gives the normalized ideal-form
construction and integral-ideal parametrization without a conductor
coprimality restriction. Only that elementary arithmetic assertion and
its proof are used, not the surrounding shifted-convolution theorem.

**Lemma FormIdealShadow.** For every positive integer \(n\),
\[
 \begin{split}
 n\in P_Q&\iff\exists I\text{ primitive integral invertible}:\quad
                      n(I)=n,\ [I]=c_Q,\\
 n\in\Lambda_Q&\iff\exists I\text{ primitive integral invertible}:\quad
                      n(I)=n,\ [I]\in c_Q+2G.
 \end{split}
 \tag{F.2}
\]

**Proof.** For \(\alpha\in A\), the ideal \(I=(\alpha)A^{-1}\) is
integral, has norm \(N(\alpha)/a\), and class \(-[A]\).
Invertibility gives \(I\subset p\mathcal O\) exactly when
\(\alpha\in pA\), which is precisely coordinate nonprimitivity in the
basis \(a,\beta\). Conversely, \([I]=-[A]\) gives
\(IA=(\alpha)\subset A\); the same identities prove the first line.

Write \(A_p=z_p\mathcal O_p\) and
\(\lambda_p=N(z_p)/a\in\mathbb Z_p^\times\), using (O.4).
The basis change to this local generator is in
\(\operatorname{GL}_2(\mathbb Z_p)\), so
\(\mathscr V_{p,Q}=\lambda_p\mathscr V_p\).
Moreover \(\Theta_p(c_Q)=\lambda_p\pmod{U_p}\).
If \(n\in\Lambda_Q\), take primitive local elements
\(\alpha_p=z_pw_p\in A_p\) with norm \(an\). Their primitive integral
ideals \(w_p\mathcal O_p\) glue as in Section 1 to an invertible
integral ideal \(I\) of norm \(n\). They are the unit lattice outside
finitely many primes. Locally
\(n/N(w_p)=\lambda_p\), hence \(\Theta([I])=\Theta(c_Q)\).
The kernel in (O.5) is \(2G\), giving the coset in (F.2).
Conversely, for such an \(I\), write \(I_p=i_p\mathcal O_p\).
The equality of the two \(\Theta\)-classes says that
\(n/(\lambda_pN(i_p))\) is the norm of a unit \(u_p\).
Then \(z_pi_pu_p\) is primitive in \(A_p\), with norm \(an\).
This proves local primitive representation, including at conductor
primes. \(\square\)

The ideal in (F.2) is existential. A conductor-divisible numerical value
need not determine a unique genus, as the earlier discriminant \(-100\)
example shows. The target class is changed, not this quantifier.

**Theorem FormClassError.** Fix two forms \(Q_1,Q_2\) of discriminant
\(\Delta\). If \(2G\ne0\), then for every fixed
\(0<\eta<\delta_*(2G)\), uniformly for \(N\ge3\),
\[
 0\le R_{\Lambda_{Q_1},\Lambda_{Q_2}}(N)-R_{Q_1,Q_2}(N)
 \ll_{\Delta,\eta}\frac{NJ_{\mathcal O}(N)}{(\log N)^{1+\eta}}.
 \tag{F.3}
\]
If \(2G=0\), then \(P_Q=\Lambda_Q\) for every such form and the
difference is zero. The exponent \(1/5\) is available in every fixed
order, with order-dependent constants and no endpoint claim.

**Proof.** Start with an ideal in \(c_Q+2G\) from (F.2). The reserved
orientations in Section 2 change its class through \(B_n\le2G\).
If \(B_n=2G\), they cancel the difference from \(c_Q\); all
conductor-local ideals remain unchanged. Thus a hidden failure has
exactly the same maximal-subgroup cover (O.7), with the same
\(\Omega_M,K_M\) and density \(d_M\).
Charge a mixed pair with an exceptional summand to that summand and a
subgroup in its cover. Both shadows forbid every good inert prime.
Mark only the exceptional summand at \(\Omega_M\). The forbidden
densities are still \(2-\mathbf1_{p\mid N}\) roots modulo an inert
prime and \(p-1\) roots modulo \(p^2\) at a marked split prime.
The exact CRT, Selberg and finite-mark argument proving (O.3) therefore
proves (F.3) unchanged, uniformly in all target valuations. No zero-free
or weighted-dilation input is used. \(\square\)

**Lemma FormGenusCoreFilters.** Let \(S'\) be any fixed finite prime set
containing all primes dividing \(2\Delta\). Set
\[
 L'=8\prod_{p\in S',\ p\text{ odd}}p,\quad H'=\ker(\chi\bmod L'),
 \quad\mathcal C'=\langle p:p\notin S',\ \chi(p)=1\rangle.
\]
For \(S'\)-supported \(r\), define the actual filter
\[
 A_{r,Q}=\{u\in(\mathbb Z/L')^\times:
                      ru\in\mathscr V_{p,Q}\text{ for all }p\in S'\}.
\]
It is a union of unit-square states, lies in \(H'\), and gives the
exact disjoint decomposition
\[
 \Lambda_Q=\coprod_{r\in\langle S'\rangle}
       r\{x\in\mathcal C':x\bmod L'\in A_{r,Q}\}.
 \tag{F.4}
\]
The separate local genus biases are retained in \(A_{r,Q}\).

**Proof.** Homogeneity and multiplication of a primitive vector by a
unit show unit-square stability, so the filter is well defined.
The local multipliers above are units and satisfy
\[
 \prod_{p\in S'}(\lambda_p,D_0)_p
          =\prod_{p\in S'}(a,D_0)_p=1.
 \tag{F.5}
\]
Their numerators are local field norms. Outside \(S'\), an invertible
ideal norm has even valuation at an inert prime, and no obstruction at
a split prime; the positive real norm has trivial infinite symbol.
Hilbert reciprocity gives (F.5). If \(u\in A_{r,Q}\), then
\((ru,D_0)_p=(\lambda_p,D_0)_p\) for every \(p\in S'\).
The positive \(S'\)-supported integer \(r\) has product of bad symbols
one, whereas a positive unit representative \(u\) has product
\(\prod_{p\in S'}(u,D_0)_p=\chi(u)\). Multiplying gives
\(\chi(u)=1\). Individual \(\lambda_p\) may have nontrivial genus
bias; their product, rather than each bias separately, is used.
Outside \(S'\), the primitive local set is \(\mathbb Z_p\) at split
primes and \(\mathbb Z_p^\times\) at inert primes. This proves (F.4).
\(\square\)

**Theorem GenusBinaryInvariance.** Fix four forms of discriminant
\(\Delta\), with \(Q_i,Q'_i\) in the same genus for each \(i=1,2\).
The genera in the two slots may be different. For every fixed
\(0<\eta<\delta_*(2G)\), uniformly on their common sufficiently large
locally admissible targets,
\[
 R_{Q_1,Q_2}(N)=R_{Q'_1,Q'_2}(N)
                         (1+O((\log N)^{-\eta})).
 \tag{F.6}
\]
Both counts have order \(NJ_{\mathcal O}(N)/\log N\); outside the local
set both counts are zero. When \(2G=0\), the counts are equal for every
target. In that exact case one may use any \(\eta>0\) in the displayed
comparison; it makes no assertion about an arbitrarily strong remainder
for a one-variable density asymptotic. All constants and onsets depend
on the fixed forms and \(\eta\). This result uses only the classical
fixed-filter counting theorem and the class-error sieve, not ZF/G.

**Proof.** Local integral equivalences preserve primitive vectors, so
\(\Lambda_{Q_i}=\Lambda_{Q'_i}\) in each slot. Apply the preceding
filter lemma with \(S'=S\). For any primitive quadratic form the
identity
\(\operatorname{adj}(\operatorname{Hess}Q)\nabla Q(v)=-\Delta v\)
holds. Thus the full-coset/Hensel proof of (O.9) applies with
\(e_p=2v_p(\Delta)+1\). Each local value set is a union of full
cosets modulo \(p^{e_p}\). Replacing local summands by
\(a+p^{e_p}t,b-p^{e_p}t\) gives both valuations at most \(e_p+1\):
at most two residues of \(t\bmod p^2\) are excluded. CRT therefore
selects prefixes \(r,s\mid\prod_{p\in S}p^{e_p+1}\) and a positive
root in the two separate filters \(A_{r,Q_1},A_{s,Q_2}\).

The classical finite-character theorem (O.15) allows distinct unit-square
filters in \(H\). It gives a lower bound of order \(NJ/\log N\) in
such a positive profile. There are finitely many profiles and local root
states, so their positive constants have a common minimum. The unmarked
inert sieve gives the matching upper bound for the entire mixed shadow.
Both uses are independent of the pointwise binary asymptotic.
Apply (F.3) to each of the two pairs. They differ from the same mixed
shadow by \(O(NJ/\log^{1+\eta}N)\). Subtract and divide by the
classical positive lower bound; the denominator actual count is at least
half that lower bound for large \(N\). This proves (F.6).
For \(2G=0\), (F.2) makes each global set equal its shadow, proving
exact equality. \(\square\)

## 5. Imported sources and evidence limits

The conductor algebra is classical principal-genus theory. Its role
here is to preserve primitive local ideal types and provide the precise
interface (O.6), rather than claim a new genus theorem.

- I. Khayutin, *Joint equidistribution of CM points*, Annals 189 (2019),
  145--276, DOI 10.4007/annals.2019.189.1.4.
  [Primary PDF](https://annals.math.princeton.edu/wp-content/uploads/annals-v189-n1-p04-s.pdf):
  Definition 2.8, printed p.167; Appendix A, printed pp.256--261,
  especially Proposition A.1 and Lemmas A.3--A.4; Lemma 10.2,
  printed p.240, for the elementary normalized ideal-form construction.
  The explicit ideal construction and primitivity equivalence are proved
  in Section 4.1, including conductor-divisible values.
  The appendix explicitly handles nonmaximal orders. Its algebraic
  statements have no exceptional-zero hypothesis; the main article's
  equidistribution hypotheses are not used here.
- J. Bourgain and E. Fuchs,
  [arXiv:1105.3698v2](https://arxiv.org/abs/1105.3698v2),
  Lemma 2, (5.2), and Section 6: fixed negative form discriminants,
  including nonfundamental ones, yield (O.8).
- J. S. Milne, *Class Field Theory*, v4.03, VIII Theorem 3.1:
  the cyclic Hasse norm theorem; III Theorem 4.4(d) and VIII:
  Hilbert symbols and reciprocity, in the exact retained edition.
- The maximal-order note, Sections 3 and 5: complete mixed-modulus
  Selberg proof and finite-character positive restoration, with
  their retained primary analytic dependencies. These analytical
  inputs are reused, with all changes of local hypotheses exposed above.

The [source ledger](primitive_norm_all_order_source_ledger.json) pins
the inspected editions and inherited analytic proof bytes.
The [finite checker](verify_primitive_norm_order_profiles.py) and its
[receipt](primitive_norm_order_profile_checks.json) compare complete
finite coordinate images, local tables, target profiles and deliberate
corruptions on the specified orders and moduli. Those computations do
not prove the unrestricted ideal, prime-counting or dispersion inputs.
The separate [maximal-subgroup checker](verify_class_group_error_cover.py)
and [receipt](class_group_error_cover_checks.json) check only the declared
twelve finite groups, 199 saturated-family subsets, exact union densities
and three rejected corruptions. They do not prove the infinite group,
prime-distribution, conductor-local or analytic assertions.
The theorem has internal written/source review and exact finite replay;
no external human review, formal verification or historical-priority
certificate is asserted.

A source caution is confined to the unused dyadic table in Khayutin's
Lemma A.6. Its displayed first dyadic case would retain the maximal-order
unit norm group whenever \(2\Vert f\). For
\(K=\mathbb Q(\sqrt{-2})\), \(f=2\), the order-unit norms
\(x^2+8y^2\), \(x\) odd, are instead exactly \(1\bmod8\);
the maximal order has the unit norm \(3=1^2+2\cdot1^2\).
The conductor table here is proved directly, rather than imported from
that line. Definition 2.8 and Proposition A.1 precede and do not depend
on that local table.
