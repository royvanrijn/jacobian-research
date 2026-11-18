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
\(0<\delta<1/h\), uniformly for all \(N\ge3\),
\[
 0\le R_\Lambda(N)-R_P(N)
 \ll_{\mathcal O,\delta}
        \frac{N J_{\mathcal O}(N)}{(\log N)^{1+\delta}}.
 \tag{O.3}
\]
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

For each inverse pair \(\{g,-g\}\subset G\) with \(2g\ne0\), let
\(\Omega_g\) be the rational split primes outside \(S\) whose proper
prime-ideal classes are \(g,-g\). Put
\[
 m_g=\operatorname{ord}(2g)-1,\qquad
 k_g(n)=\#\{q\in\Omega_g:v_q(n)=1\}.
\]
Prime ideals here are prime to the conductor and invertible.

**Lemma OrderClassSaturation.** If \(n\in\Lambda_{\mathcal O}\) and
\(k_g(n)\ge m_g\) for every such inverse pair, then
\(n\in P_{\mathcal O}\). In particular,
\[
 \Lambda_{\mathcal O}\setminus P_{\mathcal O}
 \subset\bigcup_g\{n:k_g(n)<m_g\}.
 \tag{O.7}
\]
If \(2G=0\), the two individual sets coincide.

**Proof.** Start with an ideal in (O.6). At a good split prime, primitivity
places the entire prime power in exactly one conjugate prime ideal.
Reserve \(m_g\) exponent-one primes from each family. Changing their
orientations gives increments \(0,-2g,\ldots,-2m_gg\), which fill
\(\langle2g\rangle\). Over all inverse pairs these subgroups generate
\(2G\). Hence one can cancel the initial class in \(2G\).
Every conductor-local ideal is unchanged, as are norm and primitivity.
The resulting principal ideal has a primitive integral generator.
\(\square\)

The fixed-order class-prime Mertens input is
\[
 \sum_{\substack{q\le y\\q\in\Omega_g}}\frac1q
       =\frac1h\log\log y+O_{\mathcal O}(1)
       \qquad(2g\ne0).
 \tag{O.8}
\]
Bourgain--Fuchs, Section 2, numbered part (5), Lemma 2, equation
(5.2), printed p.15, applies to a fixed primitive form of negative
discriminant. The introduction also expressly includes nonfundamental
discriminants \(-4a^2\). Their \(P_C=P_{C^{-1}}\) and
\(\epsilon(C)=1\) for a nonambiguous class give the inverse-pair
coefficient \(1/h\), with \(h=|\operatorname{Pic}(\mathcal O)|\).
For each fixed discriminant the stated prime-counting error is
\(O_{\mathcal O}(x e^{-c\sqrt{\log x}})\), once \(x\) is large enough.
Partial summation proves (O.8), since
\[
 \int_{x_0}^{\infty}e^{-c\sqrt{\log t}}\frac{dt}{t}<\infty.
\]
The finite initial range and deletion of the primes in \(S\) change
only the bounded constant. A qualitative Chebotarev density alone
would not supply this bounded Mertens remainder. This prime input
chooses/marks good primes only; no distribution theorem at growing
nonabelian moduli is invoked.

To prove OrderClassError, apply exactly the direct-integer random-mark
proof of ClassErrorComparison in the maximal-order note, Sections 2--3.
Delete the fixed primes \(S\). For \(a+b=N\) in the individual shadow,
every remaining inert prime excludes \(a=0,N\), with one root when it
divides \(N\). A marked prime in \(\Omega_g\) excludes exactly
\(v_q(a)=1\), namely \(q-1\) classes modulo \(q^2\). Higher powers are
permitted. Neither rule uses maximality at a deleted prime.

For clarity, the quantitative budget of the reused sieve is unchanged:
the finite Selberg support has prime-index product \(d\le z\);
the mixed CRT moduli are at most \(d^2\), with total endpoint
\(O_\epsilon(z^{4+\epsilon})\).
The Markov Euler truncation uses \(y=z^{1/16}\).
Independently marking the \(\Omega_g\) primes with probability \(t<1\)
bounds a deficiency by
\((1-t)^{-(m_g-1)}(1-t)^{k_g(a)}\).
Taking its finite expectation gives
\[
 \ll_{\mathcal O,t}
 \frac{N J_{\mathcal O}(N)}{(\log y)^{1+t/h}}
       +O_\epsilon(z^{4+\epsilon}).
\]
Choose \(z=N^{1/10}\), \(\epsilon<1/2\), and
\(\delta<t/h<1/h\). The endpoint is \(O(N^{0.45})\), absorbed in
(O.3); sum the finitely many deficient families and both orientations
of the additive pair. All target valuations remain in the exact
inert-root product (O.1). This proves OrderClassError.

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

## 5. Imported sources and evidence limits

The conductor algebra is classical principal-genus theory. Its role
here is to preserve primitive local ideal types and provide the precise
interface (O.6), rather than claim a new genus theorem.

- I. Khayutin, *Joint equidistribution of CM points*, Annals 189 (2019),
  145--276, DOI 10.4007/annals.2019.189.1.4.
  [Primary PDF](https://annals.math.princeton.edu/wp-content/uploads/annals-v189-n1-p04-s.pdf):
  Definition 2.8, printed p.167; Appendix A, printed pp.256--261,
  especially Proposition A.1 and Lemmas A.3--A.4.
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
