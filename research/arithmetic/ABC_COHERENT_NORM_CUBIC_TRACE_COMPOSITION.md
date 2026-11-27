# Trace-zero composition of coherent prime edges

A nonprincipal prime edge can be recovered from a trace-zero integer in
its [reflected cubic field](ABC_COHERENT_NORM_RING_CLASS_RESIDUES.md#cubic-fields-detect-the-first-color).
For two independent lines of \(3\)-classes, multiplication of trace-zero
integers has two components in the other two cubic fields. Each shared
prime of the input norms divides exactly one of those components. This
gives a composition rule on the cubic side and a uniform estimate for
single edges with a prescribed index divisor.

## Exact trace-zero representation

Let \(K=\mathbb Q(\sqrt{-d})\) with fundamental discriminant \(D<-4\), and let
\(D^\vee\) be the discriminant of \(K^\vee=\mathbb Q(\sqrt{3d})\).
For a nonzero line \(\langle t\rangle\subset\mathrm{Cl}(K)[3]\), fix one
cubic subfield \(L_t\) of its reflected splitting field \(S_t\). Put
\[
 \Lambda_t=\{\theta\in\mathcal O_{L_t}:\operatorname{Tr}(\theta)=0\},
 \qquad Q_t(\theta)=\frac{\operatorname{Tr}(\theta^2)}6.
\]

**Proposition (trace-zero representation).** Let \(n>1\) be squarefree and
coprime to \(6D\). Up to the choices of sign and quadratic conjugation of
\(\alpha\), elements generating \(I^3=(\alpha)\), where \(I\) is an integral
ideal of norm \(n\) with class in \(\langle t\rangle\setminus\{0\}\), correspond
to the elements of \(\Lambda_t\) with \(Q_t(\theta)=n\), up to sign. The
formulas are
\[
 \operatorname{Tr}_K(\alpha)=\operatorname{Nm}_{L_t}(\theta)=u,\qquad
 4n^3-u^2=d v^2,\qquad
 [\mathcal O_{L_t}:\mathbb Z[\theta]]
     =3^{1-a}\gcd(3,d)\,|v|/\epsilon_d,                           \tag{1}
\]
where \(\alpha=(u+v\sqrt{-d})/2\), \(\epsilon_d\) is \(1\) for
\(d\equiv3\bmod4\) and \(2\) otherwise, and
\(\operatorname{disc}(L_t)=D^\vee 3^{2a}\).
In particular for \(n=pq\) with distinct primes, the trace condition
recovers an oriented nonprincipal edge with free choice of labels \(p,q\).

For the forward map, the minimal polynomial of the corresponding root in
\(L_t\) is \(T^3-3nT-u\). There is exactly one such root in the chosen cubic
subfield. Two distinct roots in \(L_t\) would put the splitting field there,
which has degree six. Newton's identities give the trace, trace square and
norm formulas. For the inverse, \(\theta\ne0\) has trace zero and degree three,
so its minimal polynomial is the same trinomial for \(u=\operatorname{Nm}(\theta)\).
Its polynomial discriminant is \(27(4n^3-u^2)>0\) and has squareclass
\(D^\vee\). The quadratic polynomial \(T^2-uT+n^3\) therefore has roots
\(\alpha,\bar\alpha\) in \(K\). Both are integral. If \(p\mid n\) also divided
\(u\), then \(\theta^3\equiv0\bmod p\mathcal O_{L_t}\). The prime \(p\) is
unramified in \(L_t\), so the quotient ring is reduced and \(\theta\) would
lie in \(p\mathcal O_{L_t}\). This contradicts \(v_p(\operatorname{Tr}(\theta^2))=1\).
Thus \((u,n)=1\). Each \(p\mid n\) splits in \(K\), since it is unramified
and its odd valuation in \(\operatorname{Nm}(\alpha)=n^3\) rules out inertia.
Its entire ideal exponent three occurs at just one of the two prime
ideals, since otherwise \(p\) would divide the trace. This gives \(I\).
Its class is nonzero by irreducibility of the trinomial and belongs to the
specified line by uniqueness of the reflected field. Changing the sign
of \(\theta\) changes the sign of \(u\), and quadratic conjugation changes
the sign of \(v\). The remaining ambiguity of orientations for \(n=pq\)
is exactly the simultaneous conjugation from the two possible roots.

Finally \(D^\vee=3d\epsilon_d^2/\gcd(3,d)^2\). Dividing the polynomial
discriminant \(27dv^2\) by the field discriminant gives the square of the
index and proves (1). In the counting and projection arguments the relevant
trace pairing is the classical trace form on the trace-zero submodule.
[Mantilla-Soler, §§4–5](https://msp.org/ant/2010/4-6/ant-v4-n6-p02-p.pdf)

## The trace lattice with integral quadratic values

**Proposition (quadratic ideal lattice).** Choose an ideal \(I\) in the
nonzero class \(t\), with \(I^3=(\alpha)\), \(N I=n_0\), and a cube root
\(\beta\) for which \(\beta+\bar\beta\in L_t\). The map
\[
 \rho\longmapsto \theta(\rho)=\beta\rho+\bar\beta\,\bar\rho
\]
identifies \(I^{-1}\) with the subset of \(\Lambda_t\) on which \(Q_t\)
is integral, and takes the scaled quadratic norm \(n_0N(\rho)\) to \(Q_t\).
In particular that subset is a lattice, with covolume \(3\sqrt{|D|}\) in
the Euclidean trace metric. If the integral trace has image \(c\mathbb Z\),
its index in \(\Lambda_t\) is \(c\gcd(3,d)/3^a\). This also implies
\(\operatorname{disc}(L_t)\mid 27|D|\).

Over \(\mathbb Q\) the map is linear. Complex conjugation fixes the sum.
The involution fixing the chosen cubic root in \(S_t\) can be lifted
to the Kummer field so it acts as quadratic conjugation on \(K\) and fixes
\(\zeta_3\), as \(S_t\cap M=K^\vee\). That lift swaps \(\beta,\bar\beta\):
it permutes the roots of the quadratic with sum \(\beta+\bar\beta\)
and product \(n_0\), and cannot fix \(\beta\) since it conjugates its
nonrational cube. It also fixes \(\theta(\rho)\). The sum thus lies in
\(L_t\) and has trace zero. The map is injective, since a nonzero
\(\rho\) with vanishing sum would give \(\beta/\bar\beta\in K\),
then \(\beta^2,\beta^3\in K\) and \(\beta\in K\). It is an isomorphism of
the two-dimensional rational spaces. The trace trinomial gives
\(Q_t(\theta(\rho))=n_0 N(\rho)\).

If \(\rho\in I^{-1}\) then \(\alpha\rho^3\) and \(\beta\rho\) are
integral, and \(n_0N(\rho)\) is the integer norm of \(\rho I\).
Conversely if \(\theta(\rho)\) and its quadratic value \(n'\) are integral,
the element \(\alpha\rho^3\) satisfies
\(T^2-\operatorname{Nm}(\theta(\rho))T+n'^3=0\). It is integral, so the
fractional ideal \(\rho I\) is integral because its cube is integral.
The determinant of \(I^{-1}\) in the complex plane is
\(\sqrt{|D|}/(2n_0)\). The Euclidean trace pairing multiplies area by
\(6n_0\) under the map. The full trace-zero lattice has covolume
\(\sqrt{3\operatorname{disc}(L_t)}/c\), obtained by dividing the full
Minkowski covolume by the perpendicular trace spacing \(c/\sqrt3\).
Taking the ratio gives the claimed index. Squaring gives
\(\operatorname{disc}(L_t)=3c^2|D|/[\Lambda_t:\theta(I^{-1})]^2\)
with \(c\mid3\), which gives the divisibility.

In these ideal coordinates the two cubic projections below come from
\(\rho_1\rho_2\) and \(\rho_1\bar\rho_2\) in the corresponding product
ideal lattices. This relates the rule to quadratic ideal multiplication.
Mantilla-Soler's Theorem 5.5 relates the trace form to a Hessian class for
cubic fields of fundamental discriminant prime to \(3\). The calculation
here uses the integral-value sublattice and applies also to the
nonfundamental reflected discriminants.
[Mantilla-Soler, Theorem 5.5](https://msp.org/ant/2010/4-6/ant-v4-n6-p02-p.pdf)

## Multiplication and complementary prime contents

For a nonzero integer \(\omega\) in a number field \(F\), define its
integer content as the largest positive \(m\) with \(\omega\in m\mathcal O_F\).

**Proposition (complementary prime contents).** Take two independent lines
of \(\mathrm{Cl}(K)[3]\), and any two trace-zero integers
\(\theta_1\in L_1,\theta_2\in L_2\) with \(Q_i(\theta_i)=n_i>1\) squarefree
and prime to \(6D\). In the compositum \(E=L_1L_2\), of degree nine, there
are cubic subfields \(L_+,L_-\) for the other two lines. The components
\[
 P_\pm=\tfrac13\operatorname{Tr}_{E/L_\pm}(\theta_1\theta_2)
\]
are integral with trace zero and satisfy
\[
 \theta_1\theta_2=P_++P_-,\qquad
 Q_\pm(P_\pm)=n_1n_2,\qquad
 \operatorname{cont}(P_+)\operatorname{cont}(P_-)=\gcd(n_1,n_2).
                                                                    \tag{2}
\]
The two contents are squarefree and coprime. If \(n_1=p_1p_2\) and
\(n_2=p_2p_3\) for distinct primes, precisely one component is divisible
in its ring of integers by \(p_2\). Dividing it by \(p_2\) gives a trace-zero
integer of quadratic value \(p_1p_3\). The other component is primitive and
has quadratic value \(p_1p_2^2p_3\). The short component corresponds to the
coherent closing edge after matching the orientations at \(p_2\).

Choose the two Kummer radicals \(\beta_i\) so that
\(\beta_i+\bar\beta_i=\theta_i\) and \(\beta_i^3=\alpha_i\).
Here \(\beta_i\bar\beta_i=n_i\). The two expansions of the product are
\[
 P_+=\beta_1\beta_2+\bar\beta_1\bar\beta_2,\qquad
 P_-=\beta_1\bar\beta_2+\bar\beta_1\beta_2.                         \tag{3}
\]
They are integral roots of the trace trinomials for
\(\alpha_1\alpha_2\) and \(\alpha_1\bar\alpha_2\), respectively. Their
Kummer classes lie on the other two lines and each quadratic value is
\(n_1n_2\).

To identify these roots inside \(E\), let \(S=S_1S_2\), with group
\(C_3^2\rtimes C_2\). Each chosen \(\theta_i\) is fixed by one reflection
in the corresponding \(S_3\) quotient. Above the nontrivial automorphism
of \(K^\vee\), any pair of reflections in the two quotients lifts uniquely,
because \(C_3^2\) maps bijectively to the product of the two order-three
groups. The common stabilizer is of order two, with fixed field \(E\).
Lift its involution to \(SM\), \(M=K(\zeta_3)\), so that it acts
nontrivially on \(K\) and fixes \(\zeta_3\); this is possible since
\(S\cap M=K^\vee\). It fixes each \(\theta_i\), so permutes the roots of
\(T^2-\theta_iT+n_i\). It cannot fix \(\beta_i\), because it sends
\(\alpha_i\) to the distinct conjugate. It therefore swaps both pairs in
(3) simultaneously, fixing the two sums. Thus \(L_\pm=\mathbb Q(P_\pm)\)
are subfields of \(E\). In computing the relative trace to \(L_+\), the
three embeddings of \(E\) over \(L_+\) fix \(P_+\) and send \(P_-\) to its
three conjugates. Indeed the two line quotients on \(C_3^2\) are independent,
so the kernel of the \(+\) quotient acts transitively on the roots for the
\(-\) quotient. As \(P_-\) has trace zero, this proves the trace formula
for \(+\), and the other formula follows by interchanging the lines.

For each prime \(p\mid\gcd(n_1,n_2)\), the two ideals \((\alpha_i)\) have
valuation three at a single orientation over \(p\). Exactly one of
\(\alpha_1\alpha_2,\alpha_1\bar\alpha_2\) uses opposite orientations and
is divisible by \(p^3\) in \(\mathcal O_K\). The corresponding product
of the radicals divided by \(p\) is integral, because its cube is
integral. Its sum with the conjugate is integral and lies in \(L_\pm\).
For the other product, the ideal is supported at a single orientation
over \(p\), so its rational trace is a unit at \(p\).
This trace is the cubic norm of the corresponding component, which
therefore cannot be divisible by \(p\). There are no other content primes
and no higher content powers: for any trace-zero integer \(w\),
\(\operatorname{Tr}(w^2)/2\) is integral, so if \(p^e\) divides a component,
\(p^{2e}\mid 3n_1n_2\). This proves (2) including the assertion at \(3\),
since \(3\nmid n_1n_2\). For the single common prime, taking the product
with opposite orientations and dividing \(\alpha_1\alpha_2\) or
\(\alpha_1\bar\alpha_2\) by \(p_2^3\) is exactly the coherent ideal
cancellation at \(p_2\).

## Counting edges with a prescribed divisor

Let \(f\ge1\) with \((f,6D)=1\). Write
\[
 T_f(K)=|\operatorname{Pic}(\mathbb Z+f\mathcal O_K)[3]|-|A_f[3]|
\]
with \(A_f\) as in the [residue criterion](ABC_COHERENT_NORM_RING_CLASS_RESIDUES.md#ideal-class-lifts-and-exact-residue-patterns).
Here \(|A_f[3]|\) is the size of the group of invertible fractional ideals of
\(\mathbb Z+f\mathcal O_K\) whose cubes equal the unit ideal. To see this,
write \(R=\mathbb Z+f\mathcal O_K\). The kernel of extension from the
invertible fractional \(R\)-ideals to the maximal-order ideal group has
local representatives \(a_\ell R_\ell\) at \(\ell^e\Vert f\), with
\(a_\ell\in\mathcal O_{K,\ell}^\times\), and the unit ideal at all other primes.
Two such local representatives define the same lattice exactly when
their ratio is in \(R_\ell^\times\), namely when their residues modulo
\(\ell^e\) differ by a rational unit. Each choice glues to a unique lattice,
by taking the inverse image of the chosen local sublattices in
\(\mathcal O_K/f\mathcal O_K\), and the lattice is invertible since it is
locally principal. So the kernel is \(A_f\). The maximal-order ideal
group is torsion-free, hence its kernel contains all the \(3\)-torsion ideals.
Thus \(T_f(K)\) is the class-group minus ideal-group difference in
[Bhargava–Varma, Theorem 8](https://arxiv.org/pdf/1401.5875).
Their mean there is ordered by the discriminants \(D f^2\) of orders.
The application below instead uses the first moment on fields ordered by
\(D\), so its error does not depend on a fixed-conductor asymptotic.

**Proposition (uniform prescribed divisor bound).** The number of
nonprincipal coherent single-edge records with primes in \([Z,2Z)\),
prime to \(6D\), and \(f\mid v/\epsilon_d\), is
\[
                  \ll \left(1+\frac{Z^2}{f\sqrt d}\right)T_f(K).
                                                                    \tag{4}
\]
The constant is absolute. Count the primes with order and orientation and
the two signs of the generator, allowing these bounded factors. Uniformly
for \(X\ge2,Z\ge2,f\ge1\), summing over squarefree \(d\in[X,2X)\) with
\(D<-4\) and \((f,6D)=1\) gives, for every \(\varepsilon>0\),
\[
                    \ll_\varepsilon f^\varepsilon
                           \left(X+\frac{Z^2\sqrt X}{f}\right).     \tag{5}
\]

In an integral basis of \(\mathcal O_{L_t}\), the determinant of
\(1,\theta,\theta^2\) is the signed index, a cubic form on \(\Lambda_t\).
At \(\ell\mid f\), this form is, up to a unit after unramified base change,
the Vandermonde product of the three pair differences among the coordinates
in the split cubic algebra. The trace-zero lattice reduces onto the
trace-zero plane since \(\ell\nmid3\). Distinct difference lines on that
plane meet only at zero, since all three coordinates equal and trace zero
imply all three vanish. A root direction modulo \(\ell\) is rational exactly
when Frobenius fixes the associated unordered pair of coordinates. There
are respectively three, one or zero fixed pairs for the split,
transposition and three-cycle types. Each root direction is simple and
lifts uniquely to a projective root at every \(\ell^e\) by Hensel's lemma.
Every \(\theta\) under consideration is primitive in \(\Lambda_t\), since
a rational divisor \(m\) of \(\theta\) would have \(m^2\mid3pq\).
Consequently its index is divisible by \(\ell^e\) exactly when its vector
lies on one of those root directions modulo \(\ell^e\).

The Chinese remainder theorem gives a union of rank-two sublattices of
index \(f\) in \(\Lambda_t\), one per collection of root directions. The
number of collections is
\[
   |A_f[3]|\,\mathbf1_{\theta_{3,f}(t)=0}.
\]
In fact \(\ell\) inert in \(K^\vee\) gives the transposition type and one
root. For \(\ell\) split in \(K^\vee\) there are three roots if the first
color vanishes and zero otherwise. Each split prime of \(K^\vee\)
contributes a factor three to \(|A_f[3]|\), and the colors at powers of
\(\ell\ne3\) are determined by the first color. By (1) the index has the
same valuations as \(v/\epsilon_d\) at all these primes.

In the Euclidean norm \(\|\theta\|^2=\operatorname{Tr}(\theta^2)\), the
covolume of \(\mathcal O_{L_t}\) is \(\sqrt{\operatorname{disc}(L_t)}\).
If the image of the integral trace is \(c\mathbb Z\), then \(c=1\) or \(3\)
and the perpendicular spacing for the trace-zero plane is \(c/\sqrt3\).
It follows that each sublattice has covolume
\(f\sqrt{3\operatorname{disc}(L_t)}/c\asymp f\sqrt d\).
The disk has radius \(O(Z)\). If its lattice points span a line, at most
two nonzero points can be primitive in \(\Lambda_t\). Otherwise, take a
shortest vector of the sublattice and consider its parallel rows.
The perpendicular spacing \(h\) satisfies \(h\ge\sqrt3/2\) times the
shortest length \(l\), by subtracting the nearest multiple of the shortest
vector from a representative of an adjacent row. Rank two in the disk
forces its radius \(R\ge h\) and \(R\ge l\). There are at most
\((1+2R/h)(1+2R/l)\ll R^2/(hl)\) points, with \(hl\) the covolume.
Every recording of a prime edge from one \(\theta\) costs at most an absolute
constant by (1) and unique ideal factorization.

There are \((|\ker\theta_{3,f}|-1)/2\) eligible nonzero lines, and
\(|\mathrm{Pic}(\mathbb Z+f\mathcal O_K)[3]|=|A_f[3]|\,|\ker\theta_{3,f}|\).
This proves (4). Since \(|A_f[3]|\le3^{\omega(f)}\) and
\(|\ker\theta_{3,f}|\le|\mathrm{Cl}(K)[3]|\), the
[Davenport–Heilbronn first moment](sources/bst-2013-davenport-heilbronn.pdf)
gives (5). Its discriminant cutoff applies with \(|D|\le 8X\).
The bound is for one prescribed divisor and a single edge. Summing over
any further set of divisors costs the sum of its right-hand sides, barring
additional constraints. To improve a prime-path moment, one needs an
estimate for joining edges from different cubic fields along a prime
while retaining the rank-one and field-window restrictions. The projection
in (2) describes exactly what such a join must do to its other two
cubic fields.

## A coherent example

For the [triangle of discriminant \(-4027\)](ABC_COHERENT_NORM_RING_CLASS_RESIDUES.md#exact-fixtures),
take the two inputs at edges \(12,23\). The input quadratic values are
\(13\cdot17\) and \(17\cdot19\). Their product of quadratic generators has
trace \(16040945=17^3\cdot3265\) and coefficient \(545343=17^3\cdot111\).
The product with the conjugate has trace \(31714029\) and signed
coefficient \(-333971\). Its quadratic value is \(71383=13\cdot17^2\cdot19\),
and \(17\nmid31714029\). These data verify that exactly one component has
content \(17\), and its quotient gives the observed edge \(13\). They
follow by multiplying the two displayed quadratic pairs. The exact norm
and cancellation tests are integer identities.

For an index-three trace lattice, let \(\tau^3-4\tau-1=0\).
The polynomial is irreducible by the rational-root test and has prime
discriminant \(229\), so \(1,\tau,\tau^2\) is an integral basis. The traces
of \(\tau,\tau^2,\tau^3,\tau^4\) are \(0,8,3,32\). In the trace-zero
basis \(\tau,3\tau^2-8\), the form is
\(4x^2/3+3xy+16y^2\). It is integral exactly for \(3\mid x\),
where its integral form is \(12x_0^2+9x_0y+16y^2\) with discriminant
\(-687\). The element \(3\tau+3\tau^2-8\) has quadratic value \(37\) and
norm \(322\). The quadratic generator \(161+6\sqrt{-687}\) has norm
\(37^3\) and trace prime to \(37\), and its reflected trinomial is the
minimal polynomial of this trace-zero element. Thus the example realizes
the index-three case for \(d=687\).

An internal PARI/GP check of the \(-4027\) trinomials gave cubic field
discriminant \(108729\) in all three cases and indices \(28,111,139\),
respectively. At the index prime \(37\) for edge \(13\) the field
factor degrees are \(1,2\), consistent with the one permitted
projective root of the index form. The integer identities above give
the exact fixture independently of a field search.
