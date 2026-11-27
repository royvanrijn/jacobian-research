# Wild Dickson fields from a genus-three Pell norm

The genus-three octic with three marked elliptic quotients has
a norm-one function supported at two split rational fibers.
It is a virtual cube exactly on the specializations ramified
at three, where the local value is always a noncube.
After a power base change this gives degree-\(3^n\) Dickson
fields totally ramified at three, with explicit discriminants
and ray conductors. Actual prime ratios can prescribe the
unramified Frobenius orders at both moving primes.

## A quotient Pell norm and the local obstruction

Use the octic \(F\) of
[HG1](GENUS_THREE_THREE_MARKED_ELLIPTIC_QUOTIENTS.md).
Set
\[
 z=t-2/t,\quad \rho=R/t^2,\quad
 G(z)=z^4-6z^3+7z^2-18z+97,\qquad \rho^2=G(z).
\]
This is the quotient by \((t,R)\mapsto(-2/t,4R/t^4)\),
the simultaneous sign change on the conic of HG1.
Put
\[
 \begin{split}
 A(z)&=z^6-5z^5+2z^4-15z^3+76z^2+65z+101,\\
 B(z)&=z^4-2z^3-3z^2-14z-5,\\
 \gamma&=\frac{(A+B\rho)^2}{7776}.
 \end{split}
\]
The polynomial Pell identity and trace are
\[
 A^2-G B^2=7776,\qquad
 N(\gamma)=1,\qquad
 T=\operatorname{Tr}(\gamma)=\frac{A(z)^2}{1944}-2 .
 \tag{1}
\]
These are polynomial multiplications. The leading terms of
\(A,B\) at either infinity of the quartic give a zero of
order six on one branch and a pole of order six on the
other for \(A+B\rho\). Pullback and squaring show that
\(\gamma\) has signed orders \(\pm12\) at each of the
fibers \(t=0,\infty\). The polynomial Pell interpretation
is classical. [Barroero–Capuano–Zannier, Proposition 1.2](https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/betti-maps-pell-equations-in-polynomials-and-almostbelyi-maps/999A05A799D0BC64A35C64D4AA955166).

**Proposition WD1.** For any nonzero rational \(t\) such
that \(K=\mathbb Q(\sqrt{F(t)})\) is quadratic, \(\gamma\)
has valuations divisible by three away from \(3\). It is
virtual at \(3\) precisely when \(t\) is a \(3\)-adic unit
and
\[
                        t\equiv1,4,5,7\pmod9 .
                        \tag{2}
\]
These are exactly the specializations in which \(3\)
ramifies. In those cases \(K_3=\mathbb Q_3(\sqrt3)\)
and \(\gamma\) is a unit with
\(v_{K_3}(\gamma+1)=1\) in the uniformizer normalization.
In particular it is a local noncube.

Write \(k=v_p(t)\). The two values of \(\gamma\) are the
roots of \(X^2-TX+1\). When \(T\) is integral they are units.
When \(v_p(T)<0\) the two distinct Newton slopes give split
completions and valuations \(v_p(T),-v_p(T)\).
For \(p\ne2,3\), the trace is integral if \(k=0\).
If \(k\ne0\), \(v_p(z)=-|k|\), and the leading term of
\(A\) gives \(v_p(T)=-12|k|\).
At \(2\), \(v_2(z)=\min(k,1-k)\le0\). If this is zero,
\(z\) is odd and \(A(z)\) is odd. In both cases
\[
                       v_2(T)=12\min(k,1-k)-3 .
                       \tag{3}
\]
At \(3\) and \(k\ne0\), the leading term instead gives
\(v_3(T)=-12|k|-5\). This is nonvirtual.

Suppose \(t\) is a \(3\)-adic unit. When \(t\equiv1\pmod3\),
the identity \(z+1=(t-1)(t+2)/t\) gives \(z\equiv-1\pmod9\).
Here \(A(-1)=135\), \(A'(-1)=-171\), so Taylor expansion
gives \(v_3(A(z))=3\). Also \(F(t)\equiv3\pmod9\).
When \(t\equiv-1\pmod3\), put \(z=1+9r\), with
\(r=(t-2)(t+1)/(9t)\). Then
\[
 G(1+9r)=81(1-2r-5r^2-18r^3+81r^4).
\]
For \(t\equiv2,8\pmod9\), \(r\equiv0\pmod3\), so
\(K_3\) splits and \(v_3(A(z))=2\), by \(A(1)=225\).
Thus \(v_3(T)=-1\). For \(t\equiv5\pmod9\),
\(r\equiv4\pmod9\): \(z(t)\bmod81\) is constant
on this residue class, since
\[
 z(t)-z(5)=(t-5)\left(1+\frac2{5t}\right),\qquad
 v_3\left(1+\frac2{5t}\right)\ge2 .
\]
Evaluation at \(5\) gives \(v_3(A(z))=3\), and the
bracket in the displayed formula for \(G\) is \(3\bmod9\).
This gives squareclass \(3\) for \(G\), hence also for
\(F(t)\).

On both ramified cases \(v_3(T+2)=1\).
The identity \((\gamma+1)^2/\gamma=T+2\) and the unit
condition imply \(v_{K_3}(\gamma+1)=1\). The cube of a
unit reducing to \(-1\) is \(-1\) modulo the cube of
the maximal ideal of \(K_3\), because \(v_{K_3}(3)=2\).
This excludes a local cube.

## Power slices and the degree-\(m\) field

**Theorem WD2.** Fix \(m=3^n\), \(n\ge1\), and choose an
integer \(a\) with \(4a\equiv1\pmod m\). Use
\[
 t=2^a u^m,\qquad v_2(t)<0,\qquad t\equiv1\pmod3 .
 \tag{4}
\]
These conditions give a nonempty open rational parameter
slice. Exclude the parameters for which \(K\) is not
quadratic. Write \(D=\operatorname{disc}K>0\) and
\(D_0=D/3\). Then \(\gamma\) is \(m\)-virtual and locally
noncubic above \(3\).
Let \(\mathcal D_m\) be the Dickson polynomial defined by
\(\mathcal D_m(Y+Y^{-1},1)=Y^m+Y^{-m}\), and let \(L_m\)
be the field of a root of
\[
                         \mathcal D_m(X,1)-T .
                         \tag{5}
\]
Its normal closure \(N_m/\mathbb Q\) has group
\(C_m\rtimes(\mathbb Z/m)^\times\) in the usual affine
action, and \(L_m\) has signature \((1,(m-1)/2)\) and
absolute discriminant
\[
                   |\operatorname{disc} L_m|
                        =3^{nm}D^{(m-1)/2}.
                   \tag{6}
\]
There are \(\gg_m X^{1/(4m)}/(\log X)^2\) different \(D\le X\),
and \(\gg_m Y^{1/(2m(m-1))}/(\log Y)^2\) different \(L_m\)
of absolute discriminant at most \(Y\). After excluding
finitely many parameters, the three twists from HG1 have
positive rank and common conductor \(10D^2\).

At \(2\), the negative valuation in (3) is
\(12v_2(t)-3\), divisible by \(m\). At \(p\ne2,3\)
the corresponding valuation is zero or \(-12|v_p(t)|\),
again divisible by \(m\) on (4). WD1 handles \(3\).
The specified neighborhoods exist by weak approximation,
since \(u\) can have arbitrarily negative valuation at
two and the unit \(u\bmod3\) can be prescribed.

The real field \(K\) is disjoint from \(\mathbb Q(\zeta_m)\),
whose unique quadratic subfield is imaginary. The cubic
root of the noncube \(\gamma\) over \(K\) gives a nonnormal
degree-three extension of \(K\), since a normal cubic
extension containing one cube root would contain \(\zeta_3\).
It cannot be contained in the abelian extension
\(K(\zeta_m)/K\). Thus \(\gamma\) stays noncubic there.
Kummer theory gives a cyclic extension of degree \(m\)
by \(\theta=\gamma^{1/m}\) over \(K(\zeta_m)\).
The two involutions on \(K\) and on the cyclotomic factor
act on \(\theta,\zeta_m\) separately by inverting \(\theta\)
and inverting \(\zeta_m\), respectively. Their product
fixes all \(m\) distinct roots
\(\zeta_m^j\theta+\zeta_m^{-j}\theta^{-1}\).
On these roots translations and the cyclotomic units act
as the full affine group. Thus \(N_m\) is the fixed field
of that product involution and the root has degree \(m\).
Since \(\gamma>0\) on the real branches and \(\gamma\ne1\),
the root with \(j=0\) using positive real \(\theta\) is
the only real one.

Away from \(3\), the radical over \(K(\zeta_m)\) is
unramified by \(m\)-virtuality and the cyclotomic field
is unramified. At any \(p\ne3\) dividing \(D\) (necessarily
odd by HG3), inertia has order two and acts as an affine
involution with one fixed point, contributing \((m-1)/2\)
to the discriminant exponent. All other primes away from
three are unramified.

At \(3\), \(\mathcal D_m(X,1)\equiv X^m\pmod3\) and
\(\mathcal D_m(1,1)=-2\). The latter follows by taking
\(Y\) a primitive sixth root in its defining identity,
as \(m\equiv3\pmod6\). Thus (5) translated by \(X=1+Z\)
is Eisenstein at \(3\) by \(v_3(T+2)=1\).
It gives a totally ramified degree-\(m\) extension with
uniformizer \(Z\) and ring of integers \(\mathbb Z_3[Z]\).
For the derivative, differentiation of the defining
identity yields
\[
     (\mathcal D_m'(X,1))^2
          =m^2\frac{T^2-4}{X^2-4}
\]
at a root. In the root field \(v(X^2-4)=1\), since
\(X^2-4=(Z-1)(Z+3)\), and \(v(T^2-4)=m\).
The different exponent is \(nm+(m-1)/2\).
Multiplying by the exponents at primes dividing \(D_0\)
gives (6).

For the count, the fixed form for (4) has degree \(8m\),
is squarefree, and clearing its fixed denominator by
a square preserves the quadratic fields. Apply the
degree-two case of
[Kulkarni–Levin, Theorem 1.2](sources/kulkarni-levin-hilbert-irreducibility-2111.15582.pdf)
on fixed neighborhoods inside (4). Parameter height \(H\)
gives \(D\ll_m H^{8m}\), hence the stated count.
Formula (6) distinguishes the fields for distinct \(D\).
The degree-eight maps of HG1 are nonconstant and send
the hyperelliptic involution to negation. For a fixed
elliptic curve only finitely many torsion points have
degree at most two, by Northcott applied to the canonical
height. Each has finitely many preimages. Thus the three
twists have positive rank outside finitely many parameters.
The local squareclass at three is \(3\), and \(2,5\)
split by HG3. The conductor calculation of HG2 therefore
applies without requiring its independence congruences.

## Exact ray conductors through the tower

**Theorem WD3.** In WD2 put
\(M=\mathbb Q(\sqrt{-D_0})\) and
\(B_k=M\mathbb Q(\zeta_{3^k})^+\) for \(1\le k\le n\).
Use the same \(K,t\) for all levels. The fixed field of
the translations in \(N_{3^k}\) is \(B_k\).
The extension \(N_{3^n}/B_n\) is cyclic of degree \(3^n\),
unramified away from the unique prime \(\mathfrak q\)
of \(B_n\) above \(3\), and totally ramified there. A
character of exact order \(3^k\) has conductor
\[
                     \mathfrak q^{(3^k+1)/2}.
                     \tag{7}
\]
In particular it is a cyclic quotient of the ray class
group of modulus \(\mathfrak q^{(3^n+1)/2}\) of exact
conductor, and
\[
   \operatorname{disc}(N_{3^n}/B_n)
                  =\mathfrak q^{(3m^2+4m-7)/8},\qquad m=3^n .
                  \tag{8}
\]

The translation-fixed field follows by taking the
simultaneous involution on \(K(\zeta_{3^k})\). The two
generators remaining are the product quadratic \(M\)
and the real cyclotomic layer. The quadratic field \(M\)
has fundamental discriminant \(-D_0\), and \(M_3\)
is the unramified quadratic extension of \(\mathbb Q_3\).
Thus \(B_k\) has one prime over \(3\); the completion
has residue degree two and ramification index \(3^{k-1}\).
The Kummer argument of WD2 also works locally over
\(K_3=\mathbb Q_3(\sqrt3)\), which does not contain
\(\zeta_3\) and is disjoint from \(\mathbb Q_3(\zeta_{3^k})\).
The radical stays noncubic over their abelian compositum
by the same nonnormality argument. The local normal
closure therefore realizes the full affine group at
every level. Away from three its inertia is trivial or
of order two, so the extension over the translation-fixed
field is unramified.

For local fields write \(a(\chi)\) for the conductor
exponent. We use the induction formula
\[
 a(\operatorname{Ind}_{E/F}\chi)
 =v_F(\operatorname{disc}(E/F))+\! f(E/F)a(\chi)
 \tag{9}
\]
for a character \(\chi\) of \(E\); for the trivial character
it is the conductor of the permutation representation.
The general formula, including ramified \(E/F\), is
[Schwein, equation (4.4), p. 3723](sources/schwein-formal-degree-2024.pdf).

In the permutation representation on the \(3^k\) roots,
the translation subgroup has one copy of every character.
The units permute them transitively at each fixed exact
order; those of order less than \(3^k\) give the permutation
representation on the roots at level \(k-1\).
The primitive summand is induced from any primitive
character of the translation subgroup. Thus its conductor
exponent over \(\mathbb Q_3\) is the difference between the
different exponents of the two root fields, namely
\[
 \left(k\,3^k+\frac{3^k-1}{2}\right)
 -\left((k-1)3^{k-1}+\frac{3^{k-1}-1}{2}\right)
                  =(2k+2)3^{k-1}.
\]
At \(k=1\) the lower field is \(\mathbb Q_3\).
The full cyclotomic discriminant exponent is
\((2k-1)3^{k-1}\). Its relative quadratic extension
over the real layer is tamely ramified with different
exponent one. The unramified quadratic compositum
therefore gives discriminant exponent
\((2k-1)3^{k-1}-1\) for \(B_k/\mathbb Q\) locally at
three. Formula (9) with residue degree two gives
\(a(\chi)=(3^k+1)/2\) for each primitive order-\(3^k\)
character over \(B_k\).

At level \(n\) the translation subgroup maps onto the\ntranslation subgroup of each lower normal closure,\nso \(N_{3^k}\cap B_n=B_k\).\nIt remains to pass a lower-order character to \(B_n\).
Write \(a_k=(3^k+1)/2\) and let \(\chi\) be a primitive
character at level \(k<n\).
The first cyclic cubic step \(B_{k+1}/B_k\) has
conductor exponent \(a_k\). Indeed
\(\pi=\zeta+\zeta^{-1}-2\) is a uniformizer of the local
real layer for \(\zeta=\zeta_{3^{k+1}}\). For the
automorphism \(\sigma:\zeta\mapsto\zeta^{1+3^k}\),
\[
 \sigma(\pi)-\pi=(\zeta^{3^k}-1)(\zeta-\zeta^{-1-3^k}).
\]
Its valuation in the full cyclotomic field is \(3^k+1\).
In the real layer it is \(a_k\), unchanged by the
unramified quadratic factor. The single lower and upper
break of the cubic step is \(a_k-1\).
Any nontrivial character \(\lambda\) of \(B_n/B_k\)
has conductor at least \(a_k\), since a power of it is
a nontrivial character of the first cubic step.

The involution of the unramified quadratic factor fixes
the cyclotomic characters \(\lambda\) and inverts the
character \(\chi\), by the affine action. In particular
if \(a(\lambda)=a_k\), the two restrictions to the
last principal-unit quotient \(U^{a_k-1}/U^{a_k}\)
cannot cancel, as they have opposite eigenvalues over
\(\mathbb F_3\). Hence
\(a(\chi\lambda)=a(\lambda)\) for every nontrivial
\(\lambda\). By local class field theory, restriction of
\(\chi\) to \(B_n\) is the character pulled back by the
local norm. Inducing its restriction to \(B_k\) gives
the sum of the twists \(\chi\lambda\). The discriminant
term in (9) equals \(\sum_\lambda a(\lambda)\), and the
residue degree of \(B_n/B_k\) is one. Subtracting gives
conductor \(a_k\) over \(B_n\) as well. The degree-three
quotient is ramified, so the cyclic extension has no
nontrivial unramified quotient at three. This proves total
ramification. Summing the conductor exponents over all
nontrivial characters gives (8).

## Actual moving primes

**Corollary WD4.** For each fixed \(m=3^n\) and each pair of
divisors \(e_P,e_Q\mid m\), there is a constant \(s_m\) such
that all conclusions of WD2–WD3 apply on
\(t=s_m(P/Q)^m\) for genuine distinct primes
\(P,Q\equiv1\pmod m\) in two positive-density Chebotarev
sets. The primes split completely in \(B_n\). Their
Frobenius orders in \(N_m/B_n\) are \(e_P,e_Q\) at every
place above the respective rational prime. The resulting
different real quadratic discriminants up to \(X\) number
\(\gg_m X^{1/(4m)}/(\log X)^{C_m}\), with constants
allowed to depend on the fixed orders.

Choose \(s_m=\pm2^{a-mN}\) with \(N\) large enough
to make \(a-mN<0\) and with sign making
\(s_m\equiv1\pmod3\). The parameter satisfies (4),
absorbing the sign and \(2^{-mN}\) in \(u^m\).
At both moving primes the octic splits, by its boundary
coefficients \(F(0)=16\) and leading coefficient \(1\).
They also split in the cyclotomic field, hence in \(B_n\).
The pole branch of \(\gamma\) in (1) has leading
coefficient \(1/1944\) at infinity and \(2^{12}/1944\)
at zero. Remove the valuation, \(-12m\) at each moving
prime. The resulting residue class modulo \(m\)-th powers
has representative
\[
             A_P=\frac{2^{12}s_m^{-12}}{1944},\qquad
             A_Q=\frac{s_m^{12}}{1944},
\]
respectively. The other branch gives its inverse.
Each constant has \(3\)-adic valuation \(-5\), so is
noncubic over \(\mathbb Q\). A nonnormal rational cubic
radical field cannot be contained in the cyclotomic
field. Consequently each radical field
\(\mathbb Q(\zeta_m,\sqrt[m]{A_P})\) and
\(\mathbb Q(\zeta_m,\sqrt[m]{A_Q})\) has degree \(m\)
over the cyclotomic field. Taking Frobenius classes
of the prescribed orders in the cyclic translation
subgroups gives positive-density sets by Chebotarev.
These are separate choices at different primes.

For the count take \(P,Q\asymp H\) in disjoint fixed
intervals. Chebotarev gives \(\gg_m H^2/(\log H)^2\)
pairs. The prime-factor truncation of
[PP3](ABC_COHERENT_PRIME_PROGRESSION_NORMS.md#quantitative-count-of-genuine-fields)
applies to the fixed squarefree binary form of degree
\(8m\), giving \(\omega(D)=O_m(\log\log H)\) on all but
\(O_m(H^2/(\log H)^4)\) pairs. The same fixed-form argument
as in PP2 gives rank \(O_m(1+\omega(D))\) for the Jacobian
of the quadratic twist of this genus-\((4m-1)\) curve.
Over the fixed splitting field of the form the two-torsion
is rational, so Kummer theory bounds the rank using the
primes of bad reduction and those dividing \(D\).
The Abel–Jacobi image has fixed degree and contains
no positive-dimensional abelian translate.
[Gao–Ge–Kühne, Theorem 1.1](sources/gao-ge-kuehne-uniform-mordell-lang-2105.15085v4.pdf)
then bounds the parameter multiplicity by a power of
\(\log H\). Taking \(H\asymp_m X^{1/(8m)}\) proves
the claim.

## Verification and comparison

The [exact checker](verify_genus_three_pell_wild_dickson_fields.py)
replays the quotient, Pell, local residue, Dickson
discriminant and character conductor arithmetic.
This is the wild counterpart to the local-power
construction of
[YR4–YR5](MARKED_OCTIC_CYCLOTOMIC_THREE_POWER_REFLECTION.md).
The ray quotients are ramified at three and do not give
the asserted orders in the ordinary imaginary class
group. At each fixed \(n\) the family has a common choice
of local conditions. The construction does not assert
one fixed quadratic field with unbounded \(n\).

## Every rational Pell coordinate has the same dyadic obstruction

**Theorem WD5.** On the HG1 genus-three curve,
the conjugate difference at each of its four
rational split fibers has exact order twelve.
A function with divisor twelve times that
difference has constant norm of squareclass
\(6\). Every nonconstant norm-one function
supported on one rational split fiber has
nonzero dyadic valuation at both points of
every quadratic specialization.
Thus no rational change of the degree-two
coordinate on this curve supplies an
integral Pell unit of norm \(1\) or \(-1\)
on a quadratic specialization.

The transport between fibers is a general
parity constraint. On a hyperelliptic curve
write \(\Delta_P=P-\iota P\).
If \(P,Q\) are rational and
\(N(P-Q)\) is principal, take
\(f\) with \(\operatorname{div}(f)=N(P-Q)\).
Then \(f/\iota(f)\) has norm one and divisor
\(N(\Delta_P-\Delta_Q)\). Consequently
the functions with divisors \(N\Delta_P\)
and \(N\Delta_Q\), if they exist, can be
chosen with the very same constant norm.
For norm-one functions with divisors
\(2N\Delta_P\) and \(2N\Delta_Q\), this
transport multiplies by a square.
It therefore preserves the parity of
any split local valuation at a specialization
where the functions are defined.

For the present curve, write
\(G_0=g_0+\ell_0 R\), where
\[
 g_0=t^6-3t^5-3t^4-12t^3+6t^2-12t-8,
                   \qquad \ell_0=t^2+2 .
\]
The identity \(g_0^2-F\ell_0^2=-200t^6\)
and the leading terms give
\(\operatorname{div}(G_0)=6((0,4)-\infty_+)\).
Thus \(\gamma_0=G_0/\bar G_0\) has
divisor \(6\Delta_{(0,4)}-6\Delta_{\infty_+}\).
[The marked norm identities](GENUS_THREE_MARKED_TORSION_DEFORMATIONS.md)
give this as their \(a=0\) row.
The function \(A+B\rho\) in (1) of
this note has divisor
\(-6\Delta_{(0,4)}-6\Delta_{\infty_+}\).
Indeed \(z=t-2/t\) has a simple pole
both at zero and at infinity, and
\(\rho/z^2\) tends to \(+1\) on the
branches \(R\to4\) and \(R/t^4\to1\),
respectively. The leading terms there
give poles of order six. Its constant
norm \(7776\) gives the opposite zeros,
and there are no other poles.
It follows that
\[
 V_\infty=\frac{7776}{(A+B\rho)\gamma_0},
 \qquad
 U_\infty=\frac{V_\infty^2}{7776}
                     =\gamma^{-1}\gamma_0^{-2}
\]
have divisors \(12\Delta_{\infty_+}\)
and \(24\Delta_{\infty_+}\), respectively.
The first has norm \(7776\) of squareclass
\(6\). The second has norm one.

By WD1, every quadratic specialization
of the curve splits at two, and the
two values of \(\gamma\) have opposite
odd dyadic valuations. The square
factor \(\gamma_0^{-2}\) implies both
valuations of \(U_\infty\) are odd.
All rational fibers consist of the
pairs at \(t=-1,0,2,\infty\), by
[HG6](GENUS_THREE_THREE_MARKED_ELLIPTIC_QUOTIENTS.md).
HG6 also gives exponent twelve
for the rational Jacobian.
The transport lemma applies with
\(N=12\). At each of the rational
split fibers it gives a function
of divisor \(12\Delta_P\) with
constant norm \(7776\), and a
norm-one function of divisor
\(24\Delta_P\) with odd valuations
at two on every quadratic
specialization. The functions
are defined and nonzero at these
specializations because their
supports are rational split
fibers.

The order of \(\Delta_P\) is
divisible by four. Otherwise a
function of divisor \(12\Delta_P\)
would be an even power of one
with divisor the exact order
times \(\Delta_P\), up to
multiplication by a rational
constant. Its norm would then
have squareclass one.
The order is divisible by
three as well: at each of the
four positive branches, the
first elliptic map in HG1
gives \((0,5/2)\) or \((0,-5/2)\),
and the conjugate branch maps
to its negative. The
difference has order three.
Hence the exact order is twelve.
Both the rational primitive
divisor function and its
norm squareclass are thereby
determined up to constants, with
norm squareclass \(6\). The
least positive order of an
infinity divisor of a norm-one
function is twenty-four.
All nonconstant norm-one
functions supported at that
fiber are powers of the
indicated norm-one function,
up to rational constants of
norm one. None can have zero
dyadic valuations there after
quadratic specialization. A
norm-minus-one unit would
have a norm-one square and
is excluded too.

For a degree-two coordinate
of the same hyperelliptic
curve the fiber at infinity
is a fiber of the original
map. A ramified or nonsplit
rational fiber cannot support
a nonzero rational degree-zero
divisor, and all split
rational fibers have just
been treated. The argument
therefore excludes changing
the rational point sent to
infinity as a way to obtain
an integral Pell unit on
this fixed curve. The
obstruction is the same
dyadic squareclass throughout
the rational torsion packet.
