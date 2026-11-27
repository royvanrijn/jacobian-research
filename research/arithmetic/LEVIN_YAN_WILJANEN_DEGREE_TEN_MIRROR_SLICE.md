# Four and five virtual directions on the Levin–Yan–Wiljanen curves

The degree-ten genus-four curve in Levin–Yan–Wiljanen admits a fixed
nonempty rational parameter slice with four independent virtual cubic
classes, all cubes above three. Consequently there are
\(\gg X^{1/5}/(\log X)^2\) distinct real quadratic fields of discriminant
at most \(X\) with real class \(3\)-rank at least three and imaginary
mirror class \(3\)-rank at least four. The mirror extensions split at
three. This is a different family from
[MO3](QUADRATIC_MIRROR_CLASS_RANK_FOUR_FROM_MARKED_OCTIC.md).
Its three specified elliptic twists and their conductor laws are not
part of this result. The real unit-line question for MO3 remains open.

A cubic base change also makes the five individual norm functions
virtual and locally cubic at three. The two constructions give:

| Parameter equation | Real class 3-rank | Imaginary mirror class 3-rank | Distinct-field lower bound |
|---|---:|---:|---|
| \(y^2=f(t)\), LY2–LY3 | at least 3 | at least 4 | \(\gg X^{1/5}/(\log X)^2\) |
| \(y^2=f(-2a u^3)\), LY4–LY5 | at least 4 | at least 5 | \(\gg X^{1/15}/(\log X)^2\) |

Here \(X\) bounds the real fundamental discriminant, and the mirror
has absolute discriminant at most \(3X\). LY6 gives infinite unramified
pro-\(3\) mirror extensions with the full three-primary class group as
abelianization: two specified rational primes split completely for the
first family, and four for the second.

The mechanism is cancellation at the two rational infinities. Each
individual function \(y-y_i\) has pole orders not divisible by three;
the four ratios \((y-y_i)/(y-y_5)\) cancel those poles and tend to one
at one infinity. Approaching that point at the finitely many bad primes
controls all four ratios together without a cubic base change.

## Exact model and literature boundary

Use precisely \(f=f_4\) and \(y_i=y_{4,i}\), \(1\le i\le5\), in
[Levin–Yan–Wiljanen, Theorem 9](sources/levin-yan-wiljanen-large-three-rank.pdf),
arXiv:1910.12276v1, printed pp. 12–14. Their coefficients, in ascending
order, are retained as `F_COEFFS` and `Y_COEFFS` in the
[exact checker](verify_lyw_degree_ten_mirror_slice.py). These arrays
are the exact inputs for every statement below. Put
\[
 C:y^2=f(t),\quad a=11285328049647162,
 \quad \alpha_i=\frac{y-y_i(t)}{y-y_5(t)}\quad(1\le i\le4).
 \tag{1}
\]
The cited theorem already establishes genus four and four independent
rational Jacobian \(3\)-torsion classes, with divided divisors given by
the ratios in (1). Its proof checks independence in Magma. The argument
here uses separate residue tests to enforce independence in every
specialization on our slice; it does not infer specialization
independence from geometric independence alone.

The authors explicitly leave the enumerative application of their
degree-ten example unpursued (printed p. 13, after Corollary 10).
Their Theorem 3 requires \(6\mid\deg f\), and their Theorem 7 requires
\(\deg f\not\equiv2,4\pmod6\). Neither applies directly to this
degree-ten polynomial. The local argument below treats the ratios
instead. Theorem 5 supplies the general specialization context but
is not needed to replace the explicit residue tests or to count fields.
This comparison identifies the relevant argument, not historical
priority among later papers.
Theorem 9 also already supplies genus fourteen and five rational torsion
directions after \(t=u^3\). LY4 uses \(t=-2a u^3\) to control the
common leading term at every bad prime. The fifth geometric direction
is known; the additional work is the uniform local slice, specialized
independence and paired field count.

**Proposition LY1.** The polynomial \(f\) is squarefree of degree ten
and positive on \(\mathbb R\). All five \(y_i\) have degree five and
leading coefficient \(a\), and the leading coefficient of \(f\) is
\(a^2\). There are integer cubics \(x_i\) satisfying
\[
 y_i^2-f=x_i^3,\qquad
 N(\alpha_i)=(x_i/x_5)^3.
 \tag{2}
\]
The exact cubics, writing \(x_i=c_i h_i\), are
\[
\begin{array}{c|r|l}
i&c_i&h_i(t)\\ \hline
1&26&3327579857t^3-236377090446t^2+4162221598416t-468970154183\\
2&2&21571897004t^3-166911461493t^2-6525260234283t-6096612004379\\
3&-2&114744133t^3+712614460299t^2-4891775410473t+48838428515024\\
4&-22&6769903847t^3-125714528331t^2+962384304273t-40486600394636\\
5&22&19735990876t^3+789848985597t^2+9374196436563t+145782790908731.
\end{array}
\]

Polynomial multiplication over \(\mathbb Q\) proves (2). Squarefreeness
is the exact gcd test \(\gcd(f,f')=1\). For positivity the checker
computes the rational Sturm sequence, of degrees \(10,9,\ldots,0\).
Its signs at negative and positive infinity are respectively
\[
 (+,-,-,+,+,-,-,+,+,-,-),\qquad
 (+,+,-,-,+,+,-,-,+,+,-).
\]
Both have five variations, so \(f\) has no real root. Its positive
constant term gives positivity everywhere. This is a finite exact
Sturm certificate for the specified polynomial, not a sampled sign test.

## A uniform local slice

Let \(\Delta=\operatorname{disc}(f)\ne0\), and set
\[
 S=\{p:p\mid6a\Delta\},\quad
 k_p=2v_p(2a)+3,\quad m=\prod_{p\in S}p^{k_p}.
 \tag{3}
\]
This defines a specific finite integer modulus; factorization is not
needed to prove its existence. One may replace \(m\) throughout by
\((6a|\Delta|)^6\), whose valuations are at least the stated \(k_p\).
For the four additional primes use this table:
\[
\begin{array}{r|r|r|r|c|c}
q&t_q&r_q&\zeta_q&(\alpha_1,\alpha_2,\alpha_3,\alpha_4)\bmod q
 &\text{character exponents}\\ \hline
73&4&33&8&(33,2,38,2)&(1,2,1,2)\\
79&3&2&23&(26,21,8,1)&(2,0,0,0)\\
97&2&2&35&(27,65,31,92)&(0,2,1,1)\\
103&1&7&46&(102,95,58,19)&(0,0,2,1).
\end{array}
\tag{4}
\]
Here \(r_q^2=f(t_q)\not\equiv0\pmod q\), \(\zeta_q\) has order
three, and exponent \(e\) means that the corresponding ratio raised
to \((q-1)/3\) equals \(\zeta_q^e\). The checker verifies every entry
and that these primes lie outside \(S\).

**Theorem LY2.** Impose
\[
 v_p(t)\le-k_p\ (p\in S),\qquad
 t\equiv t_q\pmod q\ (q=73,79,97,103).
 \tag{5}
\]
These are simultaneously attainable nonempty open conditions on
\(t\in\mathbb Q\). For every such parameter for which
\(K=\mathbb Q(\sqrt{f(t)})\) is quadratic, all four classes (1) are
virtual units, independent in \(K^\times/K^{\times3}\), and cubes
in both completions above three. The field is real, two and three
split in it, and, writing its fundamental discriminant as \(D\),
\[
 \operatorname{rk}_3\operatorname{Cl}(K)\ge3,\qquad
 \operatorname{rk}_3\operatorname{Cl}(\mathbb Q(\sqrt{-3D}))\ge4.
 \tag{6}
\]
The imaginary mirror has an unramified \((\mathbb Z/3)^4\) extension
split above three, and hence at least forty complex cubic fields of
fundamental discriminant \(-3D\).

First control the primes in \(S\). Write \(z=1/t\),
\(f^*(z)=z^{10}f(1/z)\), and \(y_i^*(z)=z^5y_i(1/z)\).
All coefficients are integers, and
\[
 f^*(0)=a^2,\quad y_i^*(0)=a,\quad
 f^*(z)-a^2\in z\mathbb Z_p,
 \quad y_i^*(z)-a\in z\mathbb Z_p.
\]
Put \(e=v_p(a)\), \(b=v_p(2)\), and \(h=v_p(z)\ge2(e+b)+3\).
Then \(f^*/a^2\in1+p^{h-2e}\mathbb Z_p\). For odd \(p\), its
square root \(u\) close to one satisfies \(v_p(u-1)\ge h-2e\).
At two, \(h-2e\ge5\), and the square root close to one satisfies
\(v_2(u-1)\ge h-2e-1\). Thus the two branches are rational over
\(\mathbb Q_p\). On the branch \(y^*=-au\),
\[
 v_p(y^*-y_5^*)=e+b,
 \quad
 \alpha_i-1=\frac{y_5^*-y_i^*}{y^*-y_5^*}
       \in p^{h-e-b}\mathbb Z_p\subset p^3\mathbb Z_p.
 \tag{7}
\]
The error from \(-2a\) in the denominator has strictly higher
valuation: at odd \(p\) it is at least \(h-e\), and at two at
least \(h-e-1\). This proves the denominator equality used in (7).
For \(p\ne3\), \(1+p\mathbb Z_p\) consists of cubes (including
\(p=2\)); for \(p=3\), \(1+9\mathbb Z_3\) consists of cubes.
Equation (7) therefore makes all four evaluations cubes on this
branch. Equation (2) makes their evaluations cubes on the other
branch as well. The \(x_i(t)\) are nonzero on a quadratic fiber,
since a zero would make \(f(t)=y_i(t)^2\). This also ensures that
the ratios and their norms are defined and nonzero.

Next let \(p\notin S\). If \(t\in\mathbb Z_p\), the \(y_i(t)\)
and \(f(t)\) are integral. If \(f(t)\) is a unit, in each completion
above \(p\) at least one of \(y-y_i,y+y_i\) is a unit, because
\(2y\) is a unit. Their product is \(-x_i^3\); hence each has
valuation divisible by three. If \(f(t)\) is not a unit, every
\(y_i(t)\) is a unit. Otherwise (2) would give \(x_i(t)=0\bmod p\),
and differentiation of \(f=y_i^2-x_i^3\) would give
\(f'(t)=0\bmod p\), contrary to \(p\nmid\Delta\). Both factors
are then units. This proves virtuality of the ratios in the integral
case. All valuations here are normalized in the actual local field.

If \(v_p(t)<0\), use the reciprocal coordinate again. Now \(2a\)
is a unit, \(f^*\equiv a^2\pmod p\), and the negative branch has
all four ratios in \(1+p\mathbb Z_p\). The norm identity proves
virtuality on the other branch. This covers every remaining finite
prime, including primes moving with either coordinate of \(t\).

At each prime in (4), choose the split completion with
\(y\equiv r_q\pmod q\). Its four ratios are units with the displayed
characters. Any cube relation among the four specialized classes
would therefore be in the kernel of the displayed matrix. Its
determinant is \(1\bmod3\), so the relation is zero. In particular
no Hilbert irreducibility exceptional set is being silently removed.

The conditions are concrete coordinate congruences. Write \(t=n/d\)
and impose
\[
 n\equiv1\pmod m,\quad d\equiv0\pmod m;
 \qquad n\equiv t_q\pmod q,\quad d\equiv1\pmod q.
 \tag{8}
\]
The Chinese remainder theorem supplies fixed residue classes for
\(n,d\) modulo \(m\cdot73\cdot79\cdot97\cdot103\). A common
factor of \(n,d\) is a unit at every specified prime, so reducing
the fraction does not change (5). Positivity in LY1 makes \(K\)
real. The split calculation at two gives \(D\equiv1\pmod8\);
splitting at three gives \(3\nmid D\). Thus \(-3D\) is the
fundamental discriminant of the imaginary mirror.

For completeness, the class-group transfer uses the same argument
as [MT2](QUADRATIC_MIRROR_RANK_TWO_FROM_MARKED_MESTRE_TWISTS.md#the-common-local-conditions-and-the-mirror-classes).
The exact virtual-unit sequence has a one-dimensional unit kernel
in a real quadratic field, so four independent virtual classes give
real class rank at least three. In \(M=K(\zeta_3)\) they remain
independent by the degree-two norm. Their cubic radicals give
\(N/M\) with group \((\mathbb Z/3)^4\), unramified away from
three by virtuality, and split at three by (7).

Let \(k=\mathbb Q(\sqrt{-3D})\). The involution of \(M/k\)
inverts both \(\zeta_3\) and the four Kummer classes, the latter
by (2). Consequently it acts trivially on \(\operatorname{Gal}(N/M)\).
The resulting central extension by an order-two group splits;
the fixed field of its order-two factor is an abelian degree-81
extension of \(k\). Its inertia at any finite prime has order
dividing both a power of three and the degree two of \(M/k\),
since \(N/M\) is unramified there, so is trivial. At three the
same argument with local degrees gives complete splitting.
Class field theory proves the second inequality in (6).
There are \((3^4-1)/2=40\) cyclic cubic quotients. Conjugation
over \(\mathbb Q\) acts by inversion on the ideal class group of
\(k\); their normal closures are therefore dihedral degree-six
extensions, with distinct complex cubic subfields up to isomorphism
and fundamental discriminant \(-3D\), as in MO3.

## The degree-ten count

**Corollary LY3.** The fields in LY2 include at least
\[
 \gg\frac{X^{1/5}}{(\log X)^2}
 \tag{9}
\]
distinct real fundamental discriminants \(D\le X\). Their imaginary
mirrors are distinct and have absolute discriminant at most \(3X\).

Let \(H(n,d)=d^{10}f(n/d)\), the integral binary form of degree
ten. It is squarefree, since \(f\) is squarefree and the leading
coefficient is nonzero, so the point at infinity is not a root.
It is neither a constant times a linear power nor divisible by a
nonconstant square. Apply
[Stewart–Top, Theorem 2](sources/stewart-top-ranks-of-twists.pdf),
printed p. 954, with \(k=2,r=10\), in the two coordinate classes
(8). It gives \(\gg Y^{2/10}/(\log Y)^2\) distinct squarefree
parts of size at most \(Y\). All are positive, by LY1. Apart
from the single squarefree part one, they give quadratic fields
satisfying LY2. The possible witnesses with \(d=0\) have only
squarefree part one and are discarded. Even degree ensures that
reducing a common factor of \(n,d\), or dividing by \(d^{10}\),
does not change the squareclass. Since two splits, the positive
squarefree part equals \(D\). This proves (9), with no restriction
on the degrees of the irreducible factors: that restriction belongs
to Stewart–Top's stronger Theorem 1, not the theorem used here.

## Five virtual classes on a degree-thirty cover

Set
\[
 c=-2a=-22570656099294324,\qquad t=cu^3,
 \quad \gamma_i=y-y_i(t)\quad(1\le i\le5).
 \tag{10}
\]
The useful identity is
\[
                    c t^5=(c^2u^5)^3.
 \tag{11}
\]
Thus the common leading term of all five \(\gamma_i\) on the negative
branch is a rational cube, including its constant coefficient.

**Theorem LY4.** For the same finite set \(S\) as in (3), impose
\[
 v_p(u)\le-v_p(2a)-1\quad(p\in S),
 \tag{12}
\]
and the following five congruences on \(u\). The table also records
the exact residue certificates, with the character convention of (4).
\[
\begin{array}{r|r|r|r|r|c|c}
q&u_q&t_q&r_q&\zeta_q&(\gamma_1,\ldots,\gamma_5)\bmod q
 &\text{character exponents}\\ \hline
73&12&51&7&8&(55,12,29,45,22)&(2,1,2,1,0)\\
79&4&41&37&23&(57,32,2,64,12)&(0,2,1,0,0)\\
97&1&21&4&35&(11,48,2,21,11)&(2,2,1,2,2)\\
103&5&3&5&46&(79,100,101,100,88)&(0,0,1,0,2)\\
109&1&88&18&45&(106,84,51,61,55)&(2,1,2,2,0).
\end{array}
\tag{13}
\]
These conditions form a fixed nonempty rational parameter slice.
On every quadratic fiber \(K=\mathbb Q(\sqrt{f(cu^3)})\) of this slice,
the five \(\gamma_i\) give independent virtual cubic classes, all cubes
in both completions above three. If \(D\) is the real fundamental
discriminant, then
\[
 \operatorname{rk}_3\operatorname{Cl}(K)\ge4,\qquad
 \operatorname{rk}_3\operatorname{Cl}(\mathbb Q(\sqrt{-3D}))\ge5.
 \tag{14}
\]
The imaginary mirror has an unramified \((\mathbb Z/3)^5\) extension
split above three and at least \(121\) complex cubic fields of
fundamental discriminant \(-3D\). Every prime in \(S\), in particular
\(2,3,5,11\), splits in the real field.

Put \(n_p=v_p(2a)\). Condition (12) gives
\[
 v_p(t)=n_p+3v_p(u)\le-2n_p-3=-k_p.
\]
On the negative branch in the proof of LY2, write \(y^*=-av\) with
\(v^2=f^*/a^2\) and \(v\) close to one. Equation (11) gives
\[
 \frac{\gamma_i}{(c^2u^5)^3}
 =\frac{y^*-y_i^*}{c}
 =1+\frac{a(1-v)-(y_i^*-a)}{c}.
 \tag{15}
\]
With \(h=v_p(1/t)\), the error numerator has valuation at least
\(h-v_p(a)-v_p(2)\), including at two. Division by \(c\) leaves
valuation at least \(h-2n_p\ge3\). Thus (15) is in
\(1+p^3\mathbb Z_p\), hence is a cube. Its normalizing factor is
an actual rational cube, so each \(\gamma_i\) is a cube on this
branch. Its norm \(x_i(t)^3\) makes it a cube on the other branch.
This proves local cubeness and virtuality at every prime in \(S\).

For \(p\notin S\) with \(t\) integral, apply the integral argument
in LY2 to each individual \(y-y_i\). If \(f(t)\) is a unit, one
factor of \((y-y_i)(y+y_i)=-x_i^3\) is a unit. If it is not, the
good discriminant forces \(y_i(t)\) to be a unit. Each factor
therefore has valuation divisible by three. If \(t\) has a pole,
then \(c\) and \(2a\) are units and the negative branch gives
\[
            v_p(\gamma_i)=5v_p(t)=15v_p(u)\equiv0\pmod3.
\]
The cube norm gives the same conclusion on the other branch.
All moving primes are covered. On a quadratic fiber all \(x_i(t)\)
are nonzero, so all evaluations and norms are nonzero.

The checker verifies \(t_q=cu_q^3\), \(r_q^2=f(t_q)\ne0\), all
five unit residues, and their characters. The five-by-five matrix
in (13) has determinant \(2\bmod3\). Thus the specialized classes
are independent, including the fifth class. Equation (11) does not
change any cubic character.

For an explicit congruence description put
\(m_5=\prod_{p\in S}p^{n_p+1}\). One can replace it by
\((6a|\Delta|)^2\). Write \(u=N/E\), and require
\(N\equiv1, E\equiv0\pmod{m_5}\), and
\(N\equiv u_q,E\equiv1\pmod q\) for the five primes in (13).
The Chinese remainder theorem gives compatible fixed coordinate
classes. The checker verifies \(5,11\mid\Delta\); two and three
belong to \(S\) by definition. All these primes split by the local
calculation. The real unit kernel loses at most one of five
directions. The mirror descent in LY2, now with five classes, gives
degree \(3^5=243\) and \((243-1)/2=121\) cubic quotients. This
proves LY4. No condition on a fundamental unit is needed, and no
real class rank five is concluded.

**Corollary LY5.** The paired rank bounds (14) hold for at least
\[
                   \gg X^{1/15}/(\log X)^2
 \tag{16}
\]
distinct real fundamental discriminants \(D\le X\), with distinct
imaginary mirrors of absolute discriminant at most \(3X\).

The polynomial \(f(cu^3)\) is squarefree of degree thirty. A common
zero of it and its derivative would give a multiple root of \(f\),
or \(u=0\), contrary to \(f(0)\ne0\). It is positive on the real
line by LY1. Apply Stewart–Top Theorem 2 to its degree-thirty
homogenization in the fixed coordinate congruences just given.
The squarefree parts number at least
\(\gg X^{2/30}/(\log X)^2\). As in LY3, discard the single
squarefree part one; common coordinate factors change values by
thirtieth powers, and splitting at two makes the remaining positive
squarefree parts equal the real fundamental discriminants.

For these class-rank thresholds, the inspected Levin–Yan–Wiljanen
Theorem 1 gives \(\gg X^{1/30}/\log X\) fields, with a different
local specialization. However,
[Kulkarni–Levin, Theorem 1.5](sources/kulkarni-levin-hilbert-irreducibility-2111.15582.pdf),
arXiv:2111.15582v1, printed p. 3, already proves both the real
rank-four and imaginary rank-five counts with exactly the lower bound
in (16). Its proof on printed p. 8 replaces the earlier quantitative
Hilbert irreducibility input by their Theorem 1.2. Thus the exponent
and logarithmic denominator in (16) are not a new enumerative improvement
over that literature. The additional statements to compare are the
simultaneous mirror pairing, the five specialized local-cube classes
and the prescribed splitting structure. The cited count alone makes
no assertion of their novelty. The targeted primary-source search did
not supply an infinitude theorem for real rank five or the LY4 unit
condition; this is not a comprehensive historical-priority determination.

## Infinite towers with several completely split primes

**Corollary LY6.** Both families have infinite everywhere unramified
pro-\(3\) mirror extensions \(L/k\), Galois over \(\mathbb Q\), with
\[
         \operatorname{Gal}(L/k)^{\mathrm{ab}}
                    \simeq\operatorname{Cl}(k)_{(3)}.
 \tag{17}
\]
For the degree-ten family, the primes above \(2,3\) split completely
in \(L/k\), and the number of such imaginary quadratic fields with
absolute discriminant at most \(X\) is
\(\gg X^{1/5}/(\log X)^2\). For the degree-thirty family, all primes
above \(2,3,5,11\) split completely, and the corresponding count is
\(\gg X^{1/15}/(\log X)^2\).

Use the abelianization-preserving
[tower-cut lemma MP6](MIXED_THREE_POWER_CLASS_GROUPS_FROM_TWO_PRIME_FIBRES.md#infinite-towers-with-prescribed-moving-prime-decomposition).
Its precise classical inputs are
[Koch–Venkov, §1 and Satz 2](sources/koch-venkov-imaginary-quadratic-p-towers-1975.pdf),
printed pp. 57–58, 60–62, and
[Hajir–Maire–Ramakrishna, Theorem 1.4 and Propositions 1.6, 1.9](sources/hajir-maire-ramakrishna-cutting-towers-1901.04354v1.pdf),
pp. 4–5. These give the balanced odd-depth presentation and the
effect of killing Frobenius elements in the closed commutator subgroup.
The universal group argument is supplied in MP6; this application
requires the following prime-class and polynomial checks.

In either mirror \(k=\mathbb Q(\sqrt{-3D})\), two is inert since
\(-3D\equiv5\pmod8\), and three ramifies since \(3\nmid D\).
The unique prime above two is principal; the square of the unique
prime above three is \((3)\). Their three-primary ideal classes
are therefore zero. In the second family, five and eleven split
in the real field and \(-3\) is a nonsquare modulo both. They are
inert in \(k\), so their unique prime ideals are principal as well.
Each selected prime class has exact three-primary order one, and
each Frobenius can be killed with one relation of depth at least
two while preserving the entire abelianization (17).

Put \(d=\operatorname{rk}_3\operatorname{Cl}(k)\). With two cuts
and \(d\ge4\), MP6 bounds the polynomial by
\(1-dt+dt^3+2t^2\). At \(t=2/5\) this is
\[
                       (165-42d)/125\le-3/125<0.
\]
For the second family, \(d\ge5\) and four cuts give
\[
             1-dt+dt^3+4t^2\bigm|_{t=2/5}
                      =(205-42d)/125\le-1/25<0.
\]
Both quotients are infinite by the filtered Golod–Shafarevich
criterion. Each selected rational prime has a single prime above
it, so the cuts are stable under quadratic conjugation and \(L\)
is Galois over \(\mathbb Q\). The counts follow from LY3 and LY5
with the real discriminant bound \(X/3\). Complete splitting holds
in the entire infinite extension. The classical tower mechanism
supplies no additional virtual class or real class-rank improvement.

## The real rank-five target

**Proposition LY7 (unit separation and its finite local test).** On a
quadratic LY4 fibre let \(W_K\) be the span of its five marked classes
in the virtual-unit space \(V_3(K)\), and let
\(U_K=\mathcal O_K^\times/\mathcal O_K^{\times3}\). Then
\[
 \dim\operatorname{im}(W_K\longrightarrow\operatorname{Cl}(K)[3])
                         =5-\dim(W_K\cap U_K).
 \tag{18}
\]
Consequently unit separation \(W_K\cap U_K=0\) gives real class
\(3\)-rank at least five, retaining the rank-five imaginary mirror
and splitting statements of LY4 and LY6. A sufficient condition is
the existence of **any integral unit** \(\eta\) whose image in one
completion above three is not a cube. A fundamental-unit certificate
is unnecessary for this sufficient direction.

Indeed apply the exact sequence in MO7 to this five-dimensional
space. Localization vanishes on \(W_K\) by LY4. If an integral unit
has nonzero localization, the map on the one-dimensional \(U_K\)
is nonzero, hence injective, proving separation. Equivalently, if
\(\eta=\pm\epsilon^j\) for a fundamental unit \(\epsilon\), its
noncube localization forces both \(3\nmid j\) and noncube localization
of \(\epsilon\). The criterion on \(\epsilon\) is sufficient, not
necessary for \(W_K\cap U_K=0\).

For an explicit test write \(K=\mathbb Q(\sqrt d)\), with \(d>0\)
squarefree. On LY4, three splits, so \(d\equiv1\pmod3\) and both
completions are \(\mathbb Q_3\). For a certified integral unit
\(\eta=(A+B\sqrt d)/2\), choose either root \(r\) of \(d\bmod9\).
Its local image is a cube if and only if
\[
                  (A+Br)/2\equiv1\text{ or }-1\pmod9.
 \tag{19}
\]
To see this, decompose \(\mathbb Z_3^\times=\{\pm1\}\times
(1+3\mathbb Z_3)\). Cubing maps the second factor onto
\(1+9\mathbb Z_3\), since the convergent three-adic logarithm
identifies cubing with multiplication by three. Thus its unit cubes
are precisely \(\pm(1+9\mathbb Z_3)\). The conjugate of a global
unit is \(N(\eta)/\eta\), and \(\pm1\) are cubes, so the two
places give the same cube/noncube verdict.

The norm equation also gives convenient coefficient tests. If
\(N(\eta)=1\), (19) holds exactly when \(B\equiv0\pmod9\):
the conjugate images are inverses, so (19) implies their difference
\(Br\) is zero, and conversely \(B=0\bmod9\) and the norm
equation imply each image squares to one modulo nine. If
\(N(\eta)=-1\), (19) holds exactly when \(A\equiv0\pmod9\),
because their sum is \(x-x^{-1}\). In particular noncube
localization of a fundamental unit is equivalent to
\(v_3(\log(\epsilon^2))=1\); larger valuation gives a local cube.
This is a regulator condition on the actual field, not a consequence
of the parameter congruences (12)–(13).

The missing infinitude statement is now: prove unit separation on
infinitely many **distinct fields** of LY4. Infinitely many parameters
alone do not establish it. The stronger sufficient route seeks an
integral noncube unit in each of infinitely many such fields. Even
that infinitude would not automatically retain the count (16); that
requires a separate quantitative argument or a uniform admissible
slice. The faster degree-ten real rank-four target below remains a
separate question. Raising exponents in existing cyclic class-group
factors does not increase the five-dimensional image in (18).

## Exact bounded fibre probe

The checker provides a first bounded unit search on the **full** LY4
slice, rather than treating the local control \(u=1/9\) as a global
example. Set
\[
 E=\prod_{p\mid6a\Delta}p^{v_p(2a)+1},\qquad
 Q=73\cdot79\cdot97\cdot103\cdot109.
\]
Take the least nonnegative solution \(N_0\) of
\(N_0\equiv u_q E\pmod q\) for the five rows of (13), then take
the first three positive integers \(N=N_0+jQ\), \(j\ge0\), with
\(\gcd(N,E)=1\). Their parameters \(u=N/E\) satisfy (12) and
(13). Here \(N_0=23040989\), and \(E\) has 91 decimal digits.
The checker factors \(6a|\Delta|\) using a bounded trial step,
verifies every remaining factor is prime, and verifies the product;
it rejects an uncertified composite cofactor.

Put \(T=cN^3\), \(B=E^3\). The exact positive integer
\[
                R=\sum_{i=0}^{10}F_i T^i B^{10-i}
\]
represents \(f(cu^3)\) modulo rational squares. Remove the largest
even power of each prime in \(S\) from \(R\), obtaining \(R_0\).
No full squarefree factorization is claimed or required. Each retained
\(R_0\) is a quadratic nonresidue modulo 59, certifying a quadratic
fibre, and is \(1\bmod3\). Since \(K=\mathbb Q(\sqrt{R_0})\),
an integer Pell solution gives an integral unit of \(K\) even if
\(\mathbb Z[\sqrt{R_0}]\) is a nonmaximal order.

The internally executed probe used Python 3.12.3 and SymPy 1.14.0:

| Numerator \(N\) | Digits of \(R_0\) | Convergents tested | Pell unit found |
|---:|---:|---:|---|
| 23040989 | 2570 | 2000 | no |
| 50266021573 | 2570 | 2000 | no |
| 113069747303 | 2570 | 2000 | no |

For complete quotients \(P_j,Q_j\) of \(\sqrt{R_0}\), use
\(P_0=0,Q_0=1,a_0=\lfloor\sqrt{R_0}\rfloor\) and
\[
 P_{j+1}=a_jQ_j-P_j,\quad
 Q_{j+1}=(R_0-P_{j+1}^2)/Q_j,\quad
 a_{j+1}=\lfloor(a_0+P_{j+1})/Q_{j+1}\rfloor.
\]
For the convergent \(h_j/k_j\), the determinant identity for successive
convergents gives
\(h_j^2-R_0k_j^2=(-1)^{j+1}Q_{j+1}\). Thus \(Q_{j+1}=1\)
certifies a Pell unit. Convergent residues modulo nine suffice for
(19); the complete quotients are computed exactly. The probe outputs
the exact parameter, a SHA-256 of the decimal radicand, a nonsquare
witness and the step budget, so these inputs can be regenerated.

No maximal-order fundamental unit or local separation verdict was
obtained in these three fibres. The result only excludes a Pell unit
among the stated prefixes. It supplies neither a separation frequency,
distinctness of these three fields, nor evidence against infinitude.
A locally cubic unit of a nonmaximal order would also leave the
fundamental-unit question undecided: its exponent in the maximal-order
unit group might be divisible by three. The checker demonstrates this
failure mode with \(3+\sqrt{10}\), which is a local noncube, and its
cube \(117+37\sqrt{10}\), which is locally cubic.

The next useful step is a source of certified units on a thinner
LY4 subfamily, or a parameter model with substantially smaller field
discriminants that preserves all five virtual classes. Merely extending
these prefixes does not provide a mechanism for infinitude. Any
explicit unit construction must establish integrality, its norm and
the noncube residue in the same fibres, and then prove distinct-field
infinitude; congruences on the parameter alone do not do this.

## Replay and limits

From the repository root in its existing Linux environment:

```sh
.venv/bin/python research/arithmetic/verify_lyw_degree_ten_mirror_slice.py
.venv/bin/python research/arithmetic/verify_lyw_degree_ten_mirror_slice.py --unit-probe --steps 2000 --fibres 3
```

Internally replayed with Python 3.12.3 and SymPy 1.14.0. The checker
verifies all five cube identities, squarefreeness, the exact Sturm
signs, the four good-prime character rows and determinant, and the
local control at \(t=3^{-5}\), where the four negative-branch ratios
are \(1\bmod27\). For LY4–LY6 it also verifies squarefreeness of
the degree-thirty cover, the global cube normalization (11), all five
character rows and determinant, and the five normalized values
\(1\bmod27\) at \(u=3^{-2}\). The two tower polynomial values are
checked with exact rational arithmetic. It rejects a perturbed norm
identity and duplicated character rows. The default replay also checks
the modulo-nine unit criterion on a noncube, its locally cubic cube,
and a half-integral unit, and rejects invalid local-unit data. The
separate bounded probe has exactly the scope recorded above. The universal local argument,
Kummer descent, tower cuts and application of the counting theorem are written proofs
above, not conclusions established merely by this finite replay.
Independent replay, formal verification, external review and
historical-priority review have not occurred. A specialist audit of
the local inequalities, mirror descent and exact counting interface
remains appropriate before any manuscript or release.

The degree-ten family retains real class rank at least three; its
unit loss is not removed. The degree-thirty family reaches real
class rank four by adding a fifth direction and retains a possible
one-dimensional unit loss. Neither supplies the fundamental-unit
condition in MO7. The counts and towers concern these class-group
constructions; the three fixed elliptic-twist properties of MO3 have
not been transferred. No generic fifth direction on the exhausted
cubic covers of the marked octic is asserted.

A separate target on the faster degree-ten family is real
class rank four on infinitely many LY2 fields. If \(W\) is its
four-dimensional virtual space, \(W\cap\mathcal O_K^\times/
\mathcal O_K^{\times3}=0\) suffices. As in MO7, a fundamental unit
that is not a cube in one completion above three implies this
separation, because all of \(W\) localizes to zero there. No such
infinite subfamily is established here. Retaining the exponent
\(1/5\) needs a uniform admissible slice or a separate quantitative
argument; LY4–LY5's degree-thirty count does not settle that target.
