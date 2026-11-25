# Two fixed-phase moments from a surjective coherent Frobenius step

A coherent Frobenius step from a full de Rham character space to a
coherent line is surjective at every good prime in the prescribed
residue class. Its kernel can contain one specified input. Two moment
orders suffice to test a basis of inputs in the rank-two case. The
second moment uses both terms in the regular lift of the negative
angular monomial. This gives fixed four-branch families whose moments
are detected even when coherent Frobenius on the whole auxiliary
Jacobian is nilpotent.

## The strip and its two moments

Let \(n\) be even and let \(3\le a\le n/2\). Let \(P(z)\) be a
polynomial over a characteristic-zero field \(k\), of degree \(a+1\),
with exactly two distinct nonzero geometric roots of multiplicities
\(m_1,m_2\). Suppose \(s\) is a unit representative modulo \(n\) with
\(1\le s<n\), and that
\[
  2as<n(a+1),\qquad
  1+\left\lfloor\frac{as}{n}\right\rfloor
    -\left\lfloor\frac{m_1s}{n}\right\rfloor
    -\left\lfloor\frac{m_2s}{n}\right\rfloor=2.
\]
On the sphere \(2xy+t^2=1\), use the literal Laurent chart
\(y=(1-t^2)/(2x)\). Write \(\mathcal A\) for the functional that
takes the constant term in \(x\) and integrates it over the interval
using the normalized measure. In particular \(\mathcal A(t^j)\) is
\(1/(j+1)\) if \(j\) is even and zero if \(j\) is odd.

**Theorem.** Fix any sphere phase \(f\) satisfying
\[
           j+2i\le n,\qquad (a+1)j+2ai\ge (a+1)n-2a
\]
on the support of its literal Laurent monomials \(x^it^j\). Prescribe
the angular polynomial \(z^{-1}P(z)\) on \(j+2i=n\), where the
coefficient at angle \(i\) means the coefficient of \(x^it^{n-2i}\).
All compatible smaller-weight terms are arbitrary. Every such phase
has least sphere degree \(n+1\) and is an absolute nonpower. For all
sufficiently large rational primes \(p\equiv s^{-1}\pmod n\), put
\[
                         m=\frac{sp-1}{n}.
\]
At least one of \(\mathcal A(f^m)\) and \(\mathcal A(f^{m+p})\) is
nonzero. This set of certifying primes has natural density
\(1/\varphi(n)\). If \(-1\in\langle s\rangle\subset (\mathbb Z/n)^\times\),
the coherent return to character one of the auxiliary cover
\(W^n=z^{n-1}P(z)\) is zero in every good reduction over those primes.

For comparison, the [earlier even-conductor strip](CM_NEWTON_EVEN_CONDUCTOR_STRIP_FACES.md#general-even-conductor-criterion)
allows a weight deficit \(n-j-2i\) of at most \(2(n+i)/(n-1)\)
at a given angle \(i\). Here the allowance is \(2(a-i)/(a+1)\).
Both allowances are two at angle \(-1\); the present allowance is
smaller at every nonnegative angle. In return, the two-moment test
applies at every good sufficiently large prime in the progression,
with no endomorphism hypothesis on the auxiliary Jacobian.

### Regularity and support

Put \(b=2a/(a+1)\). The lower bound is \(j+bi\ge n-b\).
Together with the upper bound it gives \(i\le a\). At a negative
angle \(-v\) the two bounds allow a width
\(2v-b(v-1)=2(a+v)/(a+1)\) in the height \(j\). For \(v>1\)
this is smaller than \(2v\). A nonzero regular coefficient of
\(x^{-v}\) is divisible by \((1-t^2)^v\) and needs width at least
\(2v\), so all angles lie in \([-1,a]\). At angle \(-1\) the
heights lie between \(n\) and \(n+2\). Divisibility by \(1-t^2\)
forces the entire term at this angle to be
\[
              P(0)x^{-1}(t^{n+2}-t^n)=-2P(0)yt^n.
\]
At angle zero only heights \(n-1,n\) are allowed. The negative
extreme and the top terms at nonnegative angles all have regular
lifts, since \(n\ge2a\). In the normal sphere basis with no mixed
\(x,y\) monomial, the nonnegative angles have degree at most \(n\),
and the negative term has degree \(n+1\). Passing to the normal
basis by \(2xy=1-t^2\) cannot increase degree. The top angular
function has order \(-1\) at \(z=0\). With \(z=x/t^2\), taking
the leading coefficient as a rational function of \(t\) over
\(\overline{k}(z)\) therefore rules out any proper power, even a
power of a rational function.

### The second pole and the two inputs

Work at a prime of good reduction for the coefficients and branch
configuration, with \(p>s+n\). Write
\[
              U=[z^m]P(z)^m,\qquad V=[z^{m+p}]P(z)^m,
              \qquad P_j=[z^j]P(z).
\]
All congruences below are in its residue field. Because the maximum
angular-constant height in either moment is at most
\((s+n)p-1<p^2-1\), multiplying by \(p\) makes all summands
integral and retains only even heights with denominator divisible
by \(p\).

For the first moment, the support of \(f^m\) in angle zero
has maximum height \(nm=sp-1\) and width at most \(bm<p\).
Its only possible pole has height \(sp-1\), so
\[
                         p\mathcal A(f^m)\equiv U/s.
\]
The parity uses that \(s\) and \(p\) are odd.

For the second moment use \(f^{m+p}\equiv f^m f^p\), with
\(f^p\) obtained by multiplying all angles and heights of \(f\)
by \(p\) and raising its coefficients to the \(p\)-th power. If
the angle from \(f^m\) is \(i_0\) and the angle from \(f\) is
\(i_1\), an angular-constant term has \(i_0+p i_1=0\).
Since \(m<p\), the lower bound \(i_0\ge-m\) excludes
\(i_1>0\). The lower bound \(i_1\ge-1\) leaves only
\((i_0,i_1)=(0,0)\) and \((p,-1)\).

If \(j_0\) is the height from \(f^m\), an interval denominator
divisible by \(p\) requires \(j_0\equiv-1\pmod p\). The
support at angle \(i_0\in\{0,p\}\) has upper height
\(sp-1-2i_0\) and width at most
\(bm-(2-b)i_0\le bm<p\). Thus the only possible \(j_0\)
is that upper height, with coefficient \(U\) at angle zero
and \(V\) at angle \(p\). Both heights are even. At angle
zero in the factor \(f\), the height \(n-1\) has odd parity
and cannot contribute. At angle \(-1\), the two terms forced
by regularity both contribute. This gives the congruence
\[
 p\mathcal A(f^{m+p})\equiv
  \frac{UP_1^p}{s+n}
  +VP_0^p\left(\frac{1}{s+n}-\frac{1}{s+n-2}\right).
\]
The coefficient of \(V\) is nonzero since \(P_0\) is a unit.
The two residues therefore determine whether both input
coefficients are zero. Strict lower terms, including odd-height
terms, do not change either congruence.

### Surjectivity in every good reduction

On the normalized connected cover \(W^n=z^{n-1}P(z)\), use deck
characters indexed by the powers of \(W\). The character-one
coherent line is generated by \([W/z]\). Indeed the normalized
finite generator is \(W\), and its generator at infinity is
\(W/z^2\), because the total polynomial degree is \(n+a\).

If the roots of \(P\) are \(\alpha_1,\alpha_2\), put
\[
 Q_s(z)=\prod_{\nu=1}^2(z-\alpha_\nu)^{\lfloor sm_\nu/n\rfloor},
 \qquad E_s=\frac{W^s}{z^{s-1}Q_s(z)}.
\]
The finite eigensheaf has generator \(E_s\), and the condition
on the floors gives generator \(E_s/z^3\) at infinity. Thus
its two coherent basis classes are \([E_s/z]\) and
\([E_s/z^2]\). The classes \([W^s/z^s]\) and
\([W^s/z^{s+1}]\) also form a basis. In the preceding basis
they are the classes of \(E_sQ_s(z)/z\) and \(E_sQ_s(z)/z^2\),
with triangular change of basis of diagonal \(Q_s(0)\ne0\).
Frobenius takes them to character one because \(sp=nm+1\).
Their images are the classes of \(Wz^{-m-1}P(z)^m\) and
\(Wz^{-m-p-1}P(z)^m\), whose entries on \([W/z]\) are \(U,V\).
[Bouw, §§4–5](https://doi.org/10.1023/A:1017513122376)
· [Lin–Mantovan–Singhal, §2.4](https://arxiv.org/html/2406.09632v1#S2.SS4)

The de Rham spaces at both characters have dimension two. In fact
each primitive character sees all four branch points, and the two
conjugate coherent dimensions sum to two by the Kummer eigensheaf
formula and Serre duality. In characteristic \(p\), relative
Frobenius on de Rham cohomology factors as the projection onto the
Frobenius-twisted coherent source followed by the inclusion of the
conjugate subspace of de Rham cohomology at the target. The inclusion
is injective by the Cartier conjugate spectral sequence, since the
curve has coherent cohomological dimension one. On the selected
source, the Hodge quotient already has full dimension two, so the
inclusion maps it isomorphically onto the entire two-dimensional
target de Rham space. Projecting there to the coherent line gives
a surjection. Thus \(U,V\) cannot both vanish. This uses the same
structural Frobenius factorization as the Hasse–Witt triple.
[Katz, §7](https://web.math.princeton.edu/~nmk/old/nilpconn.pdf)
· [Moonen, §§2.5 and 3.1](https://www.math.ru.nl/~bmoonen/Papers/CDIPart1.pdf)

Character \(-s\) has coherent dimension zero by conjugate duality.
It lies on the orbit of character one when \(-1\in\langle s\rangle\),
so the full coherent return is zero. The arithmetic progression has
the claimed density after omitting finitely many bad rational primes.
For an arbitrary characteristic-zero coefficient field, specialize
its finitely generated coefficient algebra to a number field after
including both roots and inverting their product, difference and the
leading coefficient. This preserves the branch multiplicities and
the support and specializes each moment. A nonzero specialized
moment implies that the original moment is nonzero. This one
specialization works for all the primes used for that fixed phase.

## Fixed families and finite checks

The following rows apply to \(P(z)=A(z-1)^{m_1}(z-c)^{m_2}\)
for every \(A\ne0\) and \(c\ne0,1\) over any characteristic-zero
field. In each row the prime class generates the primitive
character group, which contains a zero coherent dimension.
The primitive Jacobian factor therefore has \(p\)-rank zero
at all the certifying primes.

| \(n\) | \(a\) | \((m_1,m_2)\) | \(s\) | Prime class | First moment order |
| ---: | ---: | --- | ---: | --- | --- |
| 10 | 5 | \((3,3)\) | 3 | \(p\equiv7\pmod{10}\) | \((3p-1)/10\) |
| 14 | 5 | \((3,3)\) | 3 | \(p\equiv5\pmod{14}\) | \((3p-1)/14\) |
| 18 | 8 | \((2,7)\) | 5 | \(p\equiv11\pmod{18}\) | \((5p-1)/18\) |

In particular, at these narrower strips the finite exceptional
parameter sets of the [one-moment argument](CM_NEWTON_FOUR_POINT_MONODROMY.md)
are unnecessary. For algebraic parameters, use any coefficient
place over each certifying rational prime, with no residue-degree-one
condition. The relative composition is a single absolute Frobenius
step and the displayed congruences are valid in every finite residue
field away from the bad primes.

At \(n=10,p=17,m=5\), the polynomial
\(P=(z-1)^3(z-2)^3\) gives \((U,V)=(14,7)\) modulo 17.
The normalized residues at orders \(5,22\) are \(16,8\).
For \(P=(z^2-1)^3\) at the same prime, the first input is zero
by parity, while the second is \(\binom{15}{11}\equiv5\pmod{17}\).
Here the two normalized residues are \(0,16\). The exact replay
computes moments by multiplying the full literal Laurent polynomial,
including odd smaller-weight terms. [Exact replay](verify_cm_newton_two_input_strip_moments.py)

## Unbounded conductor and Jacobian p-rank zero

For every \(h\ge3\), put \(N=3^{h-1}\) and choose
\[
 n=6N,\qquad a=2N,\qquad (m_1,m_2)=(N,N+1).
\]
Use all \(1\le s<N\) with \(s\equiv5\) or \(11\pmod{18}\).
Writing \(s=6q+5\) gives
\(\lfloor as/n\rfloor=2q+1\) and
\(\lfloor m_1s/n\rfloor=\lfloor m_2s/n\rfloor=q\), since
\(m_2s/n=s/6+s/n<q+1\). Thus the coherent dimension is two,
and \(s<n/6\) gives the width inequality. The four-branch cover
is connected because the exponent at zero is \(n-1\).

Each such \(s\) generates the units modulo \(2\cdot3^h\).
Its reduction modulo nine is 5 or 2, and \(s^2-1\) has
3-adic valuation one. The binomial expansion gives
\(v_3(s^{2\cdot3^j}-1)=1+j\), so \(s\) has order
\(2\cdot3^{h-1}\) modulo \(3^h\) and hence modulo \(n\).
For \(p\equiv s^{-1}\pmod n\), the coherent Frobenius orbit
of any character is thus the full set of characters of its conductor.

For conductor \(n\), the orbit contains the zero at character
\(-5\). For conductor \(3^h\), the character with index \(4\) in
\(\mathbb Z/n\) has coherent dimension two, by
\(4a/n=4/3\) and \(4m_i<n\). Its conjugate has dimension zero.
Every other nontrivial character conductor divides \(a\). The
exponent at infinity is zero modulo such a conductor, so at most
three branches contribute to that character. Its coherent dimension
and that of its conjugate sum to at most one. These conjugates lie
on the same Frobenius orbit and at least one dimension is zero.
The trivial character has coherent dimension zero since the quotient
is the projective line. Thus coherent Frobenius is nilpotent on the
entire first coherent cohomology, and the whole auxiliary Jacobian
has \(p\)-rank zero in every good reduction over these primes.
Since \(N\equiv9\pmod{18}\), exactly \(N/9\) values of \(s\)
are used. The union of their disjoint prime classes has natural
density \((N/9)/\varphi(6N)=1/18\), independently of the
conductor. Each fixed conductor gives unmarked nonzero moments of
phases of least sphere degree \(n+1\) on that prime set.

## A full coherent input space of arbitrary rank

**Theorem.** In the same even-conductor strip, let \(P\) have \(R\ge2\)
distinct nonzero roots of multiplicities \(m_1,\dots,m_R\). Replace the
coherent dimension condition by
\[
 1+\left\lfloor\frac{as}{n}\right\rfloor
      -\sum_{i=1}^R\left\lfloor\frac{sm_i}{n}\right\rfloor=R,
 \qquad 2as<n(a+1).
\]
For every fixed compatible phase and all sufficiently large primes
\(p\equiv s^{-1}\bmod n\), at least one of the \(R\) unmarked moments
with orders \(m+kp\), \(0\le k<R\), is nonzero, where \(m=(sp-1)/n\).
The degrees and absolute nonpower conclusion are as before. If
\(-1\in\langle s\rangle\), the coherent Frobenius return on the orbit
of character one is zero at these primes. If \(s\) generates the units
modulo \(n\), the primitive conductor-\(n\) factor has \(p\)-rank zero.

Use the same normalized finite eigensheaf generator \(E_s\), with
\(Q_s\) now taken over all \(R\) roots. At infinity the generator is
\(E_s/z^{R+1}\), and the \(R\) classes \([W^s/z^{s+j}]\), \(0\le j<R\),
form a basis by the triangular multiplication by \(Q_s\). The
coefficients of their coherent Frobenius images in character one are
\[
                    A_j=[z^{m+jp}]P(z)^m.
\]
Both primitive character de Rham spaces have dimension \(R\), since
all \(R+2\) branch points contribute. The source coherent quotient
has dimension \(R\) and the target coherent quotient is a line.
The Cartier factorization used for two inputs therefore gives a
surjection at every good prime in the progression. The \(A_j\) cannot
all vanish.

The \(k\)-th moment gives a triangular linear form in the \(A_j\).
In characteristic \(p\) write \(f^{m+kp}=f^m(f^k)^p\). The angles
in an angular-constant term are \((jp,-j)\), with \(0\le j\le k\),
since \(f^m\) has angles at least \(-m>-p\) and \(f^k\) has angles
at least \(-k\). At angle \(jp\), the support interval in \(f^m\)
has width at most \(bm<p\), where \(b=2a/(a+1)\). Its upper
height is \((s-2j)p-1\). The coherent dimension condition implies
\(\lfloor s/2\rfloor\ge R-1\), because \(a\le n/2\), so this height
is nonnegative for the relevant \(j\). It is the only height in the
interval congruent to \(-1\) modulo \(p\), and its coefficient is
\(A_j\). For \(p>s+(R-1)n\), all interval denominators in the moment
have valuation at most one. The contribution of an even height \(u\)
from angle \(-j\) in \(f^k\) to \(p\mathcal A(f^{m+kp})\) is
\(A_j[x^{-j}t^u](f^k)^{{[p]}}/(s-2j+u)\), where \({[p]}\) raises
coefficients alone to the \(p\)-th power. Odd heights give zero.

The diagonal comes from the unique way to achieve angle \(-k\) in
\(f^k\), using \(P_0^k t^{nk}(t^2-1)^k\). The diagonal coefficient,
with \(L=s+kn-2k\), is
\[
 P_0^{kp}\sum_{\nu=0}^k\frac{(-1)^{k-\nu}\binom{k}{\nu}}{L+2\nu}
 =\frac{(-1)^k 2^k k!P_0^{kp}}{\prod_{\nu=0}^k(L+2\nu)}.
\]
For instance, the identity follows from the partial fractions of the
product on the right. Every diagonal entry is nonzero at the primes
under consideration. This recovers the \(A_j\) from unmarked moments,
and proves the theorem over characteristic zero by the same fixed
specialization. If \(-1\in\langle s\rangle\), each primitive orbit
containing a zero-dimensional coherent character has zero return.
In a progression where \(s\) generates all the units there is only
one orbit, and the full primitive factor has \(p\)-rank zero.

For unbounded conductors at each fixed rank \(R\), choose \(L\ge2\) with
\(7R\le2\cdot3^L\). For \(u\ge2L+1\), put \(M=3^{u-L}\),
\(n=2\cdot3^u\), \(a=RM\), and use \(R-1\) multiplicities \(M\) and
one multiplicity \(M+1\). Use all \(1\le s<M\) satisfying
\(s\equiv-7\bmod 2\cdot3^L\). Write \(s/(2\cdot3^L)=v+1-7/(2\cdot3^L)\)
with \(v\ge0\) integral. Then \(\lfloor as/n\rfloor=Rv+R-1\) and
each \(\lfloor sm_i/n\rfloor=v\), since \(s/n<1/(2\cdot3^L)\).
The width condition holds since \(2s<n\). These \(s\) are all eleven
modulo eighteen and generate the units modulo \(n\). The primitive
factor has \(p\)-rank zero on all these progressions. Since
\(M\equiv3^L\bmod 2\cdot3^L\) and \(2\cdot3^L-7>3^L\), the number
of progressions is \((M-3^L)/(2\cdot3^L)\). Their total density is
\((1-3^L/M)/(4\cdot3^{2L-1})\), bounded below independently of \(u\).

For \(R=3\) one can instead take \(n=6N\), \(N=3^{u-1}\),
\(a=3N\), and multiplicities \((N,N,N+1)\). The same \(N/9\) values
\(s<N\), \(s\equiv5,11\bmod18\), used above give coherent dimension
three and generate the units. On the resulting prime set of density
\(1/18\) the entire primitive factor has \(p\)-rank zero.

For an exact third-input example, use \(n=18\), \(a=8\), \(s=5\),
\(p=83\), and \(P=(z^3-1)^3\). Here \(m=23\) and the three inputs
are \((0,0,74)\) modulo 83. The two initial indices are not multiples
of three. The third coefficient is \(\binom{69}{63}\), and its nonzero
diagonal factor is 41 modulo 83. The three normalized moment residues
are \((0,0,46)\). The [exact replay](verify_cm_newton_two_input_strip_moments.py)
includes an odd strict lower term and computes all three residues in
the literal Laurent chart.
