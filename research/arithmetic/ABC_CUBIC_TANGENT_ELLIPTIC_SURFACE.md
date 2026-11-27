# The Mordell–Weil group of the cubic tangent surface

The cubic tangent family comes from one elliptic
surface whose six geometric Mordell–Weil directions
are all explicit. The arithmetic rank depends on
the intersection of the ground field with the normal
closure of \(\mathbb Q(\sqrt[3]2)\). The extra section
also survives on the dyadic squarefree prime families.

## The geometric span and its field of definition

Let \(r=2^{-1/3}\) be real and let \(\zeta=\zeta_3\).
Use \(K=\mathbb Q(\zeta)\), \(L=\mathbb Q(r)\) and
\(M=KL\). Over \(\mathbb Q(t)\) consider
\[
 \begin{split}
 \mathcal E:\quad y^2&=x^3-\frac{3t^4-6t^2-1}{4},\\
 P_\pm&=\left(t^2\pm t,\ t^3\pm(3t^2+1)/2\right),\\
 Q&=\left(r^2(t^2-1),\ t(t^2-3)/2\right).
 \end{split}
\]
The \(P_\pm\) are the prime-progression points of
[CT2](ABC_CUBIC_TANGENT_NORM_RANK_AND_CONDUCTORS.md#rank-from-a-second-descent)
at \(t=p/h\) after rescaling by \(h^2,h^3\).
The formula for \(Q\) follows from
\[
     (t^2-1)^3-(3t^4-6t^2-1)=t^2(t^2-3)^2
\]
and \(r^6=1/4\). The automorphism
\([\zeta](x,y)=(\zeta x,y)\) acts on the Mordell–Weil group.

**Theorem TS1.** The geometric Mordell–Weil group is
torsion-free of rank six, and every section is
defined over \(M\). Put \(W=F\cap M\) for a number
field \(F\subset\overline{\mathbb Q}\). The rank over
\(F(t)\) is \(2\) if \(W=\mathbb Q\), \(3\) if \(W\)
is a cubic subfield, \(4\) if \(W=K\), and \(6\) if
\(W=M\). The field \(M\) is minimal for defining all
geometric sections.

The surface has four type II fibres at the simple
quartic zeros and type IV at infinity. In the chart
\(s=1/t\) the short model is
\[
             Y^2=X^3+s^2(-3+6s^2+s^4)/4 .
\]
This follows from Tate's algorithm, or the resolution
of the \(A_2\) surface singularity at \(s=X=Y=0\).
The fundamental line bundle is \(\mathcal O(1)\), so
\(\chi=1\) and the surface has second Betti number
\(10\). The geometric rank is at most six by
[Schütt–Shioda, Corollary 6.13](https://arxiv.org/pdf/0907.0298),
because the only reducible fibre contributes
two root-lattice directions. Torsion is zero
over \(\overline{\mathbb Q}(t)\): at a type II fibre
the component group is trivial and the
nonsingular reduction is additive, while
the characteristic-zero formal group is torsion-free.

The two \(P_\pm\) are independent over \(\mathbb Q(t)\)
by the height calculation below, or by specializing
to the fixed residue tests of CT2. They are then
independent over \(\mathbb Q(\zeta)\) in the rationalized
geometric group. Conjugating a proposed relation
with coefficients \(a+b\zeta\), and subtracting,
reduces it to two rational relations, since
\([\zeta]-[\zeta^2]\) is a nonzero isogeny.
The point \(Q\) is nonzero and torsion-free.
A generator of \(\operatorname{Gal}(M/K)\)
acts on \(Q\) by a nontrivial \([\zeta]^i\) and
fixes both \(P_\pm\). Averaging over that cyclic
group kills \(Q\) and fixes the \(P_\pm\), proving
that the three are independent over \(\mathbb Q(\zeta)\).
They supply six geometric directions, hence the
upper bound is attained.

The rationalized geometric group accordingly has
the \(S_3\) representation \(2\mathbf1+2\epsilon+
\mathrm{Std}\), where \(\epsilon\) is the sign
character and \(\mathrm{Std}\) has dimension two.
In detail, each of the two \(P_\pm\) spans a
\(\mathbb Q(\zeta)\)-line on which the cubic subgroup
acts trivially and complex conjugation acts on
the coefficients. The line from \(Q\) is the
standard representation, with nontrivial cubic
eigenvalues and complex conjugation acting on
the coefficients. Since the geometric group
has no torsion, every automorphism fixing \(M\)
therefore acts trivially on the full group. It
cannot be a smaller common field, since sections
\([\zeta]P_+\) and \(Q\) recover \(\zeta\) and \(r^2\)
from their abscissas. The fixed-space dimensions
for the subgroup fixing \(W\) give the rank table.

## Integral bases over the real fields of smallest degree

**Proposition TS2.** The points \(P_-,P_+\) are an
integral Mordell–Weil basis over \(\mathbb Q(t)\).
The three points \(P_-,P_+,Q\) are an integral basis
over \(L(t)\). The height matrices in the
Shioda pairing are \(2I_2\) and \(2I_3\).

The three sections are integral in both charts,
are disjoint from the zero section, and
their values in the chart at infinity
have nonzero \(X,Y\). Thus they meet the
identity component there. Every section over
a real field also meets the identity component
at the type IV fibre: its two nonidentity
components are conjugate over \(\mathbb Q(\sqrt{-3})\),
as the exceptional tangent cone at infinity
is \(Y^2+3s^2/4=0\).
A section meets the smooth part of the fibre,
so it cannot meet the common point of the
nonidentity components.

In the height pairing of
[Schütt–Shioda, §11.8](https://arxiv.org/pdf/0907.0298),
all corrections thus vanish for real-field sections
and each displayed section has height \(2\).
The two \(P_\pm\) meet once at infinity, with
abscissa difference \(2s\), and nowhere at
finite \(t\). The section \(Q\) meets \(P_+\)
once at \(t=-1\) and is disjoint from it at
infinity. The other root of the abscissa
difference is \(t=-r^2/(1-r^2)\), where their
ordinates are nonzero and opposite.
The intersection at \(t=-1\) is transverse
since \(2r^2-1\ne0\).
The involution \(t\mapsto-t\) takes \(P_+\)
to \(-P_-\) and \(Q\) to \(-Q\), so the
intersection with \(P_-\) also has multiplicity
one. Hence all off-diagonal height terms are zero.

The height pairing on all real-field sections
is integral, with even self-pairings. A section
in the rational span of the displayed basis
must have half-integral coefficients, by pairing
with each basis element. Modulo the generated
lattice one can take each coefficient to be
zero or \(1/2\). A nonzero such vector in rank
two or three has height less than \(2\), contrary
to the even positive self-pairing. The ranks
of TS1 now give both integral basis assertions.

## Specializing the cubic section at prime witnesses

**Theorem TS3.** On the signed squarefree prime
families of CT4, for all sufficiently large
parameter pairs, the three explicit points
are independent over \(L\), and their images
together with their \([\zeta]\)-translates
are independent over \(M\). The corresponding
rank lower bounds over \(\mathbb Q,L,K,M\) are
\(2,3,4,6\). For each sign the prime-\(k\) subclass
gives \(\gg Y^{1/5}/(\log Y)^4\) distinct rational
conductors at most \(Y\) realizing all four
rank lower bounds.

More generally the rank assertions hold for
the primitive even-\(h\) pairs of
[CT6](ABC_CUBIC_TANGENT_NORM_RANK_AND_CONDUCTORS.md#independence-from-the-full-simple-prime-support)
when \(D_1\ge p^\eta\) for fixed \(\eta>0\)
and \(p\) is sufficiently large.
In the integral model over the parameter pair,
the new point is
\[
        Q_{p,h}=(r^2(p^2-h^2),\ p(p^2-3h^2)/2).
\]
It is non-torsion once \(D_1\) is large:
take a prime \(q\) above an absolute fixed
threshold dividing \(D_1\). The cubic field
is unramified at \(q\), and at a prime over
\(q\) the equation reduces to a cusp because
\(v_q(R)=1\). Both coordinates of this point
are units, so it specializes to a nonzero
element of the additive nonsingular group
of order \(q\). Rational
torsion over the fixed field \(L\) has order
prime to all sufficiently large \(q\), uniformly
for the curves in question. For this last
fact one can reduce at a fixed prime above
\(5\) in \(L\). The prime-to-\(5\) torsion
injects into the finite quotient by the formal
group. The residue field has at most \(125\)
elements and the additive component group
has size at most four after minimizing the
Mordell model, so the quotient order is
bounded uniformly. This contradicts a torsion
specialization at \(q\).

The two rational points are independent by
CT6 or CT2. Their independence over the
endomorphism field and their independence
from \(Q_{p,h}\) use the same conjugation and
cubic averaging as in TS1. Torsion is harmless
after tensoring the specialized Mordell–Weil
groups with \(\mathbb Q\).
The rational conductor count is CT5.
For conductors measured by the norm of
the conductor ideal over \(K,L,M\), respectively,
the same subclass gives
\(\gg Y^{1/10}/(\log Y)^4\),
\(\gg Y^{1/15}/(\log Y)^4\), and
\(\gg Y^{1/30}/(\log Y)^4\) distinct norms
at most \(Y\). At all primes away from \(2,3\)
these fixed fields are unramified, and
the conductor factor from CT5 is raised
to the field degree when taking the norm.
The factors at \(2,3\) have bounded
valuations by the fixed valuations of the
model discriminant over \(\mathbb Q\).
The large-prime part still determines the
rational conductor up to bounded multiplicity.

The generic rank table gives upper bounds
on the number of independent sections of
this surface. The ranks of its arithmetic
fibres can be larger.

## The ranks of all fixed constant sextic twists

**Corollary TS4.** For \(c\in\mathbb Q^\times\), the curve
\[
        \mathcal E_c:\quad y^2=x^3-c(3t^4-6t^2-1)/4
\]
has rank over \(\mathbb Q(t)\) given by
\[
 \begin{array}{c|c}
 c \text{ modulo }\mathbb Q^{\times6} & \operatorname{rank}\\ \hline
 1,\ -27   &2\\
 4,\ -108  &1\\
 \text{all other classes}&0 .
 \end{array}
\]
For \(c=4\), the point \((t^2-1,t(t^2-3))\)
is an integral generator. It is the point of
the third coefficient colour of the coherent
prime triangles.

Fix a sixth root \(u\) of \(c\). The isomorphism
over \(\overline{\mathbb Q}\) to \(\mathcal E\)
is \((x,y)\mapsto(u^{-2}x,u^{-3}y)\).
On the rationalized geometric group with
\(\operatorname{Gal}(\overline{\mathbb Q}/K)\) action, its
descent condition is that a point transform
under the character of \(u\) via the
order-six automorphism group of \(\mathcal E\).
By TS1 the only available characters over \(K\)
are the trivial one with multiplicity two and
the character from \(Q\) with multiplicity one.
The latter is the character of \(c=4\):
its sixth root can be \(2^{1/3}\), and the
displayed point maps to \(Q\).
The automorphism group acts faithfully on a
\(\mathbb Q(\zeta)\)-line, so equality of the
characters is equality of the corresponding
sextic Kummer classes in \(K^\times/K^{\times6}\).
All other classes have rank zero already over
\(K(t)\).

The kernel of the restriction from rational
sextic classes to \(K\) consists of \(1,-27\).
Indeed if \(u\in K^\times\) and \(u^6\) is rational,
then \(\bar u/u\in\mu_6\). Multiplying by a
sixth root of unity can make \(u\) rational
when this ratio is in \(\mu_3\).
Otherwise first multiply by \(\sqrt{-3}\)
to bring the ratio into \(\mu_3\).
Taking sixth powers gives precisely the
two stated classes. Thus the possible
rational classes are \(1,-27,4,-108\).
The parameters \(c\) and \(-27c\) give
rationally \(3\)-isogenous curves, hence
have the same rank. This is the standard
\(j=0\) isogeny, written for the two
Mordell models in
[Levin–Yan–Wiljanen, equation (2)](https://arxiv.org/pdf/1910.12276).
The point for \(c=4\) is independent and
has height \(2\), as its image \(Q\) does.
The type IV component pair is again
nonreal, so every rational section has
even height. This proves saturation in
rank one. With \(t=p/h\), the same point
rescales to \((p_1p_3,p(p^2-3h^2))\)
on \(y^2=x^3-E\), using the third norm
trace and coefficient \(2h\) of CT1.
