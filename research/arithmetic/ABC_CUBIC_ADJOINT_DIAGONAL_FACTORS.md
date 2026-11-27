# Diagonal Tamagawa factors of the cubic Frey adjoint

Let \(a+b=c\) be a primitive positive triple, let
\(E:y^2=x(x-a)(x+b)\), and put \(h=\log c\) and
\(A=\log\operatorname{rad}(abc)\). Use either primitive cubic
character \(\eta\) modulo \(7\) and \(K=\mathbb Q(\zeta_3)\).
For a coefficient place \(\lambda\mid p\), write \(O=O_{K_\lambda}\),
\(v_\lambda(\mathfrak m_\lambda)=1\), and \(Q_\lambda=|O/\mathfrak
m_\lambda|\). Set
\[
 T_s=\operatorname{ad}^0(T_p E\otimes O)\otimes\eta(s),\qquad
 V_s=T_s\otimes K_\lambda,\quad s=0,1.
\]
The trace pairing identifies the adjoint lattice and its dual at
odd \(p\). Use the de Rham determinant lattice given by a Néron
differential and its cup-dual complement on \(H^1_{\rm dR}(E)\),
and the integral Artin character line for the twist.
Changing the complement by a multiple of the differential
preserves the determinant lattice. The same lattices give the
critical period normalization outside \(2,7\), as checked below.

Write \(d_{p,\lambda}(T_s)\) for the norm of the
Fontaine–Perrin-Riou local Tamagawa ideal at \(p\) relative to
this determinant lattice. It is an integral power of
\(Q_\lambda\), possibly with negative exponent. When
the Bloch–Kato exponential identifies \(t_V=D_{\rm dR}(V)/
\mathrm{Fil}^0\) with \(H_f^1(\mathbb Q_p,V)\), this is
\[
 d_{p,\lambda}(T)=
 \mu(H_f^1(\mathbb Q_p,T))\,
 \big|\det(1-\varphi\mid D_{\rm cris}(V))\big|_\lambda.
\tag{DT1}
\]
Here \(|\pi_\lambda|_\lambda=Q_\lambda^{-1}\),
the exponential image of the tangent determinant lattice has
measure one, and the volume includes the order of the torsion.
For other crystalline representations use the determinant of
the full Bloch–Kato finite sequence.
[Diamond–Flach–Guo, equation (66)](https://arxiv.org/abs/2512.02348).
[Dummigan, §4 and §7](https://intlpress.com/site/pub/files/_fulltext/journals/mrl/2003/0010/0006/MRL-2003-0010-0006-a003.pdf).

For \(p\ge5,\ p\ne7\), the diagonal factor is one at every
good reduction prime. It is also one at \(p=3\) for good
Frey reduction. At a multiplicative prime \(p\ge5,\ p\ne7\), put
\(n=v_p(q_E)=-v_p(j(E))=2v_p(abc)\). For both values of \(s\),
\[
 d_{p,\lambda}(T_s)
   = Q_\lambda^{\,\mathbf 1_{\eta(p)=1}v_\lambda(n)
                    -s\,v_\lambda(p)} .
\tag{DT2}
\]
The multiplicative computation works for any finite unramified
character of order prime to \(p\), over any finite coefficient
extension containing its values. The diagonal variable-prime
product for the critical Tate twist therefore has logarithm
\(-2A+O(h/\log(h+2))\). There is no ordinarity requirement
at the good reduction primes.

## Good reduction

Let \(M\) be the integral filtered crystalline realization of
\(H^1(E)\) over \(\mathbb Z_p\), extended to \(O\). Integral
comparison for the smooth proper model identifies its filtration
with the de Rham lattice and its Hodge line with the Néron
differential line. For \(\operatorname{Hom}(M,M(s)\otimes
D_{\rm cris}(\eta))\), the two filtered modules before taking
Hom have all jumps in the common interval \([-1,1]\). Its
length two is at most \(p-2\). The rank-one unramified character
lattice is strongly divisible, with unit Frobenius, so tensoring
with it preserves the Fontaine–Laffaille category in this
range. The local Hom argument in the proof of
Diamond–Flach–Guo's Proposition 8.12 identifies the integral
Ext group with \(H_f^1\) of the Hom lattice, and the resulting
Tamagawa ideal is a unit for the determinant lattice from the
filtered Hom. Use the argument for the two local filtered
modules just specified. It has no global modularity or residual
image hypothesis. Since \(2\) is a unit, the Hom decomposes
by trace into the adjoint and the scalar line, compatibly with
the integral determinant sequence. The scalar line is covered
by the same Hom argument with rank-one filtered modules.
Consequently the factor for the trace-zero summand is one.
[Diamond–Flach–Guo, proof of Proposition 8.12, equations (60)–(66)](https://arxiv.org/abs/2512.02348).

For \(p=3\) and good Frey reduction the same unit conclusion
holds. Full rational two-torsion and the Hasse bound give
\(\#E(\mathbb F_3)=4\), hence \(a_3=0\). The filtered module
\(M\) has no nonzero quotient supported in the top Hodge
filtration, by the nilpotent Frobenius calculation in
Diamond–Flach–Guo's Lemma 8.13. Explicitly, the characteristic
polynomial of crystalline Frobenius reduces to \(X^2\), and
strong divisibility makes the reduced Frobenius nonzero
with image equal to the Hodge line. Divided Frobenius sends
that line outside itself. No nonzero proper subobject survives
both operators. An unramified twist multiplies these maps by
a unit and preserves the argument. This is the borderline
condition used immediately after equation (66) for the
shift \(s=1\) when \(p=k+1=3\). For \(s=0\) the shorter
interval already suffices. The trace line again has the
shorter interval. These are local arguments with the
crystalline lattice of \(E\), also when \(O\) is ramified
over \(\mathbb Z_3\).
[Diamond–Flach–Guo, Lemma 8.13](https://arxiv.org/abs/2512.02348).

## Integral finite cohomology at multiplicative reduction

A quadratic twist makes the Tate curve split and disappears
in the adjoint. At an odd multiplicative prime the quadratic
twist can be unramified and changes the tangent lattice at
most by a unit after unramified comparison. Use a Tate basis
in which \(g\) acts on \(T_p E\) by
\(\begin{pmatrix}\varepsilon(g)&u(g)\\0&1\end{pmatrix}\),
where \(\varepsilon\) is cyclotomic and \(u\) represents the
Kummer class of \(q_E\). Write \(e=e_{12}\),
\(d=\mathrm{diag}(1,-1)\), \(f=e_{21}\).
Put \(B_0=Oe\otimes\eta\) inside \(T_0\), and put
\(B_1=(Oe+Od)\otimes\eta(1)\) inside \(T_1\). The resulting
pieces are
\[
 \begin{array}{c|c|c}
 s & B_s & T_s/B_s\\ \hline
 0 & O(\eta)(1)&O(\eta)\ \text{extended by the lower-weight
 quotient }O(\eta)(-1)\\
 1 & 0\to O(\eta)(2)\to B_1\to O(\eta)(1)\to0 & O(\eta).
 \end{array}
\]
In the first row \(O(\eta)\) is the submodule of the quotient,
with \(O(\eta)(-1)\) as its quotient.

Each rational \(B_s\) maps isomorphically onto \(t_{V_s}\)
on de Rham tangent spaces. For the Néron tangent lattices,
use the Tate comparison basis of Dummigan's Propositions
6.5–6.10. In \(H^1_{\rm dR}(E)\) its \(e_0\) generates the
quotient by the Hodge line and \(e_1-\log_p(q_E)e_0\) is a
Hodge generator. These two elements generate the Néron
de Rham lattice since \(\log_p(q_E)\in p\mathbb Z_p\).
The torus subrepresentation corresponds to the first line
under the Tate self-duality. In the symmetric square
description of the trace-zero adjoint, the subspace for
\(B_0\) is the span of \(e_0^2\), and for \(B_1\) it is
the span of \(e_0^2,e_0e_1\), with their respective Tate twists.
Their images are bases of the tangent determinants. For
\(B_1\) the two rank-one graded lines give the same determinant,
since replacing \(e_1\) by \(e_1-\log_p(q_E)e_0\) is
triangular of determinant one. The symmetric-square to
adjoint identification only uses division by the unit \(2\).
An unramified finite character adds the same rank-one integral
line to these identifications.
[Dummigan, Lemma 6.9 and Proposition 6.10](https://intlpress.com/site/pub/files/_fulltext/journals/mrl/2003/0010/0006/MRL-2003-0010-0006-a003.pdf).

Let \(c_{p,\lambda}=v_\lambda(n)\) when \(\eta(p)=1\), and zero
otherwise. The map \(H_f^1(\mathbb Q_p,B_s)\to
H_f^1(\mathbb Q_p,T_s)\) is injective with cokernel of order
\(Q_\lambda^{c_{p,\lambda}}\). To see this, first
\(H^1(\mathbb Q_p,T_s/B_s)\) is torsion-free. For \(s=1\)
the quotient is the unramified line \(O(\eta)\), which is
either trivial or has Frobenius minus one a unit. For \(s=0\)
the last graded character \(\eta\varepsilon^{-1}\) has no
residual inertia invariants, so the invariants of the discrete
quotient by the lattice are those from the line \(O(\eta)\).
They lift to rational invariants when \(\eta=1\) and vanish
otherwise. The coefficient sequence gives the asserted
torsion-freeness.

At \(s=0\) the Kummer class of the boundary from the diagonal
line, when \(\eta=1\), is \(-2[q_E]\) in
\(H^1(\mathbb Q_p,O(1))\). At \(s=1\) the boundary from the
trivial quotient projects to \([q_E]\) in
\(H^1(\mathbb Q_p,O(1))\). These formulas follow by conjugating
\(d\) and \(f\) with the Tate matrices. In the second row
the map \(H^1(B_1)\to H^1(O(\eta)(1))\) is onto because
\(H^2(O(\eta)(2))=0\), by local Tate duality and the nontrivial
residual cyclotomic character in the dual. Its restriction
to the finite parts is onto, with kernel \(H_f^1(O(\eta)(2))\).
Indeed over the fraction field the preimage of the finite
part has dimension two, the dimension of \(H_f^1(B_1)\),
and the integral conditions are defined by preimage.

If \(\eta=1\), the quotient of \(H^1(B_s)\) by its finite
part is identified with \(O\) by the Kummer valuation,
directly for \(s=0\) and through the last graded line of
\(B_1\) for \(s=1\). The boundary has valuation
\(-2n\) or \(n\). If \(\eta\ne1\), these quotient spaces
are zero by local Euler characteristic on the rank-one
graded lines. Over the fraction field the finite part of
\(B_s\) maps isomorphically to that of \(V_s\) by naturality
of the exponential and the tangent isomorphism. An integral
finite class of \(T_s\) maps to a torsion class in
\(H^1(T_s/B_s)\), so it comes from \(H^1(B_s)\).
Modulo the image of the finite part of \(B_s\), the only
remaining quotient is the valuation line modulo the image
of the boundary, as asserted. This calculation includes
torsion created in the finite part by dividing out the
Kummer boundary.

## The measure and its Frobenius factor

On the rank-one line \(O(\eta)(r)\), \(r=1,2\), the factor
for the integral unramified tangent lattice is a unit, by
the local Fontaine–Laffaille Hom calculation above (use
the trivial source and the rank-one target). Its crystalline
eigenvalue is a unit times \(p^{-r}\), so in the normalization
of (DT1) the finite-cohomology volume is
\(Q_\lambda^{r v_\lambda(p)}\).
For \(B_0\) this gives volume \(Q_\lambda^{v_\lambda(p)}\).
For \(B_1\) the exact sequence of integral finite parts and
the tangent determinant product give volume
\(Q_\lambda^{3v_\lambda(p)}\).
The index calculated above multiplies these volumes by
\(Q_\lambda^{c_{p,\lambda}}\) on passing to \(T_s\).
The crystalline space of \(V_s\) is the monodromy kernel, a
line with Frobenius a unit times \(p^{-1-s}\). Thus (DT1)
divides the volume by \(Q_\lambda^{(1+s)v_\lambda(p)}\).
This proves (DT2), including the factor for a nontrivial
unramified cubic character. For the untwisted case the
result agrees with Dummigan's two diagonal calculations
for the symmetric square.
[Dummigan, Proposition 7.4](https://drive.google.com/file/d/1PdkMLtp31BZGzmLUk_jiJKpigRwtvvDJ/view).

At \(p=7\) with multiplicative reduction, the diagonal
factors are bounded above and below by fixed positive
constants. Here \(\eta\) is ramified, and for \(\lambda\mid7\)
its reduction on inertia is a primitive cubic character.
The two characters in \(T_0/B_0\) have tame exponents
\(\pm2,\pm2-1\) modulo six, and \(T_1/B_1\) has exponent
\(\pm2\). None is trivial, so the quotients have no rational
invariants or torsion in \(H^1\). Also \(H^2(O(\eta)(2))=0\),
since the dual modulo \(\lambda\) has exponent
\(\mp2-1\), which is nonzero. The same integral finite-part
argument gives index one for the two maps from \(B_s\).
The three de Rham tangent identifications used above
tensor with the fixed Artin line. The rank-one factors for
\(O(\eta)(1)\), \(O(\eta)(2)\) at \(7\) are fixed. Their
crystalline spaces, as well as those for \(V_s\), are zero,
since the finite part of inertia acts by \(\eta\) on the
potentially semistable realization. Thus the two diagonal
factors are fixed products of rank-one volumes in this
normalization. Changing the Artin basis only changes them
by constants independent of \(E\).

## The determinant period for the critical value

For the rational curve, let \(u,v\) be the real and imaginary
periods on the rectangular Frey period lattice, so \(|uv|=P(E)\).
Take the de Rham basis \(\omega,\theta\) with \(\omega\) Néron and
\(\theta\) cup-dual to \(\omega\). The comparison rows for
\(H^1\) are \((u,v),(u',v')\) with
\(\Delta=uv'-vu'=\pm2\pi i\) by the cup product. On
\(\operatorname{Sym}^2 H^1(E)(2)\), the plus Betti space is
generated by the squares of the two Betti eigenvectors, and
the tangent quotient has basis \(\omega\theta,\theta^2\).
Up to the sign convention for comparison, the projected
columns are
\[
 \frac{(2\pi i)^2}{\Delta^2}
    \binom{-2v'v}{v^2},\qquad
 \frac{(2\pi i)^2}{\Delta^2}
    \binom{-2u'u}{u^2}.
\]
Their determinant is a fixed unit times \(2(2\pi i)uv\).
The adjoint Tate twist \(s=1\) is this motive, with the
lattice identification away from \(2\). For the Artin line
the comparison is a Gauss sum or its inverse, with the
choice of character convention fixed once. The Gauss sums
are units away from \(7\), since both are integral and
\(G(\eta)G(\bar\eta)=7\). In the cyclic cubic field their
a Gauss sum generates the integral character line away
from \(7\). In fact it is a unit in the extended integer
ring, and dividing any section in the same character line
by that unit gives a Galois invariant. These invariants
are exactly the coefficient ring, by flat extension of
the invariants of the integer ring of the cubic field. The even twist contributes the
square of the Artin comparison to the determinant.
Consequently outside \(2,7\) the normalized fixed cubic
critical value \(X_E\) has the valuations of the critical
value in this determinant normalization (with the chosen
character or its conjugate). On norms the character
convention has no effect.
[Fixed cubic normalization and Néron period](ABC_FIXED_CUBIC_ADJOINT_VALUE.md#a-value-in-one-fixed-quadratic-field).

## The remaining finite size estimate

Let \(\mathcal D(E)=\prod_{p\ge5,\ p\ne7}\prod_{\lambda\mid p}
d_{p,\lambda}(T_1)\). Since \(\sum_{\lambda\mid p}
[K_\lambda:\mathbb Q_p]=2\), equation (DT2) gives the exact
formula
\[
 \mathcal D(E)=
 \prod_{\substack{p\mid abc\\p\ge5,\ p\ne7}}
          p^{\,2(\mathbf 1_{\eta(p)=1}v_p(2v_p(abc))-1)} ,
 \qquad \tfrac12\log\mathcal D(E)=-A+
            O\!\left(\frac{h}{\log(h+2)}\right).
\tag{DT3}
\]
For the error, \(\sum_{p\mid abc}\log(2v_p(abc))\) is bounded by
twice the log of the divisor function of \(abc\), as in the
off-diagonal product proof. There are only the fixed missing
rational primes \(2,3,7\).

For a nonexceptional triple and \(C(E,\eta,1)\) the product
over \(p\ne\ell\), the positive rational number
\[
 Z_E=\frac{\mathrm N_{K/\mathbb Q}X_E}
               {C(E,\eta,1)\mathcal D(E)}
 \quad\text{satisfies}\quad
 \tfrac12\log Z_E=h+A+
                 O\!\left(\frac{h}{\log(h+2)}\right).
\tag{DT4}
\]
Therefore abc is equivalent to the uniform upper bound
\(\tfrac12\log Z_E\le(2+\varepsilon)A+O_\varepsilon(1)\) for
every \(\varepsilon>0\). The \(o(h)\) term is absorbed by
using a slightly smaller \(\varepsilon\) in either direction.
For example, at \(p=5\) the triple \((1,4,5)\) gives factor
\(5^{-2}\) in (DT3). At \(p=13\) with
\((1,13^{13}-1,13^{13})\) the character is trivial at \(13\),
so \(v_{13}(2v_{13}(abc))=1\); the correction cancels
the factor \(13^{-2}\) in that local product.

### Conditional Selmer formula outside the fixed primes

Put \(S=\{2,3,7\}\), and let \(V_1^D=V_1^\vee(1)\) and
\(T_1^D=T_1^\vee(1)\). If the integral special-value conjecture
holds for the twisted motive at all places \(\lambda\nmid S\),
then it gives
\[
 \sum_{p\notin S}v_p(Z_E)\log p
    =\sum_{\lambda\nmid S}
       \log\# H_f^1(\mathbb Q,V_1^D/T_1^D).
\tag{DT5}
\]
The denominator of the determinant lattice in the special-value
formula contains the Selmer Fitting ideal and the product of
the Tamagawa ideals. The fundamental-line generator is the period determinant
divided by the critical value. This gives
the positive sign of the Selmer size and the division by local
factors in (DT4). The global invariant-group terms are trivial
outside \(S\), as is the archimedean local term at odd
coefficient primes. The conjectural comparison also yields
finite support of the sum in (DT5); individual Selmer finiteness
is already unconditional.
[Diamond–Flach–Guo, equation (58) and Conjecture 8.9 for the determinant conventions](https://arxiv.org/abs/2512.02348).
[Off-diagonal factors, global torsion and Selmer rank](ABC_CUBIC_ADJOINT_TAMAGAWA_FACTORS.md).

To use (DT5) for abc still requires the integral twisted
special-value comparison, a uniform size estimate for the
finite Selmer groups, and a bound on the signed contributions
of \(Z_E\) at coefficients in \(S\). The unresolved local
diagonal factors within \(S\) are at \(2\), at \(3\) when
\(E\) is multiplicative, and at \(7\) when \(E\) is good.
If the signed contribution from \(S\) is \(o(h)\) uniformly,
then under (DT5) abc is equivalent to bounding half the
aggregate Selmer logarithm by \((2+\varepsilon)A+
O_\varepsilon(1)\) for every \(\varepsilon>0\).

The local factors can also be grouped by rational geometric
prime. At primes splitting in the cubic field, this grouping
recovers the entire Tate exponent that occurs in the classical
to Shimura modular-degree ratio. It gives an exact comparison
of the corrected critical norm with rational Hecke projector
denominators.
[Shimura projector denominators and the cubic local packages](ABC_CUBIC_ADJOINT_SHIMURA_PROJECTORS.md).
