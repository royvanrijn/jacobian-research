# Symmetric prime lines with varying cyclic outer fields

The line \(p_1+p_3=2p_2\) is one member of a family of
genuine coherent prime lines with two cubic companions.
For each fixed \(r\equiv1\bmod6\) the line
\(p_1+p_3=2r p_2\) has the same two discriminant laws.
The cyclic field governing the outer prime marks of
the wild companion varies with \(r\). One can ensure
arbitrarily many different such cyclic fields and
prescribe their outer Artin classes on each fixed line.

## Coherent norm curves and their two isogeny marks

**Proposition SL1.** For a positive integer
\(r\equiv1\bmod6\), use
\((p_1,p_2,p_3)=(rp-h,p,rp+h)\) and \(t=p/h\).
The normalized \((1,1,1)\) radicand is \(h^6F_r(t)\), with
\[
 F_r(t)=(1-4r^3)t^6+12r(3r^3-1)t^4+24r^2t^2+4.
 \tag{1}
\]
On \(C_r:y^2=F_r(t)\), the two consecutive cubic
norm functions are geometrically independent.
The two degree-two elliptic projections give the
quotient and product directions. The second elliptic
curve has two rational \(3\)-isogenies, one with
rational kernel and one with kernel over
\(\mathbb Q(\sqrt{-3})\).

Write \(a=1-4r^3,\ b=8r^3-1\).
The two coordinate pairs are
\[
 (X_1,Y_1)=(at^2-4r,ay),\qquad
 (X_2,Y_2)=(4(t^{-2}-r^2),4y/t^3).
\]
Substitution gives
\[
 Y_1^2=X_1^3+(6r^2X_1+2b)^2,\qquad
 Y_2^2=X_2^3+(6rX_2+4b)^2 .
 \tag{2}
\]
Let
\[
 \alpha_1=\tfrac12(t^3-6r^2t^2-2+y),\qquad
 \alpha_2=\tfrac12(t^3+6r^2t^2+2+y).
\]
Their norms are \((t(rt-1))^3,(t(rt+1))^3\).
The marked functions in (2) satisfy
\[
 \begin{split}
 (Y_1+6r^2X_1+2b)\bar\alpha_1\alpha_2
       &=(ry+(6r^3-1)t^2+2r)^3,\\
 (Y_2+6rX_2+4b)\alpha_1\alpha_2
       &=(-2(r^2t^2-1))^3 .
 \end{split}
\]
These are the full normalized analogues of PQ1.

The discriminant of \(F_r\) is
\(2^{16}3^6(4r^3-1)(8r^3-1)^6\), so it is
nonzero in the stated range. Each \(\alpha_i\) has a
triple zero on one branch at its outer marked fiber,
\(t=1/r\) for \(\alpha_1\) and \(t=-1/r\) for \(\alpha_2\),
and has pole order three at both infinities.
For a putative cubic relation choose exponents from
\(\{-1,0,1\}\), using the conjugate for a negative
exponent. A cube root would be regular on the affine
curve with poles of order at most two at infinity.
It would be a polynomial in \(t\), since any term
\(yQ(t)\) with \(Q\ne0\) has a pole of order at least
three on one of the infinities.
The one-sided zero at an outer marked fiber excludes
such a polynomial for either nonzero exponent.
This also proves independence at sufficiently large
characteristics.

Put
\[
 Q=4r^2+2r+1,\quad P=3r^2+2r+1,\quad
 S=6r^3+6r^2+4r+1 .
\]
On the second curve recenter at \(X_2=-4Q\), setting
\[
 Z=X_2+4Q=4(t^{-2}+P),\qquad
 L=2(r+1)Z-4Q=8(r+1)t^{-2}+4S .
\]
Then \(Y_2^2=Z^3-3L^2\).
The alternate rational subgroup is the one at
\(Z=0\), with ordinates \(\pm4Q\sqrt{-3}\).
The original rational subgroup is at \(X_2=0\), so
these are two distinct \(3\)-isogenies.

## Genuine primes and the cubic companions

**Theorem SL2.** For each fixed \(r\equiv1\bmod6\),
there are \(\gg_r X^{1/3}/(\log X)^{B_r}\)
distinct quadratic fields of either sign and absolute
fundamental discriminant at most \(X\) arising from
genuine coherent prime triangles on that line.
The two consecutive cubic classes are independent.
With \(D\) the base discriminant, the progression
has the two cubic companions of discriminants
\(-D/3\) and \(81D\), respectively split and totally
ramified at three. A polynomial for the second field is
\[
 T^3-36(h^2+Pp^2)T+144(r+1)ph^2+72Sp^3 .
 \tag{3}
\]

To verify genuineness, use \(A=(rp-h)^3,\ B=p^3,\
C=(rp+h)^3\) in the usual Heron expression, obtaining
(1). At the outer primes the opposite coefficients
give the fixed nonzero factor \(1-8r^3\) in the
incident trace. At the middle prime the two outer
powers are opposites. For primes exceeding the fixed
bad primes and with \(\gcd(p,h)=1\), the split
orientations and coherent prime-cube ideal identities
therefore hold exactly as in PP1.
The three linear forms have nonparallel linear
parts; after requiring \(p\) odd, \(p\equiv1\bmod3\),
and \(h\equiv0\bmod18\), they remain locally admissible
at \(2,3\). At every other prime simultaneous unit residues
are available. Away from
the divisors of \(r\) use \(p=1,h=0\); at an odd
prime dividing \(r\) use \(p=h=1\).
The two independent residue rows can be chosen by
Hasse–Weil on the Kummer twists of SL1, as in HP2,
at fixed sufficiently large primes \(1\bmod9\).
Then impose those rows by CRT. This keeps the
prime-form hypotheses of PP3.
The two real intervals occur on the positive parameter
line \(t>1/r\): \(F_r\) is positive at \(1/r\), and
negative for sufficiently large \(t\).
Choose compact subintervals with genuine distinct
vertices bounded away from equality.
The separable binary sextic and fixed \(2\)-torsion
splitting field give the PP2–PP3 fixed-field
multiplicity and logarithmic losses for this fixed
\(r\).

The first companion is the cubic reflection of
\(\bar\alpha_1\alpha_2\), and the three primes agree
modulo three in this slice. Thus RF3 gives its
fundamental discriminant and splitting.
For the second, use the recentered mark and set
\[
 \rho=\sqrt{-3}\,y,\qquad
 \xi=-L/3+4\rho/(9t^3).
\]
As in PC1 this is the alternate evaluation in the
mirror field and has norm \((Z/3)^3\).
At an odd prime outside \(3\) dividing \(Q\), all
three prime forms are units, and
\[
 Z/4\equiv (h/p)^2-r^2= -p_1p_3/p^2\pmod q .
\]
Thus the potentially bad constant at \(q\mid Q\)
cannot interfere with virtuality. At primes not
dividing \(2\cdot3\cdot Q\) the usual pole or unit
argument applies.
At \(2\), odd \(r\) gives \(P\equiv2\bmod4\) and
\(S\) odd. The local squareclass of the base field
is \(5\), hence the mirror splits at \(2\). As in
PC1 the two terms of \(\xi\) have valuation two
and give branch valuations \(3,6\).
At \(3\) the base field has squareclass \(-3\), since
\(1-4r^3\equiv-3\bmod9\). The mirror splits, and
\(\rho/t^3\) is \(3\) times a unit congruent to \(1\)
or \(-1\bmod3\). Also \(L\equiv2\bmod3\).
On one branch \(4\rho/t^3-3L\) has valuation one,
so \(\xi\) has valuation \(-1\); its other value
is \(1\bmod3\) by the cube norm. In particular the
marked class is noncubic.
The ramified reflection then has relative conductor
\(9\mathcal O_K\) by the Eisenstein different
calculation of PC2. This gives discriminant \(81D\)
and total ramification at three. Scaling the
cubic trace polynomial by \(3p\) gives (3).

## The varying cyclic outer fields

**Theorem SL3.** In (3) the middle prime always
splits. At either outer prime, away from fixed bad primes,
the splitting is governed by
\[
 G_r(U)=U^3-QU+(2r+1)Q/3,\qquad
 \operatorname{Disc}(G_r)=Q^2(2r-1)^2 .
 \tag{4}
\]
If a prime \(q>3\) divides \(Q\) to order one, \(G_r\)
is irreducible and defines a cyclic cubic \(M_r\) ramified
at \(q\). There are infinitely many fixed profiles
\(r\equiv1\bmod6\) with pairwise distinct \(M_r\).
For any fixed such profile one can independently
prescribe the two outer Artin classes in
\(\operatorname{Gal}(M_r/\mathbb Q)\), and the three
split or inert tests for the fundamental companion.
The lower counts of SL2 retain their exponent on
each fixed profile and prescription.

At an outer prime put \(h=\pm rp\) in (3)
and scale the root by \(6p\). Formula (4) follows;
the square discriminant comes from
\(4Q-3(2r+1)^2=(2r-1)^2\).
At the middle prime put \(p=0\) and scale by \(6h\),
giving \(U^3-U\).
The polynomial \(G_r\) is equivalent to the
[Shanks polynomial](https://doi.org/10.2307/2005372)
\(X^3-aX^2-(a+3)X-1\) for the rational parameter
\(a=6(r+1)/(2r-1)\): depress the Shanks cubic and
scale by \(-(2r-1)/3\).
At \(q>3\) dividing \(Q\) exactly once, \(2r+1\)
is a unit, so (4) is Eisenstein.
For every \(q\equiv1\bmod3,\ q>3\), one may choose
\(2r\bmod q\) as a primitive cube root of unity.
The two roots of \(Q\) are simple there, so a
lift modulo \(q^2\) can force the valuation one.
CRT then gives \(r\equiv1\bmod6\).
Choose successively a \(q\) outside the support
of the earlier field discriminants. The new outer
field ramifies at that \(q\) and is distinct.

The conductor of \(M_r\) has support dividing
\(3Q(2r-1)\), and is coprime to \(2r\) on this slice.
It has exponent one at its ramified primes other
than \(3\), and exponent at most two at \(3\),
since its Dirichlet character has order three.
At the selected ramified \(q>3\), choose the
two outer residues giving their prescribed Artin
classes, fixing the local residue \(p=1,h=0\) at
the other conductor primes. The cubic local
character at \(q\) is surjective. The residues
can be negated independently without changing
their classes, so their sum can be made a unit
at \(q\). Then \(p,h\) are recovered by taking
\(p\) equal to the sum divided by \(2r\) and
\(h\) equal to half the difference.
These congruences are compatible with the
admissibility and auxiliary tests of SL2.

For the fundamental companion, the normalized
outer character is that of \(1-8r^3\), from the
outer incident trace test of PS1, and the middle
test is that of \(2\), up to an inverse.
Both are noncube rational constants, the first
because \(8r^3-1\) lies strictly between consecutive
cubes. Prescribe their trivial or nontrivial cubic
classes in the pure cubic splitting fields over
\(\mathbb Q(\zeta_3)\).
The three prime forms are \(1\bmod3\) on our CRT
slice, so the restrictions to each maximal abelian
subfield are already met.
Frei–Tinková's theorem as in PS2 applies jointly
with the conductor congruences for \(M_r\).
Its positive local factors and the fixed-sextic
multiplicity of SL2 complete the proof.
No uniform lower-bound constant in \(r\) is used.

For example \(r=7\) gives \(G_7(U)=U^3-211U+1055\)
of field discriminant \(211^2\), in place of the
conductor-seven field at \(r=1\).
With \(p=223,h=18\), the primes are
\(1543,223,1579\). The discriminants of the two
companions are \(33136136483622217\) and
\(-8052081165520198731\).
The three marks in the second field are inert,
split, split, as predicted by \(G_7\).
The [exact checker](verify_symmetric_prime_lines_cyclic_outer_fields.py)
replays the formulae over \(\mathbb Q(r,t)\)
and this genuine prime fixture.
