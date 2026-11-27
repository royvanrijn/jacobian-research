# Primitive pencils for a repeated line in top rank two

Let \(F=(P,Q,R)\) be a complex three-variable Keller map of ordinary
degree five with independent \(P_5,Q_5\) and \(\deg R=4\).
Suppose \(C=R_4\) has a repeated line. Use the projective pencil
of the [rank-two valuation argument](QUINTIC_KELLER_RANK_TWO_QUARTIC_PENCIL.md#valuations-along-a-projective-pencil),
but take a closed generator \(p/q\) of the relative algebraic
closure of \(\mathbb C(P_5/Q_5)\) in the projective plane function
field. Here \(p,q\) are coprime homogeneous forms of degree \(m\).
Existence and the generic irreducibility of the pencil follow from
[Petravchuk–Iena, Lemma 2 and Theorem 2](sources/petravchuk-iena-math-0701588v2.pdf).

**Theorem QuinticRankTwoRepeatedLinePrimitivePencil.** If \(F\) is
noninvertible, then \(m\le3\). The following necessary restrictions
hold, up to an affine linear frame and rescaling.

- For \(m=1\), all three leading forms are binary in the same
  directions.
- For \(m=2\), the pencil contains a double line \(z^2\).
  The quintic span is \(z\) times a two-dimensional space of
  binary quadratics evaluated on \(z^2,H\), for another member
  \(H\). The leading quartic is \(z^4\), \(z^2J\), or \(J^2\),
  where \(J\ne z^2\) is a conic member. In the repeated-line scope
  the last case requires \(J\) to be a product of distinct lines.
- For \(m=3\), the quintic span is \(\langle z^5,z^2H\rangle\)
  for a cubic member \(H\). The leading quartic is \(z^4\) or
  \(zJ\), where \(J\ne z^3\) is another cubic member. In the
  repeated-line scope \(J\) in the latter case is a product
  \(L^2M\) of distinct lines.

## Member valuations

Write \(A=P_5,\ B=Q_5\) and
\(A=g A_0(p,q),\ B=g B_0(p,q)\), where \(A_0,B_0\) are coprime
binary forms of common degree \(e\), and \(\deg g+me=5\).
The determinant makes \(C^5/B^4\) algebraic over
\(\mathbb C(A/B)\), so it belongs to \(\mathbb C(p/q)\).
For every pencil member \(M\) there is a rational constant
\(\rho_M\) such that
\[
                     5\nu_f(C)-4\nu_f(g)=\rho_M\nu_f(M)
                 \qquad(f\mid M).
\tag{1}
\]
Every factor of \(g\) and \(C\) occurs in a member because
\(\deg g\le4\). Distinct members are coprime. Outside the
members that contain \(g\), a member contributes an integral
multiple of its primitive power root to \(C\). Indeed the ratios
\(\nu_f(C)/\nu_f(M)\) are constant there and their common
denominator divides the greatest common divisor of the
multiplicities in \(M\).

For \(m=5\) the common factor \(g\) is constant. Every member
contributing to \(C\) has to be a fifth power of a line, since
otherwise it contributes at least five degrees. Two pure
power line members would make the pencil binary, contradicting
primitivity. Thus \(C=z^4\), and the quintic span is
\(\langle z^5,H_5\rangle\). The polynomial \(H_5(x,y,1)\) has
degree five and generates the closed pencil field. We exclude
this case below.

For \(m=4\) we have \(e=1\) and \(g=z\) linear. Its member is
\(z^4\), since any other factor \(f\) of that member would
give, by (1), \(5\mid\nu_f(M)\). With \(z^4\) in the pencil,
no other member is a full square, since the ratio of two
squares in that case would be composite. Thus every other
contribution to \(C\) has degree divisible by four. This
leaves \(C=z^4\) or \(C\) a single member disjoint from \(z\).
In the latter case the span is \(\langle z^5,zC\rangle\)
and \(C(x,y,1)\) generates the closed pencil field.
The
[special quintic span exclusion](QUINTIC_KELLER_RANK_TWO_QUARTIC_PENCIL.md#the-special-quintic-span)
applies, using only that \(C\) is not divisible by \(z\)
and the dehomogenization is non-composite. For \(C=z^4\),
the quintic span is \(\langle z^5,zH_4\rangle\), with
\(H_4(x,y,1)\) of degree four and generating the closed
pencil field. This remaining option is excluded below.

For \(m=3\) the common factor has degree two. An irreducible
conic factor would be accompanied by a line in its member,
which contradicts (1) modulo five. If \(g\) consists of
two distinct lines in one member, any extra distinct factor
is ruled out in the same way. The member would then be
\(z^2w\) or \(zw^2\). The two numerators in (1) are congruent
to one modulo five at the two \(g\)-factors, which is
incompatible with multiplicities two and one. If they are
in distinct members, both members must be pure cubes, again
contrary to primitivity. Thus \(g=z^2\). Its member is \(z^3\)
by congruence modulo five, and every other member has
multiplicity greatest common divisor one. The possibilities
for \(C\) and the span follow. A second pure cube is
impossible, so a cubic \(J\) with a repeated factor can
only be \(L^2M\) with distinct lines \(L,M\).

For \(m=2,e=2\) we have \(g=z\); its member is \(z^2\) by
congruence modulo five. For \(m=2,e=1\) the degree of \(g\)
is three. If it has an irreducible conic factor \(H\), this
factor is its own member, and the remaining linear factor
\(z\) must lie on a double-line member by (1).
If \(g\) is a product of lines, a reduced pair of lines in
one member must have equal multiplicities in \(g\), again
by (1) modulo five. A \(g\)-line with no companion \(g\)-line
in its member must lie on a double-line member. There is
at most one such member. Thus \(g=z^3\) or \(g=zH\)
with \(H\) a reduced conic member, allowing a pair of
distinct lines. These forms all give the stated quintic
span. Outside \(z^2\), a contribution to \(C\) is a full
conic member counted with its multiplicity, since (1)
makes the \(C\)-orders equal on a reduced line pair.
Thus the \(z\)-order of \(C\) is even. Two distinct conic
members disjoint from \(z\) would give a reduced quartic,
contrary to the repeated-line hypothesis. This gives the
list for \(m=2\). When \(m=1\) all factors
are binary and the leading forms are binary as asserted.

## Closure of the two remaining large primitive degrees

The cases \(m=4,5\) still to exclude have
\[
 P_5=z^5,\quad R_4=z^4,\quad
 Q_5=z^{5-d}H_d,\quad d=4\text{ or }5,
\]
where \(r=H_d(X,Y,1)\) has degree \(d\) and generates the
closed pencil field. Use the formal slice of the quartic
component, as in the
[rank-one quartic pencil closure](QUINTIC_KELLER_RANK_ONE_QUARTIC_PENCIL_CLOSURE.md).
Set \(x=\tau X,y=\tau Y,z=\tau Z\), with
\(\epsilon=\tau^{-1}\), and choose \(Z=1+O(\epsilon)\) from
\(R=\tau^4\). Write the resulting series for
\(\tau^{-5}P,\tau^{-5}Q\) as \(\Pi,\mathcal Q\).
The determinant gives \(\{\Pi,\mathcal Q\}_{X,Y}=O(\epsilon^{11})\),
since the unscaled bracket on the slice is
\(\pm\kappa\tau^2/R_z=O(\tau^{-1})\).
Write \(r_j=R_j(X,Y,1)\), \(p_j=P_j(X,Y,1)\),
\(q_j=Q_j(X,Y,1)\).

The coefficient of \(\epsilon\) in \(\Pi\) is
\(p_4-\tfrac54 r_3\), which commutes with \(r\).
It is constant if \(d=5\), and is \(\lambda r+c_0\) if \(d=4\).
Use \(\lambda=0\) in the first case. Subtract the commuting
formal expression \(\epsilon\lambda\mathcal Q\).
The coefficient of \(\epsilon^2\) is now
\[
 p_3-\tfrac54 r_2-c_0r_3-\tfrac5{32}r_3^2
                    +\lambda\bigl(r_3 r/4-q_4\bigr).
\tag{2}
\]
For the term involving \(\lambda\), the top homogeneous
component contributing to \(p_4\) is \(\lambda H_4\), whereas
the top of \(Q\) is \(zH_4\), giving
\(\epsilon\lambda(1-Z)H_4(X,Y,Z)\) before the lower terms.
The other terms of (2) are the same square expansion as
for the rank-one pure powers. Expression (2) commutes with
\(r\), so it is a polynomial in \(r\).

For \(d=5\), (2) excludes binary degrees three and two for
\(r_3\), since the square would have degree six or four,
respectively. For \(d=4\), binary degree three would give
degree seven if \(\lambda\ne0\) and degree six otherwise,
both impossible. If \(r_3\) has degree two or one, its
product with \(r\) has degree six or five. Thus
\(\lambda=0\) for either of these degrees.

Suppose \(d=4\) and \(r_3\) has degree two. Take the binary
terms of largest possible degree \(2i\) at order \(\epsilon^i\)
in \(\Pi\) and \(Z\), as in the rank-one slice. If \(L\) is
the leading quadratic of \(r_3\) and \(A=\epsilon L\), these
terms satisfy \(Z^4+AZ=1\), and \(\Pi\) is \(Z(1+AZ/4)\).
Its \(\epsilon^2\) term of binary degree four is
\(-5L^2/32\) and must be a multiple of the leading term of \(r\)
by (2). Subtract the corresponding multiple of
\(\epsilon^2\mathcal Q\) to remove it. The highest possible
binary terms of \(\mathcal Q\) at order \(\epsilon^i\) have
degree \(4+2i\) and are those of the leading binary term of
\(r\) multiplied by \(Z\). All lower homogeneous components
are below this binary-degree bound. The \(A^3\) coefficient
of \(Z(1+AZ/4)\) is zero, while the subtraction of
\(-5A^2Z/32\) contributes \(-5A^3/128\).
Thus the next coefficient has nonzero degree six and must
commute with \(r\), which is impossible. Hence \(\deg r_3\le1\)
also when \(d=4\).

For either \(d\), if \(r_3\) has degree one we now have
\(\lambda=0\). Expression (2) is constant and determines
the quadratic binary part of \(p_3\). Write \(L\) for the linear
part of \(r_3\) and \(K\) for the quadratic part of \(r_2\).
Taking the terms of binary degree \(i\) at order
\(\epsilon^i\), put \(A=\epsilon L,\ B=\epsilon^2K,\ W=Z^2\).
Exactly as in the rank-one formal slice,
\[
 W^2+AW+B=1,\qquad
 \Pi=Z(1+AW/4+B/4+5A^2/32).
\]
The cubic term \((5/128)L(L^2-8K)\) must vanish, because
all preceding coefficients are constant and binary degree
three is below \(d\). For \(L\ne0\) this gives \(K=L^2/8\).
The resulting \(A^4\)-coefficient is \(-5/1024\), which is
already incompatible with \(d=5\). If \(d=4\), subtract
the multiple of \(\epsilon^4\mathcal Q\) canceling the
degree-four term. The top binary terms of \(\mathcal Q\)
under these bounds are again the leading term of \(r\)
multiplied by \(Z\). Thus the \(A^5\)-coefficient
\(-1/1024\) changes by \(-5/4096\) upon subtraction, using
\(Z=1-A/4+O(A^2)\). It remains nonzero, giving a forbidden
binary degree five. This excludes a nonconstant linear
part of \(r_3\).

Translate \(z\) to remove the resulting pure \(z^3\) term of
\(R_3\), so \(r_3=0\). Suppose \(r_2\) has nonzero quadratic
binary part \(K\). For \(d=4\), (2) may still have the
\(-\lambda q_4\) term. Now \(Z=1-\epsilon^2 r_2/4+O(\epsilon^3)\).
After subtracting \(\epsilon\lambda\mathcal Q\), the
\(\epsilon^3\)-coefficient has leading degree-six term
\((\lambda/4)K r_{\mathrm{top}}\), where \(r_{\mathrm{top}}\) is
the degree-four part of \(r\). The other terms have degree
at most four, including the order-three contribution of
\(-\epsilon\lambda\mathcal Q\) from its lower components.
Removing the order-two coefficient by a multiple of
\(\epsilon^2\mathcal Q\) changes the order-three coefficient
by a multiple of \(q_4\) only, since \(Z_1=0\).
Consequently degree six is forbidden and \(\lambda=0\).

We can therefore use the same binary terms of degree \(i\)
at order \(\epsilon^i\) in both cases, now with \(A=0\) and
\(B=\epsilon^2K\). They give \(\Pi=(1-B)^{1/4}(1+B/4)\).
Its \(B^2\)-coefficient \(-5/32\) at order four excludes
\(d=5\). For \(d=4\), after removing the order-four
coefficient by a multiple of \(\epsilon^4\mathcal Q\), the
order-five coefficient has degree at most four, since
the top \(Z\)-coefficient of order one vanishes. It can be
removed by a multiple of \(\epsilon^5\mathcal Q\), with no
effect on degree six at the next order. At order six, the
\(B^3\)-coefficient \(-5/64\) changes by \(-5/128\) under
the order-four subtraction, using \(Z=1-B/4+O(B^2)\).
It remains nonzero and gives a forbidden degree six.
Thus \(K=0\).

We have proved that \(R_3\) is univariate in \(z\) and \(R_2\)
has no quadratic term in \(x,y\). Consequently
\(R=f(z)+a(z)x+b(z)y\) with \(f\) monic quartic and \(a,b\)
affine. The absence of critical points forbids a common
root of \(a,b\), since a nonzero slope lets one solve
\(R_z=0\) at that root and otherwise \(R\) is univariate.
A constant Bézout combination of \(a,b\) equal to one now
gives a coordinate completion of inverse degree at most
four as in the
[rank-one proof](QUINTIC_KELLER_RANK_ONE_QUARTIC_PENCIL_CLOSURE.md).
The [bounded coordinate-completion theorem](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
proves inversion.
All arguments transfer to characteristic zero by embedding
a finitely generated coefficient field, with pencil and
factorization statements made over the algebraic closure.

[The exact series checker](verify_quintic_rank_two_repeated_line.py)
replays the three coefficient corrections in the degree-four
pencil case and the uncorrected coefficients used for degree five.
