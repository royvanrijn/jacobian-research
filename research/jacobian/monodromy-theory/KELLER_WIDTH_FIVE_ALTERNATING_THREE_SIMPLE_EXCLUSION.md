# Exclusion of three simple critical points in the width-five elliptic slice

**Theorem SixSheetWidthFiveAlternatingThreeSimpleExclusion.**
In the six-sheet \(A_6\) profile with one unramified omitted
prime, a width-five coordinate of slice genus one cannot have
three simple critical points on the immersed locus of the sole
ramified branch.

Use the
[width-five marked elliptic restriction](KELLER_WIDTH_FIVE_ALTERNATING_ELLIPTIC_MARKS.md).
The relatively minimal model is \(Y^2=X^3+D(p)\) with \(D\)
a squarefree cubic. The rational end of the resolved special
fibre carries one simple pole section. Use the unblown triple
pole as the origin \(O\). Denote the three simple sections
\(T_i=(u_i,v_i)\), indexed so that \(T_3\) is the blown
section. All four sections on the relatively minimal model
are pairwise disjoint on the finite base, since only one
section passes through either contraction centre. In
particular all \(u_i,v_i\) are polynomials.

## Equations for the marked sections

For \(f_i=(Y+v_i)/(X-u_i)\) put \(F_{ij}=f_j(T_i)\) if \(i\ne j\).
These evaluations are polynomial on the finite \(p\)-line,
since \(f_j\) has only the poles at \(O,T_j\) and no vertical
pole. They are skew in \(i,j\). In fact they equal
\((v_i+v_j)/(u_i-u_j)\) if the denominator is nonzero, and
equal \(-3u_j^2/(2v_j)\) when \(T_i=-T_j\) generically.

With \(\omega=dX/(2Y)\), the marked factor gives
\(dq=\varphi^2\omega\). As in the preceding marked restriction,
we can write
\[
                      \varphi=X+\sum_i C_i f_i+L,
             \qquad (C_1,C_2,C_3)=(c_1,c_2,p)
\tag{1}
\]
after constant rescaling of \(\varphi,q\) and an affine
change of \(p\) taking the cusp projection to zero.
The coefficient of \(X\) is a nonzero constant because \(O\)
is unblown. The nonzero constants \(c_1,c_2\) come from the
unblown simple sections. The twice blown simple section
\(T_3\) gives a coefficient with exactly one simple finite
zero, which the base change normalizes to \(p\).

Set \(r=c_1+c_2\), \(t_0=c_1^2+c_1c_2+c_2^2\),
\(S=p+r\), \(U=\sum_i C_i u_i\), and \(V=\sum_i C_i v_i\).
The residue at each \(T_i\) gives
\[
                       L+u_i+\sum_{j\ne i} C_j F_{ij}=0.
\tag{2}
\]
The Laurent pole comparison at \(O\) gives
\[
 \begin{aligned}
 L&=-p^2-rp-t_0,\\
 U&=(p+r)(p^2+rp+t_0),\\
 0&=L^2+2V+2SU-\sum_i C_i^2 u_i .
 \end{aligned}
\tag{3}
\]
To see these equations, match the highest poles in
\(\delta q=\varphi^2\), where \(\delta=d/\omega\).
The pole space with divisor \(3O+T_1+T_2+T_3\)
has basis \(1,Y,X,f_1,f_2,f_3\).
The simple poles at \(T_i\) and the two highest terms
at \(O\) give \(q=Y/3+SX-\sum_i C_i^2 f_i\) up to a
function of \(p\). In a local parameter \(\tau\) at \(O\)
with differential \(\omega\), the terms needed are
\[
 X=\tau^{-2}+O(\tau^4),\qquad
 Y=-\tau^{-3}+O(\tau^3),\qquad
 f_i=-\tau^{-1}-u_i\tau+v_i\tau^2+O(\tau^3).
\]
The coefficients of \(\tau^{-2},\tau^{-1},\tau^0\)
give (3). In particular \(\delta f_i\) has constant
coefficient \(-u_i\). The
[exact checker](verify_keller_width_five_alternating_three_simple.py)
replays these comparisons and the section coefficients
used below.

By the integral-section calculation in the marked restriction,
the possible degrees of \(u_i\) are one, two and four.
If one has degree four, every \(F_{ij}\) involving that
section has growth at infinity of degree at most two.
This follows from the quotient formula for different
abscissas, using the explicit form of the degree-four
abscissas, and also from the regular value formula when
the abscissas coincide. Equation (2) is then impossible
for that index: each \(C_jF_{ij}\) has degree at most
three and \(L\) has degree two. Thus all degrees are
at most two.

The coefficient of \(p^3\) in the middle equation of (3)
gives \(u_3=p^2+a p+b\). Its ordinate has degree three
with leading coefficient \(1\) or \(-1\). In degree four
in the last equation of (3), the first, third and
fourth terms sum to \(2p^4\) at leading order. Hence
\(v_3=-p^3-\tfrac32a p^2+O(p)\).

## The possible unblown leading terms

Call an unblown section *linear* if its abscissa has
degree one; its ordinate then has degree at most one.
For an unblown quadratic section with leading terms
\(u_i=\alpha p^2+\cdots,\ v_i=\beta p^3+\cdots\)
we have \(\beta^2=\alpha^3\).
Suppose \(\alpha\ne1\). The \(p^2\) coefficient in
(2) for this index is
\[
                        -1+\alpha+\frac{\beta-1}{\alpha-1}.
\tag{4}
\]
Here the other unblown evaluation has degree at most
one unless its quadratic abscissa has the same leading
coefficient as \(u_i\). We first treat this latter
situation to justify that it does not change (4).

If the leading abscissas agree, write their difference
as \(d p+e\) when it is nonzero. When both leading
ordinates have the same sign \(\beta\), a constant
abscissa difference yields growth degree three in the
unblown evaluation, impossible in (2). A linear
difference yields an extra term \(2\beta/d\) in the
\(p^2\) coefficient for the first index and its negative
for the other. When their leading ordinates are
opposite, a linear difference produces no extra
\(p^2\) term; a constant nonzero difference also
produces none, because the \(p^2\) terms of the
ordinates then cancel. If the abscissas are identical,
the sections are opposite and their regular evaluation
has degree one.

In the remaining same-sign situation with \(d\ne0\),
write \(c_i,c_j\) for the two unblown residues. The
two leading equations are
\[
                E+2\beta c_j/d=0,\qquad
                E-2\beta c_i/d=0,\qquad
                E=-1+\alpha+\frac{\beta-1}{\alpha-1}.
\]
They give \(r=0\). The \(p^2\) coefficient of \(U\)
in (3) then gives \(a=0\). Subtract the two
residue equations. Since \(c_i=-c_j\) and
\(F_{ji}=-F_{ij}\), the unblown interactions cancel
and we obtain \(u_i-u_j+p(F_{i3}-F_{j3})=0\).
The evaluations are polynomial, so the constant term
gives the same constant abscissa term \(s\) for both.
Write \(u_i=\alpha p^2+\ell_i p+s\) and similarly
for \(j\), with \(\ell_i\ne\ell_j\).
The ordinates have the expansions
\[
 v_i=\beta p^3+J\ell_i p^2+(Js+K\ell_i^2)p+e_i,\qquad
        J=\frac{3\beta}{2\alpha},\quad K=\frac{3\beta}{8\alpha^2}.
\]
These terms are forced by degrees five and four
in \(v_i^2-u_i^3=D(p)\). Since the two sections
are disjoint, equal abscissas at \(p=0\) force
\(e_i=-e_j\ne0\).

Put \(k=(\beta-1)/(\alpha-1)\). The evaluations
\(F_{i3},F_{j3}\) have degree at most one with
linear coefficient \(k\). Comparing degree two
in \(v_i+v_3=F_{i3}(u_i-u_3)\) gives constant
coefficient \((J-k)\ell_i/(\alpha-1)\) for
\(F_{i3}\), because \(a=0\). Their difference is
\(-(\ell_i-\ell_j)\) by the subtracted residue equation,
so \(J-k=-(\alpha-1)\) and \(F_{i3}=kp-\ell_i\).
Comparing degree one now shows that
\[
                 (K+1)\ell_i^2=(k-J)s+(3/2-k)b
\]
for both \(i,j\). If \(K=-1\), then \(\beta^2=\alpha^3\)
would give \(\alpha=9/64,\beta=-27/512\), which fails
\(J-k=-(\alpha-1)\). Thus \(\ell_i^2=\ell_j^2\)
and their distinctness gives \(\ell_j=-\ell_i\).
It follows that \(u_j(p)=u_i(-p)\) and
\(v_j(p)=-v_i(-p)\). Since both sections are on the
same curve, \(D(p)=D(-p)\), impossible for a cubic.
This excludes the exceptional same-sign situation
and establishes (4) for all remaining quadratic sections
with \(\alpha\ne1\).

## The equal leading abscissa \(1\)

Suppose an unblown quadratic has \(\alpha=1,\beta=-1\).
The evaluation \(F_{i3}\) has growth degree at least two
because its numerator has leading term \(-2p^3\) and
the nonzero abscissa difference has degree at most one.
It cannot have degree three: then \(pF_{i3}\) in (2)
would have degree four. Thus its denominator has degree
one. The unblown evaluation with the other section
must have degree three in (2). This forces the other
section to have the same leading quadratic and linear
abscissa terms, the same leading ordinate, and a
different constant abscissa term. Write their common
slope \(a_*\), constant terms \(b_1,b_2\), and
\(\Delta=b_1-b_2\ne0\). Their leading residue equations
give \(c_1=-c_2=c\) and \(\Delta=c(a_*-a)\).
The middle equation in (3) then gives \(a=0,\ b=c^2\).

For the first unblown mark, clear the denominators in
(2) by \((u_1-u_2)(u_1-u_3)\).
Its coefficients of degrees four and three are
\[
        2(a_*c-\Delta),\qquad
        c\bigl(2(b_1-c^2)+3a_*^2\bigr)
                                     -\tfrac32a_*\Delta .
\]
Indeed \(L+u_1=u_1-u_3\), and the coefficients of
\(p^3,p^2\) in \(v_1+v_2\) are \(-2,-3a_*\),
while in \(v_1+v_3\) they are \(-2,-3a_*/2\).
Thus \(b_1=c^2-3a_*^2/4\). Interchanging indices
and replacing \(c\) by \(-c\) gives the same value
for \(b_2\), a contradiction.

If \(\alpha=1,\beta=1\), then \(F_{i3}\) has leading
term \(3p/2\). When \(u_i-u_3\) is linear or constant
this follows by comparing the first differing
abscissa coefficient in the expansion of their
opposite-sign ordinates; if the abscissas coincide
it follows from the regular evaluation at opposite
sections. In (2) the \(p^2\) coefficient of \(L+u_i\)
vanishes, whereas that of \(pF_{i3}\) is \(3/2\).
Cancellation requires the other unblown section to
have the same leading coefficients \((1,1)\), by
the exclusion of \((1,-1)\) above. The two linear
abscissa coefficients \(a_1,a_2\) must be distinct
to give degree two in \(F_{12}\). Its leading
coefficient is \(2/(a_1-a_2)\). Thus \(c_1=-c_2=c\)
and \(a_1-a_2=4c/3\). The middle equation of (3)
gives \(a=0\) and \(b_1=b_2\).

For either index with \(a_i\ne0\), polynomial
division in \(v_i+v_3=F_{i3}(u_i-u_3)\) gives
\(F_{i3}=3p/2+3a_i/8\). The linear coefficient
of (2) gives \(11a_i=9(a_1+a_2)\): the contribution
of the unblown interaction is
\(-9(a_1+a_2)/8\). If both slopes are nonzero,
the identities make them equal. If one is zero,
the equation on the other forces it to be zero.
Both contradict their distinctness.

## Linear sections and the remaining quadratic type

Once (4) has been justified, an unblown quadratic
with \(\alpha\ne1\) satisfies
\(\beta=1-(\alpha-1)^2=2\alpha-\alpha^2\).
Together with \(\beta^2=\alpha^3\) and \(\alpha\ne0,1\),
this gives \(\alpha=4,\beta=-8\).
For this quadratic type with \(u_i=4p^2+\gamma_i p+\cdots\),
comparison in \(v_i+v_3=F_{i3}(u_i-u_3)\) gives
\(F_{i3}=-3p-3a/2\).
For a linear type \(u_i=\gamma_i p+\cdots\) with
\(\gamma_i\ne0\), it gives \(F_{i3}=p+\gamma_i+a/2\).

If both unblown sections are linear, (3) gives
\(a=2r\). The \(p\) coefficients in (2) give
\(2\gamma_1+c_2[p]F_{12}=0\) and
\(2\gamma_2-c_1[p]F_{12}=0\).
If their slopes differ, \(F_{12}\) is constant.
If their abscissas coincide, they are opposite,
and disjointness at finite values forces their
ordinates to be nonzero constants. Their regular
evaluation then has degree two, also impossible
in (2). Thus the slopes are equal and the
constant terms distinct. The two linear equations
now give \(r=0\), and the constant coefficient
of \(U\) in (3) forces the constant terms equal.

If one unblown section is linear and the other is
quadratic of type \((4,-8)\), label the quadratic
as index two. Equation (3) gives \(a=2(c_1-c_2)\).
The unblown evaluation \(F_{12}\) has leading
term \(2p\). The \(p\) coefficient of (2) for
the linear section then reduces to \(2\gamma_1\ne0\).

Finally if both are quadratic of type \((4,-8)\),
the two \(p^2\) coefficients of \(L+u_i+pF_{i3}\)
vanish. The evaluations between the distinct
unblown sections have growth degree at least two:
their leading ordinates agree, and the abscissa
difference has degree at most one and is nonzero.
This contradicts (2), whose other terms have
degree at most one. All possible types are excluded.

## The remaining critical-point pattern

**Corollary SixSheetWidthFiveAlternatingEllipticCriticalPattern.**
In this profile, a width-five genus-one coordinate has one
double and one simple critical point on the immersed branch
normalization, with distinct critical values.

By the
[transverse cusp result](KELLER_ALTERNATING_ELLIPTIC_TRIPLE_POLES.md#directional-and-total-degree-bounds),
the polynomial projection of the branch normalization has
local degree two at the affine cusp. Its degree is five, so
its derivative has degree four, with three remaining zeros
counted with multiplicity. Three simple zeros were excluded
above. If all three occurred at one immersed point, the local
degree there would be four. Its critical value could not be
the affine cusp value, since the two local degrees would sum
to six. The immersed index-three local power map then gives
a source fibre singularity analytically equivalent to
\(u^3=v^4\). This is a unibranch singularity of delta invariant
three. The compactified fibre is smooth elsewhere, using the
[relative compactification and polar sections](KELLER_SIX_SHEET_CUSP_SLICE_GENUS.md#fibre-genus-and-canonical-bundle).
It is connected and reduced, so this fibre is integral. Its
arithmetic genus would be at least three, contrary to genus
one. The multiplicities must therefore be two and one.

The double and simple immersed points give an ordinary triple
point and an ordinary cusp in the source fibres, respectively.
If they had the same critical value, their contributions to
\(\sum_x(\delta_x-b_x+1)\) would sum to two in a genus-one
fibre. Thus their values are distinct.
