# Quartic invertibility: written source audit and exact identities

This audit accepts only the invertibility conclusion of Wang-Yan,
Quartic Keller Maps and Tameness in Three Variables, arXiv2610.36495v1,
Theorem3.1 and its characteristic-zero collision consequence. It makes no
claim about the separate tameness theorem. The source is the27 October2026
version; the primary submission history checked on14 November lists onlyv1.
This is a written mathematical review with exact algebraic replays, not a
formal verification or identified independent external refereeing.

## Source and exact trust boundary

The [primary HTML](https://arxiv.org/html/2610.36495v1), the source archive and
its TeX are retained locally. HTML SHA256 is
04dfa22c70f97783543eb1a011d90e4d4b9f4fb8f7af394ef6969c1aa491ca29;
archive SHA256 is
3ecded584936981cf52fa53a4a6fd5c7646729e8b56dd6287b35bd76a05fb9a5.
All of the invertibility argument, Sections2-3, has now been read. The
following paragraphs record the logical audit; the two checker receipts
specify exactly which algebraic calculations were independently reconstructed.
Neither receipt is advertised as an independent proof of all geometric
arguments or of the imported Moh plane theorem.

The established external ingredients are planar invertibility through100,
already imported locally as KPAIR1, and Ax's injective-to-surjective theorem
plus the standard fact that a bijective etale map here is an isomorphism.
For our collision exclusions it suffices to use the fibre injectivity, so
no family inverse is assumed. Finite-generation embedding transfers a
hypothetical characteristic-zero collision toC. Formal inverse uniqueness
also supports the source's stronger inverse descent, but is unnecessary
for the census exclusion.

## Closed generators and leading factorization

The normal-subalgebra lemma is sound without a hidden finite-generation
assumption: a least-degree monic polynomial makes C[t] integral overB.
The minimal polynomial of t has coefficients inB. Its nonconstant
coefficient has at least the chosen degree, and total-degree additivity
forces the remaining factor of u(T)-u(t) to be constant. The tower law
then identifies the fraction field and B=C[u].

For a polynomial r, its integral closure in the ambient rational field
lies in the polynomial ring by normality. Restriction to a line where r
is nonconstant is injective: a nonzero integral b gives a nonzero c(r)
in bB, which cannot vanish identically on that line. This proves the
closed-generator statement. Clearing a leading coefficient in an algebraic
equation proves equivalence with relative algebraic closure. The differential
criterion follows from separable extension of derivations. These establish
the rational and polynomial kernel statements used later, rather than
assuming an unproved rational-to-polynomial passage.

For three independent equal-degree leading forms with singular Jacobian,
Euler gives projective differential rank one. A suitable source line
parametrizes the image curve dominantly, making its normalization rational.
Writing its homogeneous parameter p/q and the normalization coordinates
Phi gives H_i=g*Phi_i(p,q). Coprimeness of the substituted Phi follows
from powers of both parameter variables in their ideal. Hence g is a
polynomial, and degree4 has exactly the four listed degree patterns.
No compatibility of an identity linear part with a leading normal form
is assumed: the linear part remains arbitrary throughout.

## Directional and binary/common-factor reductions

The affine lemma uses absence of common zeros of the z-derivative entries,
closedness of the third component, and the kernel identity. A common factor
of p,q would have a complex zero, so p,q are polynomials in r and the
remaining determinant is a plane Jacobian. The degree bound is explicit.

For constant second directional derivative, the independent affine-derivative
case reduces to a nonzero three-function Wronskian. Its determinant is
v_y*(Omega*v+Lambda), impossible over C(x)[y]. In the quadratic case a common
closed generator for p,q must be linear unless1,p,q are dependent. Otherwise
positive distinct degrees give nonzero W,K, and the constant determinant
is u_y times a genuinely quadratic polynomial in u; the leading coefficient
K'-3K*W'/W cannot vanish by degree comparison. The surviving constant
combination straightens a coordinate with plane degrees at most7/8.
All integrations are in the polynomial ring over C(x), where the required
nonzero Wronskians are units; no exceptional coefficient stratum is dropped.

The binary-leading lemma's high z coefficients kill the cubic z^3 term,
then its z^2 linear vector. Rank zero/one uses the absence of a constant
normal annihilator; rank two is normalized to(x,y,0) and the next coefficient
contradicts independence. The checker reconstructs the universal high
coefficients for both s=0 ands=1, including arbitrary binary remainders.
For a quadratic common factor of z-degree2, the exact adjugate gives
4g*S_z=3g_z*S. Leading z-degree would require4*deg_z S=6, so S=0.
The syzygies and a source translation remove remaining cubic z^2 terms;
the next two coefficients force the z column of the linear part to vanish
on x=y=0, contradicting the Keller condition. The adjugate and these high
coefficient identities are independently reconstructed.

## All conic-pencil strata

A coprime nonsquare pencil has a smooth generic conic: an identically singular
ternary pencil either has a common linear factor or is binary with a rank-one
member. The generic conic has a constant base point and rational function
field C(lambda)(t). The invariant equation forces f^2=q^d*c(p/q).
Coprimeness forbids poles, and squarefreeness of every member forces all
binary-factor multiplicities even. Thus the invariant ring is exactlyC[p,q].
Successive odd-degree invariant equations and elementary syzygies then force
all corrections of a nonsquare Veronese leading part to factor through two
polynomials, contradicting nonzero Jacobian.

A square member reduces the pencil to(x^2,y^2),(x^2,yz) or(x^2,y^2+xz).
The first is binary. For each other case the independent checker derives
the coefficient operator directly from the actual quartic Jacobian, verifies
its complete cubic kernel(rank20,nullity10) and quadratic kernel(rank12,
nullity6), and uses every rational left-nullspace compatibility equation.
Their ideal is exactly(a*alpha,b*alpha), with no coefficient parameter
inverted. It verifies all proposed quadratic coefficients and all14 selected
linear-obstruction identities. The parabolic proof uses a*M=a*gamma*U,
never division by a, so its a=0 stratum remains covered. These give a singular
arbitrary linear Jacobian in both nonbinary cases. This closes every leading
pattern of quartic rank3 after the written geometric reductions.

## Cubic scalar classification and residue arguments

If quartic rank is smaller, a target combination R has degree at most3.
The quadratic case is the radical argument already isolated in GLDC5.
For a cubic R, the Keller identity makes its residue2-form the derivative
of the regular1-form P*dQ/kappa on every smooth fibre. At infinity a reduced
nodal cubic has a simple-pole residue with a nonzero period: at a node it
is a unit times du/u, and on a smooth compact cubic a nonzero holomorphic
form cannot have all periods zero. A normal-circle tube in the smooth
boundary locus gives this period as a limit of integrals of an exact
2-form, a contradiction. Loops avoid boundary singularities, so no smoothness
of the entire projective closure or proper affine inverse is assumed.

The second residue lemma works in C[x,x^-1,z]. A closed1-form differs
from an exact one by c*dx/x. Restriction along the d generic branches at
x=0 gives p1(z_i)=c(t)*p0'(z_i). The root sum derivative is independent of
t by Viete, so c(t) is constant and p1=lambda*p0'. This gives a critical
point. The argument also covers deg p1=deg p0; the leading root-sum quotient
still has t-independent derivative because d>=2. These are written geometric
and Laurent-de Rham arguments, not claims inferred from numerical periods.

The remaining leading types are binary, cusp and conic-tangent. In the cusp
and tangent normal forms, every eliminant leading coefficient used to obtain
a critical point is checked symbolically, including the nonzero coefficient
strata where division is used. Their zero strata are treated separately in
the source proof. The residue lemma eliminates the residual cuspidal case
and the inappropriate tangent coefficient. The surviving tangent case has
an explicit Nagata-type completion; both inverse compositions, determinant1
and inverse degree5 are independently checked. The planar noncritical cubic
classification follows by Bezout in the squarefree leading case and direct
elimination in the repeated-root cases. Binary leading forms then give a
graph or the displayed degree4 inverse completion, apart from r=x+x^2y.
Every use here depends only on degreeR<=3, absence of critical points and
exactness of its fibre form, not on degreeP or degreeQ. This observation is
used explicitly in the separate degree5 extension, rather than silently
importing a stronger stated proposition.

## Exceptional component and conclusion

For r=x+x^2y, D=-x^2*d_x+(1+2xy)*d_y becomes d_u in coordinates(r,u=1/x).
There is no polynomial slice Dg=constant: a rational correction h(r) would
have denominator dividing x, impossible for a nonconstant denominator in r.
The formal residue identity for a component monic in z produces exactly
such a slice. The affine-pair case uses the earlier closed-generator lemma.
For total degree4 the highest coefficient ratio p^n/q^m in C(r),
squarefree degree3 fibres and degree bounds force a constant ratio and then
(m,n)=(2,1), p=lambda*q^2. A polynomial target shear makes both components
affine in z; their possibly increased total degrees do not affect that
last obstruction. No total-degree-preserving claim is used for this last
nonlinear target shear.

Consequently all cases of Theorem3.1 are covered: independent quartic leading
parts by the leading-pattern argument, and dependent ones by the scalar
classification and exceptional obstruction. The coordinate completions have
inverse degree at most5, so plane degrees at most20 are within the existing
Moh100 import. Accepting the written audit gives quartic invertibility,
with independent algebraic controls specified above and the geometric
trust boundary explicit. No new smaller positive witness follows.

## Replay

```
.venv/bin/python research/jacobian/quartic-three-audit/verify_square_pencils.py
.venv/bin/python research/jacobian/quartic-three-audit/verify_structural.py
```

The source archive and HTML are retained, not externally modified. No author
code was executed. No tameness audit or global simplicity optimum is claimed.
