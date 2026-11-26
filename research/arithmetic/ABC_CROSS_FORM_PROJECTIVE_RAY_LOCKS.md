# Projective ray class counts with a divisor of another form

For a binary norm form, a divisor of a different coprime form
can determine the projective input modulo the whole divisor.
The ray class moment then counts the first form at half its
radical cost, including when that form has degree two.
For weighted squares and a cubic target, this gives a
small-budget saving with no bound on individual prime sizes.

## A divisor of a different form

Use the parameters and polynomial coefficient and input
height bounds of (PR1). Let \(F\) be the irreducible form
of degree \(d\ge2\), and let \(G\) be a homogeneous binary
form of fixed degree. Replace \(M_1\mid v,M_2\mid u\) in
(PR1) by \(M\mid G(u,v)\). At every prime dividing \(M\),
require that \(G\) have distinct projective roots, that
\(F,G\) have no common projective root, and that
\(2F(1,0)F(0,1)\operatorname{Disc}(F)\) be a unit. Retain the
radical cap on \(F(u,v)\) outside the optional prescribed
support \(P\) and all the other hypotheses of (PR1).
Then (PR2) still bounds the count, including for signed
primitive input pairs. In particular the bound is
\[
 X^{(\mu+\ell/2+\eta)L+\varepsilon},
 \qquad |\operatorname{Disc}(K)|\ll X^{g_F L},
 \tag{XL1}
\]
up to the implied constant of (PR2).

For each \(p^e\parallel M\), take a projective chart in
which the input denominator is a unit. Its coordinate is
congruent modulo \(p^e\) to one of the simple roots of
\(G\) in that chart. This follows from Hensel's lemma.
Recording the chart and root costs
\((2\deg G)^{\omega(M)}=X^{o_X(1)}\). After this record,
choose a normalized vector \((u_p,v_p)\) at each prime so
that \((u,v)\) is a unit scalar multiple of it modulo
\(p^e\). Write \(\alpha=a_0u-\theta v\) in the norm field
of \(F\), as in (PR1). By the Chinese remainder theorem
choose an integral \(\alpha_*\) congruent to
\(a_0u_p-\theta v_p\) on each entire \(p^e\mathcal O_K\).
This is a unit there, since the norm of that residue
is a unit by coprimality of \(F,G\). Then
\(\alpha/\alpha_*\) reduces to a rational unit modulo
\(M\mathcal O_K\), with one rational residue obtained by
the Chinese remainder theorem from the local scalars.
This is exactly the residue condition in (PR3) with
the full ordinary modulus. The moment, diagonal norm
comparison and unit count in (PR2) apply with this
reference. The sign of an input plays no role in those
bounds.
[Projective ray class moments (PR1–PR5)](ABC_PROJECTIVE_RAY_CLASS_NORM_MOMENTS.md)

## Cubic targets with asymmetric square kernels

Fix \(f,g,k\ge0\) with \(k\le f\). Consider positive
primitive triples \(a+b=c\le X\) written
\[
 a=A x^2,\qquad b=B y^2,\qquad c=C z^3,
\]
where \(A,B\) are squarefree and \(C\) is the canonical
cube-free defect. Fix one of the summands, denoted \(U\),
and let \(D_U=A\) or \(B\) be its squarefree kernel. Impose
the height-relative constraints
\[
 \operatorname{rad}(abc)\le c^L,\qquad
 AB\le c^{fL},\qquad \operatorname{rad}(C)\le c^{gL},
 \qquad D_U\ge c^{kL}.
\]
If
\[
                  f+g-\frac12 k<\frac16,
 \tag{XL2}
\]
there is a fixed \(\delta>0\) such that the count is
\(O_{f,g,k,L,\varepsilon}(X^{(1/3-\delta)L+\varepsilon})\)
for every sufficiently small fixed \(L>0\).
In particular one may remove the lower cap on the kernel
of a specified summand and get the same conclusion under
\(3f/4+g<1/6\). Bin the actual joint kernel mass and
choose the summand with at least half that mass.
For exact cube targets with nearly all the kernel mass
on one summand, the allowed joint cap approaches \(1/3\).

Work in height shells. By the carrier construction in
(WG4), fixing the kernels, the cubic defect and an
ideal class costs \(X^{(f+g)L+o_X(1)}\) in total. The
class count is averaged over \(D=AB\), using
Davenport–Heilbronn for \(|\mathrm{Cl}(\mathbf
Q(\sqrt{-D}))[3]|\). For each record take the norm-\(m\)
ideal in the inverse class, \(m\ll\sqrt D\), and write
\(\alpha=\beta\gamma^3\) with
\(\gamma=(r+s\sqrt{-D})/2\). The underlying projective
forms for the square roots have degree three and are
rational multiples of the real and imaginary coordinates
of \(\beta\gamma^3\). The form for the cube root has
degree two and is proportional to \(r^2+Ds^2\).
Use the reduced primitive parameter of (WG4).
[Weighted carrier and ideal-class count (WG4)](ABC_WEIGHTED_GAUSSIAN_CYCLOTOMIC_DESCENT.md#projective-recovery-with-arbitrary-ideal-classes)

Let \(J_0\) be the squarefree product of all primes dividing \(AB C m\) and the fixed primes at most seven.
Outside that support all coefficients can be made
integral without a nonunit scalar. The two linear forms
\(\gamma,\bar\gamma\) are independent there, the two
coefficients \(\beta,\bar\beta\) are units, and the
underlying forms have distinct roots and mutually
disjoint root sets in projective reduction. Indeed
their equations are the two nondegenerate degree-three
binomials in \(\gamma,\bar\gamma\) and the equation
\(\gamma\bar\gamma=0\). As the equal-degree forms for
\(a,b,c\) have no common projective zero in this reduction,
passing to the primitive triple removes no factor at
these primes. Likewise the distinct rational factors
of a single underlying form have disjoint value supports
there. Full factors supported on \(J_0\) can be fixed
at \(X^{o_X(1)}\) cost.

Choose an irreducible factor \(F\) of degree \(d=2\) or
\(3\) of the selected summand's underlying form if one
exists. Its field discriminant satisfies
\[
        |\operatorname{Disc}(K_F)|
            \ll AB\,\operatorname{rad}(C)^2 .
 \tag{XL3}
\]
In fact its splitting field is contained in
\(\mathbf Q(\sqrt{-D},\zeta_3,
(\pm\bar\beta/\beta)^{1/3})\). This follows by taking
the ratio of the two independent linear forms in either
degree-three binomial. Outside \(6DC\), the valuations of
\(\bar\beta/\beta\) at all primes of the quadratic field
are multiples of three by
\((\beta)=\mathfrak a\mathfrak c\mathfrak t^{-3}\).
A cube root of a local unit at residue characteristic
different from three lies in an unramified extension.
At a ramified prime over \(p\mid D\), \(p\nmid6\), the
same ratio is a unit because conjugation fixes the place.
Thus the full splitting field has ramification index
at most two there. In a degree-two or degree-three
factor field the tame discriminant exponent is at most
one at such a prime. At \(p\mid C\), \(p\nmid6\), the
tame discriminant exponent is at most \(d-1\le2\).
The local discriminant exponents at the fixed primes
\(2,3\) are bounded in degree at most three. These bounds
give (XL3) without including the norm \(m\) in its support.

Suppose the selected summand has radical cap
\(X^{\ell_U L}\). The radical of the factor value outside
\(J_0\) is at most \(X^{(\ell_U-k)L+o_X(1)}\), since
the entire kernel \(D_U\) is in \(J_0\). A positive-weight
divisor of the underlying form on either other letter,
supported off \(J_0\), can be used for (XL1). Its prime
support needs the additional end-coefficient condition
in (XL1). Choose a finite list of unimodular bases whose
columns are consecutive vectors \((1,j),(1,j+1)\)
for bounded \(j\). At any sufficiently large prime, one
of these bases has both columns off the at most \(d\)
forbidden lines where \(F\) vanishes. Allocate the primes
of the divisor to bases and keep a block of positive
numerical weight. This leaves the radical cost at most
\(\mu\), and the height and discriminant estimates persist
after the basis change. The other good-prime conditions
of (XL1) hold by the projective root properties above.
The resulting exponent in units of \(L\) is at most
\[
       (f+g)+\mu+\frac{\ell_U-k}{2}
                       +\eta .
 \tag{XL4}
\]

If the selected cubic splits into three rational lines,
choose two whose value radicals outside \(J_0\) have
smallest combined cap after dyadic refinement. That
combined cap is at most \(2(\ell_U-k)/3\) up to the
refinement loss. Their values are units at the external
modulus and their ratio is prescribed by its projective
branch record. Fix their factors on \(J_0\), and then
factor each remaining value into its exact valuation
groups. Cap a high-valuation tail by its radical at
arbitrarily small fixed cost, using the polynomial value
height, and fix it. The retained roots have pairwise
disjoint prime supports in each row. Prime colouring
by position costs only \(X^{o_X(1)}\) and makes this
disjointness hold across the signed ratio dictionaries.
Apply character orthogonality and the singleton moments
of (JD1) to that fixed ratio congruence. Their proof uses
only the congruence and the bounded raw root caps. An
external modulus of fullness comparable to \(1/L\) and
positive numerical weight therefore counts the two
linear values at half their combined cap, or at most
\((\ell_U-k)/3+\eta\) in units of \(L\). These two
values determine the reduced projective parameter since
the lines are distinct.
The split branch therefore also saves under
\[
                    f+g-\frac{k}{3}<\frac29                  \tag{XL5}
\]
when the selected underlying carrier form splits on
the retained rows. The exponent for the split branch
is \(f+g+(\ell_U-k)/3+\mu+\eta\).

[Signed root moments (JD1)](ABC_JOINT_FREE_LETTER_DICTIONARIES.md#a-signed-joint-dictionary-estimate)

## Obtaining the external divisor at balance

If (XL2), or (XL5) on its split branch, failed to give a uniform saving, choose
violating fixed budgets tending to zero after the inner
height limits. The limiting one-third bound forces a
near-saturating box to use \(1-o(1)\) of the ambient
budget. The balance and pair nonatomic conclusions apply,
so the selected cap is \(\ell_U=1/3+o(1)\) and the
two complementary letters have combined nonatomic mass
at least \(1/3\). Across the two complementary letters, the support
excluded by \(J_0\) has total radical mass at most
\(3f/2+g-k+o(1)\). The primes of the selected kernel
\(D_U\) occur on neither complementary letter by
primitivity. The other summand parity kernel costs
at most \(f-k\), the cubic defect at most \(g\), and
the auxiliary carrier ideal norm at most \(f/2+o(1)\).
The remaining fixed primes have vanishing cost.
Under (XL2), \(3f/2+g-k\le(3/2)(f+g-k/2)<1/4\).
Under (XL5), \(3f/2+g-k\le(3/2)(f+g-k/3)<1/3\).
Take bounded positive-scale bands in the other pair
whose nonatomic mass exceeds the excluded support cap
by a fixed margin. Remove small windows around finitely
many atoms so that each exact-degree position left has
raw mass at most an arbitrarily small fixed allowance.
The uniform first moment controls the high-scale tail.
The selected bands and deletions can be made with
continuity endpoints before using the limiting boxes.
[Balanced endpoint and pair nonatomic masses](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

Prime powers at those positions outside \(J_0\) still
have positive total radical mass. At fixed carrier
record, bin the retained sizes at subpower cost. Greedy
grouping into blocks of small raw mass yields a block
of cost at most any desired fixed \(\mu>0\), with
numerical weight bounded below along a subsequence.
Use it on one of the two complementary letters. When
passing to the underlying square or cube root form,
divide its exponents by two or three, respectively.
At these good primes the reduction is exact. The
resulting divisor has fixed positive weight, unchanged
radical, and fullness comparable to \(1/L\). Any
finite subdivision by projective basis preserves a
positive-weight block. The total choices of its full
divisor values cost at most \(X^{\mu L+o_X(1)}\).
The margins of this construction are fixed before the
inner height limits where (XL1) and (JD1) are applied.

With \(\mu,\eta\) chosen below the strict margin,
the count in (XL4) or the split count preceding
(XL5) contradicts the exponent \(1/3\) in its
respective sector.
Dyadic shells give the stated cumulative count.
The weighted square and cube defect hypotheses are part
of this sector. The argument leaves the general
balanced endpoint open.
