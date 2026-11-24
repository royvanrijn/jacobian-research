# Polynomial odd-return cores and an exact Cayley obstruction

Work over a field \(K\) of characteristic zero. A positive odd integer
\(a\) is fixed unless stated otherwise. A normalized polynomial pair
\(C,D\in K[x]\), \(C(0)=-D(0)=1\), belongs to the fibre of the
normalized rational endpoint \(u(0)=0,u'(0)=1\) when the phase

\[
 p=t^a\frac{C(x)+t^2D(x)}x,\qquad t^2+2xy=1,
\]

has every positive pure moment zero and its selected endpoint is
\(u=xX^a/C\), where \(X^2=-D/C\) and \(X(0)=1\).
The pair \((1,-1)\) is called trivial. Orders and pole statements
below are geometric, on the projective line, unless specified otherwise.

This note proves that every surviving pair without a polynomial
reciprocal direction descends at the same height to a polynomial pair
whose return is literally \(z\mapsto-z\). If that outer pair is
nontrivial, the inner map is polynomial and the outer endpoint vanishes
at infinity. It also excludes the entire family
\(z/(1-\lambda z^2)^m\) at every positive odd height. At height one,
an endpoint with one pole pair has an exact coefficient criterion;
an exact diagnostic excludes all pole orders for 44 specified rational
divisors. These statements do not decide whether a nontrivial odd-return
core exists in general.

The necessary retained inputs are the
[fixed-coordinate reciprocal-unit theorem](ODD_START_FIXED_COORDINATE_FIBRE_CLASSIFICATION.md),
the [canonical rational reflection quotient and pole box](RATIONAL_REFLECTION_QUOTIENT_AND_FINITE_FIBRES.md),
and the sharp lower bound and equality classification in
[AllStartSharpMinimumUnsafeDegree and AllStartMinimumDegreeClassification](ALL_START_HEIGHT_SHARP_UNSAFE_DEGREE.md),
Sections 1–3. Their uses are identified below. No general Orbit theorem,
integer-lattice search, or Padé existence theorem is a premise.

## 1. Same-height polynomial descent

Use the canonical maximal reflection quotient of the fixed coordinate:

\[
 \begin{gathered}
 u=U(h),\qquad h(0)=0,\quad h'(0)=1,\\
 U(0)=0,\quad U'(0)=1,\quad U(-z)=-U(z),\\
 \eta=V(h),\qquad V(0)=1,\quad V(z)V(-z)=1.
 \end{gathered}
\]

Here \(\eta\) is a surviving reciprocal unit. The original pair is
reconstructed by

\[
 X=\frac{u(1+\eta)}{2xu'},\qquad
 C=\frac{xX^a}{u},\qquad D=-\frac{xX^{a+2}}u.
\]

Assume that there is no nonconstant polynomial reciprocal unit. By the
retained polynomial-direction criterion, this is equivalent to the
absence of a finite nonzero quotient point \(r\) with
\(h^{-1}(r)=\{\infty\}\). This criterion is unchanged by extending
the coefficient field: any such exceptional value would equal
\(h(\infty)\in K\).

Define rational outer functions

\[
 \bar X(z)=\frac{U(z)(1+V(z))}{2zU'(z)},\qquad
 \bar C(z)=\frac{z\bar X(z)^a}{U(z)},\qquad
 \bar D(z)=-\frac{z\bar X(z)^{a+2}}{U(z)}.
 \tag{R1}
\]

They are regular at zero, with
\(\bar X(0)=\bar C(0)=-\bar D(0)=1\). Direct substitution gives

\[
 \begin{aligned}
 X&=\bar X(h)\frac h{xh'},\\
 C&=\bar C(h)\frac{(h/x)^{a-1}}{(h')^a},\\
 D&=\bar D(h)\frac{(h/x)^{a+1}}{(h')^{a+2}}.
 \end{aligned}
 \tag{R2}
\]

**Theorem NoDirectionOddReturnDescent.** The outer functions
\(\bar C,\bar D\) are polynomials over \(K\). They form a normalized
pure pair at the same height \(a\), with endpoint \(U\) and literal
marked return \(z\mapsto-z\). The maximal reflection field of this
outer problem is \(K(z)\), so its canonical quotient is the identity
and it too has no polynomial reciprocal direction.

**Proof.** Let \(r\ne0\) be a finite quotient point. There is a finite
source preimage \(\beta\) of \(r\): otherwise the whole fibre would
be infinity, contrary to the hypothesis. Such \(\beta\) is nonzero.
Put \(e=\operatorname{ord}_\beta(h-r)\). At \(\beta\), the factor
\(h/x\) is a unit and \(h'\) has order \(e-1\). Thus

\[
 \begin{aligned}
 \operatorname{ord}_\beta C
   &=e\operatorname{ord}_r\bar C-a(e-1),\\
 \operatorname{ord}_\beta D
   &=e\operatorname{ord}_r\bar D-(a+2)(e-1).
 \end{aligned}
 \tag{R3}
\]

Both left sides are nonnegative. Therefore both outer orders are
nonnegative. This argument includes finite nonzero zeros and poles of
\(U\), its ordinary critical points, and ramification of \(h\);
it retains the full orders, not only function values. At \(r=0\),
regularity already follows from (R1). Hence neither outer function has
a finite pole, proving polynomiality. No condition at quotient infinity
is needed for this step.

Oddness and the normalization \(U'(0)=1\) identify the marked inverse
branch as \(-z\). The identity \(V(z)V(-z)=1\), together with (R1),
is exactly the reciprocal form of the all-order endpoint criterion in
the fixed-coordinate theorem. It proves purity and the endpoint
assertion. Literal reflection preserves the full field \(K(z)\),
so the maximal reflection field is already that field. Its identity
quotient has no finite nonzero fibre supported only at infinity.
All formulas are over \(K\). \(\square\)

## 2. The inner map is polynomial unless the outer pair is trivial

Suppose \(\beta\) is a finite pole of \(h\), of order \(e\ge1\),
and put \(c=\deg\bar C\), \(d=\deg\bar D\). Here \(\beta\ne0\).
The respective orders of \(h/x\) and \(h'\) are \(-e\) and
\(-e-1\). Formula (R2) gives the exact identities

\[
 \operatorname{ord}_\beta C=a-e(c-1),\qquad
 \operatorname{ord}_\beta D=a+2-e(d-1).
 \tag{R4}
\]

In particular, polynomiality would force \(d\le1+(a+2)/e\le a+3\).

**Theorem PolynomialOrTrivialOddCore.** Exactly one of the following
possibilities holds for the descended pair.

1. It is trivial: \(\bar C=1,\bar D=-1\). Then
   \(\bar X=1,U=z,V=1\). Upstairs \(\eta=1\) and
   \(X=u/(xu')\) is polynomial, giving the already classified
   Davenport–Zannier singleton branch.
2. It is nontrivial, and \(h\) is polynomial. Thus any remaining
   nontrivial existence problem has a polynomial inner map and a
   same-height polynomial outer pair with an odd rational endpoint.

**Proof.** The retained minimum-degree theorem and its equality
classification say that a nontrivial normalized pure pair has
\(d\ge a+3\), with equality precisely when, for some
\(\lambda\in K^*\),

\[
 \bar C=(1+\lambda z)^{a+1},\qquad
 \bar D=-(1+\lambda z)^{a+3},\qquad
 U=\frac z{1+\lambda z}.
 \tag{R5}
\]

A finite pole of \(h\) would force this equality by (R4). But the
endpoint in (R5) is not odd for any nonzero \(\lambda\). Hence a
nontrivial outer pair allows no finite pole of \(h\), proving that
\(h\) is polynomial. The equality classification is essential here:
the lower bound alone would not rule out a simple finite pole.

For the trivial pair, the selected root is \(\bar X=1\); its endpoint
is \(U=z\), and (R1) forces \(V=1\). Thus \(X=b=u/(xu')\)
upstairs. This rational function has a simple zero at every finite
unmarked zero or pole of \(u\), is a unit at zero, and has possible
finite poles only where \(u\) has finite nonzero value and \(u'=0\).
At those latter points \(C=x b^a/u\) could not be polynomial if
\(b\) had a pole. Therefore \(b\) has no finite poles and is
polynomial. These alternatives are disjoint. \(\square\)

For later constructions, all finite-source conditions in (R2) are
recorded explicitly:

| Source point | Exact orders of \(C,D\) or condition |
|---|---|
| \(x=0\) | Both are units with \(C(0)=-D(0)=1\). |
| \(\beta\ne0\), \(h(\beta)=0\), local degree \(e\) | Orders \(a-e\) and \(a+2-e\); the exact condition is \(e\le a\). |
| \(h(\beta)=r\ne0,\infty\), local degree \(e\) | The two orders in (R3) must be nonnegative. |
| A finite pole of \(h\), order \(e\) | The two orders in (R4) must be nonnegative. |

These cases exhaust the finite source. Poles of \(h\) impose
additional conditions and cannot weaken outer polynomiality. In the
nontrivial branch they are absent. If \(g=\deg h\) for a polynomial
inner map, taking orders at infinity in (R2) gives

\[
 \deg C=g(\deg\bar C-1)+1,\qquad
 \deg D=g(\deg\bar D-1)+1.
 \tag{R6}
\]

These are exact identities for an admissible pullback. They do not say
that an arbitrary polynomial \(h\) satisfies the finite conditions.

## 3. Every nontrivial odd core vanishes at infinity

In this section work on the quotient line, and write \(X,C,D\) for
the outer functions. We call such a normalized polynomial pair with
odd endpoint \(U\) an **odd-return core**. Oddness and reciprocity
give

\[
 \frac1{X(z)}+\frac1{X(-z)}=\frac{2zU'(z)}{U(z)}.
 \tag{R7}
\]

**Theorem OddCoreInfinity.** For every odd-return core, \(X\) is
regular at infinity and \(\deg D\le\deg C\). Every nontrivial
odd-return core satisfies \(U(\infty)=0\). More precisely, put
\(k=\operatorname{ord}_\infty U\).

If \(X(\infty)=\ell\ne0\), then
\(\ell=-1/k\) and \(\deg C=\deg D=1+k\). The case \(k=-1\)
is exactly the trivial core; otherwise \(k\) is positive and odd.
If \(X\) has a zero at infinity of order \(t\ge1\), then

\[
 \deg C=1+k-at,\qquad
 \deg D=1+k-(a+2)t.
 \tag{R8}
\]

Here \(t\) and \(k\) are positive odd integers,
\(k\ge(a+2)t\), and the nontrivial minimum-degree bound strengthens
this to \(k\ge(a+2)(t+1)\), with \(k\) odd.

**Proof.** A nonzero odd rational function has a Laurent leading term
of odd nonzero degree at infinity. Thus \(k\) is an odd nonzero
integer, and the right side of (R7) has limit \(-2k\ne0\).
If \(X\) had a pole at infinity, the left side would tend to zero,
a contradiction. Also \(D/C=-X^2\), proving the degree comparison.

If \(X\) tends to a nonzero \(\ell\), (R7) gives
\(2/\ell=-2k\); the degrees follow from \(C=zX^a/U\) and
\(D=-zX^{a+2}/U\). Their nonnegativity gives \(k\ge-1\).
At \(k=-1\) both coefficients are constant and normalized, hence
the pair is trivial. All other allowed odd \(k\) are positive.

If \(X\) vanishes to order \(t\), the same formulas give (R8).
Polynomiality already forces \(k>0\). The factor \(U/(zU')\)
is a unit at infinity, so \(1+V\) vanishes there to order \(t\).
Write \(w=1/z\). Since \(V(w)V(-w)=1\) and \(V(0)=-1\) in
this local coordinate, comparison of the first nonconstant terms forces
\(t\) odd. From \(\deg D\ge0\) one gets
\(k\ge(a+2)t-1\), improved by parity to \(k\ge(a+2)t\).
For a nontrivial core, \(\deg D\ge a+3\) gives the stronger bound.
All these restrictions are necessary, not existence assertions.
\(\square\)

In particular, a polynomial \(X\) in an odd core is constant, since
it is regular at infinity. Its normalization makes \(X=1\); then
(R7) gives \(zU'/U=1\), hence \(U=z\) and the core is trivial.
Every hypothetical nontrivial core therefore has rational nonpolynomial
\(X\), regular at infinity.

The no-direction hypothesis in the descent cannot be dropped. For
\(u=x/(1+x)\), take \(h=u,U=z\), and
\(V(z)=(1+z)/(1-z)\). The source reconstruction gives
\(X=(1+x)^2\), \(C=(1+x)^{2a+1}\),
\(D=-(1+x)^{2a+5}\), a pure polynomial pair. But its outer functions
are \(\bar X=(1-z)^{-1}\), \(\bar C=(1-z)^{-a}\),
\(\bar D=-(1-z)^{-a-2}\), which are not polynomials. The missing
finite fibre is \(h^{-1}(1)=\{\infty\}\), exactly the polynomial
reciprocal direction excluded in Section 1.

## 4. A complete rational-residue formulation at height one

**Theorem HeightOneOddCoreResidues.** Let \(C\in K[z]\) satisfy
\(C(0)=1\). It is the first coefficient of a height-one odd-return
core if and only if the following two conditions hold:

1. The even rational function
   \(E(z)=\tfrac12(C(z)^{-1}+C(-z)^{-1})\) has zero residue at
   every finite geometric point.
2. For its rational primitive \(U\), normalized by \(U(0)=0\),
   the function \(-U^2C^3/z^2\) is polynomial.

Then the endpoint and second coefficient are exactly

\[
 U'(z)=\frac12\left(\frac1{C(z)}+\frac1{C(-z)}\right),\qquad
 D(z)=-\frac{U(z)^2C(z)^3}{z^2},\qquad X=\frac{UC}{z}.
 \tag{R9}
\]

The residue condition is equivalently equality of the residues of
\(dz/C(z)\) at every reflected pair:

\[
 \operatorname{Res}_{z=\alpha} E(z)\,dz
 =\frac12\left(
    \operatorname{Res}_{z=\alpha}\frac{dz}{C(z)}
   -\operatorname{Res}_{z=-\alpha}\frac{dz}{C(z)}\right)=0.
 \tag{R10}
\]

An absent pole contributes residue zero.

**Proof.** At height one the endpoint criterion is
\(T(z)+T(-z)=2\), where \(T=1/(CU')\). Since an odd endpoint
has even derivative, it is equivalent to the first identity in (R9).
The other identities follow from the selected endpoint and root.

For the converse, partial fractions give a rational primitive exactly
when all simple-pole coefficients vanish; polynomial terms and poles of
order at least two integrate rationally in characteristic zero.
Equation (R10) follows from the substitution \(z\mapsto-z\) in
the differential. The normalized primitive belongs to \(K(z)\):
one can integrate the partial fractions over a finite splitting field,
and uniqueness after fixing its constant makes the result Galois
invariant. Since \(E\) is even, uniqueness also gives \(U(-z)=-U(z)\).
Moreover \(U'(0)=E(0)=1\). If the stated \(D\) is polynomial,
its value at zero is \(-1\), the selected root is \(X=UC/z\),
and the endpoint is \(U\). Equation (R9) gives the exact all-order
endpoint identity. Thus the pair is pure. \(\square\)

This criterion permits zeros of \(C\) away from the poles of \(U\).
The two summands in \(E\) can cancel their poles at points where
\(U\) remains regular; in the reciprocal description these are
\(\eta=-1\) points. Consequently one endpoint pole pair is not the
same condition as one root pair of \(C\). A bare product involving
only the endpoint poles can miss additional even factors of \(C\).

## 5. An entire repeated-pole family is excluded

**Theorem RepeatedPoleOddFamilyEmpty.** For every \(\lambda\ne0\),
every integer \(m\ge1\), and every positive odd height \(a\), the
normalized endpoint

\[
                  U(z)=\frac z{(1-\lambda z^2)^m}
 \tag{E1}
\]

has no pure polynomial pair. The conclusion holds after arbitrary
extension of the coefficient field.

**Proof.** Work geometrically and choose \(\alpha^2=1/\lambda\).
The only poles of \(U\) are \(\pm\alpha\), both of order \(m\).
Since the marked return is \(-z\), every admissible reciprocal unit
has the form

\[
 \eta(z)=\phi(z)^n,\qquad
 \phi(z)=\frac{\alpha-z}{\alpha+z},\qquad
 |n|\le1+\left\lfloor\frac m{a+2}\right\rfloor.
 \tag{E2}
\]

Indeed, its divisor is anti-invariant under reflection, supported on
the endpoint poles; its orders at the fixed points zero and infinity
are zero. The bound is the exact retained pole cap, applied to both
members of the pair.

Direct differentiation gives

\[
 U'(z)=\frac{1+(2m-1)\lambda z^2}{(1-\lambda z^2)^{m+1}}.
\]

Let \(\beta^2=-1/[\lambda(2m-1)]\). This is a simple ordinary
critical point, with finite nonzero endpoint value. The necessary full
jet condition includes \(\eta(\beta)=-1\). Thus \(n\ne0\), and
\(\rho=\phi(\beta)\) is a root of unity. Its trace is exactly

\[
                   \rho+\rho^{-1}=\frac{2(m-1)}m.
 \tag{E3}
\]

Here is an elementary integrality argument, valid over any
characteristic-zero coefficient field. Set \(P_0(T)=2\),
\(P_1(T)=T\), and \(P_{j+1}(T)=TP_j(T)-P_{j-1}(T)\).
Then \(P_j(\rho+\rho^{-1})=\rho^j+\rho^{-j}\), and \(P_j\)
is monic with integer coefficients for \(j\ge1\). If \(\rho\)
has order \(L\), its trace is a root of the monic integer polynomial
\(P_L(T)-2\). A rational root of a monic integer polynomial is an
integer: writing it as a coprime fraction and clearing the denominator
shows that this denominator divides a power of the numerator.

The rational number on the right of (E3) lies in \([0,2)\), so it
is either zero or one. Therefore \(m=1\) or \(m=2\). For these
two possibilities,

\[
 \begin{array}{c|c|c|c}
 m&\text{equation for }\rho&\operatorname{ord}\rho
     &\text{necessary }|n|\text{ for }\rho^n=-1\\\hline
 1&\rho^2+1=0&4&|n|\ge2\\
 2&\rho^2-\rho+1=0&6&|n|\ge3
 \end{array}
 \tag{E4}
\]

But the exact cap in (E2) is one for both \(m=1,2\), because
\(a+2\ge3\). This contradiction excludes the entire family.
\(\square\)

This argument uses the positions of the critical points, not only
their number. It does not exclude odd rational endpoints with the same
pole pair and additional finite zero pairs.

## 6. Exact Cayley criterion for one endpoint pole pair at height one

Work at height \(a=1\). A nontrivial odd-return core whose only finite
endpoint poles form one pair has no pole at infinity, by OddCoreInfinity.
After scaling over a coefficient field containing that pair, normalize
the poles to \(x=\pm1\), of common order \(m\ge1\). Thus
“only poles at \(\pm1\)” in this section means all poles on
\(\mathbf P^1\), not only the finite poles. Scaling transports the
endpoint as \(U(\alpha x)/\alpha\) and the coefficient pair as
\(C(\alpha x),D(\alpha x)\), so normalization is preserved.

The reciprocal unit cannot be constant. If it were one, then
\(X=U/(xU')\) could have finite poles only at ordinary critical
points; polynomiality of \(C\) excludes those poles. OddCoreInfinity
then makes \(X\) constant, forcing the trivial endpoint \(U=x\).
After replacing \(x\) by \(-x\) in the pair if needed, write

\[
 r=\frac{1+x}{1-x},\qquad \eta=r^n,\qquad n\ge1,
 \qquad F(r)=U\!\left(\frac{r-1}{r+1}\right).
\]

The endpoint is unchanged by this reflection because it is odd, while
the reciprocal exponent changes sign.

**Theorem SinglePolePairCayleyCriterion.** Such a core exists with
specified positive integers \(m,n\) if and only if there is a monic
reciprocal divisor \(A(r)\mid r^n+1\), with \(A(-1)\ne0\),
satisfying the following conditions. Put
\(d=\deg A=2s\) and \(k=2m-d\). Then \(k\ge0\),

\[
                       [r^m](r+1)^kA(r)=0,
 \tag{C1}
\]

and

\[
                            m\ge3(n-1).
 \tag{C2}
\]

Divisors in this equivalence may have coefficients in the chosen
characteristic-zero field; they need not be rational-coefficient
divisors. “Reciprocal” means \(A(r)=r^dA(1/r)\).

**Construction and sufficiency.** Since the only inversion-fixed
roots are \(\pm1\), neither of which is allowed in \(A\), its
roots occur in inverse pairs. Thus \(d=2s\), \(A(0)=1\), and
\(A(1)\ne0\). Define

\[
 F'(r)=c\,r^{-m-1}(r+1)^kA(r),\qquad
 c=\frac1{2^{k+1}A(1)},\qquad F(1)=0.
\]

The differential on the right is a Laurent polynomial. Its only
possible logarithmic term in a primitive is the \(r^{-1}\) term;
its coefficient is zero precisely by (C1). Hence its normalized
primitive \(F\) is rational, in fact a Laurent polynomial, with
exact pole order \(m\) at both zero and infinity. The latter
assertion follows from the nonzero constant and leading coefficients
of \((r+1)^kA\), whose degree is \(2m\).

Reciprocity gives \(F'(1/r)=r^2F'(r)\), so
\(F(1/r)=-F(r)\) after fixing \(F(1)=0\). Therefore the
reconstructed \(U(x)=F((1+x)/(1-x))\) is odd. Since
\(dr/dx=2\) at zero, the stated scalar makes \(U'(0)=1\).
At \(r=-1\), the primitive is regular and inversion fixes the
point, so \(F(-1)=0\); its zero order there is exactly \(k+1\).
This also directly gives \(U(\infty)=0\).

The coefficients and selected root are

\[
 C(x)=\frac{1+r^n}{(r+1)^2F'(r)},\qquad
 D(x)=-\frac{U(x)^2C(x)^3}{x^2},\qquad X=\frac{UC}{x}.
\]

The finite-source checks are exhaustive:

| Cayley point | Source point | Orders of \(C,D\) |
|---|---|---|
| \(r=0\) | \(x=-1\) | \(m+1,\ m+3\) |
| \(r=\infty\) | \(x=1\) | \(m+1-n,\ m-3(n-1)\) |
| \(r=1\) | \(x=0\) | Both zero orders; \(C(0)=1,D(0)=-1\). |
| A root of \(A\) other than \(-1\) | An ordinary finite source point | The factor cancels in \(C\) because \(A\mid r^n+1\); both coefficients are regular. |
| All other \(r\ne-1\) | All remaining finite source points | Both coefficients are regular. |

The point \(r=-1\) is source infinity and imposes no polynomiality
restriction. Condition (C2) is exactly the nonnegativity of the last
potentially negative order of \(D\); it also gives
\(m+1-n\ge2(n-1)\ge0\). Thus \(C,D\) are polynomials.
Oddness and reciprocity prove the all-order endpoint identity, or
equivalently the formula for \(C\) gives (R9). The pair has the
prescribed poles and is nontrivial.

**Necessity.** Given such a core, its endpoint \(F\) is a Laurent
polynomial with poles of exact order \(m\) at zero and infinity.
Consequently
\(F'(r)=c r^{-m-1}P(r)\), where \(P\) is monic of degree
\(2m\) with nonzero constant term. Oddness of the endpoint makes
\(P\) reciprocal. The height-one endpoint identity gives exactly
the displayed formula for \(C\). Every root of \(P\) except
possibly \(-1\) must therefore be cancelled by \(r^n+1\).
This polynomial is squarefree in characteristic zero, so such roots
of \(P\) are simple and form a reciprocal divisor \(A\) of
\(r^n+1\). Writing \(P=(r+1)^kA\) gives precisely the indicated
\(k\ge0\), even degree \(d\), and \(A(-1)\ne0\).
Normalization fixes \(c\). Rational integration forces (C1),
and the exact order at \(x=1\) forces (C2). No other factor is
discarded: roots of \(r^n+1\) not used in \(A\) remain roots
of \(C\). This includes the extra factors mentioned in Section 4.
\(\square\)

## 7. All pole orders for a fixed divisor, and the exact finite diagnostic

**Theorem AllPoleOrdersForFixedDivisor.** For a fixed reciprocal
divisor \(A(r)=\sum_{j=0}^{2s}a_jr^j\) as above, the coefficient
condition (C1) for integers \(m\ge s\) is equivalent to a single
nonzero polynomial equation \(F_A(m)=0\), where

\[
 F_A(m)=\sum_{j=0}^{2s}a_j
     \prod_{i=0}^{|j-s|-1}(m-s-i)
     \prod_{i=|j-s|+1}^{s}(m-s+i).
 \tag{C3}
\]

Empty products are one. The degree is exactly \(s\), with leading
coefficient \(A(1)\ne0\). In particular, a fixed divisor has only
finitely many candidate pole orders. For rational-coefficient \(A\),
exact rational polynomial factorization and the integer roots of its
linear factors decide every integer \(m\), with no pole-order cutoff.

**Proof.** Expand (C1) as
\(\sum_j a_j\binom{2m-2s}{m-j}=0\), interpreting a binomial
coefficient outside its range as zero. Divide by the nonzero central
coefficient \(\binom{2m-2s}{m-s}\). If \(\ell=|j-s|\),
the exact factorial ratio is

\[
 \frac{\binom{2m-2s}{m-j}}{\binom{2m-2s}{m-s}}
  =\frac{\prod_{i=0}^{\ell-1}(m-s-i)}
          {\prod_{i=1}^{\ell}(m-s+i)}.
\]

This holds also for out-of-range numerators: the falling product
contains a zero. Multiplication by
\(\prod_{i=1}^{s}(m-s+i)\), a positive integer for \(m\ge s\),
gives (C3) without creating or losing a root in the allowed range.
Each summand has degree \(s\) and leading coefficient \(a_j\);
their sum has leading coefficient \(A(1)\), which is nonzero because
\(A\mid r^n+1\). The case \(s=0\) is the nonzero constant
polynomial. \(\square\)

**Computational Proposition BoundedRationalDivisorDiagnostic.** For

\[
                       n\in\{12,15,18,20,24,30\},
\]

every monic divisor over \(\mathbb Q\) of \(r^n+1\) that excludes
the factor \(r+1\) fails to give an integer \(m\ge s\) satisfying
both (C1) and (C2). There are exactly **44 divisors** in this domain.

The exact checker factors each \(r^n+1\) over \(\mathbb Q\), removes
\(r+1\) when present, and enumerates every subset product of the
remaining irreducible factors. Those factors are reciprocal of even
degree, so this enumerates precisely the claimed divisor domain.
It forms (C3), computes its rational factorization and all integer
roots exactly, and checks the two integer inequalities. The receipt
retains the actual divisors, polynomials, factorizations, and roots.
Thus the conclusion covers every pole order for these inputs, rather
than only an initial interval of pole orders.

The qualifier on the divisor domain is essential. This diagnostic
does not cover other exponents \(n\), divisors with nonrational
coefficients, multiple endpoint pole pairs, or heights other than one.
For example, the exact control

\[
 A=\Phi_{18}=r^6-r^3+1,\qquad
 F_A=m^3-21m^2+92m-120
\]

occurs at \(n=9\), outside the stated list. It has no integer root,
but its values at \(15,16\) are \(-90,72\), so it has a real root
beyond \(n-1=8\). A proposed extension based on placing every real
root below \(n-1\) would therefore be false.

## 8. Sources, verification scope, and the remaining existence question

The canonical reflection and fixed-coordinate inputs supply the unit
factorization, full local polynomiality conditions, and all-order
endpoint identity. The sharp minimum-degree equality classification is
used only to exclude finite poles of a nontrivial inner map and to
strengthen the last infinity inequality. The residue criterion,
root-of-unity trace argument, Cayley construction, and coefficient
clearing are proved here by exact identities.

The actual author-hosted primary text of Adrianov, Pakovich and
Zvonkin, *Davenport–Zannier Polynomials and Dessins d'Enfants*, was
inspected: [book PDF](https://pakovich.github.io/book.pdf), Section 6.5.6,
Propositions 6.10–6.11 (printed pages 71–72), and Section 6.6.3
(printed pages 76–77). These give related Padé forms and differential
relations for their specified Davenport–Zannier families. They do not
supply a theorem deciding (C1)–(C2). This comparison is context only;
no Padé theorem or historical-priority assertion is imported.

The companion [exact checker](verify_odd_return_cores.py),
[source ledger](odd_return_cores_source_ledger.json), and
[verification receipt](odd_return_cores_verification.json) record the
fixed rational-divisor diagnostic and scoped controls for the pullback
identities, residue criterion, repeated-pole exclusion, and Cayley
integrability versus polynomiality. The written universal results are
not deductions from sampling these controls. A successful replay is
an internal exact computation, not external review or formal verification.

The remaining question is precise: **does any nontrivial normalized
pure polynomial pair at a positive odd height have an odd rational
endpoint?** If one exists, its identity reflection quotient already
gives a counterexample to universal no-polynomial-direction rigidity.
Conversely, any hypothetical nontrivial survivor in the unresolved
no-direction branch descends to such a core at the same height, with
a polynomial inner map. The core must have \(U(\infty)=0\) and
rational nonpolynomial \(X\) regular at infinity. The repeated-pole
family above is impossible; at height one and one endpoint pole pair,
the remaining problem is exactly (C1)–(C2), allowing algebraic divisors
and unbounded \(n\). The 44 rational-divisor cases do not resolve
that question.
