# A canonical rational reflection quotient and explicit finite fibres

Fix a characteristic-zero field \(K\) and a rational coordinate
\[
 u\in K(x),\qquad u(0)=0,\quad u'(0)=1,\qquad N=\deg u.
 \tag{Q1}
\]
The odd-height fibre at \(a\ge1\) consists of the normalized polynomial
pairs \(C(0)=-D(0)=1\) for the exact phase
\[
 p=t^a(C+t^2D)/x,\qquad t^2+2xy=1,
 \tag{Q2}
\]
with all positive pure moments zero and selected endpoint coordinate
\(u\). All heights in this note are positive and odd.

The previous explicit pole-twist classification assumed that a
nonconstant polynomial reciprocal unit existed. We remove that
assumption. Every coordinate has a canonical **rational** reflection
quotient. It is extracted from a linear system with \(2N+1\) numerator
coefficients and a univariate gcd. In this quotient all pole-admitted
units form a direct product of integer intervals, with at most one
unbounded direction. This gives a smaller finite search and a uniform
bound on the size of every finite fibre.

A critical-divisor count supplies a further rigidity statement: without
a polynomial reciprocal direction, a surviving polynomial \(X\) must
come from \(\eta=1\) and polynomial \(u/(xu')\). Consequently the
eventual finite branch is exactly the already classified singleton
case; any additional solutions are confined to finitely many low heights.

The result classifies the fibre of a specified rational coordinate;
it does not give a degree atlas of all endpoint coordinates or classify
arbitrary sphere phases. Classical field and decomposition methods are
acknowledged in Section 7. No historical priority is asserted.

## 1. The canonical reflection field

Let \(v(s)=s+O(s^2)\in K[[s]]\) be the formal inverse of \(u\), and
let \(\sigma\) act on \(K((s))\) by replacing \(s\) with \(-s\).
Embed \(K(x)\) into \(K((s))\) by \(x=v(s)\), and put
\[
 \mathcal F=K(v(s)),\qquad
 E=\mathcal F\cap\sigma(\mathcal F).
 \tag{Q3}
\]
Both fields in the intersection contain \(K(s)\), since
\(u(v(s))=s\). Thus \(E\) is nonconstant and \(\sigma(E)=E\).
The involution is nontrivial on \(E\), because it sends \(s\) to
\(-s\). The identification with \(K(x)\) will be used without
changing notation.

**Theorem CanonicalReflectionQuotient.** There is a unique rational
function \(h\in K(x)\) satisfying all four conditions
\[
 E=K(h),\qquad h(0)=0,\qquad h'(0)=1,
                  \qquad \sigma(h)=-h.
 \tag{Q4}
\]
There is a unique \(U\in K(z)\) with
\[
 u=U(h),\qquad U(-z)=-U(z),\qquad U(0)=0,
                  \qquad U'(0)=1.
 \tag{Q5}
\]
Writing \(d=\deg h\) and \(D=\deg U\), one has \(N=dD\).
Every normalized rational reciprocal unit of the marked return has
the form
\[
 \eta=V(h),\qquad V(0)=1,\qquad V(z)V(-z)=1.
 \tag{Q6}
\]
Conversely every function in (Q6) is such a unit. The same \(h\)
generates the reflection field after algebraic closure of \(K\).

Here the marked return is the formal branch \(j(x)\) satisfying
\(u(j(x))=-u(x)\), \(j(0)=0\). It need not be rational on the
original \(x\)-line. Equivalently, reciprocity is imposed on the
irreducible component of \(u(x)+u(y)=0\) through \((0,0)\).

**Proof.** The field \(E\) is an intermediate field between \(K\)
and \(K(x)\), so Lüroth's theorem gives \(E=K(w)\). The marked
place is \(K\)-rational. It is unramified under \(x\mapsto w\):
the zero of \(u\) at \(x=0\) is simple, and local orders multiply
in \(u\in K(w)\subset K(x)\). Choose the generator with
\(w(0)=0,w'(0)=1\). An automorphism of \(K(w)\) is fractional
linear; this also follows from the rational-map degree formula.
The involution fixes the marked place and has derivative \(-1\)
there, so for a unique \(c\in K\)
\[
                     \sigma(w)=-\frac{w}{1+cw}.
 \tag{Q7}
\]
The generator \(h=w/(1+cw/2)\) satisfies (Q4).

Any other generator with the same involution is a fractional-linear
function of \(h\). A fractional-linear odd function preserving zero
is a scalar multiple of \(h\); its derivative normalization fixes
that scalar. The requirement that \(h\) **generate** \(E\) is
essential: \(h+h^3\) would otherwise meet the displayed local and
parity conditions as well.

Since \(u\in E\), composition gives (Q5). Its oddness and
normalization follow from those of \(u,h\). Degrees multiply.
If \(\eta\) is reciprocal, \(\sigma(\eta)=\eta^{-1}\) belongs
to \(\mathcal F\), hence \(\eta\in E\). This proves (Q6);
the converse follows by substitution and the marked formal identity.
Equality on this formal branch is equality in its irreducible
function field.

For descent, construct the normalized generator \(\bar h\) over
\(\bar K\). Every element of \(\operatorname{Gal}(\bar K/K)\)
preserves the two fields, the involution and the marked normalization.
Uniqueness fixes \(\bar h\), so its normalized rational-function
coefficients lie in \(K\). It belongs to \(E\), while
\(E\subset\bar K(\bar h)\cap K(x)=K(\bar h)\). The last equality
follows by applying Galois invariance to the unique rational expression
in the transcendental variable \(\bar h\). Hence \(h=\bar h\).
No unproved commutation of field intersections with base change is
being assumed. \(\square\)

In particular, \(E\) is the largest subfield of \(K(x)\) on which
the marked formal reflection is a rational automorphism. This is a
statement about a quotient of the original line, not a rationality
claim for \(j\).

## 2. Extraction from a bounded linear space

Write \(u=P/Q\) in coprime form, with \(Q(0)\ne0\) and
\(\max(\deg P,\deg Q)=N\). Let \(H(x,y)\) be the irreducible
factor through \((0,0)\) of
\[
                      P(x)Q(y)+P(y)Q(x).
 \tag{Q8}
\]
That factor is unique: the full polynomial has nonzero gradient at
the marked point. In characteristic zero it is geometrically
irreducible. Indeed a smooth rational point lies on one geometric
component, which Galois must preserve; irreducibility over \(K\)
then excludes further components.

**Theorem BoundedLinearReflectionExtraction.** Form the vector space
\[
 W=\left\{\frac{A}{PQ}:\deg A\le2N,\quad
 H(x,y)\mid A(x)P(y)Q(y)+A(y)P(x)Q(x)\right\}.
 \tag{Q9}
\]
It is computed by linear algebra over \(K\), and \(K(W)=E\).
For a basis \(f_i=a_i/b_i\) in reduced form, let
\[
 G(T)=\gcd_i\bigl(a_i(T)-f_i(x)b_i(T)\bigr)
                       \quad\text{in }K(x)[T],
 \tag{Q10}
\]
with the gcd made monic. It is the minimal polynomial of \(x\)
over \(E\). Any nonconstant coefficient of \(G\) generates \(E\).
Normalize such a coefficient to a generator \(w\) with
\(w(0)=0,w'(0)=1\). Then the exact output is
\[
 h=\frac{w}{1+bw},\qquad
                   b=\frac{w''(0)-u''(0)}2.
 \tag{Q11}
\]
The marked factorization, this finite linear system and the gcd give
a terminating exact extraction over presented number fields. The
dimension bound is on the linear system; no polynomial runtime or
practical bound on the entire computation is claimed.

**Proof.** Since \(U\) is odd, its value at the fixed point
\(z=\infty\) is \(0\) or \(\infty\). Its order there is a
nonzero odd integer. At a pole of \(h\), local orders in \(u=U(h)\)
therefore give
\[
          -\operatorname{ord}_\beta h
                         \le |\operatorname{ord}_\beta u|.
 \tag{Q12}
\]
Thus \(h\) lies in the space of rational functions with pole divisor
bounded by \(|\operatorname{div}u|\). The finite part of that divisor
is the zero divisor of \(PQ\); its infinity coefficient is
\(|\deg P-\deg Q|\). The space is consequently exactly
\[
 \left\{A/(PQ):\deg A\le\deg P+\deg Q+
                       |\deg P-\deg Q|=2N\right\}.
 \tag{Q13}
\]
Condition (Q9) is precisely anti-invariance on the marked curve,
after clearing nonzero rational denominators. It consists of linear
equations in the \(2N+1\) coefficients of \(A\). Every element
of \(W\) lies in \(E\), and (Q12) puts \(h\) in \(W\).
Therefore \(K(W)=E\).

To justify (Q10), the common roots are exactly the images of \(x\)
under the \(E\)-embeddings of \(K(x)\) into an algebraic closure
of \(K(x)\): equality of all \(f_i\) is equality on their generated
field. There are no denominator roots, since each \(a_i,b_i\) is
coprime. A nonconstant rational function is separable in
characteristic zero, so its polynomial at a transcendental generic
value has no repeated roots. The monic gcd thus equals the monic
minimal polynomial, with no multiplicities lost or introduced.

Write a Lüroth generator as \(r=A_0/B_0\), of degree \(d\).
The minimal polynomial is
\[
 \frac{A_0(T)-r(x)B_0(T)}{(A_0)_d-r(x)(B_0)_d}.
 \tag{Q14}
\]
Each coefficient is fractional linear in \(r\). At least one is
nonconstant, since otherwise \(x\) would be algebraic over \(K\).
Every nonconstant one therefore generates \(K(r)=E\). A pole at
zero can be removed by first replacing the chosen generator by its
reciprocal; translation and scaling then normalize \(w\).

Finally write \(w(v(s))=s+bs^2+O(s^3)\). Taylor expansion gives
\(b=(w''(0)-u''(0))/2\). Comparing \(w(v(-s))\) with (Q7)
gives \(c=2b\), proving (Q11). The use of the second jet here
determines an already proved fractional-linear involution; it is
not a finite-series inference of a general rational identity.
\(\square\)

Once \(h\) is known, \(U\) is obtained by solving the coefficient
equations in \(u B(h)=A(h)\), with degrees at most \(D=N/\deg h\).
The rational function \(U\) is unique, and direct substitution checks
the output. No complete enumeration of intermediate fields is needed.

## 3. The exact reciprocal pole box

Work geometrically over \(\bar K\). For a finite nonzero value
\(r\) on the \(h\)-line, let
\[
 F_r=\{\beta\in\mathbb A^1:h(\beta)=r\}.
 \tag{Q15}
\]
Only **finite** source preimages enter this definition. Define the
allowed pole order on the quotient by
\[
 q_r=
 \begin{cases}
 +\infty,&F_r=\varnothing,\\
 0,&F_r\ne\varnothing,\ U(r)\ne\infty,\\
 \left\lfloor m_r/(a+2)+1/e_r\right\rfloor,
       &F_r\ne\varnothing,\ U(r)=\infty,
 \end{cases}
 \tag{Q16}
\]
where \(m_r\) is the pole order of \(U\) and
\(e_r=\max_{\beta\in F_r}\operatorname{ord}_\beta(h-r)\).
There is at most one value with \(q_r=\infty\), namely a possible
nonzero finite \(r=h(\infty)\) whose entire fibre is \(\{\infty\}\).

Choose one representative \(\alpha\) of each pair \(\{r,-r\}\)
for which \(q_r>0\) or \(q_{-r}>0\). This is a finite list: apart
from the possible exceptional value, it consists of paired poles
of \(U\). Put
\[
 \phi_\alpha(z)=\frac{\alpha-z}{\alpha+z},\qquad
 V_{\mathbf n}(z)=\prod_\alpha\phi_\alpha(z)^{n_\alpha},
          \qquad -q_\alpha\le n_\alpha\le q_{-\alpha}.
 \tag{Q17}
\]
Infinite endpoints are interpreted as one-sided unbounded intervals.

**Theorem UniversalReciprocalPoleBox.** Formula (Q17), followed by
\(\eta=V_{\mathbf n}(h)\), is a bijection onto all normalized
reciprocal units satisfying the finite pole conditions of the
fixed-coordinate theorem. No integer-kernel condition remains.
The entire fibre consists exactly of those candidates satisfying
\[
               \eta\equiv-1\pmod{J_{a,u}}.
 \tag{Q18}
\]
The polynomial pairs are then
\[
 X=\frac{u}{2xu'}(1+\eta),\qquad
                  C=\frac{xX^a}{u},\quad D=-\frac{xX^{a+2}}u.
 \tag{Q19}
\]
Over \(K\), retain exactly the products \(V_{\mathbf n}\in K(z)\).
The geometric candidates need not descend individually.

Here \(J_{a,u}\) has order \(\lceil s/a\rceil-1\) at an
unmarked zero of \(u\) of order \(s\), order \(\operatorname{ord}u'\)
at a finite point with nonzero finite \(u\)-value, and order zero
elsewhere. Zero orders are omitted.

**Proof.** A rational \(V\) satisfying (Q6) has an anti-invariant
divisor under \(z\mapsto-z\). Its orders at the fixed points
\(0,\infty\) are zero, and its other zero and pole orders are
opposite on each pair. The factors in (Q17) have exactly those
divisors, are normalized to one at zero, and satisfy the reciprocal
identity. Divisor factorization and normalization prove uniqueness.

If \(V\) has a pole of order \(k\) at \(r\), its pullback has
order \(ek\) at a finite preimage of local degree \(e\). When
\(U\) has pole order \(m_r\), the exact retained allowance is
\[
 ek\le1+\left\lfloor\frac{e m_r}{a+2}\right\rfloor
 \quad\Longleftrightarrow\quad
 k\le\frac{m_r}{a+2}+\frac1e.
 \tag{Q20}
\]
The maximum over finite preimages gives (Q16). If \(U(r)\) is
finite, poles are forbidden at every finite preimage. If there is
no such preimage, no finite pole condition is imposed. The signs
in (Q17) follow because \(\phi_\alpha\) vanishes at \(\alpha\)
and has its pole at \(-\alpha\).

The remaining exact polynomiality condition is (Q18), by the retained
reciprocal-unit theorem; its roots avoid every candidate denominator.
It proves (Q19) and exhausts the fibre. Finally, since \(h\in K(x)\)
and substitution is injective, \(V(h)\) is Galois invariant exactly
when \(V\) is. This proves the descent assertion. \(\square\)

## 4. The sole possible polynomial direction

**Theorem PolynomialDirectionCriterion.** A nonconstant polynomial
reciprocal unit exists if and only if
\[
 c=h(\infty)\in K^*,\qquad h^{-1}(c)=\{\infty\}.
 \tag{Q21}
\]
When it exists, the canonical normalized generator is
\[
                         R=\frac{c+h}{c-h}.
 \tag{Q22}
\]
Every polynomial reciprocal unit is \(R^k\), \(k\ge0\).
The single unbounded interval in (Q17) is precisely this direction.
The function \(R\) is not a proper polynomial power.

**Proof.** A nonconstant polynomial unit descends to \(V(h)\).
Every pole of \(V\) must have its whole preimage at \(x=\infty\).
There is at most one such quotient point. It cannot be either fixed
point \(0,\infty\), because the divisor of \(V\) has order zero
there. It must therefore be the finite nonzero point \(c=h(\infty)\),
and its paired zero is \(-c\). Normalization gives exactly (Q22)
to a positive integer power. Conversely (Q21) makes (Q22) a polynomial
reciprocal unit with \(R(0)=1\). It generates \(E\), being fractional
linear in \(h\). If \(R=S^m\), \(m>1\), normalize \(S(0)=1\).
The marked identity \((S(x)S(y))^m=1\) forces \(S(x)S(y)=1\),
so \(S\in E\). But \(R\) generates \(E\), contradicting the
degree \(m\) extension \(K(S)/K(S^m)\). \(\square\)

For an exact finite computation, write \(h=A/B\) in coprime form.
Condition (Q21) says that \(\deg A=\deg B\), their leading-coefficient
ratio is \(c\ne0\), and \(A-cB\) is a nonzero constant. No
additional polynomial-generator search is required.

## 5. A uniform bound for every finite fibre

**Corollary UniformFiniteFibreBound.** Put \(D=\deg U=N/\deg h\).
If the geometric fibre at height \(a\) is finite, then
\[
              \#\operatorname{Fibre}_{a,u}
                       \le3^{\lfloor D/2\rfloor}
                       \le3^{\lfloor N/2\rfloor}.
 \tag{Q23}
\]
The bound is uniform in the odd height and in all coefficients of
the endpoint of the specified degree. It is an upper bound; sharpness
of the finite-fibre bound is not asserted.

**Proof.** First suppose there is no polynomial direction. All the
intervals in (Q17) are finite and arise from finite nonzero pole
pairs of \(U\). If their common pole order is \(m\), then
\[
 q_\alpha,q_{-\alpha}\le1+\lfloor m/(a+2)\rfloor\le m.
 \tag{Q24}
\]
Each pair contributes at most \(2m+1\le3^m\) choices; the last
inequality follows by expanding \((1+2)^m\). The total pole degree
of \(U\) gives \(2\sum m\le D\). Multiplication proves (Q23)
even before imposing the jet condition.

If there is a polynomial direction, remove its unbounded pair from
the box. Each remaining choice gives one forward orbit
\(\eta_0R^j\), \(j\ge0\), after shifting its exponent lower bound.
Two surviving exponents force \([R]\) to be a unit of finite order
modulo \(J\): both successful units are \(-1\) there, and their
ratio is a positive power of \(R\). A nonunit \([R]\) cannot have
two successes. Finite order would give infinitely many successes by
adding the order, so a finite fibre has at most one pair per orbit.
The same product bound for the remaining pole pairs proves (Q23).
If \(J=1\), every orbit survives and the fibre is infinite, consistent
with this argument. \(\square\)

The bound also applies to a finite fibre over \(K\) when its geometric
fibre is infinite: each twist orbit either descends completely or has
no \(K\)-points, and a descended infinite orbit stays infinite over
\(K\). Thus a finite \(K\)-fibre again has at most one survivor per
bounded choice and obeys (Q23).

### Rigidity of the stable branch without a polynomial direction

**Theorem NoPolynomialDirectionRigidity.** Assume (Q21) fails.

1. If a pair survives at any height and its \(X\) in (Q19) is
   polynomial, then \(U(z)=z\), \(\eta=1\), and \(b=u/(xu')\)
   is polynomial.
2. If \(q_\alpha,q_{-\alpha}\le1\) for every pole pair, the entire
   fibre is either empty or the singleton \(\eta=1\). That singleton
   exists exactly when \(b\) is polynomial and every finite unmarked
   zero of \(u\) has multiplicity at most \(a\).
3. Let \(M\) be the largest finite pole multiplicity of \(u\), or
   zero if there is no finite pole. Part 2 applies whenever
   \(a\ge\max(1,M-1)\). In particular, a coordinate with
   nonpolynomial \(b\) and no polynomial reciprocal direction can
   have solutions only at
   \[
                         1\le a\le M-2\le N-2.
   \tag{Q25}
   \]
   Every such hypothetical solution has nonpolynomial rational \(X\).
   The theorem does not assert that these transient solutions exist.

Thus the stable no-direction fibres are precisely the empty fibres
and the Davenport--Zannier singletons, including the separate \(u=x\)
case. The [sharp all-height stabilization theorem](ALL_ODD_HEIGHT_FIBRE_STABILIZATION.md)
supplies the uniform height after which their unit data remain fixed.

**Proof.** Let \(z_0\) be the number of nonfixed zero pairs of the
odd map \(U\), and \(t_0\) its number of nonfixed pole pairs.
The point zero is a simple zero, and infinity is either a zero or
a pole. Therefore, if \(Z,P\) count all distinct zeros and poles
on the projective line,
\[
                         Z+P-2=2(z_0+t_0).
 \tag{Q26}
\]
This is the total ramification multiplicity of the ordinary critical
points, those having finite nonzero \(U\)-value. To verify the count
directly, the differential \(dU\) has divisor degree \(-2\).
At a pole of order \(m\), its order is \(-m-1\); at a nonpole
of local degree \(f\), it is \(f-1\). The nonpole critical total
is consequently \(D+P-2\). Subtracting \(D-Z\), the contribution
at the zeros, gives \(Z+P-2\). These local formulas hold in
characteristic zero and include the point infinity.

Every ordinary critical point \(r\) has a finite preimage under
\(h\). Otherwise \(r\) would be a finite nonzero exceptional value
in (Q21), contradicting the hypothesis. If its local degree for
\(U\) is \(f\) and such a preimage has degree \(e\), the full
jet condition requires \(V(h)+1\) to vanish to order at least
\(ef-1\). Hence
\[
                   \operatorname{ord}_r(V+1)\ge f-1.
 \tag{Q27}
\]
This argument uses an actual finite source preimage and the full
jet, not just its value.

Suppose all poles of \(V\) are simple. Its divisor is supported
on the \(t_0\) nonfixed pole pairs, so \(\deg V\le t_0\).
The function \(V+1\) is not identically zero because \(V(0)=1\).
Its total zero degree equals its total pole degree, which is
\(\deg V\) when \(V\) is nonconstant and zero when \(V=1\).
Equations (Q26)--(Q27) force
\(2(z_0+t_0)\le t_0\). Thus \(z_0=t_0=0\). The simple zero
at zero is then the only zero, and infinity the only pole, so
\(U=z\). There is no remaining pole factor in (Q17), giving
\(V=1\).

For part 1, polynomial \(X\) forces every finite pole of \(\eta\)
to have order at most one: at a finite pole of \(u\), \(b\) has
a simple zero and \(\eta=2X/b-1\); elsewhere the candidate is
regular. Every nonfixed quotient pole of \(V\) has a finite
preimage, so its order is at most one as well. The preceding count
applies, and \(X=b\) is polynomial.

The caps in part 2 directly give simple poles of \(V\). For
\(\eta=1\), \(X=b\). Its only possible finite poles are ordinary
critical points of \(u\), where (Q18) would fail. At each unmarked
zero of \(u\), \(b\) has a simple zero, so polynomiality is exactly
\(a\ge\operatorname{ord}u\). This proves the stated iff condition.

Finally, if \(a\ge M-1\), every finite pole allowance for \(\eta\)
is one, since \(a+2>M\). It follows that the quotient caps are
at most one. Part 2 proves part 3. Equivalently, it is enough that
each nonfixed pole order \(m\) of \(U\) obey \(m\le a+1\),
by (Q16). No example of a transient nonconstant unit is inferred
from the possibility of larger caps. \(\square\)

The same proof gives a cheap emptiness certificate before enumerating
the box. If \(U\) is nonlinear and there is no polynomial direction,
then a necessary condition for any solution is
\[
 \sum_{\{\alpha,-\alpha\}}
         \max(q_\alpha,q_{-\alpha})\ \ge\ 2(z_0+t_0).
 \tag{Q28}
\]
Indeed the left side bounds \(\deg V\) above and (Q27) bounds it
below by the right side. A constant unit also fails the nonempty
ordinary critical divisor. This is a necessary condition, not a
converse or an existence theorem for low-height exceptions.

**Corollary SharperEffectiveFibreClassification.** Over a presented
number field, the following terminating route classifies the whole
fixed-height fibre. Extract \(h,U\) by Section 2 and form (Q17).
Without a polynomial direction, enumerate its finite Cartesian box
and test (Q18). With a polynomial direction, test each forward orbit
using the retained elementary full-jet exponent algorithm. Exact
descent gives the original-field fibre.

This replaces the old integer-kernel/Cramer-box construction and the
separate bounded polynomial-generator search. It does not import a
general Orbit theorem. The quotient extractor is implemented in the
accompanying checker; the complete arbitrary-number-field fibre
classifier is still a proved procedure, not an end-to-end implementation.

## 6. Small examples and invalid shortcuts

* For \(u=x\), the quotient is \(h=x\), \(U=z\). Infinity is a
  fixed quotient point, so there is no polynomial reciprocal direction.
  There are no finite nonzero pole pairs, and \(\eta=1\) gives the
  complete singleton at every height.
* For \(u=x/(1-x^2)\), one has \(h=x\), \(U=z/(1-z^2)\).
  The one simple pole pair gives exactly three candidates
  \(((1-x)/(1+x))^n\), \(n=-1,0,1\). All fail the jet
  congruence modulo \(1+x^2\); the fibre is empty. The full
  correspondence has another component, which must not replace the
  marked component \(x+y=0\).
* For \(u=x/(1+x)\), one has \(h=u\), \(U=z\), and
  \(c=1\) has only the preimage infinity. Formula (Q22) gives
  \(R=1+2x\), recovering the full infinite degree-one fibre.
* The polynomial coordinate \(u=x+x^2\) has a nonrational original
  return, yet \(h=u\), \(U=z\). Its sole candidate \(\eta=1\)
  fails the ordinary critical jet. Existence of a rational quotient
  must not be mistaken for rationality of the original return.
* For
  \[
     u=\frac{x(1-x)^2}{(1-3x/2)^3},
  \]
  the extraction gives \(h=u\), \(U=z\), and no polynomial
  direction. Thus the sole candidate is \(\eta=1\). It contributes
  exactly at odd heights \(a\ge3\), as direct reconstruction shows.
  This is a complete finite fibre with a nonpolynomial quotient.

The checker also tests a nontrivial polynomial direction with a
nonrational original return. Incorrect quotient normalization, a wrong
correspondence component, and omission of a denominator or ramification
condition must not pass as valid extraction or fibre evidence.

## 7. Sources and evidence scope

The essential external field input is Milne, *Fields and Galois Theory*,
version 5.10, Theorem 9.19, Lemma 9.24 and the proof on printed pages
117--120. These include Lüroth, the rational-map degree formula, and
the fact that a nonconstant coefficient of the minimal polynomial
generates the intermediate field. The retained primary PDF was
inspected for these exact assertions.

The closest algorithmic comparison is Allem, Capaverde, van Hoeij and
Szutkoski, *Functional Decomposition Using Principal Subfields*,
ISSAC 2017, DOI 10.1145/3087604.3087608,
[author PDF](https://www.math.fsu.edu/~hoeij/2017/FD_ISSAC_corrected_final.pdf),
Section 2 and Theorem 5.1 on PDF page 7. The actual statements and
proof were inspected. Intermediate-field decomposition and extracting
a Lüroth generator from a minimal-polynomial coefficient are classical.
Their full subfield-lattice algorithm and complexity bounds are not
imported here. The finite anti-invariant space (Q9), its normalization,
and its application to exact pole boxes are proved above.

The retained
[fixed-coordinate theorem](ODD_START_FIXED_COORDINATE_FIBRE_CLASSIFICATION.md)
supplies the all-order moment/reciprocity interface, full-jet
polynomiality and reconstruction. The
[elementary exponent solver](ELEMENTARY_RECIPROCAL_EXPONENT_DECISION.md)
is needed only for the effective polynomial-direction corollary.
Neither the old Cramer-box decision route nor its imported Orbit
theorem is a premise of the new quotient and finite-box proofs.

Finite exact controls test the implemented extractor and displayed
examples. They do not prove universal extraction, exhaustion or the
cardinality bound by enumeration. The written arguments establish
those statements. No formal verification, external human review or
publication is inferred from internal review and replay.
