# Quadratic-complement witnesses with smaller poles and sharp multiplier bounds

Let \(K\) have characteristic zero, let \(q\) be nondegenerate of dimension
\(N=M+2\ge3\), and write \(q=A+B\) as an orthogonal sum of ranks \(2,M\)
over \(K\). Write \(\Lambda_q\) for the normalized
\(\mathrm{SO}(q)\times\mathbb G_m\)-invariant trace on
\(K[V,q^{-1}]\). Set
\[
 S=B/q,\qquad W=1-(1+q)^2S .
\]
A simpler witness and its fixed multiplier are
\[
 F=q^{-M}(1+q)^M W^{M+1},\qquad H=q^M .
\]
For every \(m\ge1\) the pure moment \(\Lambda_q(F^m)\) is zero, while
\(\Lambda_q(HF^m)\) is one nonzero rational beta moment. Over an algebraic
closure, \(F\) is not a unit times a proper power. Its pole order along
\(q=0\) is \(2M+1\), compared with \(3(M+1)\) for the earlier
\(P^M(1+P)\) construction. In dimension three, \((1+q)W/q\) gives
the same conclusions with pole order \(2\).

The least degree of a polynomial multiplier with infinitely many nonzero
marked moments for this \(F\) is \(2M\). More precisely, for every degree
cap \(0\le D<2M\), put \(J=\lfloor D/2\rfloor\). Every polynomial
multiplier of degree at most \(D\) has zero marked moments at all orders
\[
                     m\ \ge\ 1+\left\lceil\frac{2J}{M}\right\rceil .
\]
This first order is sharp. A finite rational-function calculation also
decides eventual vanishing for every fixed regular multiplier.

## 1. Exact marked moments

The [rank-two trace calculation](QUADRATIC_COMPLEMENT_UNIFORM_NORMALIZER_COUNTEREXAMPLES.md#3-the-rank-two-beta-recurrence-on-x)
gives \(\Lambda_q(S^j)=M/(M+2j)\). In the same note, equation (7) shows
that scalar invariance extracts the Laurent constant coefficient in \(q\)
for polynomials in \(S\). The endpoint substitution is the method of
[Long, Section 5, Proposition 5.1](https://github.com/long-mathematics/mathieu-second-conjecture-spherical-pairs/blob/e9edac262b4b64eb68feb0e6e80496cefe4a1dee/mathieu_second_conjecture_spherical_pairs.tex).
Here both the angular multiplier and the exponent of the endpoint factor
can vary independently.

For integers \(\alpha,\gamma\ge1\) and \(\beta\ge0\), put
\[
 F_{\alpha,\beta,\gamma}=q^{-\alpha}(1+q)^\beta W^\gamma,\qquad
 C_{L,j}=M\int_0^1 w^{M+2j-1}(1-w^2)^L\,dw
 \quad(L,j\ge0).
\]
The integral is termwise polynomial integration over \(\mathbb Q\). It
is positive as a rational number and remains nonzero in \(K\).

**Proposition (mixed moment identity).**
For every \(m\ge1\), \(j\ge0\)
and integer \(\ell\ge(\alpha-\gamma)m\),
\[
 \Lambda_q(q^\ell S^j F_{\alpha,\beta,\gamma}^{\,m})
   =C_{\gamma m,j}\binom{\beta m-M-2j}{\alpha m-\ell}.                 \tag{1}
\]
The binomial coefficient is defined by expansion at zero, with value
zero for negative lower index. The endpoint range is sharp for this
formula.

Indeed, put \(T=1+U\), \(L=\gamma m\), \(n=\alpha m-\ell\). Scalar weights
and the angular moments express the trace as
\[
 \begin{aligned}
 [U^n]\,T^{\beta m} M\int_0^1
        t^{M+2j-1}(1-T^2t^2)^L\,dt
 &= [U^n]\,T^{\beta m-M-2j}J(T),\\
 J(T)&=M\int_0^T w^{M+2j-1}(1-w^2)^L\,dw .
 \end{aligned}
\]
This substitution can be performed termwise. The product on the right
is a polynomial because \(T^{M+2j}\) divides \(J(T)\). At \(T=1\),
\(J(T)-J(1)\) has initial term
\(M(-2)^L(T-1)^{L+1}/(L+1)\), by differentiating \(J\).
Thus it can be discarded for \(n\le L\), which gives (1).
For \(n=L+1\), the factor \(T^{\beta m-M-2j}\) has constant
coefficient \(1\) in \(U\), so the discarded coefficient is exactly
\(M(-2)^L/(L+1)\ne0\).

## 2. Reducing arbitrary multipliers

Let \(O(A)\times O(B)\) act by the two orthogonal blocks. For any of
these isometries, composing \(\Lambda_q\) with its pullback still gives
a normalized \(\mathrm{SO}(q)\times\mathbb G_m\)-invariant trace. This
uses normalization of \(\mathrm{SO}(q)\) by the isometry and uniqueness
of the trace. Each phase \(F_{\alpha,\beta,\gamma}\) is invariant under
the block group.

For a polynomial multiplier \(h\) of degree at most \(D\), take its
Reynolds projection \(\pi h\) in the degree-at-most-\(D\) polynomial
representation of the block group. The identity component is
reductive. Complete reducibility there and averaging over its finite
component group give this projection in characteristic zero. One may
work over an algebraic closure and descend the canonical projection,
since it commutes with scalar extension. The invariant linear
functional \(h\mapsto\Lambda_q(hF_{\alpha,\beta,\gamma}^{\,m})\)
vanishes on the complement of the invariant subspace. Consequently
\[
     \Lambda_q(hF_{\alpha,\beta,\gamma}^{\,m})
       =\Lambda_q\bigl((\pi h)F_{\alpha,\beta,\gamma}^{\,m}\bigr)
       \qquad(m\ge1).                                                \tag{2}
\]
For \(h=q^{-r}p\) with \(p\) polynomial, use the same projection on
\(q^{-r}K[V]_{\le\deg p}\). Equivalently, project \(p\) and multiply by
\(q^{-r}\), because \(q\) is invariant.

For clarity, the polynomial invariants of the block group are
\(K[A,B]\). Over an algebraic closure choose vectors \(v_A,v_B\) of norm
one in the two blocks. On the locus \(AB\ne0\), isometries send each
vector pair to \((s v_A,t v_B)\) with \(s^2=A,t^2=B\), by extending
the isometries between their norm lines by choosing orthogonal bases
on the nondegenerate complements. The restriction of an
invariant polynomial to this plane is even in \(s\) and \(t\), using
the sign symmetry in each block. Hence it is a polynomial in \(s^2,t^2\),
and equality on the dense locus gives the assertion. The coefficients
descend to \(K\) by the algebraic independence of \(A,B\). Thus the
averaged regular multipliers lie in \(K[q,q^{-1},S]\). A polynomial
of degree at most \(D\) averages to a linear combination of
\(q^\ell S^j=q^{\ell-j}B^j\) with
\(0\le j\le\ell\le\lfloor D/2\rfloor\). This uses homogeneity of \(q,B\)
and preservation of the polynomial degree filtration.

The projection can also be computed by linear algebra over \(K\).
Use \(E=q^{-r}K[V]_{\le D}\), where \(r\) clears the denominator
of the multiplier. The sum of the images on \(E\) of the block
rotation derivations is the complement to the invariants of the
identity component. Indeed, complete reducibility splits \(E\)
into simple modules and the Lie algebra of a connected group
acts nontrivially on every nontrivial simple module. In the quotient
by those images, the identity component acts trivially. Add the
images of \(s_A-1\) and \(s_B-1\), where \(s_A,s_B\) negate one
coordinate in the respective blocks. The two reflections generate
the component group, and the sum of their difference images is
the kernel of its averaging operator on the quotient. The image
of \(h\) in this final quotient has a unique representative in
\(q^{-r}K[A,B]_{\le D}\), using ambient polynomial degree in the
numerator. Solving for that representative gives \(\pi h\).

## 3. The witness and pole order

**Theorem (smallest pole in the indicated product family).**
Within \(F_{\alpha,\beta,\gamma}\) with \(\gamma\ge\alpha\), a phase is
pure exactly when \(\beta\ge M\) and \(\alpha\ge\beta\). Among pure
phases, a fixed regular multiplier detects infinitely many orders
exactly when \(\alpha=\beta\) or \(\alpha=\gamma\). Geometric
nonpower up to a unit is equivalent to \(\gcd(\beta,\gamma)=1\).
The least pole order among these nonpower witnesses is \(2\) for
\(M=1\) and \(2M+1\) for \(M\ge2\).

To prove purity, take \(j=\ell=0\) in (1). If \(\beta<M\), the upper
index is negative at \(m=1\) and its coefficient at \(\alpha\ge1\)
is nonzero. If \(\beta\ge M\), the upper index is always nonnegative.
The lower index exceeds it for all \(m\) exactly when
\(\alpha\ge\beta\).

For \(\alpha>\beta\) and \(\gamma>\alpha\), take any fixed term
\(q^\ell S^j\) of an averaged regular multiplier. At large \(m\) the
endpoint range holds and the lower index exceeds the nonnegative
upper index. Hence all moments eventually vanish. For
\(\alpha=\beta\ge M\), use \(j=0,\ell=M\). The lower and upper
indices in (1) are equal and nonnegative, so the marked moment is
\(C_{\gamma m,0}\) at every positive order.

The remaining case is \(\gamma=\alpha>\beta\ge M\).
Here \(\ell=-1\) selects the first coefficient past the endpoint
range. The part with \(J(1)\) vanishes since \(n=\alpha m+1\)
exceeds the nonnegative degree \(\beta m-M\). The exact omitted
coefficient computed in Section 1 gives
\[
     \Lambda_q(q^{-1}F_{\alpha,\beta,\alpha}^{\,m})
             =\frac{M(-2)^{\alpha m}}{\alpha m+1}
             \qquad(m\ge1).
\]
For polynomial multipliers in this case, \(\ell\ge0\) in the
averaged support places every term in the range of (1), and the
lower index exceeds the upper index for large \(m\).

For the power obstruction, restrict over an algebraic closure to
the punctured \(x_1\)-axis, where \(A=a_1x_1^2\) and \(B=0\).
There is a point with \(q=-1\) where its vanishing order is \(\beta\).
On the affine \(z_1\)-line with \(A=3/4\) and all other variables
fixed, choose a point with \(B=b_1z_1^2=1/4\). Both \(q\) and
\(1+q\) are nonzero there, and \(dW=-4\,dB\ne0\). The restriction
has a zero of multiplicity \(\gamma\). A unit restricts to a unit
on the open part of either line where \(q\ne0\). Hence a proper power
up to a unit would have an exponent dividing both multiplicities. If
\(\beta,\gamma\) have a common divisor \(e\ge2\), the expression is
a unit times an \(e\)-th power, since \(q\) is a unit on the complement.

Over an algebraic closure \(q\) is irreducible because a product of
linear forms has rank at most two. Modulo \(q\), the numerator
\(q-(1+q)^2B\) of \(W\) is \(-B\), which is nonzero since \(B\) is not
a scalar multiple of \(q\). Hence the exact pole order is
\(\alpha+\gamma\). For a pure phase \(\gamma\ge\alpha\ge\beta\ge
M\), so it is at least \(2M\). Equality forces all three exponents
to be \(M\), which is admissible and nonpower for \(M=1\). For
\(M\ge2\) it gives a proper power, and the choice
\((M,M,M+1)\) attains the next integer \(2M+1\). For comparison, the
earlier \(P=(1+q)W^2/q\) has pole order three, so
\(P^M(1+P)\) has pole order \(3(M+1)\).

## 4. The sharp degree bound

**Theorem (uniform polynomial multiplier cutoff).**
For the displayed \(F\) and the double-pole ternary witness, if
\(0\le D<2M\), all polynomial multipliers
of degree at most \(D\) vanish in the marked moments for
\(m\ge 1+\lceil2\lfloor D/2\rfloor/M\rceil\). This is the least
positive uniform cutoff. The minimum degree of a multiplier
detecting infinitely often is \(2M\).

In (1), take \(\alpha=\beta=M\) and either \(\gamma=M+1\) or,
when \(M=1\), \(\gamma=1\).
By (2) it suffices to use \(0\le j\le\ell\le J=\lfloor D/2\rfloor\).
When \(m\ge1+\lceil 2J/M\rceil\), the upper index
\(Mm-M-2j\) is nonnegative and the lower index \(Mm-\ell\) is
strictly larger, since \(\ell<M\). For \(J=0\) the asserted
cutoff is already \(1\). For \(J>0\), use \(h=B^J=q^J S^J\)
at \(m=\lceil 2J/M\rceil\). The upper index is then strictly
negative and the lower index is positive, so (1) is nonzero.
Finally \(q^M\) has degree \(2M\) and detects every order.
At degree at most \(2M\), the only averaged monomial which
contributes to the tail is \(q^M\), since \(\ell\ge M+2j\) and
\(\ell\le M\) force \(j=0,\ell=M\).

## 5. A test for a fixed regular multiplier

**Proposition (eventual moment criterion).**
Suppose \(\alpha=\beta=M\)
and \(\gamma>M\). For any regular multiplier \(h\), let
\(\pi h=\sum_{\ell,j}c_{\ell,j}q^\ell S^j\). Write \(a=M/2\) and
\((x)_j=x(x+1)\cdots(x+j-1)\), with \((x)_0=1\). Define
\[
 R_h(X)=
 \sum_{\ell\ge M+2j} c_{\ell,j}
 \frac{(a)_j}{(a+\gamma X+1)_j}
 \binom{MX-M-2j}{\ell-M-2j}\ \in K(X),                             \tag{3}
\]
where the binomial polynomial has a fixed nonnegative lower index.
For all sufficiently large \(m\), \(\Lambda_q(h F_{M,M,\gamma}^{\,m})
=C_{\gamma m,0}R_h(m)\). In particular eventual vanishing holds
exactly when \(R_h=0\). Otherwise the multiplier detects every
sufficiently large order.

For large \(m\), (1) applies to the finite support of \(\pi h\);
the upper and relevant lower binomial indices are then
nonnegative. Terms with \(\ell<M+2j\) vanish. For the surviving
terms, the complementary lower index is \(\ell-M-2j\).
Integrating the derivative of \(w^{M+2j}(1-w^2)^{L+1}\) on
\([0,1]\) gives
\[
  (M+2j) C_{L,j}=(M+2j+2L+2) C_{L,j+1}.
\]
Iterating gives \(C_{L,j}/C_{L,0}=(a)_j/(a+L+1)_j\), which yields
(3). A nonzero rational function has only finitely many zeros.
The criterion tests cancellations between surviving monomials
by clearing the finitely many denominators in (3).

For example, in dimension three with the uniform phase
\(F=(1+q)W^2/q\), the multiplier \(4q^4S+15q^3S-q\) has a
nonzero first marked moment but vanishes at every order \(m\ge3\).
In (3) the contribution from the first two terms is
\((4(X-3)+15)/(4X+3)=1\), and the last term contributes \(-1\).
At \(m=3\) formula (1) gives the same cancellation directly.

## 6. A pole reduction outside the product family

**Proposition (pole order \(2M\) for \(M\ge4\)).**
For \(M\ge4\) the phase
\[
  F_* = F_{M,M,M}\bigl(1+B(M-1+Mq)\bigr)
\]
has exact pole order \(2M\) and is geometrically nonpower up to
units. It is pure, \(q^M\) detects at every order, and \(2M\) is again
the least degree of a polynomial multiplier detecting infinitely
often.

In the expansion of the \(m\)-th power of the final factor, the
terms are \(B^j q^r=q^{j+r}S^j\), with \(0\le r\le j\le m\),
and coefficients \(\binom mj\binom jr (M-1)^{j-r}M^r\).
Apply (1) to \(F_{M,M,M}\) with \(\ell=j+r\). For \(m\ge2\),
the upper index \(Mm-M-2j\) is nonnegative because it is at
least \((M-2)m-M\ge M-4\). The pure lower index exceeds it by
\(M+j-r\), giving zero. With multiplier \(q^M\), the lower index
exceeds the upper one by \(j-r\). The terms \(r=j\) therefore
contribute positive rational numbers and the other terms vanish.
At \(m=1\) the unmarked terms with \(j=1\) have upper index \(-2\)
and the combination of binomials is
\[
 (M-1)\binom{-2}{M-1}+M\binom{-2}{M-2}=0.
\]
With multiplier \(q^M\), both added terms have negative lower index
at \(m=1\), leaving \(C_{M,0}\). For \(M=3\), this formula gives
\(\Lambda_q(F_*^2)=C_{6,2}\ne0\), so that dimension is outside this
construction.

The extra factor reduces to \(1+(M-1)B\) modulo \(q\), which is
nonzero since \(q\) and \(B\) are homogeneous quadratics and this
polynomial has constant term \(1\). Thus it preserves the exact
pole \(2M\) of the base phase. Over the algebraic closure, take a
point with \(q=1\) and \(B=-1/(2M-1)\). At this point the extra
factor vanishes, \(W=1-4B\ne0\), and \(dA,dB\) are independent
since \(A,B\) are nonzero in their distinct blocks. Hence \(dB\)
is nonzero on the tangent space to \(q=1\). The extra factor
has differential \((2M-1)dB\) on this tangent space, so its
zero is simple there. The base phase is nonzero at the point,
which excludes a proper power up to a unit.

For a polynomial multiplier of degree \(D<2M\), use (2) and
take an averaged term \(q^\ell S^{j_0}\), with
\(j_0\le\ell\le J=\lfloor D/2\rfloor\). In a term \(B^j q^r\)
from the power expansion, the lower index for the resulting
moment exceeds the upper index by \(M+2j_0-\ell+j-r>0\).
The upper index is nonnegative as soon as
\((M-2)m\ge M+2J\). This proves eventual vanishing for all such
multipliers. The marked identity for \(q^M\) proves the claimed
degree minimum.

The pole minimum of Section 3 ranges over the defined product family
with \(\gamma\ge\alpha\). The exact cutoff in Section 4 concerns the
uniform phase and its double-pole ternary variant. The rational-function
criterion applies to the displayed uniform phase, which has
\(\gamma=M+1\).
