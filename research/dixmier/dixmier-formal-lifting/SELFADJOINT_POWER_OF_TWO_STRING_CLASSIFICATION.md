# Selfadjoint power-of-two string operators reduce to first-order mates

Let \(F\) be a one-variable algebraic function field over its actual
characteristic-zero constant field \(C=\ker\partial\), with nonzero
derivation \(\partial\). Work in finite scalar differential operators
\[
 F[E;\partial],\qquad [E,f]=f'.
\]
The coefficient-fixing formal adjoint is \(f^*=f\), \(E^*=-E\), with
\((AB)^*=B^*A^*\). This is an algebraic antiinvolution; no Hilbert-space
adjoint, reality condition, or positivity assumption is involved.

The imported input is
[CoprimeAlgebraicStringExclusion](ALL_COPRIME_ALGEBRAIC_STRING_OBSTRUCTION.md):
two finite scalar operators over this field with a nonzero constant
commutator cannot both have orders at least two that are coprime.

## 1. Complete classification in the selfadjoint class

**Theorem PowerTwoSelfadjointStringClassification.** Let \(P\) be monic
and formally selfadjoint, of order \(n=2^e\), \(e\ge1\). The following are
equivalent:

1. Some finite \(Q\) satisfies \([P,Q]=\kappa\in C^\times\).
2. There are a monic \(f\in C[z]\), of degree \(n/2\), and \(w\in F\)
   such that
   \[
   P=f(E^2)+w,\qquad w'=\gamma\in C^\times.                 \tag{A1}
   \]

When these conditions hold,
\[
 \operatorname{Cent}_{F[E]}(P)=C[P],\qquad
 Q=A(P)-\frac{\kappa}{\gamma}E,\quad A\in C[z]             \tag{A2}
\]
is the complete set of mates with commutator \(\kappa\).
The constant term of \(f\) can be absorbed into \(w\).
The same assertions apply over every finite algebraic coefficient
extension, with its actual constant field.

Take any mate and put \(T=(Q-Q^*)/2\). Since \(P^*=P\) and constants
are fixed,
\[
 [P,Q^*]=-\kappa,\qquad [P,T]=\kappa,\qquad T^*=-T.
                                                               \tag{A3}
\]
A nonzero skewadjoint finite operator has odd order: its highest
coefficient changes by \((-1)^{\operatorname{ord}T}\) under adjoint.
An order-zero skewadjoint operator is zero. Every odd order at least
three is coprime to \(2^e\), so the imported theorem forces
\(\operatorname{ord}T=1\).

Write \(T=aE+b\). Skewadjointness gives \(b=a'/2\). The coefficient
of \(E^n\) in \([P,T]\) is \(na'\), because \(P\) is monic.
Thus \(a\in C^\times\), \(b=0\), and
\[
                  [P,aE]=-aP'=\kappa.                    \tag{A4}
\]
Every positive-order coefficient of \(P\) is constant, and its scalar
coefficient has derivative \(-\kappa/a\). Formal selfadjointness
removes all odd powers of \(E\). This proves (A1).
Conversely \(-\kappa E/\gamma\) is an actual finite mate of (A1).

### The centralizer conclusion uses the mate

Here is a direct proof of the centralizer assertion, including its
coefficient-field scope. A monic selfadjoint operator of even order has
zero \(E^{n-1}\) coefficient. Set \(W=P^{1/n}=E+O(E^{-1})\) in the formal
pseudodifferential ring and \(B_j=(W^j)_+\). The coefficient recursion
for \(W\) divides only by nonzero integers. For \(j\ge0\),
\[
                \operatorname{ord}[P,B_j]\le n-2.        \tag{A5}
\]
If a finite operator \(V\) has commutator with \(P\) of order at most
\(n-2\), its leading coefficient is constant: the coefficient in order
\(n+\operatorname{ord}V-1\) is \(n\) times its derivative.
Subtract that constant times the corresponding \(B_j\). The
commutator bound is preserved, and the order falls. Finite repetition
therefore expresses \(V\) as a finite constant linear combination of
the \(B_j\).

For (A1), all coefficients of every \(B_j\) belong to \(C[w]\):
the coefficients of \(P\) lie there, \(w'=\gamma\), and fractional-root
recursion uses only differentiation, multiplication and division by
nonzero numerical constants. In particular, any operator commuting
with \(P\), even if initially written over a finite algebraic
coefficient cover, lies in \(C[w][E]\) over the current constants.

Since \(w'\ne0\), the element \(w\) is transcendental over the actual
constant field: differentiating a separable algebraic equation over
that field would otherwise give \(w'=0\). Thus \(C[w][E]\), with
\([E,w]=\gamma\ne0\), is the indicated polynomial Weyl algebra.
Inside it the substitution
\[
                    w\longmapsto P=w+f(E^2),\qquad E\longmapsto E
                                                               \tag{A6}
\]
is an automorphism, with the displayed polynomial inverse
\(w=P-f(E^2)\). Since \([E,P]=\gamma\), expansion in powers of \(E\)
shows that an operator commuting with \(P\) has no positive \(E\)
power. Its centralizer is exactly \(C[P]\).
Subtracting the explicit first-order mate from \(Q\) proves (A2).

There is no assertion that an arbitrary selfadjoint operator has
polynomial centralizer: \(P=E^4\) already has the commuting operator
\(E\notin C[P]\). The existence of the nonzero scalar mate is essential.

## 2. Algebraic normalization and descent to the original constants

Now write the original differential ring as \(F[D;\delta]\).
Say that \(P\) admits a **selfadjoint first-order normalization** if,
after a finite algebraic coefficient extension \(\widetilde F/F\),
there are \(\rho\ne0,\beta\in\widetilde F\) for which
\[
 E=\rho D+\beta,\qquad [E,g]=\rho\delta(g),                 \tag{A7}
\]
makes \(P\) monic and fixed by the coefficient-fixing adjoint
\(E^*=-E\). This is a stated hypothesis on \(P\), not a normalization
theorem for arbitrary even-order operators. No spatially polynomial
gauge transformation is asserted.

**Corollary PowerTwoNormalizedMateDescent.** If \(P\) has order \(2^e\)
and admits (A7), then every mate \(Q\in F[D;\delta]\) has the form
\[
                 Q=A(P)+Q_1,\qquad
                 A\in C[z],\quad \operatorname{ord}_D Q_1=1,
                 \qquad C=\ker\delta.                    \tag{A8}
\]
Its centralizer in the original ring is \(C[P]\).

The theorem over \(\widetilde F\) first gives (A8) with constants in
the extended field and a first-order remainder in \(E\), hence in
\(D\). If a positive polynomial term of highest degree \(j\) survives,
its coefficient is
\[
               \operatorname{lc}_D Q/
                        (\operatorname{lc}_D P)^j.       \tag{A9}
\]
It belongs to \(F\), since the original operators do. Its
\(\rho\delta\)-derivative vanishes, so its \(\delta\)-derivative
vanishes; it therefore belongs to the original \(C\).
Subtract that coefficient times \(P^j\) and repeat. Absorb the
constant polynomial term into \(Q_1\). This gives (A8) in the
original ring. The same leading-coefficient descent proves the
centralizer claim.

Every ordinary formally selfadjoint \(P\in F[D;\delta]\) of even
order admits (A7). If its leading coefficient is \(p_n\), adjoin
\(\rho^n=p_n\) and use
\[
                        E=\rho D+\delta(\rho)/2.         \tag{A10}
\]
For the original adjoint \(D^*=-D\), this \(E\) satisfies \(E^*=-E\);
the resulting \(P\) is monic and remains selfadjoint.

The derivative and time in (A1) belong to the normalized ring:
\(w'=\rho\delta(w)\). Neither (A1) nor (A7) asserts an algebraic
primitive of one for the original \(\delta\). The descent argument
does not use regularity of \(\rho,\beta\) in a deformation parameter.

## 3. The exact quartic boundary and all completed squares

For a monic depressed quartic
\[
                  P=E^4+u_2E^2+u_3E+u_4,
\]
direct adjoint expansion gives
\[
 P^*-P=2(u_2'-u_3)E+(u_2''-u_3').
                                                               \tag{A11}
\]
Thus selfadjointness is exactly \(u_3=u_2'\). On this entire locus a
finite scalar mate exists exactly when \(u_2\) is constant,
\(u_3=0\), and \(u_4'\) is a nonzero constant.

**Corollary CompleteSquareMateCriterion.** Over the stated algebraic
one-variable differential fields, put
\[
                     P=(E^2+U)^2+W,\qquad U,W\in F.
                                                               \tag{A12}
\]
A finite mate with any prescribed \(\kappa\ne0\) exists, over \(F\)
or any finite algebraic coefficient cover, if and only if
\[
                         U'=0,\qquad W'\in C^\times.      \tag{A13}
\]
In that case every mate is
\[
                    A(P)-\frac{\kappa}{W'}E,\qquad A\in C[z],
                                                               \tag{A14}
\]
over the current constant field. This follows by expanding (A12):
\[
 u_2=2U,\quad u_3=2U',\quad u_4=U''+U^2+W.
\]

This removes both the growing-potential condition and the slow-\(W\)
condition from the earlier
[selfadjoint square obstruction](GENERAL_SELFADJOINT_SQUARE_INFINITY_OBSTRUCTION.md).
It allows all rational and algebraic coefficient shapes and arbitrary
finite mate order. It concerns the scalar selfadjoint structure, not
arbitrary quartic operators or arbitrary formal coefficient fields.

## 4. Specialization of a low-order remainder

**Lemma LowOrderRemainderSpecialization.** Let \(\mathcal O\subset F\)
be a discrete valuation ring, with maximal ideal \(\mathfrak m\),
such that
\[
          \delta\mathcal O\subseteq\mathcal O,\qquad
          \delta\mathfrak m\subseteq\mathfrak m.           \tag{A15}
\]
Suppose \(P,Q\in\mathcal O[D;\delta]\) and, in the fraction-field ring,
\[
 Q=\sum_{j\ge1}a_jP^j+V,\qquad
 a_j\in F,\qquad r=\operatorname{ord}V
       <s=\operatorname{ord}\overline P.                 \tag{A16}
\]
The sum is finite. Then all \(a_j\) and all coefficients of \(V\)
are regular. If \(t=\operatorname{ord}\overline Q>r\), then \(s\mid t\).

Choose a uniformizer \(\pi\). If some \(a_j\) had negative valuation,
multiply (A16) on the left by a positive power \(\pi^k\) that makes
their minimum valuation zero. The scaled \(V\) is regular, because
it is the difference of regular operators. By (A15), reduction is
a homomorphism of differential-operator rings. The scaled left
side reduces to zero, while the right side has a nonzero sum
\(\sum_j\overline{\pi^ka_j}\,\overline P^{\,j}\) plus an operator of
order at most \(r<s\). The summand with largest nonzero \(j\) has
order \(js\), and cannot cancel. This contradiction proves
regularity of the \(a_j\), and then of \(V\).
The same distinct-order comparison after unscaled reduction proves
the divisibility assertion.

The coefficients \(a_j\) need not be constants for this lemma. Nor
must the uniformizer be constant: it is only multiplied on the left.
A constant uniformizer together with
\(\delta\mathcal O\subseteq\mathcal O\) implies (A15), so the usual
\(h\)-adic setting is included.

The differential-ideal condition cannot simply be omitted.
In \(\mathcal O=C[h]_{(h)}\), with \(\delta h=1\), take
\(P=D^2+hD^4\) and \(Q=P^2\). Coefficientwise reduction gives
\(\overline P=D^2\), but
\(\overline Q=D^4+2D^5\), of order five. The non-differential maximal
ideal makes reduction fail to preserve products.

## 5. All selfadjoint power-of-two regular Long repairs are excluded

**Corollary PowerTwoSelfadjointLongObstruction.** In the exact-first-image
scalar Long chart, no regular finite repair can have either member of
generic order \(2^e\) admitting a selfadjoint first-order normalization
as in (A7).

The chart and classical pair are the retained ones:
\[
 [a,v]=h/v,\qquad
 S_0=4a-3ra^2+2v,\qquad
 T_0=(ra^3-a^2-va)/2.                                    \tag{A17}
\]
The parameters \(r,h\) are constants for the original derivation.
The regular scalar repair must satisfy \([S,T]=h\) and have precisely
these reductions, generically in \(r\). Its classical orders are
therefore two and three.

If \(S\) has the stated normalization and generic order, (A8) gives
\(T=A(S)+T_1\), with \(\operatorname{ord}T_1=1\), over the original
constant field. The preceding lemma, applied before any normalization,
would force \(2\mid3\). If \(T\) has the stated property, exchange
the operators and obtain \(3\mid2\). Both are impossible.

In particular all quartic repairs satisfying the normalized adjoint
criterion (A11) are excluded, with no coefficient-degree bound,
partner-order bound, unit-leading-coefficient assumption, or
regularity assumption on the generic normalization itself.
This does not exclude quartics with \(u_3\ne u_2'\), composite orders
having odd prime factors, changed first images, or matrix operators.

There is also a global rank-one consequence. If a polynomial Weyl
endomorphism has an image of power-of-two differential order admitting
(A7), (A8) supplies an actual first-order polynomial image:
\(Q_1=Q-A(P)\). The first-order generation argument in
[LowOrderDC1](LOW_ORDER_DC1_AUTOMORPHISMS.md), Section 1, then proves
that the endomorphism is an automorphism. The generic algebraic
normalization is not being treated as a polynomial Weyl automorphism;
the polynomial landing of \(Q_1\) is the explicit subtraction.

## 6. Countercontrols and the remaining composite gap

Composite-order string pairs do exist. For example
\[
 P=E^4+x,\qquad Q=P^2-E,\qquad [E,x]=1
                                                               \tag{A18}
\]
has orders four and eight and commutator one. The theorem classifies
it; it does not exclude composite pairs without the Long specialization.

Selfadjointness is also a genuine restriction. With
\[
 L=E^2+x,\qquad P=L^2+E,\qquad Q=L,
                                                               \tag{A19}
\]
one has \([P,Q]=1\), while \(P^*=L^2-E\ne P\).
Here \(u_2=2x,\ u_3=3\), so (A11) fails.
Replacing \(Q\) by \(Q+P^2\) gives a nonselfadjoint pair with both
orders composite. The top symbols of (A19) are common powers, but
\(P-L^2=E\ne0\): leading common powers do not supply an actual
operator factorization or a lower-order string pair.

For general composite orders, the imported coprime theorem leaves
the common-power leading-symbol cases open. The adjoint projection
works here because odd order is automatically coprime to \(2^e\).
If the selfadjoint order has an odd prime factor, that argument no
longer excludes every higher odd projected order. No head-chopping
or Euclidean reduction resolving that case is claimed.

## Sources and verification scope

The indispensable imported result is the retained
[all-coprime theorem](ALL_COPRIME_ALGEBRAIC_STRING_OBSTRUCTION.md),
with its inspected Beals--Sattinger forward Stokes input. Its
coefficient-field and finite-scalar hypotheses are retained.
The adjoint projection also appears in the earlier
[square proof](GENERAL_SELFADJOINT_SQUARE_INFINITY_OBSTRUCTION.md);
the present argument proves the larger power-of-two classification.

As a primary comparison, Previato--Rueda--Zurro,
[*Commuting Ordinary Differential Operators and the Dixmier Test*,
SIGMA 15 (2019), 101](https://doi.org/10.3842/SIGMA.2019.101),
Section 2.1, printed pages 4--5, distinguishes commutativity of a
centralizer from the stronger equality \(C[P]\), and Theorem 2.6
records the formal-root description. The actual statements and
the order-\((4,6)\) commuting example (1.1), printed page 2, were
inspected. No general polynomial-centralizer theorem or
operator factorization is imported from that paper; the equality
needed here follows directly from (A6).

The [source ledger](selfadjoint_power_two_string_source_ledger.json)
pins the written proof, checker and imported local sources. The
[checker](verify_selfadjoint_power_two_string.py) verifies exact
adjoint identities, positive mates, the quartic boundary,
normalization and specialization controls. These bounded checks
do not prove the theorem for all powers of two; the written
projection, coprime exclusion and descent arguments do.
No historical-priority, formal-verification or external-review claim
is made.

~~~sh
python3 research/dixmier/dixmier-formal-lifting/verify_selfadjoint_power_two_string.py \
  --repo-root . --output /tmp/selfadjoint_power_two_string_review.json
~~~
