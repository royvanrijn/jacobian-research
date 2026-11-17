# The exact scalar spectrum closes the homogeneous Euler compatibility model

**QQL1.** In the actual source Weyl division ring of QCS1, the scalar
eigenvalues of its transported Euler operator are exactly
\(5ch\mathbb Z\). In particular the eigenvalue \(-ch\) never occurs.

Consequently this homogeneous source and fixed cubic cotangent transport
admit no actual Euler factorization \(J=QP+\lambda\) with
\([P,Q]=ch\), and no fifth root of \(U\,G(J)\) for any nonzero rational
function \(G\). This closes all displayed homogeneous target orderings,
not only their finite Laurent realizations. Other quantum transports
and nonhomogeneous models remain open.

## The actual division ring and the primitive tenth-root divisor

Let \(k\) have characteristic zero and let \(h,c,c_v\in k^\times\).
In
\[
 \mathcal D=\operatorname{Frac}\bigl(k(v)[b;h\partial_v]\bigr),
 \qquad [b,v]=h,\qquad H=vb,
\]
choose an arbitrary \(F_s(H)\in k(H)^\times\), and define
\[
 \begin{split}
 \kappa&=c_vv^5F_s(H)^{-1},\\
 U&=\frac{(\kappa-3)^3}{(\kappa-2)^2},&
 R&=-\frac{(\kappa-3)(\kappa-2)}{\kappa^2},&
 J&=cRH.                                                   \tag{1}
 \end{split}
\]
In the Long application \(c_v=2/r\) and \(c=4/r^2\).
No formal parameter or regular classical specialization is needed
for the theorems below.

The relation \(Hv=v(H+h)\), together with \(b=v^{-1}H\), identifies
\(\mathcal D\) with the quotient division ring of the skew Laurent
ring \(k(H)[v,v^{-1}]\). Its least-exponent valuation extends to
the entire division ring, or equivalently one embeds in the skew
Laurent-series division ring. Every nonzero \(X\in\mathcal D\)
has a leading expression \(v^mA(H)\), with \(m\in\mathbb Z\)
and \(A(H)\in k(H)^\times\).

Here \(\nu_v(\kappa)=5\) and \(\nu_v(J)=-10\).
The coefficient of \(v^{-10}\) in \(J\) is
\[
 C(H)=-6c\,c_v^{-2}
                  F_s(H-5h)F_s(H-10h)H.                    \tag{2}
\]
Suppose \([J,X]=\rho X\) for a scalar \(\rho\in k\).
Comparing the lowest weight \(v^{m-10}\) gives
\[
                  C(H+mh)A(H)=A(H-10h)C(H).                \tag{3}
\]

Restrict divisors to the orbit \(h\mathbb Z\) in the \(H\)-line
over an algebraic closure of \(k\). Encode the signed finite
divisors of \(F_s,A,C\) on this orbit by Laurent polynomials
\(f(T),a(T),d(T)\) with integer coefficients. Formula (2) gives
\[
                  d(T)=1+(T^5+T^{10})f(T).
\]
Zeros and poles are both permitted; integral shifts of \(h\) do not
move other roots into this orbit. Equation (3) gives
\[
                  (T^{-m}-1)d(T)=(T^{10}-1)a(T).
\]
At a primitive tenth root of unity \(\xi\), one has \(d(\xi)=1\).
Thus \(\xi^{-m}=1\), proving the useful intermediate fact
\[
          [J,X]=\rho X,\ X\ne0
                    \quad\Longrightarrow\quad
                         \nu_v(X)\in10\mathbb Z.            \tag{4}
\]
This uses no finite Laurent support of \(X\), polynomial source
function, or positivity of the divisor multiplicities.

## Scalar eigenoperators lie in the quintic fixed division ring

One may first enlarge the scalar field to contain fifth roots of
unity. The original Weyl division ring embeds in the enlarged one;
no tensor product of division rings is being identified.
Every scalar eigenoperator over the original field remains one
after this extension.

Let \(\sigma(v)=\zeta_5v\), \(\sigma(b)=\zeta_5^{-1}b\).
It fixes \(H,\kappa,J\).
For a scalar eigenoperator \(X\), each deck component
\[
                  X_j=\frac15\sum_{i=0}^{4}
                                 \zeta_5^{-ji}\sigma^i(X)
\]
is again an eigenoperator with the same scalar eigenvalue.
If \(X_j\ne0\), its least \(v\)-weight is congruent to \(j\)
modulo five. Equation (4) forces \(j=0\).
Every scalar eigenoperator is therefore deck-fixed.

The fixed division ring is exactly
\[
 E=\langle\kappa,H\rangle_{\rm div}
   =\operatorname{Frac}\bigl(k(\kappa)[H;5h\kappa\partial_\kappa]\bigr).
                                                               \tag{5}
\]
For completeness, \(E\) contains \(v^5=c_v^{-1}\kappa F_s(H)\).
Conjugation by \(v\) sends \(H\) to \(H-h\) and sends \(\kappa\)
to \(\kappa F_s(H)F_s(H-h)^{-1}\); conjugation in the other direction
also preserves \(E\). Hence the span of \(1,v,\ldots,v^4\) over \(E\)
is a finite subring of \(\mathcal D\).
Multiplication by a nonzero element is injective and therefore
surjective on that finite-dimensional span, making it a division
ring. It contains \(v,b\), so it equals \(\mathcal D\).
The deck characters make the five basis terms independent and
identify the fixed ring with \(E\).

Since \(J=cR(\kappa)H\), (5) can also be written
\[
          E=\operatorname{Frac}\bigl(k(\kappa)[J;\delta]\bigr),
          \qquad
          \delta=5ch\kappa R(\kappa)\partial_\kappa.         \tag{6}
\]

## The scalar spectrum is exactly \(5ch\mathbb Z\)

Let \(0\ne X\in E\) satisfy \([J,X]=\lambda X\), \(\lambda\in k\).
Take ordinary differential order in \(J\), allowing negative orders
in the quotient division ring. Its leading coefficient is a nonzero
rational function \(f(\kappa)\). At this order the eigenoperator
equation gives
\[
                              \delta f=\lambda f.          \tag{7}
\]
This comparison applies to fractions as well: for \(X=AB^{-1}\),
the leading coefficient is \(\operatorname{lc}(A)/\operatorname{lc}(B)\);
the commutator contributes the derivative of that quotient.

At \(\kappa=\infty\), write \(f\sim a\kappa^n\), \(n\in\mathbb Z\).
Since \(R(\kappa)\to-1\), equation (7) forces
\(\lambda=-5ch n\).
Conversely the exact relation
\[
                              [J,U]=-5chU
\]
shows that \(U^n\) realizes \(-5ch n\) for every integer \(n\),
including negative integers. Thus, with the eigenvalue set defined
inside the actual division ring,
\[
 \boxed{\quad
 \{\lambda\in k:\exists\,0\ne X\in\mathcal D,\
                        [J,X]=\lambda X\}=5ch\mathbb Z .
 \quad}                                                     \tag{8}
\]
The necessary equality descends from the scalar extension used
for deck projection, and the realizers \(U^n\) already belong
to the original division ring.

## No actual Euler pair and no rational target fifth root

If actual \(P,Q\in\mathcal D\) satisfied
\([P,Q]=ch\) and \(J=QP+\lambda\), \(\lambda\in k\), then
\([J,P]=-chP\), contrary to (8). Hence no such pair exists,
irrespective of finite Laurent support, classical residues,
parameter poles, or target ordering.

More generally, let \(G(J)\in k(J)^\times\).
If \(P^5=U\,G(J)\), commutation with \(J\) gives
\[
            \sum_{i=0}^{4}
                 P^i\bigl([J,P]+chP\bigr)P^{4-i}=0.        \tag{9}
\]
This linearized power map is injective on \(\mathcal D\).
Indeed ordinary differential order in \(b\) is a multiplicative
order filtration with commutative leading symbols in
\(k(v)[\xi,\xi^{-1}]\).
For nonzero \(P,D\), all five terms \(P^iDP^{4-i}\) have the same
order and symbol \(\sigma(P)^4\sigma(D)\); their sum has nonzero
symbol \(5\sigma(P)^4\sigma(D)\).
Equation (9) would force \([J,P]=-chP\), again impossible.
Therefore
\[
                    P^5=U\,G(J)
                       \quad\hbox{has no solution in }\mathcal D
                       \quad(G\in k(J)^\times).            \tag{10}
\]
In particular \(G(J)=J(J-\beta ch)\) is excluded for every scalar
\(\beta\). Passing to degree-five or higher rational target
denominators does not evade the theorem.

## Scope and the remaining construction

These are the actual scalar division-ring obstructions for the fixed
homogeneous Euler source and scalar cotangent/Mobius cubic transport
in (1). Every rational source Euler function is covered.
The individual quantum quintic models and the nonradical cubic
model in QCS1 remain valid; their proposed compatibility is the
part excluded here. Conditional root-to-fourth-image formulas
remain algebraically correct, but their displayed premise is
ruled out in this entire model family.

There is no conflict with all-order formal lifting.
An \(h\)-adic completion can contain unbounded negative \(v\)-tails,
need not carry the integer \(v\)-adic valuation used in (4), and
can have additional scalar eigenoperators. A completed root or
eigenoperator need not belong to the actual division ring.

The next construction must alter the fixed quantum source or cubic
transport beyond this homogeneous Euler model, or use a different
route with its own finite global admission proof.
Genuinely nonhomogeneous quantum corrections, deformed first images,
and nonscalar carriers remain outside these results.
No theorem puts every bounded Long lift in the excluded family.

[Exact controls](check_quintic_ordering_laurent.py), with
[the retained receipt](quintic_ordering_laurent_controls.json), verify
ordered transported operators, rational-source lowest coefficients,
signed divisor characters, deck projections, eigenrelations and
representative positive/negative differential-order linearizations.
Earlier finite rational-norm controls remain complementary checks.
The universal conclusions use the written valuation, finite-basis,
rational leading-coefficient, and commutative-symbol arguments.

No bounded full Long lift, DC1/DC2 solution, formal verification,
independent whole-proof replay or external review is claimed.
