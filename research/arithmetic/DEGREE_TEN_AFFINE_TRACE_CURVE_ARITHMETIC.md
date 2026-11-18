# A high-genus curve controls normalized rational degree-ten pencils

There are only finitely many normalized rational even-numerator
degree-ten pencils with the paired-simple-critical passport. This is a
finiteness theorem, without an effective list. It concerns normalized
coefficients and geometric cover types; rational twists and coordinate
rescalings can still produce infinitely many coefficient models.

The nonconstant canonical-scale case lies on one explicit quartic cover
of the slope line. Its geometric monodromy is S4, and its smooth
projective normalization has genus at least twelve. Thus this curve
has no rational parametrization. The constant-scale boundary is covered
by the previously proved rational reciprocal classification.

## The exact curve and its normalization scope

Use h=2 in the [uniform odd affine-trace theorem](ODD_AFFINE_CRITICAL_TRACE_CURVES.md).
Write f(z)=N(z^2)/z with N monic of degree five and N(1)=N'(1)=0.
Its four other critical u-roots have two distinct nonzero squared values.
Their canonical value algebra B has rank two. If the canonical scale
rho is not scalar, it generates B, so the canonical trace has a unique
expression theta=s+t*rho with rational s,t.

The uniform theorem reconstructs q(R)=R^2+q1*R+q0 and N rationally
from s,t. Its only exceptional rational slope is t=-2, which admits
no good model. The last equation is N(1)=0. After taking its reduced
numerator and primitive integer part, denote it by F(s,t).
The [exact input](degree_ten_affine_trace_arithmetic_inputs.json)
contains every coefficient of F and of q0,q1. In particular

\[
 \deg_s F=4,\qquad \deg_t F=16,\qquad \deg F=16,
\]

and its primitive specialization at t=0 is

\[
 f_0(s)=9s^4-516s^3+13140s^2-75456s+230336.                 \tag{1}
\]

The coefficients of F as a polynomial in s have gcd one in Q[t].
Thus the plane equation has no additional vertical component.

The checker independently constructs S by the quadratic resultant,
inverts 2uN'-N=9(u-1)S coefficient by coefficient, checks the complete
trace-square identity, and derives F from N(1). It does not import the
uniform checker as a polynomial-generation oracle.

## The quartic cover is geometrically connected

**Theorem 1.** The splitting group of F over Q(t) and its geometric
splitting group over Qbar(t) are both S4. The smooth projective
normalization C of F=0 is geometrically integral and has genus at least12.

**Proof.** All reductions below have nonzero leading coefficient and
nonzero discriminant. At the specialization (1), the monic reductions are

\[
 \begin{array}{ll}
 p=11:&s^4+5s^3-3s^2-2s+2,\\
 p=13:&(s-5)(s^3+4s^2-2s+4).
 \end{array}                                             \tag{2}
\]

The first quartic is irreducible: the exact finite-field powers give
x^(11^4)=x modulo it and gcd(x^(11^2)-x,f)=1. The cubic in the second
row has no root in F13. To apply the integral monic form of
[Dedekind's theorem, Milne ANT3.00, Theorem8.23,p138](sources/milne-ant300.pdf),
use the polynomial 9^3*f0(X/9), which is monic integral and has the
same splitting field. Since 9 is a unit at11 and13, multiplying roots
by9 preserves the factor degrees in (2).
Hence the arithmetic group contains a four-cycle
and a three-cycle and is S4. For an elementary finite check, fixing one
four-cycle, each of the eight possible three-cycles generates S4 with it.
Good specialization at t=0 embeds this specialization group into the
generic group, so the latter is S4 as well.

The retained complete discriminant factor product has degrees and
multiplicities

\[
 \operatorname{Disc}_s F
 =c\,a_1^{12}a_3^2b_3^4a_6^4a_{30},                    \tag{3}
\]

where subscripts give degrees, c is a nonzero integer, a30 is squarefree,
and a30 is coprime to all other factors and to the leading coefficient
of F in s. These are exact polynomial product and gcd identities,
replayed from the input without a factorization oracle.

Each of the thirty distinct geometric roots of a30 has discriminant
valuation one and unit leading coefficient, hence transposition inertia.
Here is the local reason. Over the characteristic-zero local field with
algebraically closed residue field, an inertia cycle of length e
contributes at least e-1 to the root-difference discriminant valuation.
If all roots are unramified, that valuation is even. Valuation one
therefore forces precisely one two-cycle and all other sheets fixed.
This is also the usual simple discriminant branch calculation.

The geometric group is a normal subgroup of the arithmetic S4.
It contains a transposition, and all transpositions are conjugate in S4
and generate it. Thus the geometric group is S4. In particular C is
geometrically integral and C->P1_t has degree four.

Each of these thirty branch points contributes one to ramification in
that degree-four cover. Riemann--Hurwitz gives

\[
             2g(C)-2=-8+\operatorname{Ram}(C/P^1)
                       \ge-8+30=22,
\]

so g(C)>=12. Other branch points can increase this lower bound;
no exact genus is asserted. ∎

**Corollary 2.** There is no nonconstant rational parametrization of
F=0. For every fixed number field K, C(K) and the affine set F(K)=0
are finite.

**Proof.** A nonconstant rational map P1->C extends to the smooth
projective curves. Riemann--Hurwitz contradicts g(C)>=12.
Finiteness follows from [Faltings, *Endlichkeitssätze für abelsche
Varietäten über Zahlkörpern*, Satz7, printedp365](sources/faltings-1983-endlichkeitssaetze.pdf).
This is the imported finiteness theorem for smooth curves over number
fields; it gives no effective point list or uniform bound.
Every nonsingular affine K-point lifts to a K-point of C, and there
are only finitely many singular affine points. This gives the affine
finiteness assertion too. ∎

## Consequence for rational degree-ten paired pencils

**Corollary 3.** There are only finitely many monic rational N of
degree five with N(1)=N'(1)=0 giving a good paired pencil N(z^2)-wz.
In particular their geometric cover types form a finite set.

**Proof.** For nonscalar canonical rho, rank(B)=2 makes theta affine
in rho. The complete rational-slope statement of the uniform theorem
puts the model on F(Q)=0, and its rational reconstruction supplies a
unique N at each good point. Corollary2 makes these points finite.

If rho is scalar, all paired critical products uv=rho^2 are constant.
This is within the reciprocal ansatz. The
[global normalized rational reciprocal classification](MIXED_RECIPROCAL_COUPLED_ODE.md),
Theorem5, gives the unique retained model in that boundary case.
Nonconstant signed rho can also have constant square; such cases are
already covered by the nonscalar argument, so the dichotomy loses none.
This proves finiteness without determining the remaining rational points.
∎

For an unnormalized monic rational N, the unique double root u0 of the
zero-level numerator is rational: gcd(N,N') is linear. The rational
operation N_new(u)=u0^(-5)*N(u0*u) imposes the normalization above.
Over Qbar it is a source and target rescaling, so every such original
cover has one of the finitely many geometric types.

If u0 is nonsquare, this operation need not be a rational change of
the z coordinate. Consequently the theorem does not claim finitely
many Q-isomorphism classes of unnormalized covers or finitely many
raw coefficient tuples. It does not exclude families of twists or
coordinate rescalings. Over a larger fixed number field, exceptional
nonrational slopes require their own analysis before a corresponding
all-model finiteness statement can be made.

The result treats even-numerator pencils. Arbitrary degree-ten
polynomial pencils without that symmetry, rational paired degree26
existence, and the exact rational points of F remain unresolved here.

## Evidence and assurance

The [checker](verify_degree_ten_affine_trace_arithmetic.py) binds the curve
to direct norm/Euler/trace reconstruction, verifies (2) by modular
arithmetic, verifies every factor product and gcd in (3), and checks
the complete rational diagonal exceptions. Its
[receipt](degree_ten_affine_trace_arithmetic_verification.json) and
[source ledger](degree_ten_affine_trace_arithmetic_source_ledger.json)
record exact inputs and imported mathematical boundaries: the uniform
affine-trace reconstruction, Dedekind's theorem and Faltings finiteness.
No exact genus computation, formal verification, independent external
replay, external human review or novelty claim is made.
