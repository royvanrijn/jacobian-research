# Monodromy of normalized hypersurface sweeps

Work in characteristic zero, with algebraically closed geometric constants.

**KMFI24.** Let A,B,Delta be polynomial three-vectors in one variable w,
with

\[
 \det(\Delta,A',B)=\det(\Delta,B',B)=0,\qquad
 \det(\Delta,\Delta',B)=c\in k^*. \tag{1}
\]

The normalized ruled-surface sweep

\[
 S(\gamma,w,v)=A(w)+vB(w)+\gamma\Delta(w) \tag{2}
\]

has Jacobian c gamma. If its generic degree N is greater than one,
its geometric generic monodromy is the full natural S_N.
Padding and birational source/target twists preserve this group.

Thus no polynomial affine lift of this class can realize a proper
nonsymmetric subgroup in its N-sheet action. This includes every
normalized ruled member in Gao Section 4.4.1, with all free polynomials.

**KMFI25.** Gao's rational-coefficient four-dimensional Keller map F5
has geometric generic group A10, arithmetic generic group S10,
and quadratic constant field Q(sqrt(6)).
Its full polynomial Jacobian is 160/29.
This identifies an attributed actual map, not a new Keller construction.
The ruled theorem also identifies Gao F4's geometric and arithmetic
groups as S5.

Source: Gao,
[arXiv:2608.00222v2](https://arxiv.org/html/2608.00222v2),
Sections 4.3--4.4. The
[source manifest](sources/sweep-monodromy-manifest.json) pins
the retained version and the classical Jordan theorem used below.
The full Keller spectrum and the specified PSL2(11) actions remain open.

## 1. A birational coefficient curve

Put nu=B cross Delta and

\[
 R(w;Y)=\nu(w)\cdot(A(w)-Y). \tag{3}
\]

B and Delta are independent at every finite w by (1).
Each root of R therefore gives a unique pair (v,gamma) with
Y=A+vB+gamma Delta, rationally by a nonzero two-by-two minor.
The generic sweep cover equals the root cover of R.

Since nu dot Delta'=c, its entries have gcd one.
R is primitive and irreducible in k[w,Y]: it is linear in the
independent Y-coordinates, and a factor depending only on w
would divide all its coefficients.
Its degree in w is the cover degree N.

The map w -> [nu(w)] in P2 is an immersion at every finite w.
Indeed nu dot Delta=0 and differentiation gives nu' dot Delta=-c;
nu' cannot be proportional to nu.

This map is birational onto its image.
Normalize the projective image, a rational curve dominated by P1.
If P1_w mapped to that normalization with degree d>1,
Riemann--Hurwitz would give ramification 2d-2.
Infinity contributes at most d-1, leaving a finite ramification
point, contrary to the immersion.
Thus d=1, and only finitely many distinct pairs (w,u)
have proportional nu-vectors.

## 2. Double transitivity and one transposition

For distinct w,u with independent nu-vectors, the equations
R(w;Y)=R(u;Y)=0 are two independent affine linear equations
in three target coordinates.
Their incidence variety is an affine-line bundle over an
irreducible open subset of the (w,u)-plane.
The finitely many proportional pairs have target images of
dimension at most two and do not dominate target three-space.
Since N>1, the generic ordered-pair cover is nonempty and
irreducible. Monodromy is doubly transitive.

Set kappa=(nu dot A,nu1,nu2,nu3).
Its coordinate polynomials have no common factor.
If their span had dimension two, its projective image would be
a line. Projection to [nu] is birational, so that line would
have a degree-one parametrization. A primitive polynomial tuple
parametrizing a line with degree one has degree at most one,
contradicting N>1.
At least three coordinates are independent. The characteristic-zero
Wronskian criterion gives rank(kappa,kappa',kappa'')=3 at general w.

The equations R=R'=0 form an affine-line bundle over finite w:
nu and nu' are independent everywhere there.
The Wronskian rank means R'' is not identically zero on its
general line, so ordinary double roots occur.
The projective tangent-hyperplane family is not entirely outside
the affine Y-chart: that containment would make the vector
(1,0,0,0) lie in span(kappa,kappa'), forcing nu,nu' proportional.

Let D be the image of this irreducible two-dimensional tangency
incidence. At an ordinary double root, differentiating R=0
gives the target tangent plane with normal [nu(w)].
The projection has rank two there, so D is a hypersurface.
If a general Y in D had two ordinary double roots w,u,
its tangent plane would give [nu(w)]=[nu(u)].
Birationality excludes this generically.
Higher contacts and exceptional finite w occupy subsets of
the incidence of dimension at most one, so cannot dominate D.

Thus a general point of D has one double root and all other
finite roots simple. The normal [nu(w)] varies, so D is not a
hyperplane, in particular not the leading-coefficient hyperplane.
The root degree does not drop generically on D.
Tame inertia there is a single transposition.

A doubly transitive group containing a transposition contains
all transpositions by conjugation. Hence it is S_N, proving KMFI24.

Padding introduces C=gamma x, with x=C/gamma recovered rationally.
C is an independent transcendental coordinate on both sides.
Monomial maps, scalings, an invertible stage and a target twist
are birational on the nonzero chart; the field extension and
group remain unchanged.

## 3. Gao F4 and the whole ruled family

For Delta=(1,w,w^2), Gao uses

\[
 A=(w/3-4w^2/3,\ 2w^2/3-8w^3/9,\
       2w+7w^3/9-2w^4/3),\quad
 B=(w^2-1,w^3,w^4+w^2+2).
\]

These satisfy (1) with c=2, and (3) is

\[
 R=2w^5/9+2w^4/9+(8/9+Y_1)w^3
 -(4/3+2Y_2)w^2+(2Y_1+Y_3)w-2Y_2.
\]

Thus N=5. The actual F4 map is its polynomial twist.
Its affine polynomiality and fibers are in the
[canonical F4 proof](../gao-f4-fibres/PROOF.md).
The geometric group is S5; the arithmetic group over Q is
also S5, with no extra constants.

More generally, for arbitrary polynomials b,a1,a3, set

\[
 B=(b^2-b')\Delta+b\Delta'+\Delta'',\qquad
 A'=a_1\Delta+a_3b\Delta'+a_3\Delta''.
\]

Integration stays polynomial and gives exactly Gao's ruled
normalization. Equation (1) holds identically with c=2.
The theorem covers all these free choices whenever N>1.
The checker verifies an additional nonconstant free-polynomial
member as a control; the universal argument is in Sections 1--2.

## 4. The finite degree-ten sweep behind F5

Write w=w1,z=w2 and X=(X1num,X2num,X3num)/29, with

\[
\begin{split}
X1num&=160w^3-60w^2-240wz+51w+60z,\\
X2num&=60w^4-20w^3-120z^2+51z,\\
X3num&=-48w^5+15w^4+240w^3z-17w^3
       -60w^2z-240wz^2+51wz+30z^2.
\end{split}
\]

For Delta=(1,w,z), S(gamma,w,z)=X+gamma Delta satisfies

\[
 \det JS=\gamma^2. \tag{4}
\]

Eliminate gamma=Y1-X1 and define

\[
 E_2=X_2+(Y_1-X_1)w-Y_2,\quad
 E_3=X_3+(Y_1-X_1)z-Y_3,\quad
 R=-29^4\operatorname{Res}_z(E_2,E_3)/30.
\]

R has degree ten in w, with leading coefficient 20480.
E2 is quadratic in z with constant leading coefficient -120/29.
Therefore w, then z, then gamma are integral over Q[Y].
The sweep is **finite**, not just generically finite.

At Y=(0,1,1), modulo seven, its squarefree polynomial factors as

\[
 5(w+2)(w^2+2w+2)
 (w^7-3w^6-3w^5+w^4-2w^3-w^2+1). \tag{5}
\]

The exact checker proves factor irreducibility over F7.
The linear first subresultant has coefficient coprime to R
at the same good reduction.
All ten roots lift uniquely to z, giving generic degree ten.
The source is geometrically integral, so both groups are transitive.

The arithmetic group contains a Frobenius of type (7)(2)(1).
Its square is a seven-cycle fixing three letters.
An imprimitive degree-ten group lies in S2 wr S5 or S5 wr S2;
neither order is divisible by seven. Hence the group is primitive.
Jordan's prime-cycle theorem gives G_arith>=A10.
The precise attributed input is Theorem 1.1 of
[Jones, arXiv:1209.5169](https://arxiv.org/abs/1209.5169);
the stronger classification in Theorem 1.2 is not used.

## 5. Geometric sign and arithmetic constants

The finite sweep has normal source A3, so it is the full
normalization of the target in its source field.
No missing source divisor can contribute additional inertia.
In characteristic zero the codimension-one different is the
Jacobian divisor. Equation (4) makes each different exponent even.
Its pushforward, the discriminant divisor, is therefore even.

The discriminant of primitive element R differs from the field
discriminant by a square index and constant leading-coefficient
factors. Factoriality of Q[Y] yields

\[
 \operatorname{disc}_w R=d H(Y)^2,\quad d\in Q^*. \tag{6}
\]

Over algebraically closed constants the sign character is trivial.
Thus G_geom<=A10.
It is a transitive normal subgroup of
A10<=G_arith<=S10. Simplicity of A10 forces G_geom=A10.

Exact nonzero specialization at Y=(0,1,1) gives disc R=6q^2,
with the integer q retained in the model.
Thus d has square class six.
The arithmetic sign is nontrivial, also witnessed by (5).
Consequently G_arith=S10, and the constant field is Q(sqrt(6)).
This proves the monodromy assertions of KMFI25.
The arithmetic witness alone is not used as a geometric proof.

## 6. Actual affine polynomial lift for F5

Write v1=xy,v2=x^2z_src,v3=x^3t_src. Gao's triangular stage is

\[
 \gamma=1+20v_1/29,\quad
 u_1=1-49v_1/29+8v_2,\quad
 u_2=1-69v_1/29+9v_2+v_3+4197v_1^2/1682.
\]

Set w=gamma u1,z=gamma^2u2,C=gamma x and

\[
 F_5=(C,S_1/C,S_2/C^2,S_3/C^3).
\]

The expressions E_i=S_i(gamma,gamma u1,gamma^2u2)/gamma^i
are polynomials. After the stage, every monomial of source
x-weight less than i vanishes, with weights (1,2,3).
This proves global polynomiality of all three quotients.

The stage determinant is 160/29. The determinant factors are
x^6, 160/29, gamma^3, gamma^3, C^{-6};
their product is 160/29. Because the outputs are polynomials,
this proves the Jacobian globally on ordinary A4.

The inverse field change recovers x=C/gamma, u1=w/gamma,
u2=z/gamma^2, then v1,v2,v3 by the triangular inverse,
and the original variables by the indicated powers of x.
The target change is Yi=C^i F_(i+1).
Thus the affine field pair is the finite sweep with one independent
C adjoined, preserving geometric and arithmetic groups and constants.

## 7. Scope

KMFI24 excludes a proper nonsymmetric group throughout the
normalized ruled-surface route, even when its affine stage exists.
KMFI25 identifies A10 in a different branch of the method.
Neither result settles nonlinear transverse potentials, arbitrary
direction fields, the mixed-fiber Shabat completion, or the full
action spectrum. The specified PSL2(11) cases remain open.

[Checker](verify_sweep_monodromy.py) · [model](sweep_monodromy_model.json) ·
[receipt](sweep_monodromy_receipt.json).
Written proofs and attributed classical inputs are supported by exact
algebra; no formal verification, independent whole-proof review,
external review or novelty priority is claimed.
