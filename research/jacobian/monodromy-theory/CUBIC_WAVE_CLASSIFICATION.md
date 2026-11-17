# Cubic primitive pairs with constant scalar

Work over algebraically closed characteristic-zero constants. Use
\[
 P=A\begin{pmatrix}-bc&b^2\\-c^2&bc\end{pmatrix},\quad
 (Ab^2)_{ww}+2(Abc)_{wz}+(Ac^2)_{zz}=0,\quad \gcd(b,c)=1.
\]
Here max(deg b,deg c)=3 and A is a nonzero **constant**.
Normalize A=1 by absorbing its nonzero constant into the pair.

**KMFI56.** The curved family
\[
 b=w^2-\frac43z,\qquad c=w^3-wz+D                            \tag{1}
\]
is primitive and solves the wave equation for every D. Its reconstructed
sweep S=X+gamma(1,w,z) is finite of degree21 and has full geometric
monodromy A21. No affine source conjugacy of the primitive vector
puts its ratio into k(w)[z].

**KMFI57.** The member D=1 has an explicit polynomial Keller lift on
ordinary affine four-space, with generic degree21 and geometric group
A21. The construction below specifies its coefficients in a number
field of degree at most18 and proves its full constant Jacobian.
This is a new curved-pole construction for an already realized action,
not a claim of a new abstract action or of inequivalence under every
birational field rechart.

**KMFI58.** Every primitive cubic pair with constant scalar has full
geometric field group in the exact set
\[
                 \{A7,A8,A10,A12,A21\}.                      \tag{2}
\]
Every listed group has a cubic constant-scalar field witness.
Thus both PSL2(11) actions and natural A4 are excluded throughout
this coefficient-unrestricted class. The statement does not cover
nonconstant cubic scalars, primitive-pair degree>=4 or arbitrary
Keller constructions.

[Exact reconstruction and stage](cubic_wave.py) ·
[checker](verify_cubic_wave.py) · [model](cubic_wave_model.json) ·
[receipt](cubic_wave_receipt.json).

## 1. The curved field

Direct integration gives
\[
\begin{split}
 X_1={}&\frac35w^5-\frac53w^3z+\frac43wz^2+\frac43Dz,\\
 X_2={}&\frac13w^6-\frac23w^4z+\frac{16}{27}z^3
                         -\frac13Dw^3+\frac43Dwz,\\
 X_3={}&-\frac17w^7+w^5z-2w^3z^2+\frac43wz^3
                         -\frac12Dw^4+Dw^2z-D^2w.
\end{split}                                                   \tag{3}
\]
These satisfy the six reconstruction derivatives of KMFI41 and
det DS=gamma^2. Primitivity follows by substituting z=3w^2/4
into c: the remainder w^3/4+D is not the zero polynomial.

The affine reduction obstruction is elementary. If some affine
conjugacy put the primitive ratio in k(w')[z'], primitivity would
make its first vector component a polynomial in w' alone.
In the old coordinates it would be a nonzero linear combination
ell*b+m*c depending on one affine linear form. For m!=0 its cubic
top part is m*w^3, so that form must be parallel to w. But the
term -m*w*z cannot cancel. For m=0 the quadratic top part similarly
fixes the form parallel to w, whereas b has a nonzero z term.
Both cases are impossible. This excludes affine conjugacies, not
arbitrary birational equivalences of the resulting field pair.

To prove finiteness use weights(gamma,w,z)=(5,1,2).
The top weighted map is S_0, of output weights(5,6,7); every D
correction in(3) has strictly lower corresponding weight.
For w!=0 put t=z/w^2. The two equations after eliminating gamma
in S_0=0 are
\[
 P(t)=-4/15+t-4t^2/3+16t^3/27=0,\qquad
 Q(t)=-1/7+2t/5-t^2/3=0.
\]
These polynomials are coprime. For w=0 the second sweep coordinate
is 16z^3/27, forcing z=0, then gamma=0.
Thus the top weighted map has only the origin above zero.

Here is the degree implication used in this argument. A weighted
homogeneous three-variable system with zero-dimensional origin fiber
is a regular sequence. Its quotient Hilbert series is the product
of (1-t^(output weight)) divided by the product of
(1-t^(source weight)). The value at t=1 is its finite length.
Lifting a regular sequence through the weight filtration preserves
that quotient length: cancel a highest-weight term successively
against the leading sequence; its Koszul relations lift in the
same decreasing-weight induction. This also proves finiteness
for every lower-weight perturbation and every constant target fiber.
In the present case the length is
\[
                      \frac{5\cdot6\cdot7}{5\cdot1\cdot2}=21.
\]
Alternatively the retained symbolic resultant in z has exact degree21
in w and constant leading coefficient -1/156279375 for all D and
all target coordinates. Its other equation is cubic in z with
constant leading coefficient16/27, giving a second finiteness check.

For the branch-surface degree use
\[
 \rho=c/b,\quad \eta=w\rho-z,\quad \beta=-\rho.
\]
The cross product X_w cross X_z is nonzero and proportional to
(eta,beta,1). Eliminating z=-eta-beta*w from c+beta*b=0 yields
\[
 w^3+2\beta w^2+(\eta+4\beta^2/3)w+4\beta\eta/3+D=0.          \tag{4}
\]
Over k(beta), eta is therefore a degree-three rational function of w:
its denominator is w+4beta/3, and the numerator at that root is
D-16beta^3/27, nonzero in k(beta).
Consequently [k(w,z):k(eta,beta)]=3.
The Gauss field lies in the branch-image field k(X).
KMFI32 gives [k(w,z):k(X)]<=2; it also divides3 by this tower.
It is therefore1.

The finite sweep has only gamma0 as a ramification divisor.
The different order2 gives ramification index3, and its birational
map onto its image gives one3-cycle in the full affine inertia.
All remaining sheets there are unramified. The conjugates of this
inertia normally generate the geometric group, as in KMFI27:
their quotient would give a nontrivial finite etale cover of A3.
A transitive group generated by3-cycles is the full alternating
group on its support, by the connected-support lemma already used
in KMFI12. Thus G=A21. This proves KMFI56.

## 2. A regular marked zero for D=1

Define d=3P and
\[
\begin{split}
 G(t)&=3/5-5t/3+4t^2/3+(4/3)d(t)t,\\
 H(t)&=89600t^6-336000t^5+554400t^4-509040t^3\\
     &\hspace{25mm}+274050t^2-82215t+10854.
\end{split}
\]
The following identities are exact:
\[
 Q+d(t-1/2-4t^2/3)-d^2=-H/28350,\qquad
 \gcd(H,d)=\gcd(H,G)=1.                                     \tag{5}
\]
Choose any root tau of H and any w0 with w0^3=1/d(tau).
Put
\[
 z0=\tau w0^2,\qquad \gamma0=-w0^5G(\tau).
\]
Then gamma0!=0. The equations V=X2-wX1 and W=X3-zX1 at this
point reduce respectively to d=3P and(5), so
S(gamma0,w0,z0)=0. Since det DS=gamma0^2!=0, this zero is regular.
All coefficients lie in Q(tau,w0), of degree at most18.
No numerical root choice or approximate Jacobian is needed.

The homogeneous member D=0 has no regular zero over its centered
origin. The nonzero constant deformation supplies exactly the
marked point needed here.

## 3. Polynomial stage and ordinary-affine lift

Scale **w=gamma^2*u,z=gamma*v** and order the sweep outputs as
(S1,S3,S2), with exponents(e1,e2,e3)=(1,2,3). Set
\[
 E=\big(S1/\gamma,\ S3/\gamma^2,\ S2/\gamma^3\big).
\]
All three entries are polynomials. For D=1 they are
\[
\begin{split}
 E1={}&1+4v/3+(4/3)\gamma^3uv^2
                  -(5/3)\gamma^6u^3v+(3/5)\gamma^9u^5,\\
 E2={}&-u+v+(4/3)\gamma^3uv^3+\gamma^3u^2v
       -2\gamma^6u^3v^2-\gamma^6u^4/2
       +\gamma^9u^5v-\gamma^{12}u^7/7,\\
 E3={}&u+4uv/3+16v^3/27-\gamma^3u^3/3
                       -(2/3)\gamma^6u^4v+\gamma^9u^6/3.
\end{split}                                                   \tag{6}
\]
The usual opposite weight ordering would leave a negative gamma
power in X3; the scaling in(6) is essential.

Let theta=(gamma,u,v) and
theta0=(gamma0,w0/gamma0^2,z0/gamma0).
Then E(theta0)=0. At a zero of S the determinant chain gives
\[
 \det J=\det DE(\theta0)=-\gamma0^{-1},\qquad J^{-1}=L.       \tag{7}
\]
The minus sign is the output swap S2,S3.
Let c1 be the first column of L and
f20=-(D_(c1)^2 E3)(theta0)/2.
Define the polynomial stage
\[
 \theta=\theta0+
     L\begin{pmatrix}v1\\v2\\v3+f20\,v1^2\end{pmatrix}.       \tag{8}
\]
It is an affine isomorphism followed by a triangular shear,
so its inverse is explicit and its Jacobian is det L=-gamma0.

There is also a compact exact coefficient recipe avoiding an
unsolved matrix system. Replace gamma^3 by q in(6), obtaining Ehat.
In Q[tau] put
\[
 q0=-G^3/d^5,\quad u0=d^3/G^2,\quad v0=-\tau d/G.
\]
Then gamma0^3=q0, and set
\[
 B=(3q\,Ehat_q,\ Ehat_u,\ Ehat_v)|_{(q0,u0,v0)},\quad R=B^{-1}.
\]
Exact reduction modulo H gives det B=-1. Thus
L=diag(gamma0,1,1)*R. The checker retains every entry of R and
f20 as a degree-at-most-five rational polynomial in tau.
These are explicit algebraic constants.

On ordinary coordinates (x,y,z_src,t_src), set
\[
 v1=xy,\quad v2=x^2z_{\rm src},\quad v3=x^3t_{\rm src},
 \quad C=x\gamma,\qquad
 F=(C,\ E1/x,\ E2/x^2,\ E3/x^3).                             \tag{9}
\]
Every quotient is a polynomial. Give(v1,v2,v3) weights(1,2,3).
At theta0, the linear terms after L are (v1,v2,v3).
For E1 the constant term is the only forbidden weight and vanishes.
For E2 the weights below2 are its constant and v1 terms, both zero.
For E3 the complete list below3 is1,v1,v2,v1^2:
the first three vanish by E(theta0)=0 and JL=I, and the final
coefficient is canceled exactly by f20 in(8).
All higher ordinary Taylor orders have weight at least3.
This proves global polynomiality of(9).

For completeness the full determinant factors, in the source order
(x,y,z_src,t_src) and target order(9), are
\[
 x^6\det L\cdot\gamma\cdot(-\gamma^2)\gamma^3/C^6
                    =-\det L=\gamma0\ne0.                 \tag{10}
\]
The factors are respectively the stage, C=x gamma, the reordered
sweep determinant, w,z scaling, and target rescaling.
The equality holds on a dense open and hence globally, since the
coordinates in(9) are polynomials. Constant rescaling of one target
coordinate can make the determinant1 if desired.

Finally
\[
 S1=C F2,\qquad S3=C^2 F3,\qquad S2=C^3 F4.
\]
Conversely x=C/gamma,u=w/gamma^2,v=z/gamma, and the inverse
of(8) recovers(v1,v2,v3), then y,z_src,t_src by division by
x,x^2,x^3 on the generic open.
Thus the field pair of F is precisely the field pair of S1 with
an independent variable C adjoined. It has degree21 and group A21.
This proves KMFI57 on ordinary A4, including the source divisor x=0.

## 4. Exhaustion of constant-scalar cubic pairs

KMFI49 leaves constant projective top direction or a linear one
with its degree-two common factor supported on the two eigenlines.
In the latter case exchange axes if needed. There are two cases.

### 4.1 One eigenline with multiplicity two

Normalize b_top=w^3,c_top=kappa*w^2*z.
The top wave equation is kappa^2+5kappa+15=0.
Write the quadratic lower terms as
b2=a0*w^2+a1*w*z+a2*z^2,
c2=b0*w^2+b1*w*z+b2*z^2.
The ordinary degree-three equations are
\[
\begin{split}
 8(\kappa+5)a0+4(\kappa+2)b1&=0,\\
 12(\kappa+2)a1+12(\kappa+1)b2&=0,\\
 12(\kappa+1)a2&=0.
\end{split}                                                   \tag{11}
\]
Their solution space consists exactly of the two translations and
one curvature direction:
\[
 (b2,c2)=(3t w^2,\kappa s w^2+2\kappa t wz)
            +\big(-(\kappa+12)B wz/9,\ Bz^2\big).
\]
If B!=0, weights(1,2) make the last curvature terms join the
top vector(w^3,kappa*w^2*z). Its z^2 wave coefficient, after(11),
is 2B^2(291-8kappa)/81, nonzero. No lower weighted term can cancel it.
Thus B=0. Translations remove all quadratic terms.

Now write b=w^3+a*w+b*z+c,c=kappa*w^2*z+d*w+e*z+f.
The full remaining equations kill the z and constant coefficients
in b and give (6kappa+24)a+(4kappa+6)e=0,
and a^2+ae+e^2=0.
The last two equations force a=e=0:
eliminating a/e against t^2+t+1 and the kappa quadratic
has nonzero resultant48384 (and e0 gives a0 directly).
Therefore
\[
 b=w^3,\quad c=\kappa w^2z+d w+f,\quad f\ne0.                \tag{12}
\]
The last condition is precisely primitivity.
This is lambda=w^6,R=d/w+f/w^2 in KMFI37.
Its group is A8 if d0, otherwise A10.

### 4.2 Both eigenlines

The top pair is(w^2*z,kappa*w*z^2), with
2kappa^2+3kappa+2=0.
The ordinary degree-three coefficient matrix has rank4 and
kernel exactly the two source translations:
\[
 (b2,c2)=(s w^2+2t wz,\ 2\kappa s wz+\kappa t z^2).
\]
Remove them. For remaining affine terms the full wave forces
the z coefficient in the affine part of b, the w coefficient in
the affine part of c, and both affine constant terms equal zero, with
\[
 (8\kappa+12)a+(12\kappa+8)e=0,\qquad a^2+ae+e^2=0.
\]
The analogous resultant is256, so a=e=0.
Both vector components retain common factor wz; primitivity fails.
There are no primitive solutions in this case.

### 4.3 Constant leading direction

Normalize b_degree3=0,c_degree3=w^3.
The ordinary degree-three equations first force
b_z^2=0 (the coefficient of z^2 in b) and
c_z^2=-2*b_wz (the coefficient convention refers to monomials,
not derivatives). Write beta=b_wz.
If beta!=0, weights(1,2) have top vector
(beta*w*z,-2beta*z^2), whose wave is38beta^2*z^2.
So beta=0 for constant scalar.
Consequently
\[
 b=a w^2+e z+d w+f,\quad
 c=w^3+g w^2+h wz+jz+iw+l.                                  \tag{13}
\]
If e=0, this already belongs to k(w)[z].
Primitivity and KMFI33 force b constant, linear or a square of
one linear factor. Up to translations/scalings their chart data are:
- b1, polynomial cubic h, constant lambda: group A7;
- b=w, lambda=w^2,r2,s0, R cubic: group A10;
- b=w^2,lambda=w^4,r4,s0,R=w^2+L/w with L!=0: group A12.
Allowed affine R terms do not affect the group. The leading cubic
coefficient is nonzero in every case.

If e!=0, the highest (1,2)-weighted equations from(13) are
\[
 4e(a+h)=0,\quad 12a^2+6ah+6e+2h^2=0.
\]
Thus h=-a,e=-4a^2/3,a!=0. Constant diagonal scalings of source
and vector normalize a=1.
The complete lower equations are then
g=3d/2 and f=4i/3-d^2-dj-j^2.
A w translation and z translation remove d,g,i, leaving
\[
 b=w^2-4z/3-J^2,\quad c=w^3-wz+Jz+L.
\]
Apply the permitted affine covariance with
old w=w-J,old z=z-(3J/2)w and
new vector=(b,c+(3J/2)b).
It becomes exactly(1), with D=L-J^3.
Hence every curved case has group A21 by Section1.

This exhausts all top directions and all lower coefficients.
Explicit field witnesses for(2) are:
\[
\begin{array}{c|c|c}
G&b&c\\ \hline
 A7&1&w^3\\
 A8&w^3&\kappa w^2z+1,\ \kappa^2+5\kappa+15=0\\
 A10&w&w^3+\kappa z,\ \kappa^2+\kappa+1=0\\
 A12&w^2&w^3+\kappa wz+1,\ \kappa^2+3\kappa+6=0\\
 A21&w^2-4z/3&w^3-wz
\end{array}
\]
The exact checker verifies the coefficient matrices, translation
kernels, curvature obstructions, both nonzero resultants, wave
witnesses, group formulas and marked affine stage.
This proves KMFI58.

## 5. Assurance and scope

The written exhaustion and geometric group proof are universal.
Exact computations replay algebraic identities, the marked point
and its finite stage constants; no finite-field sample is used
as a monodromy proof. Formal verification and independent
whole-proof review are not claimed.

The new affine-four-space Keller map realizes A21, already in the known action
spectrum. Its importance here is an actual affine completion of a
primitive curved-pole input outside the affine reduction in KMFI55.
It does not prove a complete action spectrum or decide absolute
PSL2(11), natural A4 or binary degree8 occurrence.
The constant-scalar restriction is essential: b=1,c=w^3 with
A=z-w^4/2 is an exact wave and has group A22 by KMFI44--45,
outside(2). The remaining polynomial-wave frontier includes cubic
pairs with nonconstant scalar and primitive-pair degrees at least4.
