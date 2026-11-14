# The full first-jet cusp test retains all sixteen marked candidates

14 November 2026. PF2RULED22. This completes the noncontact pole(6,9)
**first boundary jet** test left after [PF2RULED21](UNIBRANCH_POLAR33_OBSTRUCTION.md).
All sixteen marked leading candidates lift to global sections satisfying
that test, with the required clean ordinary cusp on the boundary. This does not construct or exclude
Keller pairs or finite maps on the surface. The next new equations
are stated explicitly below; the full original goal remains active.

## The exact assertion

Use the actual ring and coordinates of [PF2WRIGHT1](WRIGHT_MARKED_MODEL.md),
with R=(u=0), v=D and F=F0(v)+u*F1(v)+u^2*F2(v)+.... Put

    J(F)=F1-(3/4)*F0', h=v-s, lambda=(3/4)*h+delta,
    L(F)=h*J(F)-delta*F0'=h*F1-lambda*F0'.

For either of the two degree-eight marked coefficient algebras in
[PF2RULED19](POLAR_CUSP_LEADING_SECTORS.md), delta and s have their fixed
values -a/48 and (24b+72-11a^2)/(12a). All eight points in each algebra
have delta!=0. The following statements hold at every one of the sixteen
complex points, not just at a sampled root:

- L:V5 -> C[v]_(<=19) and L:V8 -> C[v]_(<=31) are surjective.
- With the specified leading symbols H^2,H^3, the full solution spaces
  L(P)=L(Q)=0 are nonempty affine spaces of dimensions31 and85.
- In their product a nonempty Zariski open set has gcd(P0',Q0')=v-s,
  an ordinary cusp at s, and no other boundary parameter mapping to its
  cusp image. Its boundary parametrization is birational onto its image.
- Thus the period and complete first-jet cusp tests eliminate none of
  the sixteen leading candidates. On this affine116-dimensional space,
  the next necessary Jacobian equation is a polynomial of degree<=52
  in v, giving at most53 quadratic coefficient equations over the marked
  candidate algebra. This is a finite next gate, not a claim it is solved.

The affine dimensions include constant target translations. No quotient
by target or marked-source automorphisms is taken. Other height-nine rows,
contact cases, higher heights, global finiteness and actual compatible
normalization remain separate unresolved requirements.

## A triangular right inverse, with exact nonvanishing guards

The normalized V1 generators have restriction and corrected jet

| generator | F0 | J(F) |
| --- | --- | --- |
| g0 | -v^4+v^2 | v/2 |
| g1 | -v^3+v | (v^2+1)/4 |
| g2 | -v^2 | v/2 |
| g3 | -v | 3/4 |

Both J and L are derivations relative to restriction. For every n>=1,
deg L(Vn)<=4n-1 by the all-height generator bound of PF2RULED14.
For i=0,1,2,3 put S_(n,i)=g0^(n-1)*gi. Form the four-by-four matrix
M_n(delta,s) whose i-th column is the coefficients of
v^(4n-1),v^(4n-2),v^(4n-3),v^(4n-4) in L(S_(n,i)), and D_n=det M_n.
The certificate records these exact rational-polynomial matrices and
all eight determinants. For example

    D_1=-delta*(-192delta^3+208delta^2*s-52delta*s^2
                +20delta+3s^3-3s)/8.

For each candidate algebra K=Q[y]/(P8), substitute its rational-univariate
root x(y) and its delta,s, reducing every coefficient modulo P8. For every
n=1,...,8, the resulting D_n is a unit. Sixteen exact Bezout identities
U_n*D_n+V_n*P8=1 prove this at every candidate, including complex roots.
This uses neither a floating-point rank test nor a chosen embedding.

The nonvanishing gives a constructive induction. For a polynomial R of
degree<=4n-1, solve M_n*b for its highest four coefficients and subtract
sum b_i*L(S_(n,i)). The remainder has degree<=4n-5. Repeat at n-1, ending
with degree<=3 at n=1. This constructs a linear right inverse

    right_inverse_n : C[v]_(<=4n-1) -> Vn.

The polynomial degree falls by four at each step. The only divisions
are the recorded D_j units. This proves surjectivity and supplies exact
lifting instructions over the whole reduced candidate algebra. The
section-space coverage and dim Vn=(n+1)(3n+2)/2 are PF2COORD1 inputs.
No new global maps are inferred from this linear calculation.

## Why the four highest cusp coefficients were already complete

For a common nonic H, select any T in V3 lifting its fixed leading
section; such a lift exists by the complete boundary restriction map.
The predecessor equations for a,b,c,d,delta,s are exactly the four
coefficients of L(T) in degrees11,10,9,8. They vanish at all sixteen
points. Hence deg L(T)<=7, so

    L(T^2)=2*T0*L(T) has degree<=19,
    L(T^3)=3*T0^2*L(T) has degree<=31.

The full cusp lifts are therefore explicitly defined by

    P=T^2-right_inverse_5(L(T^2))+K5,
    Q=T^3-right_inverse_8(L(T^3))+K8,
    K5 in ker(L|V5), K8 in ker(L|V8).

All sections with those leading symbols have this form, since their
differences are in V5 and V8. The dimensions are51-20=31 and117-32=85.
Equivalently, for any chosen lift one cancels its whole residual by the
same right inverse; no preferred T or source automorphism is required.
The derivation argument explains the entire first-jet coverage rather
than extrapolating from a bounded highest-coefficient search.

## The required cusp is a nonempty open condition

The missing-prime coordinate A1=XY lies in V3 and has exact expansion
A1=-u+u^2*v. For k=2,3 define global sections

    T_k=(D-s)^k-k*A1*lambda(D)*(D-s)^(k-2).

They lie in V_(k+2), have (T_k)0=h^k and (T_k)1=k*lambda*h^(k-2), and
satisfy L(T_k)=0. Thus T2 in V4 and T3 in V5 can be added to P in V6,
and also to Q in V9, without changing either leading symbol.
After constant translations make P0(s)=Q0(s)=0, write

    P0=h^2*A(v), Q0=h^2*B(v),
    p=P0'/h=2A+h*A', q=Q0'/h=2B+h*B'.

The three independent modifications

    P <- P+alpha*T2, Q <- Q+beta*T2+gamma*T3

replace A,B by A+alpha,B+beta+gamma*h and p,q by
p+2alpha,q+2beta+3gamma*h. There is a nonempty open set on which

    Res_v(p+2alpha,q+2beta+3gamma*h)!=0,
    Res_v(A+alpha,B+beta+gamma*h)!=0,
    [(p+2alpha)*(q'+3gamma)-p'*(q+2beta)]_(v=s)!=0.

Each resultant is a nonzero polynomial: fix alpha and then avoid the
finitely many affine conditions at roots of its first polynomial using
beta,gamma. The last displayed polynomial has alpha*gamma coefficient6,
so is also nonzero. Their product cannot vanish identically in the
three-dimensional parameter space over C.

The first condition gives derivative gcd exactly h; the second excludes
any second boundary preimage of the cusp image. The last is the nonzero
wedge of the second and third derivatives, hence an ordinary cusp. A
polynomial parametrization is finite onto its image because C[v] is finite
over C[P0]. Its factorization through the normalized image is a finite
map of smooth curves. At this cusp the local degree is one, by its
nonzero second/third derivative wedge, and there is no other preimage.
The sum of local degrees is the field degree, so that degree is one.
This proves boundary birationality, not finiteness of the surface map.

## The next Jacobian equation, beyond first jets

Write p=P0'/h, q=Q0'/h, W=p*q'-p'*q. From P1=lambda*p and Q1=lambda*q,
the coefficient of u in the actual chart Jacobian is exactly

    [u]Jac_uv(P,Q)=2*h*(P2*q-p*Q2)+lambda^2*W.

Its value at s is delta^2*W(s), which is nonzero for the clean cusp.
A source Keller pair requires Jac_uv=kappa*u for one nonzero constant
kappa (the sign can be absorbed), so necessarily

    kappa=delta^2*W(s),
    G(v)=[2*h*(P2*q-p*Q2)+lambda^2*W-kappa]/h = 0.       (*)

The numerator is divisible by h, by its value at s. This is an exact
polynomial identity on the full affine cusp-lift space.

To bound its degree, use P=T^2+R5,Q=T^3+R8 with Rj in Vj. The bracket
of T^2,T^3 is zero. For F in Va its coefficient of u^j has v degree
at most4a-j, by the exact Cox chart or the generator degree bound. The
u coefficient in a bracket of Va,Vb has degree<=4(a+b)-3. The cross
terms V5,V9 and V6,V8 therefore have degree<=53, while V5,V8 gives<=49.
Subtracting kappa and dividing h gives deg G<=52. Since all coefficients
of P,Q depend affinely on31+85 parameters, (*) gives at most53 quadratic
equations, not53 asserted independent equations. Neither their solution
set nor the higher u coefficients, field degree or global finiteness is
certified here. A nonzero count of free parameters would not prove any
of those requirements.

## Eliminating the normal coefficients before a nonlinear search

The exact Cox first-two-jet map has kernel dimensions14 in V5 and56
in V8. These kernels consist precisely of functions whose coefficients
F0,F1 both vanish. Their F2 images are the entire polynomial spaces
of degrees<=10 and<=22 respectively, with further kernel dimensions3
and33. The certificate records exact rational nullspace bases, their
second-jet images, full ranks and degree bounds. These statements are
independent of delta,s and are verified on every complete Cox basis.

Inside ker L, the kernel of restriction F->F0 is this same two-jet
kernel: L(F)=0 and F0=0 force F1=0. Consequently the affine cusp-lift
space projects onto an affine boundary-restriction space of dimension

    (31-14)+(85-56)=17+29=46.

For any fixed boundary data, the second coefficients can be changed by

    P2 <- P2+A(v), deg A<=10,
    Q2 <- Q2+B(v), deg B<=22,

independently and without changing the leading symbols. Equation(*)
changes by 2*(A*q-p*B), of degree<=44. The remaining3+33=36 normal
directions have no effect through coefficient u^2 in either coordinate;
they must still be retained when higher Jacobian coefficients are imposed.

On the clean open set gcd(p,q)=1, with deg p=22 and deg q=34 fixed by
the noncontact symbols. The linear map

    (A,B) -> 2*(A*q-p*B), deg A<=10, deg B<=22

has rank11+23=34. Indeed a zero image would force p|A, impossible for
nonzero A of degree<=10; then B=0. Thus there are exactly34 independent
second-coefficient controls, without any parameter sampling.

For a reference lift let G0 denote the quotient(*). Its coefficients
of v^45,...,v^52 must vanish: these eight equations are independent of
all second-coefficient choices. The remaining polynomial of degree<=44
must lie in the image of the explicit45-by34 matrix with columns

    2*q*v^j (j=0,...,10), -2*p*v^j (j=0,...,22).

Its rank is34 on the clean open set, so membership gives eleven linear
consistency conditions on that polynomial once p,q are fixed. Equivalently,
after localizing at any nonzero34-by34 minor, the full next-normal test is
at most19 rational consistency conditions in46 boundary parameters.
Their coefficients depend on the boundary parameters; they are not
claimed to be19 quadratics or independent nonlinear equations. The
uneliminated formulation remains at most53 quadratic equations in116
parameters. All these consistency conditions, higher coefficients and
surface finiteness are unsolved. First checking the eight highest
coefficients is the next economical algebraic gate.

## An actual-ring countercontrol

Set s=0,delta=1 and use P=T2,Q=T3 themselves. Their boundary map is
(v^2,v^3), with derivative gcd v, a singleton ordinary cusp and the full
first-jet multiplier lambda=1+3v/4. Nevertheless their actual chart
Jacobian has u coefficient6*lambda^2, which is not constant. Thus even
the full clean cusp identity in actual global sections does not imply
the source Keller equation. This control is not one of the sixteen
pole69 candidate pairs and is not called a finite surface map.

## Reproduction and limits

The [checker](verify_full_cusp_lifts.py), [certificates](full_cusp_lift_certificates.json)
and [receipt](full_cusp_lift_receipt.json) replay the generator derivations,
triangular matrices, sixteen determinant-unit certificates, predecessor
seals, complete first-jet lifting conditions, affine dimensions, cusp
modifications, exact next Jacobian identity and actual-ring countercontrol.
The nonempty-open, curve-degree and truncated multiplication-map arguments
are written proofs. No
formal verification, independent full reproof, external review, novelty
priority, solved higher Jacobian equations, field-degree-four realization,
global finiteness, actual-normalization compatibility or completed goal is
claimed. All first open reduced pairs remain(5,9),(6,9),(7,9),(8,9).
