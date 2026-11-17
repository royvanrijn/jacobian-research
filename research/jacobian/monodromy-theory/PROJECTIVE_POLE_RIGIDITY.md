# Projective pole-direction rigidity

Work over algebraically closed constants of characteristic zero.
A nonzero polynomial nilpotent normalization can be written

\[
 P=A\begin{pmatrix}-bc&b^2\\-c^2&bc\end{pmatrix},
 \qquad \gcd(b,c)=1,
\]
\[
 {\cal W}(A,b,c):=(Ab^2)_{ww}
          +2(Abc)_{wz}+(Ac^2)_{zz}=0. \tag{1}
\]

The scalar A and primitive pair b,c are polynomials.
Where b!=0, lambda=Ab^2,rho=c/b give the earlier wave chart.
Equation(1) reconstructs a polynomial X uniquely up to translation,
with S=X+gamma(1,w,z),det JS=gamma^2, as in KMFI41.

**KMFI48.** If all nonzero entries of P are homogeneous of the
same ordinary degree n, its primitive pair has degree at most one.
Up to an invertible linear change of the direction coordinates
and the matching linear target change, there are exactly two types:

1. Constant direction:
   rho=0,lambda=U*z^n+V*w*z^(n-1) for n>=1
   (constant lambda for n0).
   The field is birational if V=0, otherwise natural S_(n+1).
2. Two distinct linear eigen-directions:
   rho=kappa*z/w,lambda=C*w^r*z^s,
   r>=2,s>=0,r+s=n, C!=0, and
   r(r-1)+2kappa(r-1)(s+1)+kappa^2(s+1)(s+2)=0.
   The degree is n+1, with natural A_(n+1) in odd degree and
   S_(n+1) in even degree.

In particular no homogeneous normalization tensor, even with a
rational multiplier having several apparent pole directions,
gives either PSL2(11) field action.
The obstruction to degree>=2 is an all-degree argument, not a
bounded scalar search.

**KMFI49.** For any polynomial pair b,c in (1), let D be the
largest component degree and divide its top homogeneous pair
by their common homogeneous factor g.
Its remaining projective degree is at most one.
If that degree is zero, the derivative of g along the constant
top direction is zero. Thus g is a power of a transverse linear
form. If that degree is one, the entire highest tensor is of
the monomial eigen-direction type in KMFI48.
These are necessary conditions for a scalar A of any degree.

Every primitive affine-linear pair b,c can be reduced by an
affine source change and matching linear target change to a
multiplier in k(w)[z], or to constant direction.
Consequently all its geometric groups are natural alternating
or symmetric, and both PSL2(11) actions are excluded as field pairs.

Passing the highest-direction test is not sufficient for existence
of A, polynomial X, or an affine Keller lift.
Genuinely curved pole candidates remain outside this theorem.

**KMFI50.** Suppose A is a nonzero constant and the primitive
pair b,c has maximum degree exactly two. Every solution of (1)
is affine-conjugate to one of the following forms, up to nonzero
scalings:

\[
 b=w^2,\quad c=\kappa wz+D,\quad
 \kappa^2+3\kappa+6=0,\quad D\ne0; \tag{2}
\]
or
\[
 b=aw+c_0,\quad c=w^2+dw+ez+f,\qquad
 a^2+ae+e^2=0. \tag{3}
\]

The primitive condition is retained.
Their full geometric groups are respectively A6 for (2);
A5 when a=e=0 in (3); and A8 when a!=0 in (3).
Thus the first quadratic primitive-pair class with constant scalar
has no new proper primitive action.
Nonconstant scalar A and higher primitive degree remain open.

[Algebraic helper](pole_direction.py) · [checker](verify_pole_direction.py) ·
[model](pole_direction_model.json) · [receipt](pole_direction_receipt.json).

These are constructor and leading-form results.
They do not decide absolute PSL2(11), natural A4 or the binary
degree-eight Keller action. The parent classification remains open.

## 1. Primitive factorization and coordinate covariance

For P=[[-q,p],[-r,q]], nilpotency gives pr=q^2.
Divide p,q,r by their common polynomial factor.
The remaining coprime outer entries are squares in the UFD
k[w,z], giving p=Ab^2,q=Abc,r=Ac^2 after adjusting constant units.
The cases of a zero outer entry are constant directions.

For an affine change p_old=T*p_new+p0, use the target change
(Y1,Yvec)->(Y1,T^-1(Yvec-p0*Y1)).
Then

P_new=T^-1 P_old(T*p_new+p0) T.

The primitive vector becomes T^-1(b,c)(T*p_new+p0) and the
scalar becomes det(T)*A(T*p_new+p0).
Thus polynomiality, the wave identity, generic degree and
monodromy are preserved. Linear changes preserve homogeneity.

## 2. The projective fixed-direction argument

Suppose P is homogeneous of degree n and the primitive pair
has homogeneous degree ell>=1.
Write M=deg A=n-2ell, so M>=0.
The ratio R(t)=c(1,t)/b(1,t) is a rational map of degree ell:
coprimality prevents a common projective zero.

The identity map R(t)=t has no solution.
The wave equation below would become n(n+1)L=0 with n>=2.

Otherwise change the projective coordinate so infinity is not
a fixed point. This is a linear change allowed above.
Then R(infinity) is finite and R(t)-t has ell+1 finite zeros.
Put lambda=w^n L(t),t=z/w and delta=R-t.
Direct differentiation of (1) gives

\[
 (\delta^2L)''+2(n+1)(\delta L)'
                         +n(n+1)L=0. \tag{4}
\]

All fixed points are simple.
If delta has zero order>=2 and L has order m>=0 there,
the n(n+1)L term is uniquely lowest order in (4), impossible.

At a simple fixed point t_i, write
delta=alpha_i*(t-t_i)+..., and let m_i be the zero order of L.
The denominator b does not vanish at a fixed point, because
that would also force c=0. Hence m_i is precisely the zero
multiplicity of A along that fixed line, and m_i<=M.
The leading term of (4) gives

\[
 \alpha_i^2(m_i+1)(m_i+2)
   +2(n+1)\alpha_i(m_i+1)+n(n+1)=0.
\]

For u_i=1/alpha_i the discriminant is

4(n+1)(m_i+1)(m_i+1-n)<0,

since m_i<=n-2ell<=n-2.
Thus u_i is quadratic algebraic over Q with normalized trace

\[
 {\rm tr}_{\mathbf Q}(u_i)=-{m_i+1\over n}. \tag{5}
\]

No analytic fixed-point theorem is needed for the next identity.
The rational function 1/(R(t)-t) has simple poles exactly at the
fixed points, residues u_i, and expansion -1/t+O(1/t^2) at infinity.
Partial fractions therefore give sum u_i=-1.
Take normalized traces in a common number field containing all u_i:

\[
 \sum_i(m_i+1)=n.
\]

But sum m_i<=M and there are ell+1 fixed points, so

n<=M+ell+1=n-ell+1.

Consequently ell<=1. This rejects every homogeneous primitive
pair of degree at least two for every polynomial scalar degree.

## 3. Equality classifies the homogeneous tensor

For ell1, equality requires all scalar zeros to lie on the two
fixed lines. The linear pair is an invertible matrix map:
otherwise its components would share a linear factor.
The fixed points are distinct, so that matrix diagonalizes.
Put the pair in the form (w,kappa*z), absorbing a nonzero scalar.
A is then C*w^(r-2)*z^s, giving exactly type2.
The wave equation gives its displayed quadratic residue condition.
KMFI34 supplies the exact degree and full group.

For ell0, put the constant pair in the direction (1,0).
Equation(1) becomes A_ww=0.
For n>=1, homogeneity gives
A=U*z^n+V*w*z^(n-1).
Integrating yields

X1=-V*z^n/n,
X2=U*z^(n+1)/(n+1),
X3=-V*z^(n+1)/(n+1),

up to translations. KMFI38 gives the stated symmetric or
birational cases. The constant n0 case is birational as well.
This proves KMFI48.

## 4. Highest terms and affine fractions

For a general pair, take its highest homogeneous vector v_D.
The highest tensor is A_top*v_D*v_D^perp.
If v_D=g*vbar with vbar primitive of degree ell, its scalar is
A_top*g^2. The highest wave equation is still zero, so KMFI48
applies. It gives ell<=1.
For ell0, after putting vbar=(1,0), (A_top*g^2)_ww=0.
Degree in w then forces g_w=0. Hence g is a power of z,
a transverse linear form in coordinate-free language.
For ell1 the equality classification supplies the full top
monomial and eigenvalue constraints.

For a primitive affine pair v=M*p+b0, choose a left eigenvector
of M as the first coordinate. Then its first component becomes
alpha*w+beta, independent of z; the second remains affine.
If the first component is nonzero, their ratio belongs to k(w)[z].
If it is identically zero, a swap makes the ratio zero.
Thus KMFI46 applies.

More precisely, a z-dependent ratio is fractional affine with
linear denominator in w; its R in the KMFI33 normal form is affine,
so KMFI34 gives natural alternating/symmetric group.
A polynomial ratio independent of z is constant or linear,
so KMFI38 or the finite heat classification with deg h1 applies.
The signed and binary even-quadratic exceptions cannot occur.
This proves KMFI49 without assuming existence from the leading test.

## 5. Exhausting the quadratic constant-scalar pair

Scale A to1.
The homogeneous-top result leaves a nonconstant linear
projective top direction or a constant one.

In the first case, the common top factor has degree1.
KMFI48 puts it on one eigenline; exchange the coordinates if
needed and scale to obtain b_top=w^2,c_top=kappa*w*z.
The top wave equation is kappa^2+3kappa+6=0.
Translations remove the w term from b and the w term from c,
so write b=w^2+az+b0,c=kappa*w*z+cz+d0.
The full wave is

\[
 2(\kappa^2+3\kappa+6)w^2+
 4c(\kappa+1)w+4a(\kappa+1)z+
                    2b_0(\kappa+2)+2c^2.
\]

Since kappa is nonreal, this forces c=a=b0=0.
Primitivity forces d0!=0, giving (2).
Its chart is lambda=w^4,rho=kappa*z/w+d0/w^2.
KMFI43 gives degree6 and A6.

In the constant top-direction case, put that direction vertical.
Equation (c_top^2)_zz=0 forces c_top to be a nonzero multiple
of w^2; scale to1. Write

b=aw+b1*z+c0, c=w^2+dw+ez+f.

Its full wave is

4b1*w+2(a^2+b1*d+ae+e^2).

Hence b1=0 and a^2+ae+e^2=0, giving (3).
If a0, then e0 and c0!=0 by primitivity: this is polynomial
h of degree2 with constant lambda, hence degree5 and A5.
If a!=0, e/a is a primitive cube root and a shift of w makes
lambda proportional to w^2.
The resulting rho is kappa*z/w+R(w)/w with R quadratic and
nonzero leading coefficient. Removing its affine part gives
the r2,s0,D2 case of KMFI43, hence degree8 and A8.
This exhausts every quadratic primitive pair with constant scalar.

The remaining inhomogeneous problem now starts with either
nonconstant scalar for a quadratic pair or higher pair degree,
subject to the universal highest-direction restrictions above.
All proofs are written; exact controls are not formal verification
or independent whole-proof review.
