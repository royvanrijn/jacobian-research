# No degree-five Keller completion of the exceptional cubic component

Let r=x+x^2y. There are no P,Q in C[x,y,z] of total degree at most5 with
Jac(P,Q,r)=kappa!=0. This extends the quartic exceptional-scalar obstruction
in Wang-Yan3.17 by a squarefree-fibre specialization argument. No novelty
priority or external referee is claimed. The proof below is universal;
the accompanying symbolic identities and ten degree cases are controls,
not a formalization or bounded coefficient search.

## Basic derivation and earlier obstructions

Put B=C[x,y] and D=-x^2*d_x+(1+2xy)*d_y. Then
Jac(P,Q,r)=D(P)*Q_z-P_z*D(Q). In the rational coordinates(r,u=1/x), D=d_u,
so ker in C(x,y) is C(r), and the polynomial kernel is C[r]. The last
intersection follows by Bezout for reduced univariate rational functions.
Every fibre r-c is squarefree: the two partial derivatives1+2xy andx^2
are coprime. Its total degree is3, including at c=0 where it is reducible.

There is no g in B with Dg in C*. Indeed g=kappa/x+h(r). In a reduced
h=A/B0, Bezout implies B0(r) divides x; degree forces B0 constant,
contradicting the pole1/x. If one component has positive z-degree m and
constant nonzero leading coefficient, rescale it to monic and put
S=P^(-1/m)=z^-1*(1+O(z^-1)). The exact formal identity

```
D(Q*P_z*S) - d_z(Q*D(P)*S) = -kappa*S
```

gives D(res_z Q*P_z*S)=-kappa, a polynomial slice contradiction. If either
component is independent of z, the bracket is a product of polynomials
and likewise forces a slice. Thus both z-degrees are positive and both
leading coefficients nonconstant.

If both components are affine in z, their leading coefficients p,q have
no common zero, as otherwise the z Jacobian column vanishes. The highest
coefficient gives D(p/q)=0. Write p/q=A(r)/B0(r) reduced; polynomiality
and Bezout give p=eta*A(r),q=eta*B0(r), and absence of common zeros makes
eta a nonzero constant. Therefore Dp=Dq=0 and the constant coefficient
is D(q*a-p*b)=kappa for P=p*z+a,Q=q*z+b, again impossible. A vanishing
leading entry is handled by the independent-of-z obstruction.

## New fibre-factor lemma

Suppose Q=q*z+b is affine in z, P has z-degree m>=2 and leading coefficient
p. If r-c divides q, then r-c divides p.

For each irreducible factor f of the squarefree r-c, D(f) belongs to(f),
by differentiating its factorization. Hence D descends to the domain B/(f).
Modulo f the component Q equals b and is independent of z, so the constant
bracket is -P_z*D(b). Their product is the nonzero constant kappa in
(B/(f))[z]. A unit of a polynomial ring over a domain has z-degree0.
Thus the leading z coefficient m*p of P_z vanishes modulo f. Characteristic
zero gives f|p. Applying this to every factor gives r-c|p. This covers the
reducible special fibre0 without assuming it irreducible.

## Complete leading-degree proof

Order m=deg_z P>=n=deg_z Q. Nonconstant leading coefficients imply m,n<=4.
Write them p,q, with deg p<=5-m and deg q<=5-n. The highest bracket term gives

```
n*q*D(p)-m*p*D(q)=0,       p^n/q^m=A(r)/B0(r).
```

A root c of A forces r-c|p; a root c of B0 forces r-c|q. This uses Bezout
for A,B0 and squarefreeness, not irreducibility, of each fibre.

First consider equal degrees. For m=n>=3, both coefficient degrees are
at most2, so a nonconstant ratio in C(r) is impossible. For m=n=2,
a nonconstant ratio forces one degree<=3 coefficient to be a scalar
multiple of r-c, hence invariant, and the differential relation forces
the other to lie in C[r] as well. Both are affine nonconstant polynomials
in r. A constant linear target combination cancels their r terms; its
leading coefficient is either a nonzero constant, already excluded, or
vanishes. A constant ratio also permits cancellation. Thus a constant
linear target shear lowers m+n while preserving total degree<=5.
Repeat until unequal degrees; equality1 is already the affine-pair case.

For m>n with m>=3, deg p<=2. Thus A must be constant. If B0 has a root,
then r-c|q. When n>=2, deg q<=3 gives q=scalar*(r-c), so Dq=0 and then
Dp=0, impossible for a nonconstant p of degree<=2. When n=1, the new fibre
lemma gives r-c|p, directly contradicting its degree bound. Therefore
A/B0 is constant. Unique factorization gives
p=alpha*l^(m/g),q=beta*l^(n/g), with g=gcd(m,n). For every unequal pair
with m in{3,4}, the integer m/g exceeds5-m, contradicting deg p<=5-m.
All five such pairs are explicitly listed in structural-validation.json.

It remains(m,n)=(2,1), with deg p<=3 and deg q<=4. If A/B0 is nonconstant,
a numerator root forces r-c|p and hence p=scalar*(r-c); a denominator
root first gives r-c|q and then the fibre lemma again gives r-c|p.
In either event Dp=0, the leading relation gives Dq=0, and p,q are affine
nonconstant polynomials in r. Choose the root of q. The fibre lemma
then makes p proportional to q, say p=lambda*q.

Write P=lambda*q*z^2+a*z+c,Q=q*z+b with Dq=0. The z coefficient forces
D(a-2lambda*b)=0, hence a=2lambda*b+h(r). The constant coefficient becomes

```
kappa = q*D(c)-(2lambda*b+h)*D(b)
      = D(q*c-lambda*b^2-h*b),
```

again a polynomial slice contradiction. If A/B0 is constant instead,
p=lambda*q^2. Subtract lambda*Q^2 from P to make both components affine
in z, contradicting the earlier obstruction. This last nonlinear target
shear can increase total degree; the affine-pair obstruction has no degree
bound, so the increase loses no hypothesis. All cases are closed.

## Consequence for an arbitrary small cubic output combination

Let F be any3D Keller map of degree at most5 with a nonzero output combination
R of degree<=3. For degree<=2 use GLDC5. For degree3, the scalar part of the
[written quartic audit](INVERTIBILITY_AUDIT.md) uses only absence of critical
points and exactness of the regular fibre2-form. The latter follows from
Jac(P,Q,R)=kappa in every degree. Its nodal-period, cusp, tangent and binary
classification therefore applies unchanged, with no hidden use of a degree
bound on P,Q. It gives either a coordinate completion of R with inverse
degree<=5 or the exceptional scalar above. In the former case the plane
fibres have degree<=25 and are excluded by Moh100; the latter case is now
excluded by this proof. This is an explicit extension of that scalar
classification argument, not an import of a stronger stated proposition.

After the separate quartic theorem, this closes exactly the three additional
profiles(5,3,3),(5,4,3),(5,5,3). The combined current strict census has54of74
profiles excluded and20 remaining. A degree5 counterexample would need
smallest invariant flag entry>=4. None is constructed, and ordinary degree
in3D remains in[5,7], not proved equal5. Positive witness bounds are unchanged.
