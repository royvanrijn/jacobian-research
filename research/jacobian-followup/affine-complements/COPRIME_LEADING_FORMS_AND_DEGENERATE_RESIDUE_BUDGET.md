# Coprime leading forms and a degenerate residue budget

PF2RULED107, 16 November 2026. This closes the stated local class at ALL
indices and degrees. It does not classify arbitrary degenerate pencils
or prove the full Wright goal. Retain the original finite normalization,
the WHOLE marked plane open, and EVERY omitted prime as in
[PF106](CYCLIC_QUOTIENT_ORDER_AND_MORSE_RESIDUE_BUDGET.md).

## 1. Hypothesis and theorem

Let h=(P,Q):(C2,0)->(C2,0) be finite, with sole reduced critical curve
A smooth, and Jac(h)=unit*X^(e-1) in coordinates A={X=0}.
Write r=ord(P), s=ord(Q), and U=P_r,V=Q_s for their ordinary homogeneous
leading forms. Assume U,V have NO COMMON ZERO in P1. This is a substantive
hypothesis on a chosen target coordinate pair; finiteness of h alone
does not imply it. Set g=gcd(r,s), a=r/g, b=s/g. Let q be the local
degree of A onto the normalization of its image.

Then

    e=r+s-1,   deg_0(h)=r*s,   q divides g,
    e <= r*s/g.                                             (1)

If q>1, necessarily

    a,b>=2,  gcd(a,b)=1,  a+b>=5,  a*b>=6,
    q divides e+1,  e>=5q-1,  deg_0(h)>=6q^2.                (2)

For an ACTUAL cyclic quotient boundary point of order n>1 on a Keller
normalization, the retained power-map LOCATION gives q=n*f, where f
is the global omitted-prime residue degree. Consequently

    marked quotient type 1/n(1,1),
    n*f divides e+1,  e>=5*n*f-1,
    ell=r*s/n>=6*n*f^2,
    d>=6*n*f^2+1>=12*f^2+1.                                (3)

Here ell is the DOWNSTAIRS ambient local degree and d the GLOBAL field
degree. All other omitted primes and finite excesses remain in the
budget. In this class f>1 requires d>=49 and e>=19; even f=1 requires
d>=13. These statements do not exclude every map of smaller degree:
they force a COMMON PROJECTIVE ZERO of the leading forms at any
remaining singular boundary residue point in those ranges.

## 2. Proof using a projective pencil

First J(U,V) is nonzero. Otherwise the homogeneous rational function
U^(s/g)/V^(r/g) has zero derivative and is constant, hence U^s=c*V^r. Its
positive-degree factors would give a common projective zero. The first
nonzero homogeneous part of Jac(P,Q) is therefore J(U,V), of degree
r+s-2. The assumed smooth critical support implies

    J(U,V)=C*X^(r+s-2), C!=0,   e=r+s-1.                   (4)

Both U(0,1) and V(0,1) are nonzero when r,s>=2. Indeed, if U vanishes
there with order t>=1, coprimality makes V nonzero there. In the chart
Y=1 Euler's identity gives

    J(U,V)(X,1)=s*U_X(X,1)*V(X,1)-r*U(X,1)*V_X(X,1).

Its order at X=0 is t-1<=r-1<r+s-2, contradicting (4). Interchange U,V
for the other case. If r=1 and s>=2, write U=cX+dY. Were d=0,
coprimality would force V(0,1)!=0, and J(U,V)(0,1)!=0 would contradict
(4). Thus d!=0, so P|A has order1 and q=1. The s=1 case is symmetric.
The r=s=1 case is unramified and outside a nonempty critical divisor.
Thus q>1 forces r,s>=2. This also agrees with the stronger all-index
corank-one primitive restriction lemma in PF105. An actual quotient
case has no invariant linear forms in the first place.

Thus ord(P|A)=r and ord(Q|A)=s. A primitive parametrization of the image
normalization followed by its degree-q map from A shows q divides BOTH
r and s. Only divisibility is asserted; leading orders alone are not
used to claim equality for a nonhomogeneous germ.

The morphism

    R:P1->P1, [X:Y] |-> [U(X,Y)^b : V(X,Y)^a]

has degree N=r*s/g: its two homogeneous coordinates have this degree
and no common projective zero. At [0:1] both are nonzero. Differentiating
R in the X/Y coordinate and using (4) shows R' has order e-1 there.
The local degree of R is therefore e. A point's local degree cannot
exceed the whole degree N. This proves the last inequality in (1);
no classification of arbitrary discriminants is imported.

For completeness the local intersection multiplicity of P=Q=0 is r*s:
coprime tangent cones leave no intersection on the exceptional line of
the first blowup, and the blowup intersection formula contributes r*s
at the origin. Equivalently their homogeneous leading forms form a
regular sequence with quotient Hilbert series
(1-T^r)(1-T^s)/(1-T)^2, of length r*s, and the filtered local complete
intersection has that length. This is the local degree of the finite h.

Now q>1 implies g>=2. If a=1, inequality (1) becomes

    g*(1+b)-1 <= g*b,

forcing g<=1, a contradiction; similarly b!=1. Coprimality then gives
a+b>=5 and a*b>=6. Since q|g, e+1=g*(a+b) is divisible by q, and the
remaining inequalities in (2) follow. These are unbounded arguments,
not consequences of a search through degrees.

At an actual singular point PF106 gives marked type1/n(1,-e). Now
n|q|e+1, so -e=1 mod n and the type is 1/n(1,1). The quasi-etale
uniformization has degree n, giving ell=r*s/n. The original whole-plane
Euler budget is d-1>=e*f+eta=ell at this point. Substituting q=n*f gives
(3). No ordinary unramified companion is relabeled an omitted prime.

In particular this singular class has r,s>=4, so BOTH quadratic parts
vanish. It supplies a separate complete sector of the degenerate case;
it does not establish a Morse pullback or extend nf|(e-1) outside PF106.

## 3. Genuine finite non-Gorenstein controls at every scale

There are examples attaining BOTH local bounds in (2). The existence
input is the rational-function passport theorem of Zannier, stated as
Theorem1 on full PDF page2 of Orevkov's author manuscript
[On germs of mappings C2->C2](https://www.math.univ-toulouse.fr/~orevkov/k-en.pdf).
Apply it to 2q entries equal3 and 3q entries equal2. Their sums are6q,
their common gcd is1, and the special local degree is5q-1<=6q.
After permuting target values and putting that special source point at
infinity, it gives coprime polynomials p(t),z(t) with

    deg p=2q,  deg z=3q,
    deg(p^3-z^2)=q+1.                                    (5)

This application uses the passport theorem, NOT the germ classification:
the latter assumes bijectivity on EVERY curve over the discriminant,
an assumption unavailable for our nontrivial residues. Also Remark2's
coordinate description has an apparent value-label inconsistency:
its equal-degree P/Q tends to1 at infinity, despite the preceding
assignment to infinity. We use Theorem1 and make the above permutation
explicit; no argument relies on that sentence of Remark2.

Homogenize p,z to P=X^(2q)*p(Y/X), Q=X^(3q)*z(Y/X).
Their leading coefficients are nonzero and they have no common projective
zero. Hence (P,Q) is a finite homogeneous polynomial map of degree6q^2.
For R=p^3/z^2, the passport uses all critical multiplicity: its zeros
contribute 4q, its poles 3q, and the special point 5q-2, totaling12q-2
=2*deg(R)-2. Thus there are no other critical points. Equivalently the
polynomial 3*p'*z-2*p*z' is a NONZERO CONSTANT (it has no finite zeros).
The homogeneous Jacobian is consequently

    Jac(P,Q)=C*X^(5q-2), C!=0.

Its single critical line has transverse index5q-1 and restriction
(constant*Y^(2q),constant*Y^(3q)), of EXACT normalized residue q.
The normalized image is the cusp with primitive pair(2,3).
This proves the claimed examples for EVERY q, conditional only on the
declared published existence theorem; it is not a finite sample proof.

For each n>1 dividing q the scalar action (X,Y)->(zeta*X,zeta*Y)
preserves both outputs. The descended normal finite map has

    source quotient 1/n(1,1), e=5q-1, f=q/n,
    degree=ell=6q^2/n, eta=q*(q+1)/n.                      (6)

For n>2 these sources are NON-GORENSTEIN. For example q=6,n=3 gives
e=29,f=2,ell=72. Thus neither A-type singularities nor nf|(e-1) follow
automatically for general finite germs outside the scalar-Morse sector.
Whether the actual WHOLE plane open forces stronger conclusions remains
open. For q>=2 the examples have zero quadratic pencil; q=1 is the
separate unquotiented endpoint control. The case e=5,n=3 left
open in PF106 is NOT constructed or excluded by this higher-index family.

These are NOT KELLER COMPLETIONS. Removing X=0 and taking the quotient
gives the ring C[X^n,X^(-n),Y/X], hence complement C* times C. The local
degree equals the global degree in (6), while the actual whole-plane
budget would demand d>=ell+1. This one-sheet discrepancy is explained
by the wrong open source; it does not disprove the finite maps.

A fully explicit q=2 control is

    p=t^4+t, z=t^6+(3/2)*t^3+3/8,
    p^3-z^2=-(1/8)*t^3-9/64,
    P=Y^4+X^3*Y,
    Q=Y^6+(3/2)*X^3*Y^3+(3/8)*X^6,
    Jac(P,Q)=-(9/4)*X^8.

Its order-two quotient has e9,f1,degree12 and zero quadratic pencil.
The q=1 control P=Y^2+X^2,Q=Y^3+(3/2)*X^2*Y has Jacobian3X^3;
it checks the unquotiented endpoint and is not a new quartic Keller map.

## 4. Exact boundary of the result

The PF105 Morse counterfamily has leading forms X^2/2 and X*Y with
COMMON factor X. Its unbounded residues and much smaller local degrees
therefore do not satisfy this theorem's hypothesis. Finiteness alone
does not make leading forms coprime. Neither does a zero quadratic
pencil establish that hypothesis. Higher-order common-tangent cases,
including common normal squares, remain OPEN. The general inequality
q<=e-1 outside the two proved sectors is still UNPROVED.

The [checker](verify_coprime_leading_forms.py) and
[receipt](coprime_leading_forms_receipt.json) retain exact polynomial,
projective-degree, gcd, quotient, budget and predecessor controls.
The [source manifest](sources/coprime-leading-forms/manifest.json)
distinguishes the written elementary projective-pencil proof from the
imported passport existence theorem. No whole Orevkov/Kulikov germ
classification is imported. Formal/external verification remains false;
closure precedes review expense. All f1 packets and actual Laurent finite
algebra remain OPEN, and published Domrina2000 quartic retirement persists.
