# Cusp defect, finite recognition and rigidity between two cusp filtrations

14 November 2026. Claims JC2CFR1 (written reproduction), JC2CFF1
(finite recognition and two-filtration extension), and JC2CFF2 (F2 contact gate). Work over C; the same
algebra works over an algebraically closed characteristic-zero field.

The external input is the public JC2 campaign's coprime cusp-filtration
argument, attributed to swarmHQ and pinned at commit
15dcad7b1e65e8cb3d39e9aaac1d4c44063ec4fc in
[sources/manifest.json](sources/manifest.json).
Its source calls the result MANUAL / PROVISIONAL, with novelty unassessed.
Sections 1--3 below reproduce its algebra in full. Sections 4--6 are the
local follow-up: a finite recognition criterion, a strict cusp-defect
inequality, and compatibility between different binomial cusp filtrations.
No literature-priority claim is made.

All-degree assertions have written proofs. The [checker](verify.py) replays
finite identities and controls; it does not certify those proofs, classify
arbitrary valuations, or supply an arbitrary-Keller reduction.

## 1. A genuine source-ring filtration

For coprime a,b>=2 put h=ab-a-b>0 and w=x^a-y^b.
Every polynomial has a unique normal expression

    sum_{0<=i<a, j,k>=0} c_ijk x^i y^j w^k.

Define delta=delta_(a,b) as its maximum weight bi+aj+hk, with
delta(0)=-infinity. This is a positive semidegree: delta(fg)=delta(f)+delta(g),
delta(f+g)<=max(delta(f),delta(g)), and only constants have degree zero.
Indeed the associated graded ring is

    A_s=C[X,Y,W]/(X^a-Y^b), weights (b,a,h).

The unique normal basis proves there are no other initial relations;
coprimality makes the displayed binomial prime. Its fraction field is
C(z,W), with X=z^b,Y=z^a and z of degree one.

For the ordinary Jacobian bracket {f,g}=f_x*g_y-f_y*g_x,

    {x,y}=1, {y,w}=-a*x^(a-1), {x,w}=-b*y^(b-1).

The normal basis and Leibniz rule show
delta({f,g})<=delta(f)+delta(g). The induced degree-zero bracket has

    {X,Y}=0, {Y,W}=-a*X^(a-1), {X,W}=-b*Y^(b-1),
    {z,W}=-z^(h+1).

Its fraction-field Poisson centre is C: commuting with z first removes
W-dependence, then commuting with W removes z-dependence. Characteristic
zero is used for both derivative conclusions.

The polynomial Euler derivation E_s=b*x*partial_x+a*y*partial_y preserves
delta and has divergence a+b. Its eigenvalue on w is ab, whereas delta(w)=h.
Those are different weights.

## 2. Degree-preserving derivations form one line

Here is the general lemma reproduced from the source. Suppose a positive
integer semidegree on C[x,y] has a domain graded ring, the Jacobian bracket
has degree at most zero, and its graded fraction field has Poisson centre C.
Divergence is injective from the vector space of filtration-preserving
polynomial derivations to C.

Choose u,v whose leading degree-zero bracket is nonzero. For a preserving
derivation L, the identity

    L{u,v}={Lu,v}+{u,Lv}-(div L){u,v}

gives delta(div L)<=0 if divergence is nonzero, so divergence is constant.
If div L=0, termwise polynomial integration gives H with L={H,-}.
For L nonzero, H is nonconstant and delta(H)>0. Its nonconstant leading
form cannot be Poisson central, so for some homogeneous leading form fbar,

    delta(Lf)=delta(H)+delta(f)>delta(f),

contrary to preservation. Thus L=0. Since E_s is a preserving derivation
with nonzero divergence, every preserving derivation is a scalar multiple
of E_s. This uses the actual polynomial source ring for Hamiltonian
integration; a rational or completed-local source would not justify it.

## 3. Transport between two binomial cusp filtrations

Let the target cusp have coprime p,q>=2, s=pq-p-q and target coordinates u,v,
with w_t=u^p-v^q, weights delta_t(u,v,w_t)=(q,p,s).
Let F^*u=P,F^*v=Q be a polynomial Keller pair, J={P,Q} in C*.
Assume for the moment that

    delta_s(F^*f)=e*delta_t(f) for every nonzero f, e>0.       (3.1)

Since delta_s(P)=eq and delta_s(Q)=ep are integers and gcd(p,q)=1,
e is a positive integer. Leading forms give an injective graded map
phi:A_t->A_s multiplying degrees by e.

The Keller condition supplies polynomial derivations on the whole source:

    D_P=(Q_y*partial_x-Q_x*partial_y)/J,
    D_Q=(-P_y*partial_x+P_x*partial_y)/J,
    D=q*P*D_P+p*Q*D_Q.

They satisfy D(F^*f)=F^*(E_t f), where E_t=q*u*partial_u+p*v*partial_v,
and div D=p+q.

Every polynomial derivation L has finite filtered order

    r=max(delta_s(Lx)-b, delta_s(Ly)-a, delta_s(Lw)-h).

This follows by differentiating the normal basis; a nonzero L attains
its order on one of the three generators. If D had r>0, its nonzero
leading derivation would kill phi(A_t), by (3.1) and preservation by E_t.
But Frac(A_s) is finite algebraic over Frac(phi(A_t)): both are finitely
generated of transcendence degree two. It is separable. Differentiating
minimal polynomials shows that a derivation zero on the latter field is
zero on the former. This contradicts the nonzero leading derivation.

Hence D preserves delta_s. Section 2 gives D=lambda*E_s, where
lambda=(p+q)/(a+b). Consequently P,Q are E_s eigenpolynomials with positive
weights q/lambda,p/lambda. They have no constant terms. J(0)!=0 implies
their linear parts are independent, so these weights must be b,a in one
order. The only monomial of E_s weight b is x, and the only one of weight a
is y, because a,b>=2 are coprime. Primitive coprime exponent pairs then
force lambda=1 and precisely one of:

    (p,q)=(a,b), P=alpha*x, Q=beta*y, alpha^a=beta^b;       (D)
    (p,q)=(b,a), P=alpha*y, Q=beta*x, alpha^b=beta^a.       (S)

The coefficient conditions follow by applying (3.1) to w_t; an unmatched
coefficient gives degree ab rather than h. In either case e=1.
Conversely (D) and (S) give (3.1). In (S), J=-alpha*beta; in (D),
J=alpha*beta.

For identical source and target cusps this reproduces the external
rigidity theorem JC2CFR1. The two-filtration conclusion is part of
JC2CFF1: exact scaling between different coprime binomial cusp types is
impossible, even allowing source/target interchange. This is still a
special class of source and target filtrations.

## 4. A sharp cusp-Jacobian bound and finite recognition

Keep both cusp filtrations. In this section P,Q may be ANY algebraically
independent polynomials: J={P,Q} is nonzero but need not be constant.
Put

    A=delta_s(P), B=delta_s(Q), C=delta_s(P^p-Q^q),
    m=max(A/q,B/p), j=delta_s(J).

A,B,C are positive, and j>=0. Then the cusp-Jacobian inequality is

    C >= s*m+j.                                           (4.1)

If p*A!=q*B, the inequality is strict. It is sharp on the equal-slope
stratum, including non-Keller examples.

Indeed, without equal slopes there is no leading cancellation, so
C=max(p*A,q*B)=p*q*m. The degree-zero bracket bound gives
j<=A+B<(p+q)*m. Since p*q-s=p+q, this proves the strict inequality.
With equal slopes A=e*q,B=e*p, coprimality gives an integer e>=1.
The exact identity

    {Q,P^p-Q^q}=-p*J*P^(p-1)

gives B+C>=j+(p-1)*A, hence C>=e*s+j, proving (4.1).

More strongly, the ZERO cusp defect is completely classified:

    C=s*m  if and only if (D) or (S) of Section 3 holds.    (4.2)

To prove this, (4.1) first forces j=0. Positivity of the source filtration
makes J a nonzero constant: the finite degree condition itself implies
the Keller condition. Strictness off equal slopes gives A=e*q,B=e*p,C=e*s.
Leading forms U,V,H define a graded map A_t->A_s:
U^p=V^q since e*s<e*p*q. The top-degree bracket identity is

    {V,H}=-p*J*U^(p-1)!=0.                                (4.3)

In a characteristic-zero Poisson field, algebraically dependent elements
have zero bracket: one element is algebraic separable over the rational
function field of the other, and differentiation of its minimal polynomial
applies. Thus V,H are algebraically independent. The graded map has image
of transcendence degree two. Since A_t is a two-dimensional finitely
generated domain, its prime kernel is zero.

Injectivity means that the leading form of every normal polynomial remains
nonzero under substitution. This proves universal scaling (3.1) from the
three degrees. Section 3 proves (D) or (S), with e=1.
Conversely these classified maps give zero defect.

An equivalent two-polynomial recognition certificate is

    delta_s(P)=e*q and delta_s(P^p-Q^q)<=e*s,               (4.4)

for any real e>0, provided P,Q are algebraically independent. The second
bound forces delta_s(Q)=e*p by leading cancellation, hence e is an integer.
Inequality (4.1) forces j=0 and equality in the bound. Thus (4.4) forces
(D) or (S), e=1, and scaling of EVERY polynomial.

The finite bridge applies to two different binomial cusp types as well as
to a single type. It needs only a nonzero Jacobian polynomial; the degree
bound then forces that Jacobian to be constant. It supplies no general
reason why a candidate should attain zero defect.

The all-order power pairs calibrate sharpness. For P=x^n,Q=y^n, n>=1,
and identical source/target cusp (a,b),

    m=n,
    delta_s(P^a-Q^b)=n*h+(a+b)*(n-1),
    delta_s(J)=(a+b)*(n-1).

Thus (4.1) is an equality for every n. For n>1 its positive defect records
the nonconstant Jacobian; it is NOT in the zero-defect rigidity class.
The checker replays representative n, while the displayed binomial
expansion proves all n.

## 5. Exact graded symbols do not supply polynomial lifts

The following injective Poisson maps have exact graded degree scaling,
but Section 3 excludes a full polynomial Keller lift with that scaling.
Here the source always has (a,b)=(2,5), A_s=C[X,Y,W]/(X^2-Y^5).

* Same cusp (p,q)=(2,5), e=3:

      U=X^3, V=Y^3, H=Y^3*W/3.

* Different cusp (p,q)=(2,3), e=5:

      U=X^3, V=Y^5, H=Y*W/5.

* The (3,5) output-degree ratio, (p,q)=(5,3), e=65:

      U=X^39, V=X^65, H=Y^226*W/65.

In all three cases U^p=V^q, {U,V}=0,
{V,H}=-p*U^(p-1), {U,H}=-q*V^(q-1).
A nonzero bracket proves injectivity exactly as in Section 4.
The last example is a symbol calibration for the F2 carrier ratio:
for C_R=x*(x*y^5-1)^2*R(x*y^5) with R of degree two,
delta_(2,5)(C_R)=65, so C_R^3,C_R^5 have degrees 195,325.
This DOES NOT identify the actual F2 source valuation with this cusp
semidegree, prove the bounds for its missing lower bands, or exclude F2.

Translations also show why P,Q degrees alone are insufficient:
P=x+1,Q=y has J=1 and delta_(2,5)(P,Q)=(5,2), but
delta_(2,5)(P^2-Q^5)=5>3. The polynomial pair x^3,y^3 is a further control:
its component degrees have scaling e=3 but its Jacobian is nonconstant
and its cusp-difference degree is 23>9.

## 6. A concrete contact gate for the retained F2 row

The retained PF2L1 B0 envelope has actual original source monomials

    P: i+j<=75, 5*i-j<=15, i,j>=0;
    Q: i+j<=125, 5*i-j<=25, i,j>=0.

This is read from [the retained band derivation](../../archive/non-elliptic/plane-jc/F2_75_125_DERIVATION.md)
and the complete intervals in its
[classifier source](../../archive/non-elliptic/plane-jc/cas/classify_f2_75_125_layers.py).
Those necessary halfspaces are retained without claiming that the envelope
is an exhaustive normal form or that every allowed monomial occurs.
Their respective source supports have 706 and 1901 lattice points.

For the deliberately chosen auxiliary GLOBAL source filtration
delta_s=delta_(2,5), each source monomial has degree 5*i+2*j.
The two inequalities for P give i<=15 and

    5*i+2*j=2*(i+j)+3*i<=195,

with equality only at (15,60). That is a forced nonzero vertex.
For Q the same argument gives the unique maximum 325 at (25,100).
Since no other source monomial can have those weights, cancellation after
cusp normalization is impossible. Therefore every actual pair in this
declared F2 envelope has delta_s(P)=195, delta_s(Q)=325.

Use the retained PF2CW1 normalization with carrier coefficients

    P_lead=C_R^3, Q_lead=-(9/5)*C_R^5,

and put Qtilde=-(5/9)*Q. The target cusp (p,q)=(5,3) has s=7 and e=65.
The precise polynomial to inspect is

    W_F2=P^5-Qtilde^3=P^5+(125/729)*Q^3.

For any actual Keller pair with these supports, Section 4 first gives
delta_s(W_F2)>=455. Equality would give exact scaling between source cusp
(2,5) and target cusp (5,3), contradicting Section 3.
Thus the adopted F2 contact gate, JC2CFF2, is

    delta_(2,5)(P^5+(125/729)*Q^3) >= 456.                 (6.1)

This is a universal necessary inequality in the declared support stratum.
It is NOT an exclusion of F2. Any rigorous upper bound <=455 for the
complete polynomial W_F2 would now exclude that candidate. If its actual
degree is >=456 the new test leaves it alive. Leading carrier cancellation
and finite jets do not establish such an upper bound.

The auxiliary cusp filtration is defined on the actual polynomial source
coordinates. It need not be identified with an existing dicritical or
endpoint valuation to apply this algebraic result. The substantive open
task is to bound the FULL cusp difference using actual Keller equations,
not to postulate that a local target cusp supplies source scaling.

The ordinary A2 cusp in the r7 target curve and its normalization degrees
(6,10) concern a curve of transcendence degree one. They do not themselves
supply a source-ring filtration of transcendence degree two.
An r7 application still needs actual source and target polynomial
coordinates and a whole-polynomial contact estimate.

For other surviving configurations, identify the actual global source
cusp filtration, target cusp type and scalar-normalized cusp difference;
derive complete-support degree bounds and compare the defect to delta(J).
Zero defect gives the rigid classification; positive defect remains open.
A central graded map, a local cusp germ or finite jets cannot replace
that whole-polynomial estimate.

The external source's full Rees-lift corollary is also valid: a polynomial
homogeneous lift with injective central Poisson map has unit Jacobian
after specializing its Rees parameter to one, and that injectivity gives
(3.1). The source proves the Jacobian assertion by divisibility by the
prime Rees parameter and a weight-zero quotient. This does not hold for
a merely formal or rational lift. Our explicit graded controls show why
an algebraization premise is substantial.

No JC2, HC4, normal-r7, geometric-degree-four, or (75,125) closure,
new degree frontier, independent external review, formal verification,
or novelty priority is claimed.
