# A marked gcd criterion and a finite purity prefix for every odd start

Date: 20 November 2026. Work over ANY characteristic-zero field K.
Fix an odd integer a>=1 and a regular normalized unit-center phase

 p=t^a(C(x)+t²D(x))/x,
 C,D in K[x], C(0)=-D(0)=1,

in K[x,y,t]/(t²+2xy-1). Purity means vanishing of EVERY positive
pure angular moment. This theorem concerns exactly this two-height
profile, not arbitrary phases.

The main test uses only rational functions built from C,D. It does
NOT import all-order rationality of the selected square root, any
Luroth theorem, rational return on the original line, or a polynomial
quotient classification.

## 1. Rational data and the marked correspondence

Define

 f=x²(-D)^a/C^(a+2)=x²+O(x³),
 W=(-D/C)^((a-1)/2)/C=1+O(x),
 H=2W/f'=1/x+O(1).                                  (G1)

These functions are rational regardless of whether X=sqrt(-D/C)
or u=sqrt(f) is rational. Let u=x+O(x²) be the selected FORMAL square
root of f, and v its formal inverse. Since f'=2u u',

 T=W/u'=uH, Phi(s)=u(v(s))*H(v(s)).

The formal endpoint identity, already proved for every odd height,
is

 Lambda(p^m)=kappa_(a,m)[s^m]Phi(s) for even m>=2,

with nonzero rational kappa; odd moments vanish by t parity.
Thus purity is exactly Phi(s)+Phi(-s)=2.

Write reduced polynomial fractions

 f=A/B, H=M/N,

where B(0)!=0, A has exact order2 and A_2/B(0)=1.
The denominator N has an exact SIMPLE zero at0 and M(0)/N'(0)=1.
Define polynomials in K[x,Y]:

 P=[A(Y)B(x)-A(x)B(Y)]/(Y-x),
 Q=A(x)[M(x)N(Y)-M(Y)N(x)]²
                         -4B(x)N(x)²N(Y)².          (G2)

The quotient in P is a literal polynomial. At the marked point,

 P(0,0)=0, P_Y(0,0)=A_2 B(0)!=0.                     (G3)

Let G be a primitive-in-Y representative in K[x,Y] of the monic
gcd of P,Q in K(x)[Y]: clear coefficient denominators and REMOVE
their K[x]-content. It is determined up to a scalar in K*.
The evaluation G(0,0)=0 is therefore well-defined. A constant gcd
is represented by1.

**Theorem OddStartMarkedGcdPurity.** The phase is all-order pure if
and only if

 G(0,0)=0.                                           (G4)

The criterion requires the MARKED evaluation, not just a nonconstant
global common factor.

## 2. Proof of the marked criterion

By(G3) the formal implicit equation P(x,j(x))=0 has a unique branch
j(x)=-x+O(x²) through0. The other local factor Y=x was removed in(G2).
It satisfies f(j)=f(x). Thus u(j)=-u(x), by leading terms and the
unique selected square root, and v(-u(x))=j(x).

Put

 S(x)=u(x)(H(x)-H(j(x))).
 
Then S=Phi(u(x))+Phi(-u(x)) has constant term2. Direct substitution
in(G2) gives the exact formal identity

 Q(x,j)=B(x)N(x)²N(j(x))²*(S(x)²-4).                  (G5)

The polynomial Q is not identically zero: its diagonal value is

 Q(x,x)=-4B(x)N(x)^4 !=0.                             (G6)

Since S+2 has nonzero constant4, Q(x,j)=0 iff S=2, which is precisely
purity.

It remains to justify that Q(x,j)=0 iff(G4), including specialization
atx0. The branch j is algebraic over K(x), since P(x,j)=0.
If Q also vanishes, its minimal polynomial divides both P,Q over
K(x), and hence divides their gcd. Clearing denominators gives
G(x,j)=0; taking the constant term gives G(0,0)=0.

Conversely Gauss's lemma gives EXACT polynomial divisibilities

 P=G P1, Q=G Q1 in K[x,Y].                            (G7)

If G(0,0)=0, equation(G3) implies

 G_Y(0,0)P1(0,0)=P_Y(0,0)!=0.

Hence G has a unique simple formal branch through0, and that branch
is the already unique j of P. Since G divides Q, Q(x,j)=0.
No issue from a degree drop in Y atx0 occurs.

Primitivity in Y is essential. Multiplication of a representative
by x would force its value at(0,0) to vanish artificially. Removing
content prevents this error. Q may have content; Gauss's lemma is
applied to the primitive G and does not require Q itself primitive.

## 3. An explicit finite prefix

Put

 n=degree of the rational map f=max(degA,degB)>=2,
 h=degree of the rational map H=max(degM,degN)>=1,
 D0=(n-1)(n+4h).                                     (G8)

**Theorem OddStartFinitePurityPrefix.** If the phase is not pure,
some EVEN exponent m satisfies

 2<=m<=D0-4, Lambda(p^m)!=0.                          (G9)

Consequently vanishing of all positive moments through D0 proves
all-order purity. In fact the even moments through D0-4 suffice.
This is a PROVED degree-dependent cutoff in the stated class,
not a finite-prefix experiment or an assumption of rational u.

Proof. In the nonpure case(G4) gives G(0,0)!=0. Divide BOTH P,Q
by G as in(G7). Their residual polynomials P1,Q1 are coprime over
K(x); no unrelated off-base common component is left. Moreover
G(x,j(x)) is a formal unit, so Q1(x,j) has exactly the order of Q(x,j).
The branch of P1 at the origin remains simple.

The bidegree bounds are

 deg_x P,deg_Y P <= n-1,
 deg_x Q <=n+2h, deg_Y Q <=2h.                       (G10)

The first follows by dividing the antisymmetric fiber polynomial
by Y-x in either variable. Polynomial division by G cannot increase
either degree, so the same bounds hold for P1,Q1.

Their resultant

 R(x)=Res_Y(P1,Q1)

is a NONZERO polynomial because P1,Q1 are coprime over K(x).
The Sylvester degree estimate gives

 deg_x R <= (deg_Y P1)(deg_x Q1)
               +(deg_Y Q1)(deg_x P1)
          <=(n-1)(n+4h)=D0.                         (G11)

The polynomial Sylvester/adjugate identity provides F1,F2 in K[x,Y]
with

 R=F1 P1+F2 Q1.

Substitute Y=j(x). The coefficients F2(x,j) have nonnegative x-order;
therefore

 ord_x Q(x,j)=ord_x Q1(x,j)
             <=ord_x R<=deg_x R<=D0.                (G12)

There is no division by a leading coefficient that might vanish
at the origin. If Q1 has Y-degree0 the same argument uses
R=Q1^(deg_Y P1) and its direct polynomial Bezout identity.

Suppose m is the first positive even index with [s^m]Phi(s)!=0.
Then S-2 has exact x-order m, since u=x+O(x²). The factor
B N(x)²N(j)² in(G5) has exact order4; S+2 is a unit.
Consequently

 ord_x Q(x,j)=m+4.

Together with(G12) this proves(G9). Conversely absence of such a
first nonzero coefficient means all-order purity. The nonzero
endpoint factors identify the same first index for angular moments.
Since n>=2,h>=1, D0>=6 and the bound is meaningful. QED.

## 4. Exact recognition and the zero center

For a unit-center input, compute the reduced rational functions(G1),
form(G2), compute the gcd over K(x), clear denominators, remove
content and evaluate(G4). These are finitely many exact rational
polynomial operations. Over a field with effective arithmetic and
equality they are a finite algorithm. Over an abstract field they
are finite algebraic conditions.

No rational-square test for X is needed. In particular the test
also rejects nonpure inputs for which -D/C is not a rational square.

For a general regular pair C(0)+D(0)=0 with nonzero C(0), first
divide C,D by C(0), which preserves all vanishing questions.
If C(0)=D(0)=0, the retained all-height zero-center theorem already
gives the separate exact criterion: purity iff x² divides both C,D.
Every such zero-center pure phase is safe. This completes purity
recognition for this two-height class at every fixed odd a>=1;
it does not decide purity of wider height profiles.

## 5. Off-base components and finite-control scope

A common component can be completely unrelated to the marked branch.
For a rational-data control, take

 u=x(x-1)(x-2)/2,
 T=1+u(x-1)²,
 f=u², H=T/u.

Here f=x²+O(x³),H=1/x+O(1), so the marked correspondence theorem
applies. The involution x->2-x sends u->-u and T->2-T, producing
a nonconstant common factor through Y=2-x. This component misses
(0,0). In fact the selected endpoint has

 Phi(s)=1+s-2s²+O(s³),

so it is not pure. The exact gcd certificate verifies that the
primitive G has G(0,0)!=0. This example is labeled rational
correspondence data; it is not asserted to come from a polynomial
C,D phase.

The finite inputs also include actual polynomial phases with
nonrational selected X, pure T=1 cases for which Q remains nonzero,
an input with second moment0 but fourth moment nonzero, and a
separate constant-in-Y resultant control. Every discovered gcd is
retained as polynomial factors and a polynomial Bezout identity
certifying coprimality of the residuals. A standard-library replay
multiplies all exact identities. Symbolic algebra may be used to
discover those certificates, not as an unrecorded oracle.

## 6. Classical context and proof dependencies

The essential imported formula is the retained formal endpoint
identity and nonvanishing of its even coefficients' rational factors.
The marked implicit branch is proved from(G3); Gauss's lemma and
the Sylvester polynomial identity are the standard algebraic steps.

Milne, *Fields and Galois Theory*, v5.10, Gauss's lemma and its
primitive-divisibility consequence (Section9, printedp118) were
inspected in the retained primary source. No Luroth theorem is used
in this proof.

A related classical result is Batenkov--Binyamini, *Uniform upper
bounds for the cyclicity of the zero solution of the Abel differential
equation* (2015), Introduction, Definition1, Remark2 and Theorem1,
printedpp2--3. It treats polynomial P,Q over C with fixed complex
endpoints and the moments integral P^k Q P' or integral P^k q.
Their moment Bautin-index bounds are d_Q+3(d_P-1)^2 and
2+d_q+3(d_P-1)^2. Those polynomial-path hypotheses are not identified
with this angular constant-coefficient functional, so their bound is
context, not an imported premise. No priority claim about finite
moment criteria in general is made.

All source editions, exact input fractions, literal gcd/Bezout
certificates and checker bytes are pinned in the receipt. Internal
exact replay does not constitute proof-assistant verification,
external human review or a historical-priority claim.


The [checker](verify_odd_start_marked_gcd_finite_purity.py),
[literal certificates](odd_start_marked_gcd_certificates.json),
[source ledger](odd_start_marked_gcd_source_ledger.json) and
[receipt](odd_start_marked_gcd_verification.json) retain the exact replay.
