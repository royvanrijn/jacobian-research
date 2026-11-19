# Canonical polynomial descent for pure odd-cubic phases

Date: 19 November 2026. This concerns regular unit-center winding-one odd-cubic-height phases. It is a reduction to a rational-return quotient, not a classification of all such quotients and not an assertion about arbitrary Image phases.

## 1. Statement and exact normalization

Let K be any characteristic-zero field and

    p=t/x*(C+t²D), C,D in K[x], C(0)=-D(0)=1,

have all positive pure moments zero. Let u=x*sqrt(-D/C)/C be the selected rational inverse coordinate, u=x+O(x²), let R be the rational primitive part of dx/C normalized R(0)=0, and set T=1/(C u'). Define

    E=K(u,R,T) subset K(x), g=[K(x):E].

**Theorem OddCubicCanonicalPolynomialDescent.** There is a unique polynomial w in K[x] generating E with

    degw=g, w(0)=0, w'(0)=1.

There exist Cbar,Dbar in K[z] and ubar in K(z), with Cbar(0)=-Dbar(0)=1 and ubar=z+O(z²), such that

    u=ubar(w),
    C=Cbar(w)/w',
    D=Dbar(w)*(w/x)²/w'^3.                             (P1)

The quotient phase pbar=t/z*(Cbar+t²Dbar) is all-pure. If Rbar is its normalized rational primitive part and Tbar=1/(Cbar ubar'), then

    R=Rbar(w), T=Tbar(w),
    K(ubar,Rbar,Tbar)=K(z).                            (P2)

Thus the quotient has a GENUINE rational return involution on z. The original return in x is not asserted rational: it solves w(xminus)=sigma(w(x)) and may be algebraic or rational.

For nonconstant C, writing c=degC,e=degD,H=number of distinct geometric roots of C and barred quantities for the quotient, one has

    c=g(cbar-1)+1, e=g(ebar-1)+1,
    H=g(Hbar-1)+1.                                    (P3)

In particular g divides c-1,e-1,H-1 and (c-e)/2. The last divisibility uses the same-parity consequence of rational sqrt(-Dbar/Cbar). If D dominates, d=g(dbar-3)+3 for exact least sphere degrees; if C dominates, d=g(dbar-1)+1. These are exact degree identities, not existence criteria for arbitrary polynomial pullbacks.

Constant normalized C gives only C1,D-1, with w=x and a trivial degree1 descent. For completeness, C1 gives rational X with X²=-D polynomial, so X and u=xX are polynomial. The exact-part identity makes the return j(x)=x-2u(x) a polynomial involution. Rational-map degree multiplication makes j linear; j(0)=0,j'(0)=-1 give j=-x,u=x,D=-1. All other pure phases have c>=2 and R nonconstant.

## 2. Stable field and the unique zero of the differential

R belongs to K(x). One can see this without a split-coefficient assumption: the rational primitive part is unique after normalization R(0)=0, its denominator divides gcd(C,C'), and its numerator degree is at most its denominator degree (regularity at infinity). Its coefficients solve linear equations over K imposing only simple poles on dx/C-dR. Existence over the algebraic closure and uniqueness imply a K-rational solution. Equivalently the canonical partial-fraction construction is fixed by every coefficient-field automorphism.

Retained rational-u descent and exact/log splitting imply

    sigma(u)=-u, sigma(R)=R-2u, sigma(T)=2-T.            (P4)

Formal substitution u->-u is a genuine involution of E. It need not extend to K(x). Over the algebraic closure the exact rational-map degree of R is c-H, because its pole order at each repeated C-root of multiplicity m is m-1, with no infinity pole. The field tower therefore gives

    g<=c-H<=c-1.                                      (P5)

Use Luroth's theorem over K to write E=K(w0). Its rational-map degree is g, unchanged by scalar extension. The differential

    omega=dx/C=T du

is the pullback of a rational differential eta on the w0-line. It has no finite zeros: at C-roots it has poles, and at every other finite point it is regular and nonzero. Its only possible zero is infinity, of order c-2.

Let b be the ramification index at x=infinity. The local pullback formula gives

    c-2=b(ord_(w0(infinity)) eta+1)-1,
    b divides c-1, b<=g.                              (P6)

If b<c-1, eta has a POSITIVE zero at w0(infinity). Every point in that fiber would pull back to a zero of omega. There is only one possible source point, infinity, so the fiber consists of infinity alone and b=g. If b=c-1, (P5),(P6) already force b=g. Hence the cover is ALWAYS totally ramified at infinity:

    b=g.                                               (P7)

This equality is proved from the differential; no Galois-cover assumption is used.

The value w0(infinity) is K-rational. Apply a K-Mobius change sending it to infinity. The resulting rational generator has its only pole at x=infinity, hence is a polynomial of degree g. The simple zero of u at x=0 makes the cover unramified there, so its derivative at0 is nonzero. Subtract its constant term and divide by that derivative to obtain w(0)=0,w'(0)=1.

Any two polynomial generators of E differ by a fractional-linear transformation. Since both have only a pole at infinity, that transformation is affine; the two normalizations make it the identity. Thus w is unique and lies in K[x].

## 3. Polynomial quotient coefficients and purity

After the coordinate normalization, eta has no finite zeros. A finite zero would have a nonempty finite fiber under the polynomial w, giving finite zeros of omega. Consequently

    eta=dz/Cbar(z), Cbar in K[z],

with Cbar(0)=1 from w'(0)=C(0)=1. Pullback gives C=Cbar(w)/w'. Since u belongs to E, write u=ubar(w), where ubar(0)=0,ubar'(0)=1. Define the rational function

    Dbar(z)=-Cbar(z)^3 ubar(z)^2/z².                   (P8)

It is regular at0 with value-1. Substitution of(P8) into u²=-Dx²/C³ gives precisely the D identity in(P1).

Suppose Dbar had a pole at a finite beta!=0. Every preimage alpha under w is finite and nonzero. If its ramification index is j>=1, the local order of the right side in(P1) is

    j*ord_beta Dbar-3(j-1)<0,

because w/x is a unit there. This would be a pole of polynomial D. Hence Dbar has no finite pole and is polynomial. In particular no assumed polynomial composition for D is being used; the factors w/x and w' are essential.

The quotient inverse coordinate is genuinely ubar. Its chosen square root is Xbar=ubar Cbar/z, with Xbar(0)=1; (P8) gives Xbar²=-Dbar/Cbar. Differentiating(P1) yields

    T=1/(C u')=1/(Cbar(w) ubar'(w))=Tbar(w).

Under the local inverse branches, w(v(u)) is the inverse branch of ubar. Therefore the quotient endpoint Jacobian is exactly the original one. Its positive even coefficients vanish, and odd pure powers vanish by height parity. Thus pbar is all-pure over K.

## 4. Admissible converse and exact marker equality

**Theorem AdmissibleOddCubicPolynomialPullback.** Conversely, let w in K[x] be any polynomial with w(0)=0,w'(0)=1. Let Cbar,Dbar in K[z] be unit-center coefficients with rational selected coordinate ubar=z*sqrt(-Dbar/Cbar)/Cbar. Define C,D by(P1), and ASSUME both are polynomials. Then the selected original coordinate and Jacobian satisfy

    u=ubar(w), X=Xbar(w)*(w/x)/w', T=Tbar(w).

The normalized original phase is all-pure if and only if the quotient phase is all-pure. For EVERY polynomial marker Gbar in K[z], one has the exact angular identity

    Lambda_x(Gbar(w(x))*p(x,t)^m)
       =Lambda_z(Gbar(z)*pbar(z,t)^m)                  (P4a)

for all integers m>=0. Consequently every fixed x-only polynomial marker witness of unsafe behavior transfers from the quotient by polynomial composition. This is not an equivalence of safety for arbitrary original or arbitrary sphere markers.

Proof. The square-root identity follows directly from(P1); all factors have center value1, so it is the selected branch. It gives u=ubar(w), and differentiation gives T=Tbar(w). The inverse branches satisfy w(v(u))=vbar(u). The endpoint/residue formula for even m is

    Lambda_x(G(x)*p^m)=kappa_m [u^m]G(v(u))*T(v(u)).

For G=Gbar(w), its right side becomes exactly the corresponding quotient expression. Odd m vanish on both sides by t parity; m0 gives Gbar(0) on both sides. This proves(P4a), purity equivalence, and fixed-marker transfer. Polynomiality of C,D is essential: no arbitrary polynomial w is asserted admissible. Nor is a rationality theorem outside the selected rational-coordinate hypothesis being imported into the converse.

If the quotient already has full stable field K(z), the equality of rational primitive parts proved next shows the full original stable field is K(w). Thus the chosen w is exactly its canonical normalized descent polynomial. If the quotient itself descends further, the canonical descent is the composition with that further normalized polynomial.

## 5. Exact field equality and rational return downstairs

Decompose the quotient differential as dz/Cbar=dRbar+Omegabar, normalized Rbar(0)=0. Pulling back gives

    omega=d(Rbar(w))+w*Omegabar.

The second term is still a logarithmic differential: each dlog(z-alpha) pulls back to dlog(w-alpha). Compare this with omega=dR+Omega. Exact/log separation implies R-Rbar(w) is constant; its value at0 is0. Hence R=Rbar(w).

By construction

    E=K(u,R,T)
     =K(ubar(w),Rbar(w),Tbar(w))=K(w).

Injectivity of rational composition therefore gives(P2). The formal return is a Mobius involution of the quotient rational-function field K(z), fixing z0. This proves the promised rational-return reduction. It gives no automatic extension to the original x-line; when g>1 that extension is a separate question.

## 6. Degree and root-count identities

Let cbar=degCbar,ebar=degDbar. Taking exact leading degrees in(P1) gives

    c=g*cbar-(g-1)=g(cbar-1)+1,
    e=g*ebar+2(g-1)-3(g-1)=g(ebar-1)+1.

For c>=2 this forces cbar>=2. The rational primitive degrees satisfy

    c-H=degR=g*degRbar=g(cbar-Hbar),

and subtraction gives H=g(Hbar-1)+1. Thus g divides H-1 as well as c-1,e-1. The normalized rational Xbar implies cbar and ebar have the same parity, so g divides(c-e)/2 too.

The exact least degree is max(c,e+2). If cbar<=ebar then D dominates upstairs and downstairs, giving d=g(dbar-3)+3. Otherwise cbar>=ebar+2 (same parity), C dominates upstairs and downstairs, giving d=g(dbar-1)+1. No unsupported tie case is omitted.

## 7. Existing examples and explicit limits

For the retained hyperbolic odd-cubic family with parameter N, set U=1+lambda*x and

    w_N=(U^N-1)/(N*lambda), w_N'=U^(N-1).

The quotient is the N=1 family with parameter lambda_bar=N*lambda, since 1+lambda_bar*w_N=U^N. This applies also to the retained Q-deformation with odd ell: Q=(1+U^(2Nell))/(1+U^(2N)). Its only repeated C-root is U0 with multiplicityN+1, and its normalized rational primitive is R=(1-U^(-N))/(Nlambda). The full stable field is exactly K(U^N), so g=N. The existing prime13 witness is thus a degree2 polynomial descent from the degree8 quotient; prime37 is a degree2 descent from the degree20 quotient. These identities illustrate the theorem, not new witness constructions.

The main indispensable classical input is Milne Fields and Galois Theory v5.10, Theorem9.19 (printedp117), inspected in the retained author PDF. Rational automorphisms of K(z) are fractional-linear by rational-map degree multiplication. The local pullback order formula follows by differentiating a uniformizer power times a unit. The inverse-coordinate and exact/log inputs are the retained odd-cubic descent and mixed-center separation theorems. No universal classification of rational-return quotients, no arbitrary-pullback preservation claim, no outside novelty assertion, and no finite-prefix purity inference is made.

The local differential order was also checked against the actual characteristic-zero calculation in [Stacks Project Section53.12](https://stacks.math.columbia.edu/tag/0C1B), preceding Lemma53.12.4. The proof above supplies the needed uniformizer computation directly. [Pakovich, Theorem1.2](https://www.cs.bgu.ac.il/~pakovich/Publications/tri.pdf) describes polynomial moment solutions by polynomial compositions; its complex polynomial segment-integral hypotheses differ from the rational differential and formal inverse identities here, so it is context rather than an imported solution of this problem.

The [source ledger](odd_cubic_canonical_polynomial_descent_source_ledger.json) pins six retained inputs. The [checker](verify_odd_cubic_canonical_polynomial_descent.py) and [receipt](odd_cubic_canonical_polynomial_descent_verification.json) bind exact finite pullback and marker controls to these proof bytes.

```bash
python3 research/foundations/extended-geometry/verify_odd_cubic_canonical_polynomial_descent.py \
  --write-receipt research/foundations/extended-geometry/odd_cubic_canonical_polynomial_descent_verification.json
```
