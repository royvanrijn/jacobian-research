# A gap-free quarter-critical profile and a one-sixth coefficient-field obstruction

This proof gives a NECESSARY condition for an actual six-zero binary quadratic block at the specific constant radius3/4. It neither constructs such a block nor closes unrestricted quadratic SFC. No quarter-integer spacing or quantitative lower bound on positive coefficient valuations is assumed in the main theorem.

Let E be a characteristic-zero real-valued nonarchimedean coefficient field, with v(2)=1, and let

 q=alpha*X^2/2+beta*XY+gamma*Y^2/2+dX+eY

have integral natural coefficients. Write R_r=Lambda(q^r)/r!,R0=1 and suppose its native residue series is nonpolynomial. Let k belong to a real-valued extension, with v(k)=3/4. All polynomial and series congruences below are finite or formal; no analytic EGF evaluation is used.

The mathematical inputs are the installed native integer coefficient/residue proof [native integer coefficient and residue proof](QUADRATIC_CONSTANT_NONCRITICAL_DYADIC_WINDOW_TRANSFER.md), and the critical detector/six-window proof [critical detector and six-window proof](BINARY_QUADRATIC_CRITICAL_DYADIC_SIX_MOMENT_WINDOWS.md). The last ramification consequence also uses the existing necessary grid obstruction in [necessary grid obstruction](BINARY_QUADRATIC_TWO_POLE_CONSTANT_SHIFTS_STRONG_FACTORIAL_WINDOWS.md). The new native profile and invariant are proved here without importing the earlier quarter-gap theorem.

## 1. Statement and the field-residue necessity

**Theorem.** Suppose six consecutive positive factorial moments of f=q+k vanish. After swapping X,Y and dividing the WHOLE phase by an E-unit linear coefficient, one can normalize d=1. Put

 C=R1=alpha+beta+gamma+1+e, t=C-1,
 lambda=k^2/2, U=k^4/8.

Then C has residue one and

                         v(C-1)=1/6.                     (1)

In particular the value group of the NONCONSTANT coefficient field E contains1/6, regardless of the field initially containing k. If v(E*) is contained in(1/e0)Z for an integer e0, then6 divides e0. Combining this with the existing nonpolynomial six-zero obstruction4 divides e0 gives12 divides e0. For a finite Q2 coefficient field this applies to its MINIMAL ramification index; it is not a claim that a counterexample exists at index12.

Here is the residue/start reduction, recalled with its native degree argument. The installed critical theorem shows that an actual six-zero block at v(k)=3/4 requires

 B(T)=(1+mu*T^4)Rbar(T)

to be a polynomial, where mu is the residue of U. Its degree is at most5. The necessary eight-periodic detector is the coefficient prefix B0,...,B7, with B0=1. A six-zero circular block can start only at residue1 or2. At1 it would make B=1; then D=(1+mu*T^4)N, forcing degD4 and degN0, contradicting the native implication degD4 =>N2=dbar*ebar!=0. At2, degB<=1 and the constant case is already excluded. Native degrees force degD3,N=1 and D=(1+rho*T)^3, rho^4=mu. Factoring the two native quadratics, up to swap one obtains

 alphabar=ebar=0, gammabar=dbar=betabar=rho!=0,
                         n=2 modulo8.

Divide by the E-unit d. Its residue is rho, so the normalized data obey

 alphabar=ebar=0, betabar=gammabar=dbar=1,
 Cbar=1, Ubar=1, n=8H+2.                                 (2)

All subsequent natural coefficients and C still belong to E. This normalization changes neither zero moments nor the constant valuation. No lower bound such as1/4 is inferred for a positive E-valuation.

## 2. A full modulo-two prefix and a monic-factor argument without a gap

For this section work in the ring A=O/2O. Its characteristic is two, but it can have nilpotents. Set lambda=k^2/2 as above. It has valuation1/2, k has valuation3/4 and U is a unit; lambda^2,k^2,k*lambda vanish in A.

The exact expansion is M_m=sum_r(m)_r*k^(m-r)*R_r. Divide each of the six zero equations m=8H+s,s2,...,7 by its nonzero r=s coefficient. Every r>s costs at least one valuation unit more, by the critical digit-cost split. For r<=s the odd parts of the retained rational integer ratios are one in A. The six equations, successively simplified, give

 R2=lambda+kC, R3=lambda*C,
 R4=U, R5=UC, R6=UR2, R7=UR3 in A.                        (3)

For completeness the unsimplified rows are

 0=R2+lambda+kC,
 0=R3+lambda*C+kR2,
 0=R4+U+lambda*R2+kR3,
 0=R5+UC+kU+lambda*R3+kR4,
 0=R6+UR2+lambda*U+kUC+lambda*R4+kR5,
 0=R7+UR3+lambda*UC+kUR2+lambda*R5+kR6.

The first rows make v(R2)=v(R3)=1/2. Hence lambda*R2,lambda*R3,kR2,kR3 vanish in A; the remaining U-pairs cancel. This proves(3) without any coefficient gap.

The native integer discarded-factor-two argument is valid over A, not only over the residue field. Thus R=N/D with

 D=(1+alpha*T+T^2)(1+gamma*T+e^2*T^2),
 N=1+N1*T+eT^2, N1=1+e+beta.

Let

 P(T)=1+CT+(lambda+kC)T^2+lambda*C*T^3.

The prefix(3) and degree bounds degD<=4,degP<=3,degN<=2 prove the EXACT polynomial identity over A

                         D P=N(1+UT^4).                  (4)

Indeed the two sides have degree at most7 and their difference has zero coefficients through7.

Set Z=1+CT,h=k/C,l=lambda/C^2. Then h^2=l^2=hl=0 in A and

 P=h+(1+l)Z+hZ^2+lZ^3=(Z+h)K,
                         K=1+l+hZ+lZ^2.                 (5)

In particular Z=h is an ACTUAL polynomial root of P in A. It corresponds to T0=(1+h)/C, an integral unit. The value N(T0) is a unit because its residue is one. Evaluating(4) there, and using T0^4=C^(-4), proves

                              U=C^4 in A.                (6)

This is the replacement for the earlier gap-dependent evaluation at T=1/C. No estimate on v(C-1) was used.

Now(4) becomes DP=N Z^4. Since h^2=0 one has Z^4=(Z+h)Z^2(Z+h). A monic polynomial is a non-zero-divisor in a polynomial ring over ANY coefficient ring: its highest nonzero coefficient cannot disappear on multiplication. Cancel the monic Z+h to obtain

                              D K=N Z^2(Z+h).             (7)

Scale D by C^4 after T=(1+Z)/C and write F=A_x A_y, where

 A_x=C^2+alpha*C+1+alpha*C*Z+Z^2,
 A_y=C^2+gamma*C+e^2+gamma*C*Z+e^2Z^2,
 F=sum_(i=0)^4 F_i Z^i.

Equation(7) and the unit constant coefficient1+l of K show F0=F1=0 in A. Modulo the maximal ideal, F=Z^3, so F3 is a unit. Coefficient three in FP=C^4 N Z^4 gives F2=hF3. Coefficient six gives

                         C^2 e=lF3+h e^2.

The first term has valuation1/2, whereas the second has valuation greater than3/4, because e initially has positive valuation. Thus v(e)=1/2. In particular e^2=F4=0 in A and F3=gamma*C. Coefficients five, six and four now give respectively

 C^3N1=hF3, C^2e=lF3,
 C^4+C^3N1+C^2e=(1+l)F3.

Eliminating them yields (1+h)F3=C^4. Therefore the GAP-FREE full profile is

 gamma=C^3+kC^2, e=lambda, beta=1+lambda+k,
                         alpha=(C+k)(C-1)^2 in A.        (8)

The last equality uses the exact native first moment C. These identities require no ordinary specialized pivot except the displayed units. They reduce to the earlier inverse-C formulas if a separate gap makes C^4=1 in A; that equality has not been assumed here.

## 3. The necessary E6,E7 modulo four does not require U=1 modulo two

Write a_j=k^j/j! and E_s=sum_(j<=s) a_j R_(s-j) for s<=7. Only these finite coefficients are used. Equation(4) also proves

 R(T)=P(T)/(1+UT^4) in A,
                         R_(r+4)=U R_r for every r>=0.    (9)

In addition E0=1,E1=C+k,E2=E3=0 in A. None of these statements requires U=1 in A; it is only a unit with residue one.

Here is the complete arbitrary-H next-precision argument. The normalized coefficient of R_r in M_(8H+s), after dividing by its r=s coefficient, is

                  c_(s,r)=(8H)!/(8H+s-r)! * k^(s-r).

For r=s-j,j0,...,7, its ratio to a_j is binom(8H+j,j)^(-1): one modulo4 for j<4 and (-1)^H modulo4 for j>=4. The factor i4 is the only nontrivial odd ratio after removing powers of two. For r>s+8 the high digit cost is at least two and all terms vanish modulo4. For the first tail r=s+t,t1,...,8, even H again gives cost at least two. For odd H, the integer (8H)_t*(8-t)!/128 is odd, and k^8=64U^2, so

                         c_(s,s+t)=2U^(-2)a_(8-t) mod4.

The whole first tail cancels modulo4: pair w and w+4, using a_(w+4)=U a_w and R_(r+4)=U R_r in A. Only the head sign correction -2U E_(s-4) remains. For s6,7 it vanishes modulo4 because E2=E3=0 in A. Hence the hypothetical original block implies

                              E6=E7=0 modulo4.            (10)

This locally rederives the needed part of the previous parity proof under its actual weaker hypotheses; the older source need not be edited.

## 4. Removing unknown 2O coefficient errors

Let W=E6+E7, holding k fixed. It is an integral polynomial in the four natural parameters alpha,beta,gamma,e because every R_r has integer coefficients and all a_j through7 are integral. Its residue at (2) is that of R6+R7+R2+R3: only a0 and a4 have nonzero residues and both residues are one.

At the cube point the derivatives of the universal native residue R with respect to alpha,e,beta,gamma are respectively

 T/(1+T)^5, T/(1+T)^2, T/(1+T)^3, T/(1+T)^4.

Their coefficients at r2,3,6,7 are the following explicit rows:

| Parameter | r2 | r3 | r6 | r7 | Sum in characteristic two |
|---|---:|---:|---:|---:|---:|
| alpha | 1 | 1 | 0 | 0 | 0 |
| e | 0 | 1 | 0 | 1 | 0 |
| beta | 1 | 0 | 1 | 0 | 0 |
| gamma | 0 | 0 | 0 | 0 | 0 |

Thus every first partial derivative of W has positive valuation at any tuple with the cube residue. Changing any parameter by an element of2O changes W only by an element of valuation STRICTLY greater than one: the linear terms have this property and higher polynomial Taylor terms contain at least two such increments. This is integral polynomial/Hasse Taylor expansion, so no factorial division is introduced.

Set t=C-1. By(8) one may therefore replace the natural parameters, for the purpose of computing W modulo the ideal I={x:v(x)>1}, by the exact representatives

 alpha_0=(C+k)t^2, beta_0=1+lambda+k,
 gamma_0=C^2(C+k), e_0=lambda, d=1.                      (11)

Only the SUM W is asserted invariant; E6,E7 separately need not be invariant under these unknown2O corrections. Also replace U by C^4 in a4,a5,a6,a7, using(6). In a4 the error is multiplied by R2+R3 of positive valuation; all other errors already contain k or lambda. Hence these changes likewise affect W only in I.

## 5. A finite native polynomial invariant modulo STRICT valuation greater than one

In the quotient O/I, all the following relations are valid:

 4=0, 2t=2lambda=2k=0,
                         lambda^2=2, lambda*k=k^2=0.      (12)

Indeed t has positive valuation, lambda has valuation1/2 and k has valuation3/4. Moreover lambda^2=2U and U has residue one, so lambda^2-2 has valuation greater than one. In particular a finite calculation may be made in the universal quotient

 (Z/4Z)[t,k,lambda]/(2t,2k,2lambda,lambda^2-2,k*lambda,k^2).

No bound on v(t) beyond positivity is used in this calculation; no infinite series in t is evaluated. Put C=1+t, use(11), and for each r0,...,7 use the integer native formula

 R_r=sum_(i+j+h+l+m=r)
 [(2i+j+l)!*(2h+j+m)! /
   (2^(i+h)*i!*j!*h!*l!*m!)] alpha^i beta^j gamma^h e^m.   (13)

The finite exponential coefficients become

 a0=1, a1=k, a2=lambda, a3=0,
 a4=3+t^4, a5=k(1+t^4), a6=lambda(1+t^4), a7=0.

Here 1/3=3,1/15=3,1/45=1 modulo4, and the terms2 times a positive variable vanish. Direct multiplication and(12) give these TWO complete rows for the representative small moments:

 E6_tilde=2+t^6+t^10+t^12+t^14+t^18,
 E7_tilde=t^12+t^13+t^14+t^15+k*t^12+k*t^14.

These are calculations for the representatives(11), not individual invariance claims about the original E6,E7. Their sum, however, equals the original W modulo I by Section4. Consequently

 W=2+t^6+t^10+t^13+t^15+t^18+k*t^12+k*t^14 modulo I.     (14)

Equations(12),(13), the stated coefficients a0,...,a7 and the two rows are a finite universal polynomial certificate for(14). They can be checked using only integer factorial coefficients and ordinary multiplication in the displayed quotient; no experimental range in a moment index is extrapolated. The number of native monomials in R0,...,R7 before combination is792.

If v(t) is infinite then the right side of(14) has valuation exactly one, contradicting(10). Otherwise set g=v(t)>0. Every nonconstant term of(14) except t^6 has valuation strictly greater than6g. If6g<1, t^6 is uniquely least, contradicting W in4O plus an error of valuation greater than one. If6g>1, the constant2 is uniquely least, with the same contradiction. Therefore6g=1, proving(1).

## 6. Exact consequence and limits

The normalized first moment C belongs to E even when k does not. Thus v(C-1)=1/6 is an obstruction in the NONCONSTANT coefficient field itself. If a discrete grid bound e0 is available, 6 divides e0. The already proved nonpolynomial six-zero obstruction4 divides e0 then gives12 divides e0. This assertion applies to every admissible bound and in particular to the least ramification index of a finite Q2 coefficient field. It does not say such a phase exists at ramification12.

The profile also then forces v(alpha)=1/3, v(e)=1/2, v(beta-1)=1/2 and v(gamma-1)=1/6: in(8) the respective least terms are unique and below valuation one. In particular these four coefficients are nonzero. These are additional NECESSARY local data, not constructions or smaller actual-support SFC claims.

Only the radius v(k)=3/4 and nonpolynomial native residue are considered. Polynomial residues, other critical radii on more general grids, and unrestricted quadratic SFC/SFC(2)/SFC(3) remain separate. The earlier quarter-gap exclusion and all-constant quarter-grid theorem stay unchanged.
