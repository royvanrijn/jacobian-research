# A lifted native quarter-critical profile on a quarter-integer coefficient grid

This note gives necessary congruences for a hypothetical six-zero binary-quadratic moment block. They are compatible congruences, not a contradiction or a construction of such a block. It refines the one-axis critical residue profile and keeps the actual native generating identity, not an arbitrary recurrence state.

Let E be a characteristic-zero nonarchimedean valued field with v(2)=1 and v(E*) contained in (1/e0)Z, where v_2(e0)=2. Work in a valued extension containing k, with v(k)=3/4. Normalize the WHOLE phase by its unit linear-X coefficient, so

 q=alpha X^2/2+beta XY+gamma Y^2/2+X+eY,
 R_j=Lambda(q^j)/j!, C=R_1,
 R(T)=sum R_j T^j,
 f=q+k.

Suppose the critical residue profile is

 alphabar=ebar=0, gammabar=betabar=1,
 n=2 mod8, M_n=M_(n+1)=...=M_(n+5)=0,
 M_m=Lambda(f^m).

Thus Cbar=1, and every positive valuation of a nonzero E-element is at least1/e0, which need NOT be1/4 if e0 has an odd factor. For this note impose the stronger explicit quarter-gap hypothesis

 min{v(alpha),v(e),v(beta-1),v(gamma-1)}>=1/4,              (1)

with v(0)=infinity. In particular v(C-1)>=1/4. Hypothesis(1) holds when the coefficient value group is contained in (1/4)Z. The conclusion below is NOT claimed for a grid with denominator4 times a nontrivial odd integer without this extra gap.

Put

 lambda=k^2/2, kappa=k, U=k^4/8.

All congruences below are in the valuation ring modulo the ideal2O, not merely modulo its maximal ideal. Hence lambda has valuation1/2 and kappa valuation3/4 even though both have zero residue.

**Proposition.** Under these hypotheses,

 U=1 mod2O,
 e=lambda mod2O,
 beta=1+lambda+kappa mod2O,
 gamma=C^(-1)+kappa/C^2 mod2O,
 alpha=(C-1)^2/C mod2O.                                  (2)

Consequently v(e)=1/2 and v(alpha)>=1/2. These are necessary native conditions only. They do not prove the six-zero block exists, and do not exclude the remaining quarter-grid case.

## 1. Six actual moment equations give a finite generating identity

Legendre's formula gives the exact coefficient cost, after removing k^m,

 v((m)_r)-3r/4=r/4+s_2(m-r)-s_2(m).

For m=8H+s with2<=s<=7, divide the moment equation by its r=s coefficient. This coefficient has the minimum base cost; when s>=4 it ties r=s-4. Every r>s has cost at least ONE more than that minimum. Indeed j=m-r=8u+w has u<=H-1; the high digit cost s_2(u)-2u increases by at least1 compared with u=H, while the low cost s_2(w)-w/4 is nonnegative. For H=0 those tails are empty.

For r<=s, the odd parts of every retained integer falling-factorial ratio reduce to1 modulo2. Therefore the six equations, successively simplified using integrality and preceding equations, give

 R_2=lambda+kappa C,
 R_3=lambda C,
 R_4=U,
 R_5=U C,
 R_6=U R_2,
 R_7=U R_3                              modulo2O.         (3)

Here the first equation also proves v(R_2)=1/2 and the second proves v(R_3)=1/2. For clarity, before simplification the retained equations from the r=s divisions are

 0=R_2+lambda+kappa C,
 0=R_3+lambda C+kappa R_2,
 0=R_4+U+lambda R_2+kappa R_3,
 0=R_5+U C+kappa U+lambda R_3+kappa R_4,
 0=R_6+U R_2+lambda U+kappa U C+lambda R_4+kappa R_5,
 0=R_7+U R_3+lambda U C+kappa U R_2+lambda R_5+kappa R_6.

All other coefficients already have valuation at least1. Products lambda R_2,lambda R_3 and kappa R_2,kappa R_3 have valuation at least1; the displayed U-pairs cancel using preceding equations. This proves(3) without assuming only the six small initial moment orders. It is valid for EVERY n=8H+2.

The native coefficient identity is valid coefficientwise modulo2O:

 R(T)=N(T)/D(T),
 N=1+(1+e+beta)T+eT^2,
 D=(1+alpha T+T^2)(1+gamma T+e^2T^2).                    (4)

It follows from integer coefficient identities before reduction, so remains valid in the possibly nonreduced ring O/2O. Merely knowing the identity over the residue FIELD would be insufficient here. A complete native derivation is Sections1 and5 of QUADRATIC_CONSTANT_NONCRITICAL_DYADIC_WINDOW_TRANSFER.md; its factorial coefficient brackets are integers and every discarded term contains the rational factor2.

Let

 P(T)=1+C T+(lambda+kappa C)T^2+lambda C T^3.

Equations(3) say that (1+U T^4)R(T) agrees with P(T) through degree7. Since degD<=4,degP<=3 and degN<=2, this implies the EXACT polynomial identity

                   D(T)P(T)=N(T)(1+U T^4) mod2O.          (5)

Both sides have degree at most7. Thus the finite prefix determines every coefficient of their difference.

## 2. Evaluation determines the critical unit

At T=1/C one has P(1/C)=kappa/C. Each quadratic factor of D has zero residue there, and by(1) its valuation is at least1/4. Hence D(1/C)P(1/C) has valuation at least1/4+1/4+3/4>1. On the other hand N(1/C) is a unit. Equation(5) gives

                        U=C^4 mod2O.

Since Cbar=1 and v(C-1)>=1/4, characteristic2 in O/2O gives C^4=1+(C-1)^4=1 mod2O. Thus U=1 mod2O. This step is where the explicit quarter-gap hypothesis is used. It must not be silently inferred from a denominator4-times-odd grid.

## 3. Coefficient comparison at the surviving pole

Put Z=1+C T and set l=lambda/C^2,h=kappa/C, so v(l)=1/2,v(h)=3/4. In this coordinate

 P=h+(1+l)Z+hZ^2+lZ^3.

Multiply the two factors of D by C^2 and write

 A(Z)=A_0+A_1Z+Z^2,
 A_0=C^2+alpha C+1, A_1=alpha C,
 B(Z)=B_0+B_1Z+B_2Z^2,
 B_0=C^2+gamma C+e^2, B_1=gamma C, B_2=e^2.

Let F=A B=sum_(j=0)^4 F_j Z^j. Then F_0 has valuation at least1/2,F_1 at least1/4,F_3 is a unit,F_4=e^2. Since U=C^4 modulo2O, equation(5) becomes

 F(Z)P(Z)=C^4 N((1+Z)/C) Z^4 mod2O.                      (6)

The coefficient of Z^3 gives

                       (1+l)F_2=hF_3 mod2O,

because lF_0 and hF_1 have valuation at least1. Thus v(F_2)=3/4. The coefficients of Z^2 and Z then give F_1=F_0=0 mod2O.

The coefficient of Z^6 in(6) is

                      e C^2=lF_3+h e^2 mod2O.

The final term has valuation at least3/4+1/2>1 by(1), whereas lF_3 has valuation1/2. Hence v(e)=1/2. In particular e^2=0 in O/2O, and F_3=gamma C there.

Put N_1=1+e+beta. The coefficients of Z^5,Z^6 and Z^4 now respectively give

 C^3 N_1=hF_3,
 C^2 e=lF_3,
 C^4+C^3 N_1+C^2e=(1+l)F_3                    mod2O.

Eliminating the first two equations yields (1+h)F_3=C^4. Since h^2=0 modulo2O,

 F_3=C^4(1+h), gamma=C^3(1+h),
 e=lambda(1+h)=lambda,
 N_1=kappa                                  mod2O.

Using C^4=1 and h=kappa/C gives gamma=C^(-1)+kappa/C^2 and beta=1+lambda+kappa. Finally the exact first moment C=alpha+beta+gamma+1+e gives

 alpha=C+C^(-1)+kappa(1+C^(-2)).

The last term has valuation at least3/4+2v(C-1)>=5/4 and vanishes modulo2O. Therefore alpha=(C-1)^2/C mod2O, proving(2). In particular v(alpha)>=1/2.

## 4. What remains

The conditions(2) are consistent as congruences and do not force a contradiction. They show that in this quarter-gap sector the linear-Y edge coefficient has valuation exactly1/2, rather than the7/8 required by the singular H-edge pivot e^2=2gamma k. They also constrain the small diagonal coefficient and both unit cross/diagonal perturbations. A native argument beyond this first full2-adic coefficient layer is still needed.

No HGM rank or arbitrary-jet argument is used. No numerical prefix, canonical/DAG change or external/formal assurance is a premise. The coefficient-gap hypothesis(1) is explicit and essential to the presently written proof.