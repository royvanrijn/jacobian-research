# Excluding the quarter-gap critical six-zero block by a native weighted coefficient identity

This proof closes the explicit quarter-gap branch. It does not close arbitrary binary quadratic SFC, arbitrary coefficient ramification, or the higher critical radii. All valuations below are real-valued extensions of v2, with v(2)=1 and v(0)=infinity. Congruence modulo valuation at least a means the difference has valuation >=a.

The proof uses these already written inputs, with their precise scopes:

* The critical residue profile in `BINARY_QUADRATIC_QUARTER_CRITICAL_RESIDUE_PROFILE.md`: a nonpolynomial-residue six-zero block at v(k)=3/4 must start at n=2 modulo8 and has the one-axis residue profile after swapping variables and whole-phase unit scaling.
* The explicit-gap lift `BINARY_QUADRATIC_QUARTER_GRID_LIFTED_NATIVE_PROFILE.md`: under the additional gap below, it proves the full O/2O congruences in(2), not merely residue-field identities.
* The uniform start reduction in `BINARY_QUADRATIC_NATIVE_MOD_FOUR_AND_QUARTER_PARITY_REDUCTION.md`, Section3: every actual six-zero block n=8H+2 implies E6=E7=0 modulo4, independently of the parity of H. The possible parity offsets are confined to E4,E5.

No computational experiment or Groebner output is a premise of the proof below. The universal finite coefficient calculation is displayed explicitly in Section2.

## 1. Exact local statement and consequences of a hypothetical block

Put

 q=alpha*X^2/2+beta*XY+gamma*Y^2/2+X+eY,
 f=q+k, R_j=Lambda(q^j)/j!, C=R1, lambda=k^2/2, U=k^4/8.

Assume all five natural nonconstant coefficients are integral, v(k)=3/4, the native residue generating series Rbar is nonpolynomial, and

 min{v(alpha),v(e),v(beta-1),v(gamma-1)}>=1/4.             (1)

**Theorem.** Under these hypotheses there are no six consecutive zero factorial moments of f at any positive start.

Suppose otherwise. The critical residue profile forces n=8H+2. The lifted profile then gives, modulo2O,

 U=1, e=lambda, beta=1+lambda+k,
 gamma=C^(-1)+k/C^2, alpha=(C-1)^2/C,
 R2=lambda+kC, R3=lambda*C,
 R4=U, R5=UC, R6=UR2, R7=UR3.                            (2)

Here C is a unit and v(C-1)>=1/4 by(1) and the exact native first moment. In particular,

 v(alpha)>=1/2, v(e)=1/2, v(beta-1)>=1/2,
 v(gamma-1)>=1/4, v(R2)=v(R3)=1/2.                       (3)

The last equality follows because lambda has valuation1/2, kC has valuation3/4, and C is a unit. Set z=gamma-1. The same congruences imply

                       z^2+alpha*gamma in 2O.             (4)

Indeed in O/2O one has z=(C-1)/C+k/C^2, so z^2=(C-1)^2/C^2 because v(k^2)=3/2. Also alpha*gamma=(C-1)^2/C^2: the additional term k(C-1)^2/C^3 has valuation at least3/4+1/2>1. This proves(4). No specialized coefficient or matrix minor is divided by.

Write E_s=Lambda((q+k)^s)/s!, for s<=7. The uniform start reduction gives

                         E6=E7=0 modulo4.                 (5)

We will show their sum has valuation exactly one, contradicting(5).

## 2. A universal finite weighted native coefficient calculation

Let

                         Hnat=R6+R7-R2-R3.

As a polynomial in the five natural coefficients, each R_n has INTEGER coefficients. With d=1, its exact formula is

 R_n=sum_(i+j+h+l+m=n)
   [(2i+j+l)!*(2h+j+m)! /
       (2^(i+h)*i!*j!*h!*l!*m!)] alpha^i beta^j gamma^h e^m. (6)

The bracket is an integer, for example by the two vertex partitions with the cross multiplicity factorial. Set y=beta-1 and z=gamma-1. Give alpha,y,e weight2 and z weight1. Expansion of(6) for n2,3,6,7 gives the following COMPLETE list of monomials of weight at most four in Hnat. Every monomial not in this table has weight at least five.

| Monomial | Integer coefficient |
|---|---:|
| 1 | 865998 |
| z | 3841512 |
| z^2 | 8086080 |
| z^3 | 10134825 |
| z^4 | 8004360 |
| alpha | 537636 |
| alpha*z | 1236227 |
| alpha*z^2 | 1399887 |
| y | 1952758 |
| y*z | 6076200 |
| y*z^2 | 8858910 |
| e | 656716 |
| e*z | 2601096 |
| e*z^2 | 4695090 |
| e^2 | 232152 |
| y*e | 1167400 |
| y^2 | 2472454 |
| alpha*e | 288038 |
| alpha*y | 1703160 |
| alpha^2 | 465522 |

The finite calculation can be performed directly from(6): substitute beta=1+y,gamma=1+z and use the binomial theorem. Equivalently its odd coefficients follow by expanding the universal formal residue

 [1+(y+e)T+eT^2] /
     [(1+alpha*T+T^2)(1+(1+z)T+e^2*T^2)]

through T7 and weighted degree four. The constant is independently R2=12,R3=68,R6=58766,R7=807312 at alpha=e=0,beta=gamma=1, giving 865998. Thus no characteristic-positive coefficient lifting is assumed.

By(3), a weight-w monomial has valuation at least w/4. All higher-weight terms therefore have valuation at least5/4. Every even-coefficient nonconstant term in the table also has valuation at least1+1/4=5/4. The constant is2 modulo4; the three odd coefficients can be replaced by one, since their differences from one are even and their monomials have positive valuation. Consequently

              Hnat=2+z^3+alpha*z+alpha*z^2
                              modulo valuation >=5/4.    (7)

But the nonconstant part factors EXACTLY as

 z^3+alpha*z+alpha*z^2=z*(z^2+alpha*gamma).

Equations(3),(4) give this factor valuation at least1/4+1=5/4. Hence

                  Hnat=2 modulo valuation >=5/4,
                         v(Hnat)=1.                      (8)

## 3. Comparing the actual E6,E7 with the native difference

The first eight exponential coefficients are

 a0=1, a1=k, a2=lambda, a3=k*lambda/3,
 a4=U/3, a5=k*U/15, a6=lambda*U/45,
 a7=k*lambda*U/315.

Their valuations are respectively0,3/4,1/2,5/4,0,3/4,1/2,5/4. This statement concerns only coefficients through seven; no integral infinite exponential is evaluated.

Expand E6,E7 as sum_(j<=s) a_j R_(s-j) and subtract Hnat. Terms with a3 or a7 have valuation at least5/4. For the remaining terms use(2), replacing R4 by U,R5 by UC,R6 by UR2; the replacement errors, multiplied by k or lambda, have valuation at least3/2. The result is

 E6+E7-Hnat=
   kU(1+1/15)(C+R2)
  +lambda*U(1+1/45)(1+C)
  +(U/3+1)(R2+R3)
                         modulo valuation >=5/4.          (9)

Every displayed term also has valuation at least5/4. The first scalar contains16/15, the second46/45 has valuation one, and U/3+1 lies in2O since U=1 modulo2O. The remaining factors are integral and R2,R3 have valuation1/2. Therefore

                   E6+E7=Hnat modulo valuation >=5/4.     (10)

By(8), the left side has valuation exactly one. This contradicts(5), which puts both summands in4O. The theorem follows for both parities of H, at every positive start allowed by the critical residue necessity.

## 4. Coefficient-grid corollary and exact limit

If the natural nonconstant coefficients lie in a valued field E with v(E*) contained in (1/4)Z, normalize the whole phase by the unit linear coefficient selected by the quarter-critical residue profile. For a hypothetical six-zero block with nonpolynomial Rbar at v(k)=3/4, that profile makes alpha,e,beta-1,gamma-1 have strictly positive E-valuations, hence at least1/4. Thus(1) holds and the theorem excludes the block, even if the constant initially lies in a larger valued extension. No assumption that the constant belongs to E is needed for this corollary.

This is a result at the SPECIFIC constant radius3/4 and with nonpolynomial native residue. It does not exclude all larger critical radii for coefficient grid(1/4)Z, and it does not infer the gap(1) merely from a grid with denominator four times a nontrivial odd integer. Polynomial residue series and more finely ramified one-pole cancellations remain outside this proof. The unrestricted all-start binary-quadratic, SFC(2), and SFC(3) goals are not claimed closed.

