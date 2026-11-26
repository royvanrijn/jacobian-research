# Quadratic Strong Factorial windows from matchings

Fix m>=1 and a characteristic-zero field K. Let
Lambda_m(X_1^a1...X_m^am)=a1!...am!, and consider a nonzero polynomial

 f=sum_i A_i X_i^2+sum_(i<j) B_ij X_i X_j+sum_i d_i X_i,

with zero constant coefficient. For every i<j, fix a rational number c_ij
with odd denominator, allowing c_ij=0, and assume

                    B_ij^2=4 c_ij A_i A_j for every i<j.  (1)

Equivalently v_2(c_ij)>=0, with v_2(0)=infinity.
Put k=#{i:A_i!=0}, ell=#{i:d_i!=0}.

**Theorem 1.** For every n>=1, at least one factorial moment among

 Lambda_m(f^n),...,Lambda_m(f^(n+k+ell-1))

is nonzero. Consequently f satisfies its actual-support Strong Factorial
window. The theorem includes every rank-one quadratic part (c_ij=1),
all purely diagonal quadratic parts (c_ij=0), and all linear coefficient
patterns, in every dimension.

The actual number of monomials is simply

                  s=k+#{ij:B_ij!=0}+ell >= k+ell.

If k=0 then (1) makes every cross coefficient zero; f!=0 ensures ell>=1.
In every case k+ell>=1. The claimed Strong Factorial window follows from
the stronger displayed window, by taking its initial k+ell indices.

For m=2 and c_12=1, (1) is the singular-Hessian condition, and
k+ell<=4 gives a uniform four-window. Even with all c_ij=1,
condition(1) need NOT force Hessian rank1 in dimensions three and larger.
For example the symmetric matrix with diagonal1 and all off-diagonal -1
has all principal two-by-two minors zero but determinant -4 in dimension3.
Thus the theorem also covers these sign patterns, with arbitrary linear
coefficients satisfying no additional relation.

## 1. Natural coefficients and exact integrality

Work first in a finite valued extension L of Q_2, with v(2)=1 and
valuation ring O. Write

 alpha_i=2A_i, beta_ij=B_ij,
 P_n=Lambda_m(f^n)/n!, E(T)=sum_(n>=0) P_n T^n.

Assume temporarily all natural coefficients alpha_i,beta_ij,d_i are integral.
No principal-minor relation is needed for this section.

In a monomial of the coefficient expansion of P_n, let p_i count the
occurrences of alpha_i X_i^2/2, l_i those of d_i X_i, and j_e those of
a cross monomial beta_e X_i X_j. Their total is n. The exponent of X_i is

                    N_i=2p_i+l_i+sum_(e incident to i) j_e.

The exact coefficient in P_n is

 [product_i N_i!]/[product_i(2^p_i p_i! l_i!) product_e j_e!]
       =[product_i C_i] [product_e j_e!],                (2)

where

 C_i=N_i!/[2^p_i p_i! l_i! product_(e incident to i) j_e!].

Each C_i is an integer: it is (2p_i-1)!! times the multinomial coefficient
partitioning N_i into parts 2p_i,l_i and the incident j_e.
The convention (-1)!!=1 applies. Equation(2) therefore proves P_n is
an integer polynomial in all natural coefficients.

Reduce into the residue field, denoting reduction by bars.
If some cross count j_e>=2, the corresponding external factor j_e!
in (2) is even, so the whole coefficient vanishes.

Suppose all cross counts are0or1. At a vertex i let r_i be its number
of incident chosen edges. Then

 C_i=(2p_i-1)!! binom(N_i,2p_i) (l_i+r_i)!/l_i!.

If r_i>=2, the last product contains two consecutive integers and is even.
Thus any surviving set of cross edges has no repeated edge and no vertex
of degree at least two: it is a MATCHING. This is a coefficient argument;
no analytic Gaussian or Wick formula is an imported premise.

## 2. The universal matching generating function

For a single variable let q(X)=alpha X^2/2+dX, and define the formal series

 Q(T)=sum_(r>=0) Lambda_1(q^r)T^r/r!,
 I(T)=sum_(r>=0) Lambda_1(Xq^r)T^r/r!.

Their coefficients are integral by Section1. Literal coefficient comparison
gives, for q_r=Lambda_1(q^r)/r!,

 q_0=1,q_1=alpha+d,
 q_r=alpha(2r-1)q_(r-1)+d^2 q_(r-2), r>=2.

Consequently, over characteristic two, with D=1+alphabar T+dbar^2 T^2,

 Qbar=(1+dbar T)/D, Ibar=1/D.                            (3)

For completeness the second identity is also polynomially native.
The coefficient of alpha^i d^(r-i) in the r-th coefficient of I is
r+i+1 times that in q_r. Hence, as a universal parameter-series identity
before evaluating the residue parameters,

 Ibar=Qbar+T partial_T Qbar+alphabar partial_alphabar Qbar.

Substitute the first expression in(3) and multiply by D^2; its numerator
is exactly D, proving the second identity. No numerical parameter is
being differentiated.

By Section1 only matching cross sets survive in the full normalized series.
An unmatched vertex contributes Qbar_i, a matched vertex contributes Ibar_i,
and each matched edge contributes betabar_ij T. Therefore, putting
D_i=1+alphabar_i T+dbar_i^2 T^2,

 Ebar(T)=
 [sum_(M a matching) T^|M| product_(ij in M) betabar_ij
                         product_(i unmatched by M)(1+dbar_i T)]
 / product_i D_i.                                       (4)

The sum is over matchings of the complete graph on the m indexed variables;
zero edge coefficients make the relevant summand zero. Every coefficient
arises from a finite expansion of exp(Tf), so (4) is a purely formal identity.
All denominators have constant coefficient1. Formula(4) holds for ANY
quadratic with integral natural coefficients and zero constant.

## 3. Integral fixed weights force a proper rational series

Now impose beta_ij^2=c_ij alpha_i alpha_j, with all c_ij integral in O,
assume all natural coefficients integral, and assume at least one is a unit.
The rational odd-denominator condition ensures c_ij is integral in every
embedding used below. No square root or sign consistency is assumed.
Let U={i:alphabar_i!=0}, q=|U|. Put
ell_in=#{i in U:dbar_i!=0}, ell_out=#{i not in U:dbar_i!=0}.

If i is outside U, then betabar_ij^2=cbar_ij alphabar_i alphabar_j=0,
so every incident betabar_ij is zero. Such a vertex is unmatched in every
nonzero term of(4). Its factor simplifies EXACTLY as

 (1+dbar_i T)/(1+dbar_i^2 T^2)=1/(1+dbar_i T).            (5)

This is legitimate as a formal identity with unit constant denominators.

After these cancellations (4) has numerator

 N_U=sum_(M a matching on U) T^|M| product_(ij in M) betabar_ij
                             product_(i in U unmatched)(1+dbar_i T),

and denominator

 D_U=[product_(i in U) D_i] [product_(i not in U)(1+dbar_i T)].

Every i in U contributes degree2 when dbar_i!=0 and degree1 otherwise.
The outside factors contribute degree1 precisely at the ell_out nonzero
linear residues. Hence

                  deg(D_U)=q+ell_in+ell_out.              (6)

If q>=1, every numerator term has degree at most |M|+ell_in,
which is at most floor(q/2)+ell_in. Thus

 deg(N_U)<=floor(q/2)+ell_in < q+ell_in+ell_out=deg(D_U).   (7)

The numerator has constant coefficient1, so it is not zero.
If q=0, every edge residue is zero, N_U=1 and
D_U=product_i(1+dbar_i T). At least one dbar_i is nonzero:
otherwise the fixed-weight relations(1) would make all natural
coefficients nonunits, contrary to normalization.
Thus in this case too deg(N_U)=0<deg(D_U).

We have proved Ebar=N_U/D_U with D_U(0)=1 and positive degree r, with
nonzero leading coefficient and deg(N_U)<r.
No assumption about coprimeness of numerator and denominator is needed.

## 4. Exact nonvanishing window and phase normalization

A proper rational series p(T), with p_0=1 and denominator
1+a_1T+...+a_r T^r, a_r!=0, satisfies

 p_j+a_1p_(j-1)+...+a_r p_(j-r)=0, j>=r.

If r consecutive positive p_j values vanished, this relation propagates
backward until it forces p_0=0 at index r. This contradicts p_0=1.
Thus r consecutive zero coefficients are impossible.

In the situation of Section3,

 r=q+ell_in+ell_out <= k+ell,

since every unit diagonal or linear coefficient was already nonzero.
Consequently k+ell consecutive actual moments cannot all vanish:
their normalized residues would contain a forbidden initial r-subblock.

For arbitrary admissible coefficients over L, choose a nonzero natural
coefficient c of minimum valuation and divide the ENTIRE phase by c.
Every normalized natural coefficient is integral and at least one is a
unit. Relations(1) survive, as each is homogeneous of degree two.
The sets counted by k and ell do not change. The j-th factorial moment
is multiplied by c^(-j), so vanishing is unchanged. This proves Theorem1
over every finite extension of Q_2, for every support pattern.

## 5. Every characteristic-zero coefficient field

The finitely many coefficients generate a finitely generated field F over Q.
Choose a finite transcendence basis and map it to algebraically independent
elements of Q_2, then extend the embedding through the finite algebraic
extension into Q_2bar. Such a finite independent set exists because a
countable subfield has only countably many algebraic elements in Q_2,
whereas Q_2 is uncountable.
The finitely many images lie in one finite extension of Q_2.
All polynomial relations and nonzero coefficients are preserved.
A hypothetical counterexample over the original field would therefore
give one in the local setting excluded in Section4.
This proves Theorem1 over every characteristic-zero field.

## 6. Scope and classical relations

The result concerns zero-constant quadratics satisfying(1) with FIXED
rational weights of odd denominator. It does not prove Strong Factorial
for arbitrary independent quadratic coefficients or arbitrary nonzero
constants, and does not close unrestricted SFC(m).
For rank-one quadratic Hessians take every c_ij=1; in dimensions at
least three this choice also admits sign patterns of higher Hessian rank.
Taking every c_ij=0 gives arbitrary diagonal quadratics with arbitrary
linear coefficients in the SAME variables. Taking c_ij in{0,1} includes
sparse interaction graphs and disjoint rank-one blocks. No condition is
placed on the characteristic-zero linear coefficients d_i.

Formula(4) is a weighted matching or monomer-dimer polynomial divided
by the vertex factors D_i. Its use here is justified entirely by the
native factorial coefficient calculation, including the factors2
which remove repeated edges and vertices of cross-degree at least two.
General matching-polynomial theorems are not proof premises.

Closest classical Strong Factorial sectors include powers of linear forms
and independent-variable linear adjunction in Rocks, On the incompatibility
of diophantine equations arising from the Strong Factorial Conjecture,
author2015 preprint, Theorem3.3 and Theorem3.12/Corollary3.13,
https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf.
Their exact statements have been inspected. They supply context; every
mathematical input to the theorem above is supplied locally.
The linear coefficients here may use the SAME variables as the quadratic
part, so the general result is not being inferred from disjoint adjunction.
