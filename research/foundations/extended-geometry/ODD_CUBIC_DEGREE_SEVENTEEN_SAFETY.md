# Degree-seventeen safety in the full odd-cubic class

Date: 19 November 2026. This is a statement about regular winding-one odd-cubic-height phases t(C+t²D)/x, not arbitrary degree17 Image or sphere polynomials.

**Theorem OddCubicDegreeSeventeenSafety.** Every pure phase of exact least degree17 in this entire class is eventually safe for every fixed sphere and original polynomial marker. There is no pure unit-center phase of exact degree17.

The zero-center case is safe in all degrees by the retained classification. Constant normalized C gives C1,D-1 and degree2. Normalize the remaining unit center C0=-D0=1. The canonical polynomial-descent theorem gives a normalized polynomial w of degree g, a pure quotient Cbar,Dbar with full stable field, and exact degrees

    c=g(cbar-1)+1, e=g(ebar-1)+1,
    u=ubar(w), C=Cbar(w)/w', D=Dbar(w)(w/x)²/w'^3.

The quotient has a genuine rational Mobius return j fixing0. Write Rbar for its rational primitive part, Omega for its logarithmic part, and Tbar=1/(Cbar ubar'). Then

    ubar(jz)=-ubar(z), Rbar(jz)=Rbar(z)-2ubar(z),
    j*Omega=Omega, j*(dz/Cbar)=dz/Cbar-2dubar.          (Q1)

No extension of the original field involution to the original x-line is assumed. All return arguments below take place on the proved rational quotient.

## 1. Rational-return infinity lemmas

For a nonconstant quotient center Cbar, Rbar is regular at infinity. If ubar has a pole there of order r>=1, j cannot fix infinity: otherwise both Rbar and Rbar(jz) would be regular there. Thus A=j(infinity) is finite and a pole of Rbar of order r. It is a Cbar-root of multiplicity r+1.

If cbar>=3 and ubar instead vanishes at infinity with order s>=2, j MUST fix infinity. Otherwise j*(dz/Cbar) would have a positive zero at infinity by(Q1): both dz/Cbar and dubar have positive order there. This would make dz/Cbar have a zero at the finite point A=j(infinity), impossible for a reciprocal polynomial differential. Hence a moving-infinity zero has order exactly1, except that cbar2 is already covered by the retained one/two-root classification.

When j fixes0 and infinity it is -z. Then ubar is odd, and the geometric Cbar-roots are paired under negation: repeated roots by the pole support of Rbar(jz), simple roots by their nonzero logarithmic residues. There is no root0. Root pairing makes Hbar even, and the retained two-root exclusion then gives Hbar>=4. The zero order of ubar at infinity is odd.

For each paired root, with multiplicities m<=n and common ubar pole order k>=0, one has k<=n-1 and2k<=3m. Summing gives the retained reflection capacity inequality

    cbar >= (5/6)deg(ubar)+Hbar/2.                     (Q2)

If ubar has infinity zero order s, deg(ubar)>=s+1. Since s=(3cbar-ebar-2)/2, this implies

    5ebar>=3cbar+6Hbar.                               (Q3)

## 2. Two elementary quotient obstructions

A pure unit-center quotient with cbar4,ebar8 (least degree10) is impossible. A single root gives degrees12 or18, two roots are excluded, and four roots are squarefree. Three roots have pattern2+1+1. The retained one-double-root proof gives degree3(q+2)=12 in its finite-alpha branch, while its inversion branch requires q=degCbar-2 divisible by4; here q2. Thus no root pattern remains. This uses the degree alternatives in that proof, not its merely odd-prime theorem statement.

A rational logarithmic differential with exactly two nonzero finite residues summing0 also cannot occur in a pure quotient. It equals a*dlog((z-r1)/(z-r2)), a!=0. Log invariance makes that Mobius function equal at z and jz after normalization at the fixed basepoint, forcing jz=z. This is the same injective-log argument as the retained two-root exclusion and remains valid when other Cbar-roots have zero residue.

Rationality of the selected square root makes c,e have the same parity. Exact degree17 therefore splits exhaustively into e15,c<=15 or c17,e<=15, with quotient dominance preserved by polynomial descent.

## 3. Original D-dominant case: e15,c<=15

The quotient also has cbar<=ebar and

    14=g(ebar-1).

Thus g is a positive divisor of14. A finite nonzero value of u at infinity would force the original degree to be composite, by the retained finite-limit lemma; this is impossible for17.

### 3.1 A pole of ubar at infinity

For a pole order r>=1, ebar=3cbar+2r-2, so

    14=g(3cbar+2r-3).

With cbar>=2 the only integer possibilities are

    (g,cbar,r)=(1,3,4),(1,5,1),(2,2,2).

The first and last would require a Cbar-root of multiplicity r+1 exceeding its total degree. Only cbar5,r1 remains. The point A=j(infinity) must then be a DOUBLE Cbar-root.

The possible multiplicity patterns containing a double root are3+2,2+2+1,2+1+1+1. The first is excluded by the two-root theorem; the last by the retained one-double-root odd-prime exclusion, at exact degree17. In the middle case, Omega is regular at A because it is regular at infinity and j-invariant. Hence the residue at that double root is zero. The remaining double root and simple root have nonzero opposite residues (the simple residue is nonzero and all finite residues sum0). The two-log obstruction of Section2 excludes this final case.

### 3.2 A moving-infinity zero of ubar

This has zero order1 by Section1. Then ebar=3cbar-4 and

    14=g(3cbar-5).

The only integer possibilities are (g,cbar,ebar)=(2,4,8) and(14,2,2). The first is excluded by the specific degree10 quotient obstruction. The second is excluded by the retained degree<=5 unit-center classification (or one-root classification). No moving-zero case remains.

### 3.3 A fixed-infinity zero of ubar

Now j=-z. Tbar has growth order(cbar-ebar)/2. If cbar<ebar, its value at the fixed infinity would be0, whereas sigma(Tbar)=2-Tbar requires value1. Thus cbar=ebar=1+14/g. The four possibilities are

    (g,cbar)=(1,15),(2,8),(7,3),(14,2).

For cbar15 or3, the infinity zero order s=cbar-1 is even, contrary to anti-invariance. The cbar2 case is already excluded by the low-degree/one-root classification. For cbar8, the root count is even and at least4, while deg(ubar)>=s+1=8. Inequality(Q2) would give8>=26/3, a contradiction. This explicitly covers the equal-degree quotient branch; it is not omitted from D dominance.

## 4. Original C-dominant case: c17,e<=15

The quotient has cbar>=ebar+2 and

    16=g(cbar-1), cbar<=17.

The cbar2 case is excluded by the one/two-root classification. For cbar>=3, ubar has an infinity zero of order s>=(cbar), so j fixes infinity by Section1. Its root count is even and at least4. Combine(Q3) with ebar<=cbar-2:

    cbar>=3Hbar+5>=17.

Since cbar<=17, equality holds throughout:

    g1,cbar17,ebar15,Hbar4, s17, deg(ubar)18.

Equality in(Q2) forces equality in both paired-root bounds: for each pair m=2k/3,n=k+1. Thus each k is a positive multiple of3. There are two pairs and the sum of their k's is9, so, up to exchange and independent sign choices,

    pole orders(3,6), Cbar multiplicities(2,4),(4,7).

Write the signed roots as r,-r,s,-s with the smaller multiplicities at r,s. They are nonzero and r²!=s². The simple zero at0 and the order17 zero at infinity exhaust all18 zeros of ubar. Consequently

    ubar=K*z/[(z²-r²)^3(z²-s²)^6], K!=0.

Put B=(z²-r²)^4(z²-s²)^7. For a nonzero scalar A,

    1/Cbar=A*(z-r)^2(z-s)^3/B.

The necessary exact derivative identity is ubar'=(1/Cbar(z)+1/Cbar(-z))/2. Its derivative and even numerator are

    ubar'=K[-17z4+(11r²+5s²)z²+r²s²]/B,
    even(1/Cbar)=-A[(2r+3s)z4
          +s(3r²+6rs+s²)z²+r²s³]/B.

The constant coefficients give K=-As. The leading coefficients then give the literal equation r+10s=0. The middle coefficient gives4r²-3rs+2s²=0. Substituting r=-10s leaves432s²=0, impossible in characteristic zero. Thus the last C-dominant case is excluded.

## 5. Exact scope

Both degree dominances are exhausted. No unit-center pure phase of degree17 exists; zero-center pure phases are safe. This proves the theorem. No claim is made that degree17 is impossible for wider height profiles or arbitrary Image polynomials, and no finite moment prefix is used for this exclusion. The separate degree19 moving-infinity analysis is outside this theorem and checker.


## 6. Inputs and exact controls

The canonical polynomial-descent proof is a separately reviewed input, currently
at ODD_CUBIC_CANONICAL_POLYNOMIAL_DESCENT.md in the same proof packet. Its
Luroth statement over K and genuine rational quotient are indispensable here;
no original rational return is presumed. The retained one-pole, mixed-center,
degree-seven and odd-start proofs provide exactly the antecedents identified
in the source ledger. No new classical classification is imported.

Run `python3 verify_odd_cubic_degree_seventeen_safety.py`. The standard-library
checker exhausts the displayed integer degree branches, reconstructs the
whole forced product-rule and reciprocal-center numerators at the sole
reflection residual, and verifies the nonzero integer certificate432. The
checker rejects disabled assertions and a changed certificate coefficient.
These controls verify the algebra of the written reduction; they do not
substitute for its all-order field and differential arguments.

The [source ledger](odd_cubic_degree_seventeen_safety_source_ledger.json) pins the exact retained inputs, including the [canonical polynomial quotient theorem](ODD_CUBIC_CANONICAL_POLYNOMIAL_DESCENT.md). The [checker](verify_odd_cubic_degree_seventeen_safety.py) and [receipt](odd_cubic_degree_seventeen_safety_verification.json) bind the finite arithmetic and full rational-identity certificate to this proof.

```bash
python3 research/foundations/extended-geometry/verify_odd_cubic_degree_seventeen_safety.py \
  --output research/foundations/extended-geometry/odd_cubic_degree_seventeen_safety_verification.json
```
