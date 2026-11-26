# Fresh radicals and power defects of primitive pair comparisons

Put \(n=X^L\). Work with actual positive primitive rows
\(A+B=C\) in \(X/2<C\le X\), keeping the total cap
\(\operatorname{rad}(ABC)\le n\) and the three individual caps
\(\operatorname{rad}(T)\le n^{1/3+\tau}\).
For two rows and two distinct original roles \(U,V\), write
\(W\) for the remaining role and put
\[
 g_U=\gcd(U_1,U_2),\qquad g_V=\gcd(V_1,V_2),\qquad g=g_Ug_V,
 \qquad
 P=\frac{U_1V_2}{g},\qquad Q=\frac{U_2V_1}{g}.
\]
A prescribed prime pool \(S\), possibly depending on \(X\), has
product at most \(X^{C_0}\) for fixed \(C_0\). The fresh radical
of a nonzero value on this pair keeps only the primes outside
\(S\) and outside the six original letters. Its fresh
\(m\)-defect radical keeps those fresh primes for which the
valuation is not divisible by \(m\). Each result below permits
permuting the original roles.

## Reduced cross products and removal of old support

Primitivity on each row gives
\[
               \gcd(U_1V_2,U_2V_1)=g_Ug_V.
\]
This follows prime by prime, since each prime occurs on at most
one of \(U_i,V_i\) on each row. Thus \(P,Q\) are positive coprime
integers. After fixing \(P,Q\), the four divided factors are
obtained by factoring \(P,Q\), at height-subpower cost.
Consequently if \(\operatorname{rad}_{S^c}(g)\le n^\sigma\),
the original ordered pairs over fixed \(P,Q\) number
\(O_\varepsilon(n^\sigma X^\varepsilon)\). Record the outside
support of \(g\), assign the primes and the exact supported
powers to \(g_U,g_V\), and recover four complete letters.
Each row then follows from the signed additive equation.

Fix a finite list of nonzero homogeneous binary forms with
integer coefficients, positive degrees, and nonzero values at
\((1,0),(0,1)\). Enlarge \(S\) by the first and last coefficient
primes. For every fixed \(\sigma>0\), on a family of
\(N\ge n^{1/3-\tau}\) rows, all but
\(O_{\varepsilon}(N^2 n^{-c}X^\varepsilon)\) ordered pairs
have each of the following properties for some \(c>0\), at
sufficiently small fixed \(\tau,L\).

- Their shared original prime support outside \(S\) has
  radical at most \(n^\sigma\).
- The reduced comparison is tall, in the sense that
  \(\max(P,Q)>X^{h_0}\) for a fixed numerical \(h_0>0\).
- Each specified \(F(P,Q)\) is nonzero. Its support on the
  original letters outside \(S\) has radical at most \(n^{3\sigma}\).

For the first property, fix the first original row and use
its complete prime support as a prescribed zero-price pool.
By EB5 the number of second rows sharing outside mass at
least \(\sigma\) is \(O(n^{1/3-\sigma/3+\eta}X^\varepsilon)\)
for any fixed small \(\eta>0\) at sufficiently small fixed \(L\).
The bound is uniform in the first row. For the height property,
fix the first row and choose \(g_V\mid V_1\) at divisor cost.
If both \(P,Q\le X^{h_0}\), then \(V_2=g_V z\) for an integer
\(1\le z\le X^{h_0}\). CF1 on the complete letter \(V_2\),
with the original three caps, gives a fixed power saving
per anchor when \(h_0\) is sufficiently small. These two
exceptions remain power-small relative to \(N^2\) when
\(\tau\) is sufficiently small.

To check the third property, a prime on one of the four
crossed letters can divide \(F(P,Q)\) outside the fixed
coefficient pool only when it is shared by the two rows.
For an omitted-role prime \(p\mid W_1\), there is a fixed
sign \(s\in\{1,-1\}\), determined by the orientation, with
\(U_1\equiv sV_1\bmod p\). The factors \(g,V_1\) are units
at \(p\), and
\[
          g^{\deg F}F(P,Q)\equiv
               V_1^{\deg F}F(sV_2,U_2)\pmod p.
\]
The symmetric congruence for \(p\mid W_2\) uses
\(F(U_1,sV_1)\) and a unit factor \(V_2^{\deg F}\).
Fix the other row. Unless the displayed anchor value is
zero, its support is a prescribed polynomial-height pool.
Apply EB5 to the varying original row to remove the
power-small fraction where its \(W\)-overlap has mass
above \(\sigma\). A zero anchor value forces one of finitely
many exact original ratios, each determining the row by
primitivity. Discarding such anchors costs \(O(N)\) pairs
per form. The equation \(F(P,Q)=0\) itself gives at most
\(\deg F\) pairs per first row by the same projective
uniqueness for the second row. The shared support and the
two omitted-role overlaps give \(3\sigma\).

[Fixed-pool radical pricing (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Complete-letter affine features (CF1)](ABC_COMPLETE_LETTER_SMALL_FEATURES.md)

## An ordinary radical floor for every fixed comparison

**GP1.** For each fixed \(\delta>0\), each fixed
polynomial-height pool, and any finite list of fixed squarefree
binary forms as above, every \(F(P,Q)\) has fresh radical
larger than \(n^{2/3-\delta}\) on all but
\(O_\varepsilon(N^2 n^{-c}X^\varepsilon)\) pairs, for some
\(c>0\), when \(\tau\) and the fixed positive \(L\) are
sufficiently small. If the forms have no common geometric
root in pairs, their fresh supports are pairwise disjoint
after adding the primes of their pairwise resultants to \(S\).
The orders of choices are the fixed forms, \(\delta\) and
the cap on pool height, then the fixed tolerances and \(L\),
then the height limit.

For a form of degree at least three, prescribe its exact
nonzero value using its outside-\(S\) radical. The possible
supports and then their exact powers and sign cost at most
\(n^{2/3-\delta+3\sigma}X^\varepsilon\) on the pairs where
fresh mass is at most \(2/3-\delta\). This uses the preceding
old-support bound. Evertse's Thue–Mahler estimate gives
height-subpower primitive completions per exact value, since
the squarefree form has at least three distinct geometric
linear factors. For a squarefree quadratic, write
\(F(P,Q)=aP^2+bPQ+cQ^2\), where \(a\ne0\) and
\(b^2-4ac\ne0\). Complete the square and apply the integral
quadratic norm estimate to
\[
              (2aP+bQ)^2+(4ac-b^2)Q^2=4aF(P,Q).
\]
This also gives a uniform height-subpower count. In either case
the final reconstruction of \(g\) costs at most \(n^\sigma\)
up to height-subpower factors.

For \(F(P,Q)=aP+bQ\) with \(ab\ne0\), instead use the
primitive positive additive auxiliary obtained from
\((|aP|,|bQ|,|F(P,Q)|)\) in its additive ordering. Its
normalizing gcd divides \(|ab|\), which has a bounded number
of possibilities. Its numerical height is at least a fixed
multiple of \(X^{h_0}\) on the retained pairs. Enlarge \(S\)
by the coefficient primes and fix all supported powers of the
auxiliary there at height-subpower cost. The remaining caps
in base \(n\) are \(2/3+2\tau,\ 2/3+2\tau,\
2/3-\delta+3\sigma\), allowing the roles to be permuted.
MCB2 for arbitrary fixed coefficients gives exponent at most
\[
               2/3-\delta/2+\tau+3\sigma/2+\eta .
\]
Record \(g\) at price \(\sigma\) afterwards. In both degree
ranges choose \(\sigma,\tau,\eta\) small in terms of \(\delta\)
and the universal CF1 gap. Comparing the counts with
\(N^2\ge n^{2/3-2\tau}\) proves the claim. Disjointness of fresh
supports follows by reducing the two forms at a putative
shared prime and using the coprime coordinates \(P,Q\).

[Evertse, Theorem 6](https://ir.cwi.nl/pub/2463/2463D.pdf#page=14)
· [Quadratic norm count, Browning–Verzobio, Lemma 2.3](https://arxiv.org/html/2608.24512v1#S2.SS1)
· [Two cheapest radical caps (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)

Here is a simultaneous consequence using all three
orientations. Set
\[
 \Delta=A_1B_2-A_2B_1,\qquad K_{AB}=A_1B_2+A_2B_1,\qquad
 K_{AC}=A_1C_2+A_2C_1,\qquad K_{BC}=B_1C_2+B_2C_1 .
\]
The cross difference agrees with \(\pm\Delta\) in every
orientation. The other three satisfy
\[
 K_{AC}-K_{AB}=2A_1A_2,\qquad
 K_{BC}-K_{AB}=2B_1B_2,\qquad
 K_{AC}+K_{BC}=2C_1C_2.
\]
In each orientation the sum plus or minus its corresponding
difference is twice an original product. Hence, away from
the six original letters and 2, all four fresh supports are
disjoint. GP1 applied to \(s\pm t\) shows that the product
of these four fresh radicals is larger than
\(n^{8/3-4\delta}\) on the same generic set of pairs.
This argument keeps the positive-multiplicity and
occurrence regimes, and needs no common modulus or target.

## Fresh power defects and independent multiplicative classes

**GP2.** In GP1, take each \(F\) squarefree of fixed degree
at least five. For every fixed integer \(m\ge2\), the same
\(2/3-\delta\) floor holds for the fresh \(m\)-defect
radical of \(F(P,Q)\). The same estimates hold simultaneously
for any finite list of orders chosen after fixing \(L\) and
before the height limit, with inner constants depending on
the list.

Indeed it is enough to use a prime \(p\mid m\).
After the old-support removal proved above, a fresh
\(p\)-defect radical at most \(n^{2/3-\delta}\) supplies
the exact signed \(p\)-power-free coefficient of \(F(P,Q)\)
at price \(2/3-\delta+3\sigma\), up to height-subpower.
Here \(p\) is fixed during enumeration of the nonzero
valuation residues. BT1 bounds the primitive points per
exact coefficient by \(X^{o_X(1)}\). Its splitting field is
fixed at a chosen \(p\), its class group contribution is
constant and the coefficient has polynomial height. Add
\(\sigma\) for the two original gcds. This gives the
same strict pair-count deficit as the nonlinear part of GP1.
The only small-budget thresholds here come from original
support removal, independent of the chosen order. Thus
one may also allow a finite order list at each fixed \(L\).
[Fixed-form power twists (BT1)](ABC_BINARY_POWER_TWISTS.md#primitive-power-and-powerful-values-of-binary-forms)

For pairwise geometrically coprime \(F_1,\ldots,F_J\), it
follows that the values \(F_j(P,Q)\) are independent in
\(\mathbb Q^\times/\mathbb Q^{\times p}\) modulo the subgroup
generated by the six original letters and primes in \(S\).
More precisely, for any exponent vector \(e\in\mathbb Z^J\),
the fresh \(p\)-defect radical of
\(\prod_j F_j(P,Q)^{e_j}\) is greater than
\(n^{(2/3-\delta)h}\) whenever \(h>0\), where
\(h=\#\{j:p\nmid e_j\}\). Each active component supplies
its own disjoint fresh factor, so this holds simultaneously
for all vectors. For composite \(m\), the same quantitative conclusion
holds with \(h=\#\{j:m\nmid e_j\}\). For each active
coefficient choose \(p\mid m\) with \(v_p(e_j)<v_p(m)\).
The fresh \(p\)-defect of \(F_j\) forces an equal lower
bound on the fresh \(m\)-defect of \(F_j^{e_j}\).
The disjoint supports make this simultaneous for all
vectors and prove independence modulo \(m\).
For example, one may take any finite list of distinct
homogeneous cyclotomic polynomials of distinct prime
indices at least seven.

## A one-row budget for genuine power defects

**GP3.** Use the same kind of fixed pairwise geometrically
coprime forms of degrees at least five, now evaluated at
one primitive pair \((A,B)\). Assume each form also has no
geometric factor in common with \(s,t,s+t\). Fix an order
\(m_j\ge2\) for each form and a prescribed polynomial-height
pool \(S\). Write \(q_j\) for the normalized fresh
\(m_j\)-defect radical mass outside \(S\) and \(ABC\), and
\(r=\log_n\operatorname{rad}(ABC)\). Under the original
cap \(r\le1\), an additional joint constraint
\(r+\sum_j q_j\le B_*\) gives the upper count
\[
            O_{\eta,\varepsilon,L}
               \!\left(n^{\min(1/3,B_* /(3+J))+\eta}
                                 X^\varepsilon\right)
\]
for every \(\eta>0\) and sufficiently small fixed \(L\).
The constants may depend on the fixed forms and orders
and the prescribed pool height bound. For a fixed \(m_j\)
one uses a prime divisor in BT1 as above.

Add the resultant primes to \(S\), so the defect supports
of the outputs are disjoint from the original letters
outside \(S\) and can be counted at their stated prices.
Partition \(r,q_j\) into intervals of a sufficiently small
fixed width, avoiding any exchange of limits. On each box
the ordinary effective-budget bound under the original
height shell gives exponent at most \(r/3\) up to the
mesh and count allowances. The BT1 count for each output
gives exponent at most \(q_j\) up to the same allowances.
Thus if a box has count exponent \(d\), its caps must
satisfy \(r\ge3d\) and all \(q_j\ge d\) up to those errors.
Summing and optimizing gives the claimed coefficient.
A fixed budget \(B_*<1+J/3\) therefore gives a strict
saving on this sector even if the ordinary output radicals
are large. An unrestricted improvement by this route needs
a uniform joint-defect budget on sufficiently many actual
rows, with its saving exceeding the retained population loss.
[Ordinary cyclotomic output budgets for one row (BFR3–BFR4)](ABC_BINOMIAL_FACTOR_RADICAL_PRICES.md#fresh-cyclotomic-radicals-on-original-rows)

## Weak degree limits leave the exact defects undetermined

There is a formal obstruction to extracting the required cheap
defect support from the scaled degree laws alone. Choose integers
\(R(L)\to\infty\) sufficiently slowly that
\(L M(L)\to0\) for \(M(L)=4R(L)!\). On any formal exact degree
law replace each degree \(e\ge1\) by
\[
                 e'=1+M(L)\left\lfloor (e-1)/M(L)\right\rfloor .
\]
Every scaled coordinate moves by at most \(LM(L)\). The
normalized masses, their limiting laws, and all their first
moments are preserved up to \(o(1)\), including for measures
with an escaping first-moment part. Yet all retained degrees
are odd and eventually nondivisible by every fixed \(m\ge2\).
Thus a profile with the same limiting degree information can
have full mass in every fixed-order defect, including the
mixed fourth-power defect of APC9. This is an exact statement
about degree laws. It supplies no additive realization of
those laws by actual primitive triples.
[Mixed fourth-power defect (APC9)](ABC_ADAPTIVE_POWER_COEFFICIENT_PRICES.md#a-mixed-fourth-power-defect-replaces-exact-power-records)
