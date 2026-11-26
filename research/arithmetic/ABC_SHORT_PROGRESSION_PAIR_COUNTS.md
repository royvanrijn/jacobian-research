# Aggregate small-radical pairs in short progression blocks

A single count controls pairs across every unit residue and
indexed width for the same step. Reducing by their common
gcd makes the remaining index a third letter. The short-letter
dictionary bound charges its radical together with the two
original radicals, and the original caps also bound the price
of restoring the gcd. The step can have arbitrary radical.

## An aggregate pair bound at a fixed real gap

Fix \(A_0,s_0,\nu,\eta>0\) and \(0<\gamma<1\). For all
sufficiently small fixed \(L>0\) put \(n=X^L\). Fix
\(a,b\in[0,A_0]\), \(s\in[0,s_0]\), and a positive integer
step \(M\), which may depend on \(X\). Count ordered pairs
\[
  \gamma X\le U<V\le X,\qquad
  V-U=M j,\quad 1\le j\le n^s,\qquad
  \gcd(U,M)=1,\qquad
  \operatorname{rad}(U)\le n^a,\quad
  \operatorname{rad}(V)\le n^b .
\]
The residue class modulo \(M\) is allowed to vary with the
pair. Suppose \(M n^s\le X^{1-\nu}\).

**PCP1.** Uniformly in these data the pair count is
\[
  \ll_{A_0,s_0,\nu,\eta,\gamma,L,\varepsilon}
       n^{\frac12\min\{a+b,\ \min(a,b)+s\}+\eta}
       X^\varepsilon .
                                                               \tag{PCP1}
\]
Neither original integer needs a lower valuation floor.
For a common cap \(a=b=d\) the exponent is
\(\frac12\min\{2d,d+s\}\). A fixed \(s<d\) therefore
gives a strict aggregate loss from the exponent \(d\).

Put \(g=\gcd(U,V)\). The unit class and the equation
give \(g\mid j\). Write \(U'=U/g\), \(V'=V/g\),
and \(J=j/g\), so
\[
                         U'+M J=V',\qquad \gcd(U',V')=1 .
\]
Bin \(\log_n g\) at \(\chi\le s\) and the two normalized
quotient radicals at \(r_1,r_2\), with a fixed arbitrarily
small mesh and hence \(J\le n^{s-\chi+O(\mathrm{mesh})}\).
The bin caps imply \(r_1\le a\) and \(r_2\le b\) up to
mesh errors. The two quotients and their difference form
a positive primitive triple. At a fixed dyadic size of
\(g\), the quotient height scale is comparable to \(X/g\).
Its short-to-tall ratio is \(M j/V=O_\gamma(X^{-\nu})\).
Moreover \(g\le n^s\), so the quotient budget is still
sufficiently small and its heights tend to infinity
at every fixed original budget. Equivalently use a
fixed finite normalized mesh and dyadic subpower
refinements, keeping all height comparisons uniform.
Numerical constants below can be chosen using these
bounds and \(\nu\).

Apply the short-coefficient estimate (SHF4). Leave
\(M\) as the fixed coefficient on the short term
and use \(J\) as a degree-one root with dictionary
cap \(s-\chi\). For \(U',V'\) use the ordinary
anatomic factorization with arbitrarily small
coefficient cost and vanishing numerical coefficient
height. The raw root caps can be taken as \(r_1,r_2\)
up to the prescribed allowances. In base \(n\) this
counts primitive quotients at exponent
\[
           \frac12\min\{r_1+r_2,\ r_1+s-\chi,\ r_2+s-\chi\}
                                                                  \tag{PCP2}
\]
before the small errors. It is not required to choose
a residue class, a real interval or another fixed
coefficient on either tall side.
[Short-coefficient estimate (SHF4)](ABC_FIXED_SHIFT_STRICT_HALF_SAVING.md#grouped-dictionaries-with-a-large-coefficient-on-the-short-letter)
· [Anatomic reduction, Proposition 2.1](https://arxiv.org/html/2410.12234v2#S2)

After a quotient pair is fixed, restore \(g\). Its
size gives the cost \(\chi\). Outside the support of
\(U'\) its radical has cost at most \(a-r_1\), and
outside the support of \(V'\) it costs at most \(b-r_2\).
The part supported on either prescribed quotient has
height-subpower enumeration cost. Use whichever of
these three bounds is smallest. In particular the
restoration cost is at most
\[
                              y=\min(\chi,a-r_1,b-r_2)
\]
up to allowances. Even when the restored \(g\) and the
quotients share primes, no copying cost for a shared
prime is charged here. In the three arms of (PCP2),
respectively use
\[
 \begin{split}
    2y&\le (a-r_1)+(b-r_2),\\
    2y&\le\chi+(a-r_1),\\
    2y&\le\chi+(b-r_2).
 \end{split}
\]
They give \(a+b,a+s,b+s\) after doubling the total
exponent. Absorb the bin, factorization, and fixed
support allowances within \(\eta\) and
\(X^\varepsilon\). This proves (PCP1).

## Long steps with a bound on low-valuation mass

Here allow any step \(M\) without the real-span
condition. Fix a small \(\kappa>0\) with
\(2\kappa A_0<1/4\), and suppose that the \(U\) and
\(V\) radicals carried by primes of scaled valuation
below \(2\kappa\) are at most \(n^{\ell_1},n^{\ell_2}\)
respectively, where \(0\le\ell_i\le A_0\). Keep
\(j\le n^s\) and the same unit and height hypotheses.
To distinguish the two conventions, a prime on
the original term is called low here at cutoff
\(2\kappa/L\), and on a primitive quotient at
\(\kappa/L\).

**PCP3.** The pair exponent, up to arbitrary
losses as in (PCP1), is at most
\[
 \frac12\min\{a+b+\ell_1+\ell_2,\
                  a+s+\ell_1+2\ell_2,\
                  b+s+2\ell_1+\ell_2\}.
                                                               \tag{PCP3}
\]
The threshold on the fixed budget may depend on
\(\kappa\), in addition to the stated uniform
bounds. In particular both original profiles with
no atom at scaled degree zero give (PCP1) in
the height-first limit without the real-span
condition.

Normalize by \(g\) as above. A newly low prime on
a quotient had original scaled valuation at least
\(2\kappa\), and then consumes at least
\(\kappa/L\) of its valuation in \(g\). All such
primes have normalized logarithmic radical mass at
most \(Ls/\kappa\), since \(\log_n g\le s\). Thus the
low quotient masses are at most the stated \(\ell_i\)
plus a vanishing error at fixed \(\kappa\). Bin those
actual masses at \(q_1,q_2\). Fix the low full factors
at normalized cost \(q_1+q_2\) up to the mesh and
inner allowances. Their numerical heights are at
most \(\kappa q_i+o(1)\). Both tall quotients have
numerical height \(1-O(L)\), so their remaining
high variable terms still have positive numerical
height. The ordinary anatomic truncation can retain
the high floors by taking its compression index
above the cutoff and its exceptional coefficients
arbitrarily small. The remaining quotient radical
caps are \(r_1-q_1,r_2-q_2\). They use disjoint
supports from the frozen low factors on each
quotient.

Use the extension of (SHF4) with an arbitrary
fixed coefficient and a single low root. Both
tall terms now have the required high floor, and
the third term is \(M J\) with \(J\) the only
low root and normalized size at most \(s-\chi\).
The fixed low factors on the two tall terms
are allowed coefficients. Restore the gcd with
the same price \(y\), where the quotient caps
\(r_i\) include their low factors. The count
exponent is bounded by
\[
 q_1+q_2+y+
 \frac12\min\{r_1-q_1+r_2-q_2,\
               r_1-q_1+s-\chi,\
               r_2-q_2+s-\chi\}.
\]
Use the three inequalities on \(2y\) separately,
then \(q_i\le \ell_i+o(1)\), to obtain (PCP3).
When a limiting degree measure gives no mass
at zero, choose \(\kappa\) at a continuity cutoff
after fixing the desired loss, then pass to
sufficiently small fixed budgets. This keeps the
loss from \(Ls/\kappa\) arbitrarily small.
[Arbitrary fixed coefficient and one low root](ABC_FIXED_SHIFT_STRICT_HALF_SAVING.md#a-prescribed-coefficient-of-arbitrary-height-with-one-low-root)

## Arbitrary-span pair counts and fixed full divisors

By (MCB2), the bound (PCP1) holds for arbitrary
\(M\) and physical span, with the same fixed normalized
index width \(s\), even if the two quotient profiles have
low scaled valuation mass. The calculation of the
primitive quotient cap and gcd cost in (PCP1) is unchanged.
The fixed step of any height is a coefficient on the small
quotient index; its radical does not enter that term's
dictionary.
[Arbitrary-span progression count](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#arbitrary-span-progression-occupancies)

**PCP4.** More generally, in an actual primitive abc
shell prescribe a joint dictionary of full divisors
\((D_A,D_B,D_C)\), of size at most \(n^f\).
After fixing a record, let \(b_T\) bound either the
summed conditional bounded-group costs on \(T/D_T\),
or the numerical size exponent \(\log_n(T/D_T)\).
The smaller bound can be used for each role. Then the
retained rows have count at most
\[
 O\!\left(n^{\,f+(\sum_T b_T-\max_T b_T)/2+\eta}
                       X^\varepsilon\right).
\]
Here the numerical option is needed only when the
chosen bound on it is fixed in base \(n\). Use the
recorded divisors as arbitrary fixed coefficients
in (MCB2), and use a single degree-one variable for
each numerical option. The original positive primitive
sum still has maximum at least \(X/2\). The original
radical caps and exact repair restriction can all be
imposed, and full-factor removal only lowers the
available raw radical caps on the quotients.
For example, if the free caps are \(1/3\) and fixing
a full divisor on the third role leaves a quotient
of normalized radical or numerical size exponent
\(b<1/3-2f\) by a fixed margin, this saves uniformly
below one third across the whole dictionary range.
The condition is on the joint record price, not just
the existence of a large numerical divisor.

## Occupancy moments and the marked abc count

Partition a shell into any disjoint intervals of
length at most \(M n^s\). Let \(k_{\rho,I}\) be
the numbers of elements of a set in unit class
\(\rho\bmod M\) and interval \(I\), using a common
cap \(d\), with the fixed index width required in
(PCP1). Then the corresponding pair estimate
bounds the sum of \(\binom{k_{\rho,I}}2\) over
all cells at once. In particular,
\[
          \sum_{\rho,I}\binom{k_{\rho,I}}2
                     \ll n^{\frac12\min\{2d,d+s\}+\eta}
                                      X^\varepsilon .
\]
The same bound holds for any subfamily of the pairs
in those cells. When \(s<d\) by a fixed margin, a
dictionary of exponent \(d\) loses a negligible
fraction by keeping one value per occupied cell.
In one cell the pointwise upper count is at most
\[
         1+O\!\left(n^{\min\{d/2,(d+s)/4\}+\eta}
                                       X^\varepsilon\right),
\]
as well as the trivial index bound \(1+n^s\).
This transfers the numerical-gap use of (SHF6)
from highly full inputs to arbitrary radical caps.
For highly full inputs (PCP3) recovers the
long-step range and gives its aggregate version.
[Progression occupancies with highly full inputs (SHF6)](ABC_FIXED_SHIFT_STRICT_HALF_SAVING.md#progression-occupancies-with-highly-full-inputs)

For an application to actual abc rows, fix a
dictionary of divisors \(M\mid A\) of size at
most \(n^f\) and impose \(A/M\le n^s\). The two
original letters \(B,C\) are the unit pair, since
\(B\equiv C\bmod M\) and primitivity makes them
coprime to \(M\). Every pair recovers \(A=C-B\);
the input radical caps and its least-repair
condition can remain in the count. If both
\(B,C\) have normalized radical cap \(d\), the
count exponent in a fixed height shell is at
most \(f+(d+s)/2+\eta\) when \(s<d\), with
no physical span restriction by (PCP4).
A fixed margin \(2f+s<d\)
gives an absolute saving on that sector.
For positive complete multiplicity this can
apply even when \(f+s\) exceeds the required
complete dictionary exponent of \(A\).

These estimates control collisions inside a
residue and short index block. In the occurrence
regime the below-diagonal moment is consistent
with one row in each of many different cells.
The original step dictionary and index bounds
are genuine hypotheses for the marked abc
saving. The remaining aggregate incidence
problem has all three supports and can spread
over separated cells.
[Induced occurrence matchings](ABC_OCCURRENCE_INDUCED_ADDITIVE_MATCHINGS.md)

Bajpai, Bennett and Chan give constructions and
conditional restrictions for long progressions of
powerful numbers. Van Doorn gives infinite
three-term powerful progressions with specified
step of square-root size. The bounds here count
pairs and cells at small radical budget and
allow both the residue and the starting index
to vary for a fixed step. The growing-index
and modulus ranges are stated separately.
[Bajpai–Bennett–Chan, Theorem 1.1](https://arxiv.org/abs/2302.03113)
· [Van Doorn, Theorem 1.1](https://arxiv.org/abs/2605.06697)
