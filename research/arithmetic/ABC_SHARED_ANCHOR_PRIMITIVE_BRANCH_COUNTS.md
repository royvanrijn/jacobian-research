# Recovering a shared complete anchor from primitive branches

A family of exact additive relations can share one original
complete value before dividing out the separate gcds.
Recover the anchor once from the union of the normalized
supports. Prime powers that disappear on every branch are
recorded once, while each branch receives its own radical
discount. A second bound counts one primitive branch
globally and the others in fixed-anchor fibres.

## Counting with the union of primitive supports

Fix an integer \(m\ge1\), and put \(n=X^L\) with the
height limit at each sufficiently small fixed \(L\).
Consider \(m\) additive relations with signs \(\pm1\) whose first
entry in absolute value is the same positive integer \(W\).
Let \(Y_i,V_i>0\) be the two other absolute values on
branch \(i\). All values are of fixed polynomial height.
The signed relation on a branch always rearranges its
three nonzero absolute values as a positive sum. Write
\[
     G_i=\gcd(W,Y_i,V_i),\qquad
     (w_i,y_i,v_i)=(W,Y_i,V_i)/G_i .
\]
The sum after reordering is primitive, so its three
coordinates are pairwise coprime. Assume each of these
primitive triples has height at least \(X^\theta\) for
one fixed \(\theta>0\).

Use a common fixed prime pool \(P\) with bounded normalized
radical mass and polynomial height. A varying pool may
instead be specified by a fixed side record, with the
cost of that record added separately. Write
\[
 \rho(A)=\log_n\operatorname{rad}_{P^c}(A).
\]
Let \(R\) be an upper bound for \(\rho(W)\),
\(B_i\) an upper bound for \(\rho(Y_i)+\rho(V_i)\), and
\(e_i=\rho(G_i)\). Write \(t_i\) for the normalized raw
mass of the primes outside \(P\) that divide \(G_i\)
and are absent from the normalized product \(w_i y_i v_i\).
Let \(t_0\) be the mass of the primes outside \(P\) on
\(W\) that are absent from all \(m\) normalized products.
The parameters can be binned at any fixed small mesh
before the height limit, and all the original raw
radical budgets are assumed bounded in base \(n\).

**SAB1.** The number of records consisting of the
original anchor and all original branches, with the
indicated parameter bounds, is
\[
 \ll n^{\,mR/3+\sum_i B_i/3-\mathcal D+\eta}X^\varepsilon,
 \qquad
 \mathcal D=\frac13\sum_i(2e_i+t_i)-t_0 .
\]
Here \(\eta>0\) is arbitrary at sufficiently small fixed
budgets, using lower endpoints for the discounts and
upper endpoints for the prices with the mesh allowance
absorbed in \(\eta\). The bound includes the original
gcd values and the signs, whose latter choices are
bounded in terms of \(m\).

For a prime dividing \(G_i\) outside \(P\), its ordinate
support is counted three times before division. It is
counted once after division if it survives somewhere
on the primitive triple, and zero times otherwise.
For primes not dividing \(G_i\) the count is unchanged.
Consequently the primitive effective radical budget on
branch \(i\) is at most \(R+B_i-2e_i-t_i\).
The effective one-third estimate counts the primitive
branches at the sum of these prices divided by three.
Here one bins their positive physical heights on a fixed
mesh between \(X^\theta\) and the prescribed polynomial
height upper bound, and translates the raw and effective
caps to each shell. The total raw cap on each branch is
bounded in base \(n\), so the shell budgets are small
at sufficiently small fixed \(L\).
[Effective one-third estimate with prescribed supports (EB4–EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

Given the normalized products, enumerate the missing
outside part of \(\operatorname{rad}(W)\) of mass \(t_0\)
at cost \(n^{t_0+o_X(1)}\). All primes of \(W\) are
now in a pool of known polynomial height consisting of
\(P\), that record and the supports of the normalized
products. Supported powers give height-subpower many
choices for \(W\). Then every \(G_i=W/w_i\) is fixed,
and so are all original branch values. This proves
the count. The number of branches is fixed, so all the
height and exponent bins cost subpower.

In particular each prime in the last support record
is erased on every branch, so \(t_i\ge t_0\) and
\(e_i\ge t_i\). The saving has the exact expression
\[
 \mathcal D=(m-1)t_0+
        \frac23\sum_i(e_i-t_0)+\frac13\sum_i(t_i-t_0).
\]
For example, a fixed positive raw mass erased on all
branches yields the credit \((m-1)t_0\) even if none
of it survives in any normalized coordinate. A prime
erased on one branch and visible on another needs no
support record for restoration. For \(m=1\), the
saving reduces to two thirds of the mass of common
primes that survive somewhere on the primitive output.

### The two cheapest normalized roles

There is a sharper version of SAB1 retaining separate
outside caps on the original three coordinates.
For each role \(T\) on branch \(i\), let \(e_{i,T}\)
be the raw mass of outside primes of \(G_i\) remaining
in its normalized value. These masses sum to
\(e_i-t_i\). The exact normalized cost of that role is
\(\rho(T)-e_i+e_{i,T}\). Equivalently put
\[
               k_{i,T}=\rho(T)-(e_i-t_i)+e_{i,T},
\]
so that the normalized cost is \(k_{i,T}-t_i\).
Then an upper exponent for recovering all branches is
\[
 \frac12\sum_{i=1}^m\min_{T\ne U}(k_{i,T}+k_{i,U})
                          -\sum_{i=1}^m t_i+t_0+\eta.
\]
In each primitive output, first fix the \(P\)-supported
parts uniformly at subpower cost. Apply MCB2 with the
resulting arbitrary fixed coefficients and the singleton
exact valuation roots outside \(P\). The physical height
floor ensures the required \(n^C\) height condition.
Its two-cheapest-role price is
\(\tfrac12\min(k_{i,T}+k_{i,U})-t_i\).
Restoring \(W\) from all supports costs \(t_0\) as in
SAB1. Taking the minimum over pairs can only improve
the one-third-of-sum price. This proof also improves
the global branch of SAB2 by using
\(\tfrac12\min(k_{j,T}+k_{j,U})\) instead of
\((R+B_j-2(e_j-t_j))/3\).
[Two cheapest costs at arbitrary fixed coefficients (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)

## Counting one branch globally and the rest in fixed fibres

For each branch put
\[
                 c_i=\min\bigl\{\rho(y_i),\rho(v_i)\bigr\}.
\]
Again the parameters can be upper-binned at a fixed
mesh. Use the same raw and height conditions as in
SAB1.

**SAB2.** Another upper exponent for the same records is
\[
        \min_{1\le j\le m}
          \left(\frac{R+B_j-2(e_j-t_j)}3
                          +\frac12\sum_{i\ne j}c_i\right)
                    +\eta .
\]
Indeed count the primitive triple in branch \(j\),
and then recover \(W\) by recording just the primes
of \(W\) missing from its normalized supports outside
\(P\). Their mass is at most \(t_j\). As before,
enumerating the powers on known support costs subpower.
This gives exponent at most
\((R+B_j-2e_j-t_j)/3+t_j\).
For every other branch enumerate \(G_i\mid W\) at
subpower cost. Its complete anchor \(w_i=W/G_i\)
is now fixed, with original raw radical budget bounded
in base \(n\). On the two free letters apply the
fixed-letter estimate with the outside caps and the
original pool. Its cost is \(c_i/2\), up to allowances,
also when the sign orientation or the ordering of
that primitive branch differs from branch \(j\).
The positive physical height assumption makes the
same shell conversion valid. Sum the finitely many
profile and role choices, and choose the best \(j\).
[Fixed complete letters with prescribed pools (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

Either bound applies when the costs vary by branch
and gives a bound on any event with prescribed
budgets and common-support requirements. The two
counts use different information: the first shares
the missing support over all branches, while the
second uses a smaller fixed-anchor price on the
remaining primitive triples.

## Actual cross-row and same-target stars

For example, take genuine primitive rows \(A+B=C\)
in the original height and radical caps. Let the anchor
be a row \(j\), with marked complete role \(W_j\)
and opposite compared value \(Z_j\).
Select \(m\) other rows \(i_\ell\) having the same
shallow residue and a prescribed numerical neighborhood
on the compared role. For a sum splice use
\(Y_\ell=W_j+(Z_{i_\ell}-Z_j)\) and
\(V_\ell=|Z_{i_\ell}-Z_j|\); in another sign orientation
insert the appropriate fixed signs. Keep the positive
splices and distinct compared values. The original
caps are on the source rows. The budgets on \(Y_\ell\)
and \(V_\ell\) are hypotheses for using the two counts.

In the fine-layer shallow charts of HPS3, first
retain the typical cells of size at most
\(n^{d-s+o(1)}\), using the JRM2 collision bound.
A bad near-splice pair from HPS3 extends to at most
\(n^{(m-1)(d-s)+o(1)}\) such stars. Thus discarding
the bad pairs on each of the fixed number of
branches loses only a power-small fraction of the
full-exponent typical star population. On the
retained pairs the common gcd before normalization
is at most \(M n^{r_W+\xi}\), using
HPS3 with fixed tolerance \(\xi>0\) and the original
marked radical bound \(n^{r_W}\). Thus if \(W_j\)
has positive numerical height, the primitive heights
are bounded below by a fixed power at small budgets.
More generally that height lower bound can be checked
directly. In an actual complete-label matching, the
restored \(W_j\) gives the original anchor row, and
each restored signed \(V_\ell\) then gives the
original compared complete value \(Z_{i_\ell}\).
This reconstructs the original rows. If complete
values are repeated, restore the corresponding
complete-degree factors \(n^{u_W+m u_Z}\).
In a same-complete-target ripple with all
\(W_{i_\ell}=W_j\), each \((W_j,Z_{i_\ell})\)
determines its row, so only the anchor's degree
factor is needed.
[Shallow near-splices (SRS1–SRS2)](ABC_SHALLOW_REGISTER_SPLICES_AND_CARRIES.md#near-splices-in-both-residue-and-numerical-scales)
· [Powerful excess control for near-splices (HPS3)](ABC_HYBRID_CELL_PRIVATE_SUPPORT_AND_SPLICE_GCD.md#gcds-of-a-splice-with-its-compared-true-target)

For a complete matching and mixed depth \(s<d\) at
\(d=1/3\), the full-exponent cell scale has
\(n^{(m+1)d-ms-o(1)}\) ordered anchored stars:
sum the \((m+1)\)-st powers of the cell occupancies,
using the alphabet bound \(n^{s+o(1)}\), and remove
the repeated choices. The same calculation applies
to prescribed degree bins with a positive cell exponent.
Comparing this population with SAB1, a region with
fixed total branch budget \(\sum_i B_i\) and discount
\(\mathcal D\) is power-small whenever
\[
       \frac{mR+\sum_i B_i}{3}-\mathcal D
                          < (m+1)d-ms
\]
by a fixed margin. The corresponding criterion
using SAB2 replaces the left side by its mixed
fixed-anchor exponent. These are restrictions on
actual capped source stars which allow large common
multipliers after the outputs are formed. Their
application to the unrestricted bound still
requires an arithmetic extraction of a favorable
budget or support pattern on these stars.
