# A finite decision for bounded-gap angular and twisted factorial families

Fix an integer \(L\ge1\). Consider all ordinary two-monomial
profiles with distinct columns and
\[
                       \sum_i |b_i-a_i|\le L .
\]
Write \(L_C\) for this sum on any profile \(C\) and \(\Delta_C\)
for its signed sum. Omit equal rows. The number of unequal rows is then at most
\(L\), while their heights and the original dimension are
unrestricted. Write \(Q_j\) for the original factorial moment of
order \(j\), as a polynomial in the ratio of monomial coefficients.

**Theorem.** Each of the following questions admits a terminating
exact decision using real polynomial quantifier elimination.

1. Do all these profiles have disjoint consecutive root arguments
   in every degree?
2. For all these profiles and every integer \(h\), are
   \(Q_n(z)\) and \(Q_{n+1}(((n+1)/n)^h z)\) coprime in every
   degree?
3. Are the same twisted pairs coprime whenever \(L_C+|h|\le L\)?
   One may also fix an integer \(\Delta_0\) and restrict to pairs
   with \(h+\Delta_C=\Delta_0\).
4. For profiles in this family satisfying
   \[
      \sum_{i:\,\max(a_i,b_i)=r}(b_i-a_i)=0
                       \quad\text{for every integer }r,
   \]
   are the original consecutive moments coprime in every degree?

For each question the procedure computes a finite height bound.
A negative decision produces a genuine ordinary profile and a
terminal order, with an integer twist in the second and third questions.
A positive decision produces a finite height cutoff whose
coefficient boxes certify the whole family. A contact with
nonzero twist concerns the corresponding twisted property. Only the contacts with \(h=0\) furnish counterexamples to the
original physical property.

## Low orders and small curvature

Let \(J\ge4\) be a computable integer at least the
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
For \(N=n+1\ge J\) the angular theorem excludes all the contacts
in these questions, including any positive argument rescaling.
For \(n=1\), the linear root is negative real and the quadratic
factorial moment has no real root. The latter is the expected
square of a real binomial with distinct monomials and positive
exponential density, and that binomial is nonzero almost surely
at any fixed coefficient ratio. This also excludes every
positive rescaling at the first pair.

For the finitely many \(2\le n<N<J\), take their computable
rational constants \(\lambda_*(n)>0\) from the
[linked Hermite angular bound](MULTIVARIATE_FIXED_TERMINAL_HERMITE_SEPARATION.md).
With \(d_i=|b_i-a_i|\), \(r_i=\max(a_i,b_i)\), and
\(S=\sum_i d_i^2/r_i\), put
\[
           s_0=\min\left(\frac12,\,
                \min_{2\le n<N<J}\frac{n\lambda_*(n)^2}
                                           {1+n\lambda_*(n)^2}\right).
\]
If \(0<S\le s_0\), every unequal row is shared since a private
row contributes at least one to \(S\). The Hermite bound applies
to both independently normalized arrays, so it excludes all
four kinds of contact. The bound is effective also at \(n=2\),
where its even Hermite anchor has a linear square polynomial.
Thus only \(S>s_0\) and the indicated finite set of terminals
need further tests.

## Exact algebraic boxes at a height cutoff

At an integer cutoff \(B\ge2L\), retain the rows of height
\(r_i\le B\). Call the resulting profile \(C\). Every omitted
row is shared and its smaller exponent is at least \(r_i/2\).
The omitted curvature satisfies
\[
           S_{\rm tail}\le c_B:=\frac{L^2}{B+1}.
\]
Enumerate all finite ordinary cores of height at most \(B\),
gap sum at most \(L\), and curvature at least \(s_0-c_B\).
For the fourth question impose zero signed sum at each height
on these cores. These conditions contain all possible
truncations from profiles still under consideration. Only
consider \(B\) large enough that \(c_B<s_0/2\); every enumerated
core is then nonempty.

Normalize each core moment at both endpoints to obtain the
monic polynomial \(P_j^C\) with constant coefficient one.
Its coefficients are positive algebraic numbers since the
normalizing scale
\[
                 G_j^C=\bigl([z^j]Q_j^C/Q_j^C(0)\bigr)^{1/j}
\]
is a positive algebraic root. The actual normalized shape
changes the \(k\)-th coefficient of \(P_j^C\) by a multiplier
in \([1-\epsilon_B,1]\), for \(\epsilon_B=Nc_B/4\).
Indeed the omitted log-Gamma chord at \(k\) on a row is between
zero and the negative of
\(d_i^2 k(j-k)/(2j\min(a_i,b_i))\). This follows from
\(\psi'(1+x)\le1/x\) and the zero-endpoint quadratic chord
comparison. Summing gives absolute chord size at most
\(jc_B/4\), and \(e^{-x}\ge1-x\) for \(x\ge0\).
Endpoint multipliers are fixed at one. Take \(B\) large enough
that \(\epsilon_B<1/2\).

All normalized polynomials in these boxes have coefficients
at most \(\binom jk\) and positive endpoints equal to one.
Their roots satisfy \(2^{-j}<|z|<2^j\), by the monic root
bound and reversal. For the first question, test for a complex
\(z\) and a real \(v\in[2^{-(n+N)},2^{n+N}]\) at which the
two box polynomials vanish at \(z\) and \(v z\), respectively.
This is exactly the form of a common argument, with positive
\(v\).

For the other questions let \(\Lambda=N/n\) and write
\(\Delta=\sum_i(b_i-a_i)\) for an actual profile. Its physical
argument scale ratio satisfies
\[
          \log(G_N/G_n)=\Delta\log\Lambda+\sum_i E_i .
\]
For each row, \(E_i\) is the logarithm of its individual
scale ratio after removing its signed gap times \(\log\Lambda\).
The [fixed-terminal center estimate](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#normalized-polynomials-and-row-data)
gives \(|E_i|\le 2 S_i\), with \(S_i=d_i^2/r_i\).
This estimate also holds at \(n=2\). Directly, the row's
logarithmic scale at real order \(x\) is
\(\int_{a_i}^{b_i}\psi(1+x t)\,dt\).
The bounds \(0\le1-y\psi'(1+y)\le1/(1+y)\) for \(y\ge0\)
then bound its derivative after subtracting
\((b_i-a_i)/x\) by
\[
           \frac1x\int_{\min(a_i,b_i)}^{r_i}\frac{dt}{1+xt}
                     \le \frac{2 S_i}{x}.
\]
For \(d_i\ge r_i/2\), use \(d_i\le2S_i\) inside the
integral. Otherwise use
\(d_i/(x\min(a_i,b_i))\le2S_i\) for \(x\ge1\).
Integration over \([n,N]\) yields the stated bound.

At an actual integer-twisted contact, the normalized root annuli
and \(S\le L\) therefore give
\[
   |h+\Delta|\log\Lambda\le(n+N)\log2+2L .
\]
Since \(\log\Lambda\ge1/N\) and \(|\Delta|\le L\), it suffices
to test twists with \(|h|\le H_n=L+N(n+N+2L)\) on exact
profiles. In a coefficient box for \(C\), let
\(\Delta_{\rm tail}\) be the signed gap of the omitted rows.
The physical scale on the upper polynomial is
\[
    \frac{G_N^C}{G_n^C}\,
        \Lambda^{h+\Delta_{\rm tail}}\,
                 \exp\left(\sum_{\rm tail}E_i\right).
\]
For the second question enumerate the integer
\(g=h+\Delta_{\rm tail}\) over \([-H_n-L,H_n+L]\).
For the third question it suffices to use \(|g|\le L-L_C\):
indeed
\[
 |g|\le |h|+L_{\rm tail}\le L-L_C,
 \qquad g+\Delta_C=h+\Delta .
\]
If a signed target is specified, also require
\(g+\Delta_C=\Delta_0\).
For the fourth question take \(g=0\), since each discarded
height contributes zero signed gap. In each case test the
upper polynomial at
\((G_N^C/G_n^C)\Lambda^g v z\), where
\(v\in[1-2c_B,1+4c_B]\). This interval contains the omitted
exponential factor by \(|\sum_{\rm tail}E_i|\le2c_B\) and
the elementary exponential bounds for \(2c_B<1/2\).

In all these tests the real and imaginary parts of the root
equations are polynomial in the two root coordinates and in
the finitely many real error and scale variables. All interval
bounds are rational and all fixed coefficients are algebraic.
The existence tests are decidable by real quantifier elimination.
[Perrucci and Roy, Theorem 1](https://arxiv.org/pdf/1609.02879)
If all boxes at a cutoff are infeasible, the corresponding
answer is positive. A rational separation margin on these boxes
can also be computed. Use two root variables and replace the
exact contact equation by a squared distance bound between the
unit root projections, or between the two roots after the
prescribed scaling. Decrease the trial margin through dyadic
values. Compactness and the infeasibility of an exact contact
guarantee termination.

### The bounded-gap physical limit

The third question describes exactly the possible positive-curvature
limits of physical coefficient triples of the bounded-gap family.
To see necessity, any convergent subsequence has a finite ordinary
subprofile of constant rows, and the other rows escape in height
with bounded gaps. Their curvature and normalized coefficient
contributions vanish. Their signed total gap tends along a
subsequence to an integer \(h\) of modulus at most the remaining
gap budget, giving the affordable twist of the subprofile.
Conversely, for an affordable twist \(h\ne0\) append to \(C\) one
shared row of signed gap \(h\) and send its smaller exponent to
infinity. The row's curvature tends to zero, and its normalized
log-coefficients tend to zero, while it contributes \(h\) to the
signed total gap. The relative scale tends to
\((G_N^C/G_n^C)\Lambda^h\) by the estimate on \(E_i\).
The gap sum of each approximating profile is \(L_C+|h|\le L\).
For \(h=0\) use the core itself. The same realization preserves
the optional fixed signed target.

Thus if all nonzero affordable twists are excluded, the bounded-gap
original physical property itself has a two-sided decision by the
same procedure. In that case a limiting contact can only have
\(h=0\), which furnishes a genuine original witness.

## Termination and a finite height bound

Each exact profile of bounded height can also be checked by the
same root equations with zero errors. Enumerate these profiles
by increasing height. If any admissible exact profile has a
contact, this finds it, using the finite exact range of twists
above in the second question and the budget directly in the third.

Conversely suppose that the property in question is true and
there are feasible boxes along cutoffs \(B\to\infty\).
Take a subsequence with the same terminal and the same number
of retained rows. Each retained row can in turn be made constant
or divergent in height along a subsequence. Each of the finitely
many row gaps and orientations can also be made constant.
Every divergent row has curvature tending to zero, since its
gap is bounded by \(L\). Its normalized coefficient factor
tends to one by the chord estimate. The constant rows have
curvature at least \(s_0\) in the limit, by the filter
\(S_C\ge s_0-c_B\), so they form a nonempty ordinary subprofile
of gap sum at most \(L\). The box errors tend to zero.
The roots and, in the first question, their positive argument
scale stay in the fixed annuli. Thus the first question's
feasible boxes give an actual angular contact on this
subprofile.

For the second question make the chosen integer \(g\) constant
along the subsequence. The diverging rows' combined center
correction tends to zero, but their signed gap is some integer
\(g_{\rm esc}\). The scale in the limit is the constant
subprofile's scale ratio times \(\Lambda^{g+g_{\rm esc}}\).
A limit of feasible boxes is therefore an exact integer-twist
contact in the same bounded-gap family.
For the third question the limiting twist still obeys the budget:
if \(C_*\) is the constant subprofile, then
\(|g+g_{\rm esc}|\le L-L_C+(L_C-L_{C_*})=L-L_{C_*}\).
Any fixed target \(\Delta_0\) is also preserved, since the constant
signed gap is \(\Delta_C-g_{\rm esc}\).
For the fourth question the constant rows contain all retained rows at every fixed
height eventually, and each such height is balanced. Both
the constant and the diverging rows consequently have zero
signed gap. The limit is an original physical contact on an
admissible constant subprofile. Each outcome contradicts the
assumed property, so some finite cutoff has only infeasible
boxes.

Proceed by increasing height, performing the exact witness tests
and the outer box tests at each cutoff. The argument shows
that one branch terminates. In the positive case its cutoff
certifies every profile, including those of arbitrarily large
height. In the negative case the returned height contains a
finite counterexample. Thus in either case the algorithm
computes a finite height at which checking all admissible
exact profiles is equivalent to answering that question.
The algorithm takes a certified cutoff \(J\) from the cited
absolute theorem. No numerical \(J\) or height bound is evaluated
here. Physical coprimality transfers
to any characteristic-zero field by rational Bézout identities.
