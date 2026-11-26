# Four unit-gap rows suffice for high-prime exterior matching

Fix \(N=n+1\ge3\).
For every finite set of primes \(p>nN\), every choice of positive
precisions \(L_p\), and every lower bound on the exponents, there
is an ordinary mirrored profile of four unit-gap rows above the
lower bound for which, at each selected prime, the following hold.

1. The normalized polynomials at both adjacent orders have
   Newton slopes \(1,-1\) of width one. The horizontal width
   is \(j-2\) at order \(j\).
2. The valuation-one roots can be written \(p w_n,p w_N\) with
   \(w_n,w_N\in\mathbb Z_p^\times\) and
   \(w_n\equiv w_N\pmod{p^{L_p}}\). The normalized reciprocal
   units also agree modulo \(p^{L_p}\).

The profile can therefore be chosen in the
[all-order weak-curvature range](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md),
so these exterior roots can agree to prescribed precision
at many primes while the ordinary moments are coprime at every
order. This improves the dimension and prime threshold in the
[general gap exterior construction](MULTIVARIATE_SIX_ROW_PADIC_EXTERIOR_OBSTRUCTION.md).
For a single unit-gap mirror pair and \(p>nN\), the same two
exterior faces cannot match even at the first residue.

## Two mirrored unit gaps

Use smaller bases \(A,B\) for two pairs of opposite rows:
\((A,A+1),(A+1,A)\) and \((B,B+1),(B+1,B)\).
The scalar-free moment at order \(j\) is
\[
 F_j(z)=\sum_{k=0}^j\binom jk
       (jA+1)_k(jA+1)_{j-k}(jB+1)_k(jB+1)_{j-k}z^k.
\]
At each selected prime set \(A\equiv p-1\pmod{p^2}\);
these classes combine by CRT. The only factor \(jA+h\)
divisible by \(p\) for \(1\le h\le j\) is the one at
\(h=j\). Its valuation is one and it occurs only at the
endpoints. Choose the second base initially modulo \(p\) by
\[
                                B=-1/(nN).
\]
All its sampled factors are units there: for \(j=n\) they are
\(h-1/N\), with \(1\le Nh-1\le nN-1<p\), and for \(j=N\)
the roles of \(n,N\) interchange.

After substituting \(z=pw\) and dividing by \(p\), the residue
in \(w\) is affine with unit coefficients. Omitting the pair
at \(B\), the ratio of constant to linear coefficients is
\[
            \frac{j(A+1)}{p\,j(jA+1)}\equiv\frac1{1-j}\pmod p.
\]
Adding the pair at \(B\) multiplies this ratio by
\((jB+j)/(jB+1)\). Thus its residue root is
\[
           r_j(B)=\frac{jB+j}{(j-1)(jB+1)} .
\]
At the indicated \(B=-1/(nN)\), both residue roots are the same
unit \((nN-1)/(n(n-1))\). Their quotient as a function of \(B\)
is
\[
             \frac{n}{n-1}\cdot
                  \frac nN\frac{NB+1}{nB+1},
\]
whose derivative is a unit at that residue. Exactly as in the
[six-row lift](MULTIVARIATE_SIX_ROW_PADIC_EXTERIOR_OBSTRUCTION.md#matching-both-exterior-roots),
the two equations \(p^{-1}F_j(pw)=0\) in \((B,w)\) have an
invertible Jacobian modulo \(p\). They admit a common lift
\((B_p,w_p)\) in \(\mathbb Z_p^2\) by successive linear digit
lifts. An ordinary \(B\) approximating all the \(B_p\) modulo
their prescribed \(p^{L_p}\) exists by CRT and can be increased
freely by its modulus. The two simple residue roots then lift
to units congruent modulo \(p^{L_p}\). Reciprocity gives the
asserted reciprocal match.

## Why a single unit mirror does not give this match

For one unit-gap mirror pair, a nonhorizontal face can occur at
a prime \(p>nN\) only if some \(jA+h_j\) is divisible by \(p\)
with \(1\le h_j\le j\). Otherwise all coefficients are units
because \(p>j\). If this occurs at both orders, we have
\(Nh_n=nh_N\pmod p\). The integers on the two sides have
difference in absolute value less than \(nN<p\). Hence they
are equal as integers and \((h_n,h_N)=(n,N)\). Thus
\(A+1=p^e v\) for some \(e\ge1\) and \(p\)-adic unit \(v\).
All interior coefficients are units, and both endpoints have
valuation \(e\). Replacing \(z=p^e w\), the two exterior root
residues are \(v/(n-1)\) and \(v/(N-1)\), respectively. These
are distinct modulo \(p>nN\). The reciprocal faces are distinct
as well.
