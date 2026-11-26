# Curves and directions in the actual least-repaired plane

For a genuine primitive row \(a+b=c\), use
\[
        Z=Z_M(c),\qquad U=aZ,\qquad V=bZ,\qquad H=cZ .
\]
Put \(n=X^L\), \(M=\lceil K/L\rceil\), with the height limit
taken at fixed \(K,L\). Unless specified otherwise we work
in \(X/2<c\le X\) and retain the original individual caps
and total radical cap. The primitive row determines its
unique least multiplier. Conversely the repaired point
\((U,V)\) recovers the original row, since \(\gcd(U,V)=Z\).
[Actual least-target repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

## Uniform bounded-degree curves in the repaired plane

**LP1.** Fix \(D\ge1\). For any nonzero rational polynomial
\(F(U,V)\) of total degree at most \(D\) with coefficients of
uniform polynomial height in \(X\), the number of genuine
repaired points on \(F=0\), apart from its line components
\(U+V=h\), is \(O(X^{D/M+\varepsilon})\).
If the lowest nonzero homogeneous part of \(F\) has nonzero
value at \((1,-1)\), the exponent improves to \(1/M\).

Write \(F=\sum_d F_d\) as a sum of its nonzero homogeneous
parts, clear denominators, and write each \(F_d=(U+V)^{m_d}G_d\)
with \(G_d(1,-1)\ne0\) and \(0\le m_d\le d\). Use the finite
prescribed pool supporting the nonzero \(G_d(1,-1)\). Its
product has polynomial height. For any deficient target prime
\(p\mid c\) outside this pool, put \(e=v_p(c)<M\).
Primitivity makes \(a\) a unit at \(p\) and \(b\equiv-a\bmod p\).
Consequently \(G_d(a,b)\) is a unit for each \(d\), including
the constant degree when present. The exact valuation of the
\(d\)-term at the repaired point is therefore
\[
                          d(M-e)+m_d e.
\]
If there is a vanishing sum, at least two degrees attain its
minimum at this prime. For a tie between \(d\) and \(d+r\)
with \(r>0\), we have
\[
                   r M=(r+m_d-m_{d+r})e.
\]
Thus \(e/M\ge r/(r+d)\ge 1/D\), since \(d+r\le D\) and
\(m_d\le d\). If there is just one nonzero homogeneous part,
no such deficient prime outside the pool can occur.
When the lowest part has \(m_d=0\), every positive change
of degree strictly raises its valuation at a deficient prime,
and again there is no deficient prime outside the pool.

The nondeficient primes have \(e\ge M\). This bounds the
outside radical of \(c\) by \(X^{D/M}\), or by \(X^{1/M}\)
under the stronger hypothesis. Count its possible exact
values by fixed-support enumeration. For each \(c\), \(Z\)
is fixed and points are on \(U+V=cZ\). This line meets
\(F=0\) in at most \(D\) points unless it is a component,
by restricting the polynomial to the line. The excluded
components have prescribed repaired targets. Each of their
original targets lies among the divisors of \(h\), and their
primitive row fibres are controlled by FL1 using both
original summand caps.

This count is sharp for the integer repair map before
imposing the other original radical caps. When \(D\mid M\),
take \(c=R^{M/D}\) for squarefree \(R\) and \(a=1\).
Then \(Z=c^{D-1}\) and
\(F(U,V)=(U+V)^{D-1}-U^D=0\). With \(X/2<c\le X\)
this realizes \(X^{D/M-o(1)}\) points. Thus the branch
\(e=M/D\) is attainable. For degree one, each line other
than \(U+V=h\) has fibre at most \(n^{1/K}X^\varepsilon\).
[Original fixed-target count (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

## Two normalized displacement records for arbitrary pairs

**LP2.** For two distinct original rows put \(x_i=c_i/a_i\).
When \(Z_1\ne Z_2\), define the exact rational records
\[
             \delta=x_1-x_2,\qquad
             k=\frac{x_1 Z_1-x_2 Z_2}{Z_1-Z_2}.
\]
Then each exact \((\delta,k)\) record has at most
\(O(n^{3/K}X^\varepsilon)\) original ordered pairs.
If \(k\ne\delta,-\delta\), the exponent is at most \(2/K\).

The original odds determine the primitive row, so
\(\delta\ne0\) for distinct rows. The two records give
\((x_1-k)Z_1=(x_2-k)Z_2\). At a deficient prime \(p\mid c_1\)
outside the prescribed pool of primes of the nonzero
numerators and denominators in \(\delta,k\), \(x_1\equiv0\), \(x_2\equiv-\delta\) is a unit,
and hence \(p\nmid c_2\). Thus \(Z_2\) is a unit and the
equation forces \(p\mid\delta+k\). Symmetrically a deficient
prime on \(c_2\) outside that pool must divide \(\delta-k\).
If neither quantity is zero, prescribing their rational
supports puts all deficient primes into fixed polynomial-height
pools. Enumerating nondeficient supports costs \(n^{1/K}\)
per target. After \(c_i,Z_i\) are fixed, recover the two
odds from \(\delta=x_1-x_2\) and
\(k=x_2+\delta Z_1/(Z_1-Z_2)\). This fixes \(a_i=c_i/x_i\)
and the original row.

If \(k=-\delta\), first recover \(c_2\) by the nonzero
condition \(\delta-k=2\delta\). The equality becomes
\((2x_1-x_2)Z_1=x_1 Z_2\). At a prime of \(c_1\) outside
the fixed support of \(c_2,\delta\), both \(x_2\) and
\(2x_1-x_2\) are units. Thus \(v_p(Z_1)=v_p(c_1)\);
the only possible positive original valuation is \(M/2\).
This costs \(n^{2/K}\) on \(c_1\). The other case follows
by exchanging the rows, and the same recovery applies.
There are no exceptional primes from the variable summand
denominators at the unit tests, since \(\gcd(a_i,c_i)=1\)
and \(x_i\) is required to have the stated valuation.

The exponent cost is uniform in every actual polynomial-height
record. The data can be written as the difference of the
primitive target odds, together with the slope of the
repaired displacement after dividing each original row by
its own \(a_i\). No shared original \(a\) or exact \(a_1/a_2\)
is needed. In a near-full balanced family, take \(K\) sufficiently
large and then \(L\) sufficiently small and fixed. Equal
\(Z\) pairs are power-small by LM1 target recovery and
FL1, as in LG1.
The joint exact \((\delta,k)\) dictionary therefore has
exponent at least \(2/3-3/K-o(1)\), across the entire
balanced dictionary range. On a selected pair sector one
can instead compare with its actual retained exponent.

## A large outside radical on primitive repaired directions

**LP3.** Suppose \(N\ge n^{1/3-\tau}\) original rows have
balanced individual caps and target masses at least
\(1/3-\tau\). Fix any polynomial-height prime pool \(S\).
For every \(\rho>0\), at sufficiently large fixed \(K\)
and sufficiently small fixed \(\tau,L\), on all but a
power-small fraction of pairs the reduced integer direction
\[
              (f,g)=\frac{(U_1-U_2,V_1-V_2)}
                               {\gcd(U_1-U_2,V_1-V_2)}
\]
has \(\operatorname{rad}_{S^c}(|fg(f+g)|)>n^{1-\rho}\).
The same floor holds after also omitting both original
target supports, by taking \(K\) larger and the other
allowances smaller.

Here is the numerical-height input. Write \(F=U_1-U_2,\ G=V_1-V_2\),
\(h=\gcd(F,G)\), \(z=\gcd(Z_1,Z_2)\), and
\(D_0=a_1b_2-a_2b_1\). The two reduced multipliers
\(Z_i/z\) are coprime, and \(h/z\) is prime to their
product. A prime shared with one of them would divide
both original summands of the other row. Eliminating
one multiplier therefore shows \(h/z\mid D_0\).
For distinct primitive rows \(D_0\ne0\), while
\[
               Fb_2-Ga_2=Z_1D_0.
\]
Thus \(\max(|F|,|G|)\ge Z_1/c_2\). Outside a power-small
fraction EB5 bounds the shared original target radical
by \(n^\sigma\) for any small fixed \(\sigma>0\).
Then \(z\le n^{M\sigma}\), and the original lower target
cap gives \(Z_1=H_M(c_1)/c_1\ge n^{M(1/3-\tau)}/X\).
Combining these inequalities gives
\[
          \max(|f|,|g|)\ge
            X^{ML(1/3-\tau-\sigma)-4}
\]
up to fixed constants. This is at least a fixed positive
power of \(X\) when \(K\) is large. Directions with
\(F=0,\ G=0\), or \(F+G=0\) are power-small. In the first
two cases the other point is on a fixed line counted by
LP1. In the last case the other repaired target is fixed,
and the original fixed-target count applies.

The triple of absolute values of \(f,g,f+g\), put in
additive order, is a genuine primitive triple of fixed
positive numerical height and polynomial upper height.
Directions with outside cap \(n^\beta\) number at most
\(n^{\beta/3+\eta}X^\varepsilon\) by MCB2, after fixing
the exact supported powers and using that the sum of the
two cheapest outside masses is at most \(2\beta/3\).
Given the direction and the first repaired point, the
second point lies on a fixed non-anti-diagonal line
and has count \(n^{1/K}X^\varepsilon\) by LP1. Thus the
bad pair count is at most
\(O(N n^{\beta/3+1/K+\eta}X^\varepsilon)\) on the tall
subset. Taking \(\beta=1-\rho\) and choosing the allowances
inside the strict margin proves the assertion.

At a deficient prime private to one original target and
absent from all three letters of the other row, all three
of \(F,G,F+G\) are units. In particular the reduced direction
has no prime there. Omitting both original target supports
therefore loses at most the shared original radical mass
plus \(2/K\) from nondeficient primes, outside the same
power-small overlap exception. Start with a smaller
tolerance to get the same final floor.
[Two-cheapest fixed-coefficient bound (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)
· [Original prime support pricing (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
