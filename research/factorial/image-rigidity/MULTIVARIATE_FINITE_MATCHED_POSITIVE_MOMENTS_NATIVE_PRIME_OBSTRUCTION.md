# Finitely many ordinary matched positive moments do not determine native prime chords

Fix adjacent orders \(n,N=n+1\) with \(n\ge2\), and any finite set
of nonnegative integer pairs \((r,s)\). There are two ordinary
profiles with the following properties.

- Both use two private unit rows and only increasing unit-gap
  rows on their shared variables. They have exactly the same
  number of shared rows and the same specified sums
  \(\sum T_n^r T_N^s\).
- Their total shared curvature can be made arbitrarily small.
  The integer row multiplicities may be held fixed while the
  smaller exponents tend to infinity in a common progression.
- At one fixed prime \(p>N\), the degree-\(N\) Newton polygon of
  one profile is a single edge of strictly positive slope.
  The corresponding polygon of the other profile is horizontal.
  Here only scalar factorial factors have been removed.

The same equalities hold for the signed moments, since the shared
orientations coincide, and for complete moments if the two fixed
private rows are included. Thus these finite positive moment data
do not determine the original factorial array's native prime
behavior, even on ordinary profiles of arbitrarily small shared
curvature. The result concerns a particular arithmetic route
to angular separation. It gives no common-root example.

## A two-pole identity with fixed multiplicities

For an increasing unit-gap row of smaller exponent \(a\), the
[exact Gauss grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)
gives
\[
                  T_j=\frac{j-1}{2ja+j+1}.
\]
Let \(R,S\) be the largest exponents \(r,s\) in the specified
list, including \((0,0)\). Choose \(I=S+2\) distinct ordered
integers \(u_i\) and \(K=R+2\) distinct ordered integers \(v_j\),
with \(u_i<v_j\) for every cell. For an integer parameter \(B\)
set
\[
    a_{ij}(B)=\frac{B-u_i}{2nN(v_j-u_i)}-\frac{n+1}{2n}.
                                                               \tag{1}
\]
The integrality and prime conditions will be imposed below.
The two tanh values at this cell are
\[
       T_n=\frac{(n-1)N(v_j-u_i)}{B-u_i},\qquad
       T_N=\frac{(N-1)n(v_j-u_i)}{B-v_j}.                     \tag{2}
\]
In particular, their mixed \((r,s)\) monomial, apart from a
constant independent of the cell, is the rational function
\[
                    \frac{(v_j-u_i)^{r+s}}
                              {(B-u_i)^r(B-v_j)^s}.          \tag{3}
\]

Take signed rational weights
\[
        A_i=\frac{1}{\prod_{i'\ne i}(u_i-u_{i'})},
        \qquad
        C_j=\frac{1}{\prod_{j'\ne j}(v_j-v_{j'})}.
\]
The first list annihilates polynomials in \(u_i\) of degree
at most \(S\) and the second annihilates polynomials in \(v_j\)
of degree at most \(R\), by the leading coefficient formula
in Lagrange interpolation. Both lists have alternating signs.

The sum of (3) weighted by \(A_iC_j\) is zero as a rational
function of \(B\). Indeed at \(B=u_i\) the Laurent coefficient
of \((B-u_i)^{-r+h}\), \(0\le h<r\), is a constant times
\((v_j-u_i)^{r-h}\). Summing in \(j\) therefore kills the
principal part at \(u_i\). At \(B=v_j\) the corresponding
principal part has coefficients polynomial in \(u_i\) of
degree at most \(s\), so the sum in \(i\) kills it. With no
finite poles the sum is a polynomial. For \(r+s>0\) it
tends to zero at infinity and is identically zero. For
\((r,s)=(0,0)\) it vanishes because both lists sum to zero.

Multiply all \(A_iC_j\) by one positive integer to make
them integral. Give the first profile this many rows at
each positive cell, and the second profile the absolute
multiplicity at each negative cell. Equations (2)–(3)
prove exact equality of the specified moment sums and
of the numbers of rows. The weights depend only on the
pole grid and are independent of \(B\).

## Grid residues and the ordinary progression

Put \(L_0=2nN\) and \(L_1=n^2+1\).
Choose a prime \(p>\max\{L_0,L_1,I,2N+5\}\).
All target \(p\)-adic residues below are read modulo
\(p^2\) unless stated otherwise. Put \(\rho=n^2/L_1\)
in that ring. Choose \(a_*\) so that its residue modulo
\(p\) differs from \(0,1,\rho,\rho^{-1}\), and set
\[
   u_+=1,\quad u_-=a_*,\qquad
   v_+=\rho,\quad v_-=\rho a_* .
\]
In addition require that for the two opposite-sign cells
the functions
\[
       \frac{n^2(1-a_*)}{2n(\rho a_*-1)},\qquad
       \frac{n^2(a_*-1)}{2n(\rho-a_*)}
\]
avoid \(0,-1,\ldots,-(N-1)\) modulo \(p\).
These are nonconstant fractional linear functions over
\(\mathbb F_p\), since \(\rho\ne0,1\). Each of the \(N\)
forbidden values excludes at most one choice per function.
Together with the first four exclusions this leaves a
choice by the bound on \(p\).
In particular all four \(v_\pm-u_\pm\) and cross differences
are nonzero modulo \(p\).

Choose increasing integers \(u_i\) with residue \(u_+\)
or \(u_-\) according to the known sign of \(A_i\).
Also take all \(u_i\equiv0\bmod L_0\) and congruent to
one another modulo every prime at most \(I\).
These choices are compatible with the \(p^2\) conditions
because \(p>I\) and \(p\nmid L_0\).

Next choose increasing integers \(v_j>\max_i u_i\),
with residue \(v_+\) or \(v_-\) according to the known
sign of \(C_j\), and \(v_j\equiv1\bmod L_0\).
They can be chosen so that all differences
\[
                              D_{ij}=v_j-u_i
\]
are pairwise coprime. At each step avoid the residues of
all \(u_i\) modulo every prime dividing some nonzero
\(u_i-u_{i'}\) or some already chosen \(D_{i j'}\).
For a prime at most \(I\) the \(u_i\) have one residue,
and for a larger prime they occupy at most \(I\)
residues. There is always an allowed residue. At the
prescribed \(p\) and the primes dividing \(L_0\) the
fixed \(v_j\) residues already avoid all \(u_i\).
The Chinese remainder theorem permits a sufficiently
large \(v_j\) satisfying all constraints. A common prime
of two new differences would divide a difference between
two \(u_i\), so the avoidance proves the assertion.
The differences are also prime to \(pL_0\).

Formula (1) is integral precisely when
\[
                  B\equiv u_i+N^2 D_{ij}\pmod {L_0D_{ij}}.
\]
All these congruences reduce to \(B\equiv N^2\bmod L_0\)
and \(B\equiv u_i\bmod D_{ij}\) in each cell.
The pairwise coprimality makes them compatible, and we
may additionally impose \(B\equiv p\bmod p^2\).
This gives one arithmetic progression. For all sufficiently
large positive \(B\) every \(a_{ij}(B)\) is an ordinary
positive integer. All grow linearly in \(B\) and the fixed
multiplicities make the total shared curvature tend to zero.

## The native prime comparison

By (1),
\[
              N a_{ij}+1=
             \frac{B+n^2u_i-(n^2+1)v_j}{2n(v_j-u_i)}.
\]
When the two signs agree, the numerator is \(p\) modulo
\(p^2\), since \(L_1v_\pm=n^2u_\pm\bmod p^2\).
The denominator is a unit, so \(N a_{ij}+1\) has
valuation one. The other \(N-1\) factors \(N a_{ij}+h\)
are units. For opposite signs, the specified exclusions
ensure that all \(N\) factors are units, by the two
fractional linear computations.

In the first profile every shared row supplies its
valuation-one factor to precisely the nonconstant
coefficients at order \(N\), after its scalar
\((N a_{ij})!\) is removed. In the second profile
all normalized row factors are units. The two private
units cancel the binomial coefficient up to \(N!\),
which is a \(p\)-adic unit. Thus the first coefficient
valuation array is \((0,E,\ldots,E)\), where \(E\) is
the positive number of its shared rows, and the second
array is \((0,\ldots,0)\). Their Newton polygons have
the claimed different endpoint slopes.
