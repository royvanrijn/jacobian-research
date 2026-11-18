# A local-square obstruction in the balanced reciprocal degree-26 case

A normalized rational reciprocal degree-26 model in the balanced
\((6,6)\) sign class cannot have a rational-square reciprocal parameter.
More strongly, its reciprocal parameter must be a nonsquare in
\(\mathbb Q_{37}\). This follows from a complete finite-field exclusion
and the uniform denominator-control theorem.

This local-square argument is one ingredient of the
[complete rational reciprocal obstruction](DEGREE_TWENTY_SIX_RECIPROCAL_GLOBAL_OBSTRUCTION.md),
which also excludes the nonsquare balanced class by eighteen further
unit-ideal computations. General nonreciprocal pencils remain outside
these reciprocal obstructions.

## Exact labelled strong equations

Use the normalization of the
[mixed coupled-ODE note](MIXED_RECIPROCAL_COUPLED_ODE.md), with
\(h=6,L=25,r=1/\rho,\rho^2=\lambda\). Write

\[
 P=w^3+p_2w^2+p_1w+p_0,\qquad
 Q=w^3+q_2w^2+q_1w+q_0,\qquad s=P Q.
\]

The unique coefficient reconstruction gives \(F,X,Y\) with

\[
 \begin{split}
 X+2(w+2)X'&=25(D_6-rD_5)P Q,\\
 Y+2(w-2)Y'&=25(C_6-rC_5)P Q.
 \end{split}
\]

Every good model satisfies

\[
                  P^2\mid X,\qquad Q^2\mid Y,\qquad
                  G:=F(r)/(2r)=0.                         \tag{1}
\]

The two square divisibilities follow from the simple roots of \(P,Q\),
the vanishing sign factors, and the operator equations; the good locus
excludes the fixed roots \(w=-2,2\). These are necessary equations at
every coefficient-field point. Their polynomial closure, used below,
also includes degenerate points.

Equations (1) have thirteen scalar equations: six coefficients from
each remainder and one equation \(G=0\), in seven variables
\((p_2,p_1,p_0,q_2,q_1,q_0,r)\). The apparent division in \(G\)
cancels polynomially, and \(G|_{r=0}=25\).
Their exact primitive integer forms are generated symbolically and
retained in the
[receipt](degree_twenty_six_balanced_local_square_verification.json).

For fixed \(P,r\), the first six equations are affine-linear in the
three coefficients of the monic \(Q\). This makes a complete residue
census feasible without sampling seven independent coordinates.

## Theorem 1. The labelled strong closure has no \(\mathbb Q_{37}\) point

First, every point of (1) over \(\mathbb Q_{37}\) has integral
\(P,Q,r\). Square divisibility implies the raw divisibility equations
\(P\mid X,Q\mid Y,G=0\), so
[uniform integrality](RECIPROCAL_PARAMETER_FINITENESS_AND_INTEGRALITY.md)
applies at \(37>25\), including every degenerate point.
It also makes \(r\) a unit.

A complete census in \(\mathbb F_{37}\) gives no solution of (1).
Two distinct candidate-generation methods retain all singular branches:

| Method | Outer tuples \((P,r)\) | Candidate monic \(Q\) | Full solutions |
|---|---:|---:|---:|
| Six strong rows, Gaussian elimination | 1823508 | 1925562 | 0 |
| Three weak rows, Cramer rule and enumeration | 1823508 | 12275266 | 0 |

The outer domain is all \(37^3(37-1)\) tuples, with no
discriminant, root-position or good-reduction filter. The omitted
\(r=0\) case is impossible because \(G=25\).

The first method solves all six linear equations for \(P^2\mid X\).
Its numbers of consistent outer tuples at coefficient-matrix ranks
\(0,1,2,3\) are respectively

\[
                       (0,5,4021,1769940).
\]

Every free coefficient is enumerated; thus the candidate count is

\[
                  5\cdot37^2+4021\cdot37+1769940=1925562.
\]

The second method starts with the weaker \(P\mid X\), obtaining
three rows. For a nonzero determinant it uses Cramer's rule.
For a singular matrix it chooses one nonzero coefficient, enumerates
the other two variables, solves that one equation, and tests all rows.
If every coefficient vanishes, it tests the constant rows and
enumerates all three variables when necessary. Every resulting
candidate is then tested against both full square divisibilities
and \(G=0\).

Thus either enumeration covers every possible residue of an integral
strong point. The empty residue solution set proves that (1) has
no \(\mathbb Q_{37}\) point. In particular arbitrary denominators
and bad reduction have been covered.

## Corollary 2. Necessary arithmetic for rational balanced models

The reciprocal parameter of a rational good model lies in
\(\mathbb Q\), since its unique equal-value critical pairs are
Galois invariant. Suppose it were a square in \(\mathbb Q_{37}\).
Choose \(\rho\in\mathbb Q_{37}\). The aggregate plus and minus
factors are then defined over \(\mathbb Q_{37}\), since each is
the monic greatest common divisor of the critical polynomial
with the corresponding sign expression. They would give a point
of (1), contradicting Theorem 1.

Therefore every remaining rational balanced model must have

\[
                       \lambda\notin\mathbb Q_{37}^{\,2}.     \tag{2}
\]

The uniform integrality theorem applies over
\(\mathbb Q_{37}(\sqrt\lambda)\) as well, so \(\lambda\) is a unit
at 37. Since 37 is odd, a unit is square in \(\mathbb Q_{37}\)
exactly when its residue is square. Thus (2) is the precise
quadratic-residue condition

\[
                         \left(\frac{\lambda}{37}\right)=-1.
\]

In particular \(\lambda\) cannot be a rational square. A global
nonsquare whose completion splits at 37 is excluded too; the
remaining nonsquare case is unramified quadratic at 37.

The same uniform theorem restricts the global square class of
\(\lambda\) to primes at most 25 and the sign. There are at most
\(2^{10}=1024\) such signed square classes. The present single
residue condition does not eliminate all of them.

## Exact replay and scope

~~~text
.venv/bin/python research/arithmetic/verify_degree_twenty_six_balanced_local_square.py --write-receipt research/arithmetic/degree_twenty_six_balanced_local_square_verification.json
~~~

Default replay derives all thirteen equations from the actual \(F\),
checks their linear structure and coefficient denominators, compiles
the
[C++ census](verify_degree_twenty_six_balanced_mod37.cpp), and runs
both methods. It also compares complete \(X,Y\) coefficient vectors
and \(G\) against exact symbolic reconstruction at four retained
controls, including \(r=0\) and a degenerate sign factor. The two
constructions are also compared on a complete unisolvent tensor grid:
each cubic coefficient group is zero or one of its three unit vectors,
and \(r=0,\ldots,12\), giving 208 controls. Both constructions are
affine-linear in each cubic group; \(X,Y\) have \(r\)-degree at most
one and \(G\) at most twelve. The checker verifies these bounds
symbolically; the C++ multiplication and operator formulas have the
same bounds directly. Thus these controls bind every coefficient
polynomial over \(\mathbb F_{37}\), rather than just selected examples.
The two
candidate methods share the transparent low-level polynomial
arithmetic and final equation tests; they use different linear
systems and different solvers.

This is an exhaustive finite computation combined with a written
uniform integrality theorem. It is not a good-reduction exclusion,
a numerical search, or a proof that the nonsquare balanced class
has no rational model. No novelty, formal-verification or external
human-review claim is made.
