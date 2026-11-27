# Spectral monodromy at the regular finite pole of the degree-forty-two residual

**Theorem PolynomialWeylFortyTwoRegularSpectralMonodromy.** Take
the centered actual Fourier pair of the
[finite-pole theorem](POLYNOMIAL_WEYL_FORTY_TWO_FINITE_POLE.md)
in its regular singular case \(k=9\). For every spectral value
\(\lambda\in\mathbb C\) the following hold.

- The conjugacy class of the monodromy of \(P-\lambda\) at zero is
  constant in \(\lambda\), and each Jordan block size for each
  eigenvalue occurs with even multiplicity.
- The corresponding monodromy of \(Q-\lambda\) is constant in
  \(\lambda\), and each Jordan block size for each eigenvalue
  occurs with multiplicity divisible by seven.
- The entire solution spaces of \(P-\lambda\) and \(Q-\lambda\)
  have dimensions exactly \(18\) and \(63\), respectively. Each
  entire solution space is the fixed space of its monodromy
  acting on solution germs away from zero.

Here monodromy is analytic continuation of solutions once
counterclockwise around zero. In particular, every noninteger
indicial class consisting of a single pair of simple roots of
\(P\) has scalar monodromy on its two-dimensional generalized
eigenspace.

## A flat spectral equation

For \(P-\lambda\), take the usual jet column of a scalar
solution \(y\), so the companion matrix gives
\(\partial_X v=A(X,\lambda)v\) on \(\mathbb C^\times\).
Since \([Q,P]=1\),
\[
             (P-\lambda)(\partial_\lambda+Q)y=0
\]
along any analytic family of solutions of \((P-\lambda)y=0\).
On the jet column, \(Q\) induces a matrix \(B_Q\) by
reducing the first \(30\) derivatives of \(Qy\) by the
equation. The entries of \(B_Q\) are rational in \(X\),
holomorphic for \(X\ne0\), and polynomial in \(\lambda\).
One can obtain them by right division by the monic
operator \(p_{30}(X)^{-1}(P-\lambda)\), where
\(p_{30}(X)=\alpha X^{12}\).
The companion equations together with
\(\partial_\lambda v=-B_Q v\) are flat. Indeed
\(\partial_\lambda v+B_Q v\) solves the companion
equation for every analytic family \(v\). Applying
\(\partial_X-A\) to this expression gives the matrix
\(\partial_\lambda A+\partial_X B_Q+B_QA-AB_Q\)
applied to \(v\). These families span the fiber at
each ordinary point, so the matrix vanishes.
For \(Q-\lambda\) the same construction uses
\(\partial_\lambda-P\), since \([P,Q]=-1\).

The spectral transport at any fixed ordinary point
is analytic and invertible along every spectral path.
Flatness intertwines it with analytic continuation
around zero, so the monodromy conjugacy classes are
constant. This uses the exact pair and its single
finite pole.

## Scaling at large spectral value

Let
\[
 P_0=X^{-18}f(XY),\qquad Q_0=X^{-63}g(XY)
\]
be the lowest Euler components from the finite-pole
theorem. They commute and
\[
               Q_0^2=\kappa P_0^7,\qquad
               \kappa=\beta^2/\alpha^7\ne0.                  \tag{1}
\]
In fact both products have Euler grade \(-126\) and
commute with \(P_0\). Their coefficient polynomials in
\(C=XY\) satisfy the same equation
\(f(C-126)H(C)=f(C)H(C-18)\).
The quotient of the two is rational of period \(18\),
hence constant. Their leading differential coefficients
give the indicated constant.

For the \(P\)-spectral equation set \(\lambda=s^{18}\)
and \(X=z/s\). Put
\[
   P^{(s)}=s^{-18}P(z/s,sD_z),\qquad
   Q^{(s)}=s^{-63}Q(z/s,sD_z).
\]
Their coefficients are polynomial in \(s^{-1}\).
Their limits at \(s^{-1}=0\) are \(P_0,Q_0\) in the
coordinate \(z\). The leading \(D_z\)-coefficient of
\(P^{(s)}-1\) is \(\alpha z^{12}\) for every \(s\),
so the scaled solution spaces form a rank-thirty
analytic companion system in \(s^{-1}\) at each
\(z\ne0\). Let \(M(s)\) be its \(z\)-monodromy
around zero, based at \(z=1\). It is analytic in
\(s^{-1}\) at zero by analytic dependence of the
companion equation along a fixed \(z\)-loop.

On a horizontal family the \(s\) derivative at fixed
\(z\) is given by
\[
        \partial_s y=-18s^{80}Q^{(s)}y-\frac{z}{s}D_zy .
\]
Thus at the fixed \(z\)-basepoint the spectral
connection matrix \(B(s)\) has expansion
\[
             B(s)=s^{80}(B_0+O(s^{-1})),
        \qquad B_0=-18[Q_0]_{P_0-1}.                      \tag{2}
\]
The brackets denote the induced action on the ordinary
solution germs of \(P_0-1\). This action is an endomorphism
since \([P_0,Q_0]=0\). The expansion follows by reducing
the differentiated jets in the scaled companion equation.
Its leading coefficient is independent of \(s\), so
the reduction introduces no pole in \(s^{-1}\) apart
from the displayed power. By (1), \(B_0\) is diagonalizable,
with possible eigenvalues the two roots of \(18^2\kappa\).
Flatness gives \(M'(s)=[B(s),M(s)]\).

For the \(Q\)-spectral equation use \(\lambda=s^{63}\)
and the same \(X=z/s\), \(P^{(s)},Q^{(s)}\). The
leading \(D_z\)-coefficient of \(Q^{(s)}-1\) is
\(\beta z^{42}\), again independent of \(s\).
Now \(\partial_\lambda=P\) is spectral transport, so
the \(s\)-connection has leading coefficient
\(s^{80}B_0\) with
\(B_0=63[P_0]_{Q_0-1}\). This is diagonalizable
with possible eigenvalues \(63\) times the seven
roots of \(\kappa^{-1}\), by (1).

## Formal blocks and the spectral deck transformation

For either connection, take a constant basis that
diagonalizes \(B_0\). Over \(\mathbb C((s^{-1}))\)
the spectral equation is equivalent to a system
block diagonal by the distinct leading eigenvalues.
For completeness, use a gauge \(I+O(s^{-1})\)
and solve for its off-diagonal coefficients in
successive powers of \(s^{-1}\).
At each order the commutator with \(B_0\)
is invertible on those coefficients because the
two block scalars differ.
The derivative of the current gauge coefficient
enters only at a later order, so the recursion
constructs the formal gauge.

In the resulting frame, a horizontal endomorphism
has zero off-diagonal blocks. Indeed a block
\(Z\) between two distinct leading scalars satisfies
\(Z'=B_i Z-Z B_j\). If \(Z\ne0\), the right side
has leading power of \(s\) exactly \(80\) greater
than that of \(Z\), and the derivative has smaller
power than \(Z\). The same argument applies to
horizontal intertwiners between two such
connections. In particular the monodromy \(M(s)\)
preserves the formal blocks.

Write \(a=18\) or \(63\), respectively, for the
spectral exponent \(\lambda=s^a\). Let \(\zeta\)
be a primitive \(a\)-th root of unity.
The two scaled basepoints at \(s,\zeta s\) refer
to the same spectral value. Their jets are
identified by the \(z\)-transport from \(1\) to
\(\zeta^{-1}\) for the equation at \(s\), followed
by the diagonal derivative rescaling
\(\operatorname{diag}(1,\zeta^{-1},\ldots)\).
Indeed the scaled coordinate at \(\zeta s\) is
\(\zeta z\). This gives a matrix \(R(s)\), analytic
and invertible near \(s^{-1}=0\), that intertwines
the connections and monodromy:
\[
 R'=\zeta B(\zeta s)R-R B(s),\qquad
 M(\zeta s)=R M(s)R^{-1}.                                   \tag{3}
\]
The leading spectrum of the pulled-back connection
is multiplied by \(\zeta^{81}\).
For \(a=18\) this is \(-1\); for \(a=63\) it
is a primitive seventh root of unity. Thus the
deck identification permutes the possible leading
eigenvalues in a single cycle of length two or
seven. The intertwiner in (3) carries each formal
block to its corresponding permuted block, by the
same off-diagonal valuation argument.
In particular every possible block occurs and the
blocks have matching monodromy Jordan partitions,
since formal pullback by \(s\mapsto\zeta s\)
preserves those partitions.
The characteristic polynomial of total monodromy
is constant. The ranks of its constant-eigenvalue
Jordan tests over the formal field agree with the
generic analytic ranks, because \(M(s)\) is analytic
in \(s^{-1}\). They then agree with the ranks
at each finite spectral value by the flat transport.
This proves the divisibility assertions.

## Fixed spaces and the homogeneous limit

In the \(P\) limit the equation is
\(f(C_z)y=z^{18}y\) with \(C_z=zD_z\). For \(Q\)
the equation is \(g(C_z)y=z^{63}y\).
Write \(H=f\) or \(g\) and \(a=18\) or \(63\).
At zero these are regular singular equations.
Any single-valued local solution is meromorphic at
zero. This follows, for example, by writing the
Euler jet equation with bounded holomorphic
coefficients for \(z\partial_z\).
Radial integration bounds the norm by a fixed
power of \(1/|z|\), uniformly over the circle
for a single-valued solution. Its Laurent series
therefore has a finite principal part.

In a Laurent series \(\sum u_n z^n\) the limit
equation has \(H(n)u_n=u_{n-a}\).
Split by residues of \(n\bmod a\). Each residue
has at least one root of \(H(n)\) because of the
factor \((C)_a\). A nonzero series on the residue
must start at the largest integer root there.
If it starts below another root in that residue,
the recurrence propagates nonzero coefficients
up to the next root, which would require the
previous coefficient to vanish. After the largest
root the recurrence instead determines the series
uniquely from its first coefficient. It is
convergent, since \(\lvert H(n)\rvert\) grows
polynomially on the nonzero tail. Thus the
homogeneous limit has exactly \(a\) independent
fixed solutions, one for each residue.

Since the scaled monodromy \(M(s)\) is analytic
at \(s^{-1}=0\), its fixed dimension at a
generic finite \(s\) is at most the limit
fixed dimension \(a\). The finite-pole theorem
constructs at least \(a\) entire solutions
at each spectral value in this regular case.
Those are fixed by the monodromy, so equality
holds throughout. Their germs exhaust the fixed
space, proving the remaining assertions.

**Corollary PolynomialWeylFortyTwoIsolatedFiniteResonance.** Let
\(\rho,\rho+9m\) be the only two indicial roots of \(f\) in
\(\rho+\mathbb Z\), each simple, where \(\rho\notin\mathbb Z\) and
\(m\in\{1,3,5,7\}\). The polynomial Frobenius resonance
obstruction for the lower exponent vanishes identically in
the spectral value. The following calculation gives a
finite test on the coefficients of the actual pair. The
[normalized relation](POLYNOMIAL_WEYL_FORTY_TWO_NORMALIZED_CUSP_RELATION.md)
later forces at least six roots in an occupied noninteger
class, excluding the isolated pair in the degree-forty-two residual.

Write
\[
     P=X^{-18}\sum_{d\ge0}X^d f_d(C),\qquad
     Q=X^{-63}\sum_{d\ge0}X^d g_d(C),\qquad f_0=f,\ g_0=g.
\]
For the lower exponent write \(X^\rho\sum a_nX^n\)
with \(a_0=1\). For \(1\le n<9m\) the recurrence is
\[
 f(\rho+n)a_n=\lambda a_{n-18}
       -\sum_{d=1}^n f_d(\rho+n-d)a_{n-d},                  \tag{4}
\]
where \(a_j=0\) for negative indices.
The obstruction at \(n=9m\) is the right side of (4)
at that index. It is a polynomial of degree at most
\((m-1)/2\) in \(\lambda\), since only a loss of \(18\)
indices can contribute a factor \(\lambda\).

To see directly that it vanishes, start from the higher
Frobenius solution \(\phi_\lambda\) with exponent
\(h=\rho+9m\) and leading coefficient one.
Its recurrence has no subsequent zero indicial
coefficient, and all its coefficients \(b_n\)
are polynomials of degree at most \(\lfloor n/18\rfloor\).
At index \(18j\) the leading coefficient in \(\lambda\)
is \(\prod_{i=1}^j f(h+18i)^{-1}\).
The coefficient at exponent \(\rho\) in
\((\partial_\lambda+Q)\phi_\lambda\) has degree
\((7-m)/2\), with leading coefficient
\[
 \frac{g(\rho+63)}
      {\displaystyle\prod_{i=1}^{(7-m)/2} f(h+18i)} .      \tag{5}
\]
Indeed \(\partial_\lambda\phi_\lambda\) contributes
nothing at \(\rho\). For the \(Q\) term the sum of
the indices from \(g_d\) and \(b_n\) must be
\(63-9m\), so the largest degree in \(\lambda\)
can come only from \(d=0\).
In the pair formula of the finite-pole theorem,
the pair \(\rho,\rho+9m\) contributes no root of
\(g\) at index seven in its \(9\)-progression.
No other excess pair contributes there because
\(\rho+\mathbb Z\) is isolated, and the falling
factorial part of \(g\) has integer roots only.
Thus \(g(\rho+63)\ne0\).

Since \((\partial_\lambda+Q)\phi_\lambda\) is an
actual Frobenius series solution of \(P-\lambda\),
it has no exponent below \(\rho\): its lowest
exponent, if nonzero, must be in the same class
and be an indicial root there. Formula (5) gives
a nonzero coefficient at \(\rho\) for generic
\(\lambda\), so the obstruction polynomial in
(4) vanishes at those spectral values.
It therefore vanishes identically.

## Occupied indicial classes and Jordan partitions

**Corollary PolynomialWeylFortyTwoIndicialJordanDomination.**
In an indicial class \(\rho+\mathbb Z\), count the pairs of
excess \(f\)-roots in each of its \(9\mathbb Z\)-suborbits,
with multiplicity. Add one to each of the nine counts in
the integer class to include the base falling factorial.
Write \(r=(r_1,\ldots,r_d)\) for the resulting positive
counts in decreasing order, omitting zero counts. At
\(\chi=e^{2\pi i\rho}\), the homogeneous model \(P_0-1\)
has Jordan block sizes \(r_i\), each repeated twice.
The model \(Q_0-1\) has sizes \(r_i\), each repeated seven
times. The corresponding actual partitions for any
spectral value, divided in multiplicity by two and seven
respectively, both dominate \(r\).

In particular, if a noninteger class occupies only one
\(9\mathbb Z\)-suborbit containing \(r_1\) pairs, the
actual \(P\)-monodromy has exactly two blocks of size
\(r_1\) and the actual \(Q\)-monodromy has exactly seven
blocks of that size at \(\chi\). In the integer class
both normalized partitions have exactly nine parts.
If all \(e\) excess integer pairs occupy one
\(9\mathbb Z\)-suborbit, both normalized partitions are
exactly \((1+e,1,1,1,1,1,1,1,1)\).

To compute the model Jordan sizes, use the coordinate
\(w=z^a\), with \(a=18\) for \(P\) and \(a=63\) for \(Q\).
The model equation is the pullback of
\[
                     H(a wD_w)\widetilde y=w\widetilde y,
\]
where \(H=f\) or \(g\). Its indicial roots in \(w\)
are the roots of \(H\) divided by \(a\).
For a regular singular Euler companion equation the
characteristic polynomial of monodromy has these indicial
roots exponentiated by \(2\pi i\). Indeed, the monodromy
computed along a circle of radius \(\epsilon\) tends to
the exponential of the residue matrix as \(\epsilon\to0\),
and the characteristic polynomial is unchanged under
basepoint transport. For each indicial class modulo
integer shifts in \(w\), there is only one eigenvector
solution up to scale. After factoring out the chosen
complex power of \(w\), the solution is single-valued
and has finite principal part by the regular singular
growth bound. Its one-step recursion, as in the fixed
space calculation above, permits only the largest
indicial exponent in that class as a starting term.
The geometric multiplicity is therefore one, so its
Jordan block size is the number of roots in that class,
with multiplicity. Pullback takes the \(a\)-th power
of the monodromy. This preserves these block sizes
because the eigenvalues are nonzero.

Each excess pair of \(f\) has odd index difference in
its \(9\)-progression and therefore adds one root to
each of its two classes modulo \(18\). The associated
\(g\) record from the finite-pole theorem adds one root
to each of the seven residue classes modulo \(63\)
on that progression. To see the latter, for a matched
pair its record \(G(z)\) has \(G(1)=7\) and vanishes
at every nontrivial seventh root of unity by (4) of
the finite-pole theorem. Thus its index counts modulo
seven are all one. The base falling factorials give
the same occupancy in each of the nine integer
progressions. This proves the model partitions.

For an analytic family of monodromy matrices, the
rank of \((M-\chi I)^j\) at a generic parameter is
at least its rank in the homogeneous limit. The
characteristic polynomials are the same in the limit
and at all finite parameters, by the analytic
specialization and spectral conjugacy already proved.
After dividing the Jordan repetitions by their
respective factors two and seven, the rank inequalities
on the generalized \(\chi\)-eigenspace say that an actual
partition \(q\) satisfies
\[
              \sum_i\max(q_i-j,0)\ \ge\
                 \sum_i\max(r_i-j,0)
              \qquad(j=0,1,2,\ldots).                      \tag{6}
\]
These are the partition dominance inequalities.
If \(r\) has only one part, (6) at \(j=1\) and
equality of total size give at most one part in \(q\),
so \(q=r\). In the integer class, the exact
fixed dimensions \(18,63\) give nine parts after the
division by repetitions. With all excess integer
pairs in one progression, (6) at \(j=e>0\) forces
a part at least \(1+e\), and the other eight positive
parts force the remaining sizes to be one. If \(e=0\),
all nine parts are already forced by the total size.
