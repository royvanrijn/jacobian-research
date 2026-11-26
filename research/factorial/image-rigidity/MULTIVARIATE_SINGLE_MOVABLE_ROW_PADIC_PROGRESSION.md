# An explicit all-degree angular progression for one movable ordinary row

Fix any ordinary two-monomial core and choose a positive gap
\(d\) for one additional shared row. Either orientation is
allowed. There is an explicit arithmetic progression of positive
smaller exponents \(A\) on which the augmented profile has
disjoint consecutive root arguments in every degree. The
progression requires no resultant evaluation or terminal
root test. In particular this applies when the movable row
is one existing row of a profile and every other row is held
fixed. The added curvature tends to zero along the progression.

Take \(J\ge4\) at least the
[computable absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Let \(a_*\) be the largest retained exponent, or zero for
an empty core. Independently for each \(3\le N<J\) choose a
distinct prime
\[
                     p_N>\max\{J^2d,\ J(a_*+1)\}.
\]
The progression is specified by
\[
                      N A+1\equiv p_N\pmod {p_N^2}
                              \qquad (3\le N<J).              \tag{1}
\]
Since \(N\) is prime to \(p_N\), these congruences determine
one class modulo \(\prod_{3\le N<J}p_N^2\).
The absolute cutoff is unevaluated. The construction gives
arbitrarily large members of the cofinite family of the
[one-row angular theorem](MULTIVARIATE_SINGLE_ROW_ANGULAR_EXCEPTIONAL_FIBERS.md),
with the baseline determined by integer congruences.

## The two factorial arrays at one prime

Exchange the two columns, if needed, to orient the movable
row as \((A,A+d)\). Such an exchange inverts nonzero roots
and preserves equality of arguments. After removing the
row scalar \((mA)!\), put
\[
 F_m(A,z)=\sum_{k=0}^m b_{m,k}(mA+1)_{dk} z^k,\qquad
 (x)_t=x(x+1)\cdots(x+t-1),
\]
where \(b_{m,k}>0\) are the ordinary coefficients of the
retained core, including the binomial coefficient. Fix
\(N=n+1\) with \(3\le N<J\), and put \(p=p_N\).
All \(b_{m,k}\) at \(m=n,N\) are \(p\)-adic units by the
bound on \(a_*\). At \(m=N\), the first rising factor \(NA+1\)
has valuation one by (1). No other factor can be divisible by
\(p\) because the remaining offsets \(h-1\), \(2\le h\le Nd\),
lie strictly between zero and \(p\). Thus every nonconstant
coefficient of \(F_N\) has valuation one and its constant
coefficient is a unit.

At \(m=n\), the factor \(nA+h\) is congruent modulo \(p\)
to \((Nh-n)/N\). For \(1\le h\le nd\), the numerator lies
strictly between zero and \(p\). Every coefficient of
\(F_n\) is a unit. In particular its roots all have valuation
zero: at positive valuation its constant term would dominate,
and at negative valuation its leading term would dominate.

## The upper ratios have an exact fractional displacement

Work over the completion \(K_0\) of the maximal unramified
extension of \(\mathbb Q_p\), with valuation \(v(p)=1\)
and algebraically closed residue field. Put \(\pi=p^{1/N}\).
Writing the nonconstant coefficients of \(F_N\) as \(p C_k\),
with each \(C_k\) a unit, gives
\[
 F_N(A,u/\pi)=b_{N,0}+C_N u^N+\pi C_{N-1}u^{N-1}+O(\pi^2).
                                                                    \tag{2}
\]
The error has coefficient valuations at least \(2/N\).
The binomial \(b_{N,0}+C_N u^N\) has \(N\) unit roots
\(\lambda_i\in K_0\) with distinct residues, by Hensel's
lemma and \(p>N\). In \(K_0(\pi)\), Hensel's lemma applied
to (2) lifts them to all \(N\) simple roots \(u_i\).
Solving the first-order equation at each \(\lambda_i\) gives
\[
             u_i=\lambda_i+c\pi+O(\pi^2),\qquad
                       c=-\frac{C_{N-1}}{N C_N}.
\]
The number \(c\) is a unit and is the same for every root.
For distinct \(i,j\), set \(\omega=\lambda_i/\lambda_j\),
a nontrivial \(N\)-th root of unity in \(K_0\). Then
\[
              \frac{u_i}{u_j}=\omega+
                    \frac{c(1-\omega)}{\lambda_j}\pi+O(\pi^2).
\]
Since \(1-\omega\) is a unit, the distance from the quotient
to \(\omega\) has valuation exactly \(1/N\). The same statement
holds for the corresponding physical roots \(u_i/\pi,u_j/\pi\),
whose quotient does not change on rescaling.

## The lower ratios cannot have that displacement

Let \(\alpha,\beta\) be roots at order \(n\).
Each has degree over \(K_0\) at most \(n<p\).
Write the respective degrees as \(e_\alpha,e_\beta\).
Such extensions of \(K_0\) are totally tamely ramified,
because its residue field is algebraically closed.
An extension of degree \(e\) is generated over \(K_0\)
by \(p^{1/e}\). Indeed, for a uniformizer \(\varpi\),
the unit \(\varpi^e/p\) has an \(e\)-th root in the
extension by Hensel's lemma. Dividing \(\varpi\) by that
root gives an \(e\)-th root of \(p\), whose Eisenstein
polynomial has degree \(e\). The choices of \(p^{1/e}\)
differ by roots of unity already in \(K_0\).

Put \(E=\operatorname{lcm}(e_\alpha,e_\beta)\) and
\(\eta=p^{1/E}\). In the power basis of \(K_0(\eta)\),
the expansion of \(\alpha\) uses only \(\eta^j\) with
\(j\) a multiple of \(E/e_\alpha\), and similarly for
\(\beta\). Since \(\omega\in K_0\), the expansion of
\(\alpha-\omega\beta\) uses the union of these two
index sets. Nonzero terms have distinct valuations
modulo one if their indices \(0\le j<E\) are distinct.
Thus a nonzero difference has valuation in
\[
                    \frac{1}{e_\alpha}\mathbb Z
                       \ \cup\ \frac{1}{e_\beta}\mathbb Z .
\]
Neither lattice contains \(1/N\), since both degrees
are smaller than \(N\). As \(\beta\) is a unit, the
lower quotient \(\alpha/\beta\) cannot have displacement
of valuation \(1/N\) from \(\omega\).
This excludes every coincidence of an upper off-diagonal
root quotient with a lower quotient.

If the two complex root sets share a nonreal argument,
each root and its conjugate are distinct. The quotient
of the root by its conjugate is then the same algebraic
number in both degrees. Embedding the common splitting
field in \(\overline{\mathbb Q}_p\) contradicts the
quotient comparison. A shared real argument is impossible
because the even ordinary moment is strictly positive on
the real line. This also handles the first pair. For
\(N\ge J\), the absolute angular theorem applies. This
proves the all-degree assertion for every member of (1).

## The one-variable exceptional bound

The progression also certifies a nonzero polynomial
obstruction for the entire baseline fiber. Write
\(f_k(A)=[z^k]F_m(A,z)\), and let \(z_i\) be the roots
over the algebraic closure of \(\mathbb Q(A)\), with
multiplicities. Use the ordered quotient polynomial
\[
 D_m(A,u)=f_m(A)^{2(m-1)}
               \prod_{i\ne j}(u z_i-z_j)
    =\frac{(-1)^m\operatorname{Res}_z(F_m(A,z),F_m(A,uz))}
                  {(u-1)^m f_0(A) f_m(A)}.                   \tag{3}
\]
The factor \((u-1)^m\) removes the diagonal pairs.
The endpoint factors are independent of \(u\) and
nonzero at every ordinary \(A\).

For (3) to be a polynomial in \(A,u\), use the elementary
symmetric functions in the roots. The product before the
prefactor is homogeneous of root degree \(m(m-1)\), with
degree at most \(2(m-1)\) in any individual root. Its
expression in elementary symmetric functions uses at most
\(2(m-1)\) factors per term. This follows by eliminating
successive leading symmetric monomials: the leading
exponent partition \(\lambda\) is eliminated by a
product with \(\lambda_1\) elementary factors, one for
each column of the partition. Since that product also has
degree at most \(\lambda_1\) in each root, the degree bound
persists. The prefactor in (3) therefore clears all
denominators. Substituting \(e_i=(-1)^i f_{m-i}/f_m\) and
using \(\deg_A f_k=dk\) gives
\(\deg_A D_m\le d\,m(m-1)\): a term of elementary weight
\(m(m-1)\) subtracts \(d\,m(m-1)\) from the prefactor
bound \(2d\,m(m-1)\). Also \(\deg_u D_m=m(m-1)\), with
leading coefficient a nonzero product of endpoint powers.

At a single \(A\) from (1), all orders \(3\le m<J\)
are simple by their own upper Hensel mark, and degree two
is simple because its roots are nonreal conjugates.
The quotients of distinct roots are disjoint at each
adjacent pair by the proof above. Therefore for \(2\le n<N<J\)
the polynomial \(\operatorname{Res}_u(D_n,D_N)\) in \(A\)
is nonzero and has degree at most
\(2d\,n(n-1)N(N-1)\).

At an ordinary angular failure the two conjugate pairs
give a common quotient and its inverse. Both are distinct
unless the quotient is \(-1\), in which case the two
ordered pairs give multiplicity at least two at \(-1\)
in both quotient polynomials. Real common arguments
were already excluded. Hence the specialized Sylvester
matrix has corank at least two. A determinant of a
polynomial matrix with corank at least two at \(A_0\)
has a zero of order at least two there, by constant
row operations making two rows zero at \(A_0\).
Thus there are at most \(d\,n(n-1)N(N-1)\) ordinary bad
baselines at this pair. Summing gives the same cofinite
bound as the [unordered quotient calculation](MULTIVARIATE_SINGLE_ROW_ANGULAR_EXCEPTIONAL_FIBERS.md#an-integer-polynomial-for-the-ratio-obstruction),
with the additional explicit progression of good baselines.
