# All-degree angular separation along fixed affine rays of exponent bases

Fix the positive integer gaps \(d_i\) and orientations of the unequal
exponent rows of an ordinary two-monomial profile. Let the smaller
exponent in each row be
\[
                        \ell_i(T)=a_i+s_iT ,
\]
where \(a_i,s_i\) are nonnegative integers, at least one \(s_i>0\),
and \(T\) ranges over nonnegative integers. Rows with \(s_i=0\) remain
fixed and may include arbitrary private rows. Let \(J\ge4\) be at least
the [effective absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).

For each \(3\le N<J\) define
\[
 \tau_N=\max_{s_i>0}-\frac{Na_i+1}{N s_i}.
\]
Among the moving rows attaining this maximum, let \(u_N\) count the
increasing orientations and \(v_N\) the decreasing orientations.

**Theorem.** Suppose
\[
                   \gcd(N,|u_N-v_N|)=1\qquad(3\le N<J).
 \tag{1}
\]
Then the actual ordinary consecutive moments have disjoint root
arguments in every degree for all but finitely many \(T\ge0\).
The proof also gives a full explicit congruence class of such \(T\)'s,
using distinct sufficiently large primes for the relevant pairs.
No rows are appended and the signed total gap is fixed.
More generally, at each terminal one can use any rational parameter
where no lower factor vanishes and every vanishing upper factor is
at an offset between \(1\) and its row's gap. Count the orientations
of those upper factors in (1). The rightmost zero above is an
explicit choice with these properties.
In particular this holds for every ray with one moving row.
For a prescribed set of moving rows, it holds for positive integer
directions outside finitely many explicit rational hyperplanes.
For common translation of the moving rows, it holds if their
minimum-baseline cohort has one more row of one orientation than
of the other.

If \(G=\sum_{s_i>0}d_i\), the number of exceptional nonnegative
integer parameters is at most
\[
                   G\sum_{n=2}^{J-2} n^2(n^2-1).             \tag{2}
\]
The bound for a one-row ray depends only on its gap and the
absolute cutoff, independently of the retained exponents.

## The rightmost row-factor zero

At order \(m\), strip \((m\ell_i(T))!\) from each row
contribution. For the coefficient at \(k\), the remaining
factor from row \(i\) is
\[
 \prod_{h=1}^{d_i k}(m\ell_i(T)+h)\quad\text{or}\quad
 \prod_{h=1}^{d_i(m-k)}(m\ell_i(T)+h),
\]
depending on the orientation. These products with the binomial
coefficient give \(F_m(T,z)\) with polynomial coefficients in \(T\).
Their degrees in \(T\) are at most \(mG\).

For a fixed terminal \(N<J\), put \(n=N-1\) and \(\tau=\tau_N\).
For each moving row both
\[
  -\frac{ma_i+h}{m s_i}\le -\frac{ma_i+1}{m s_i}
       \quad(h\ge1),\qquad
  -\frac{na_i+1}{n s_i}<-\frac{Na_i+1}{N s_i}
\]
hold. Thus \(\tau\) is a zero only of the upper-order factors at
\(h=1\) in the stated attaining cohort, and no lower-order factor
vanishes there. Fixed-row factors are positive constants.
Choose a prime \(p>N\) outside the denominators of \(\tau\),
the primes dividing the nonzero slopes in the cohort, and the
prime divisors of the numerators and denominators of the
nonzero values of all factors at \(\tau\). These are finite
restrictions explicitly given by the exponent data.
Impose the congruence \(T-\tau\equiv p\pmod{p^2}\).
For an alternative rational parameter as in the theorem, the
same construction applies, and an offset \(1\le h\le d_i\)
affects exactly the same coefficient indices as \(h=1\).

At such a parameter every lower coefficient is a \(p\)-adic unit.
At the upper order, each of the selected factors at \(h=1\)
has valuation one. It contributes to all positive indices for an
increasing row, and to all indices below \(N\) for a decreasing row.
In particular the valuations at \(0,N\) are \(v_N,u_N\), and
each interior valuation is \(u_N+v_N\).
This uses the same two adjacent coefficient arrays from the given
ray. No independent choices of their coefficients are made.

## The first nontrivial term of an upper root quotient

Write \(u=u_N,\ v=v_N\), and take \(\pi^N=p\) in a valued
algebraic closure of \(\mathbb Q_p\). In \(F_N\) substitute
\(z=\pi^{v-u}x\) and divide by \(p^v\). The rescaled polynomial has
unit endpoints \(c_0,c_N\) in \(\mathbb Q_p\).
Its interior term at \(k\) has the form
\(b_k\pi^{Nu+(v-u)k}x^k\), with \(b_k\in\mathbb Z_p^\times\).
Under (1) we have \(u\ne v\). The minimum interior exponent is
\[
                   \nu=N\min(u,v)+|u-v|
\]
and it occurs at a unique \(k_*\), either \(1\) or \(N-1\).

The exact binomial \(c_0+c_Nx^N\) has distinct unit roots
\(\lambda_1,\ldots,\lambda_N\) in an unramified extension, by
Hensel's lemma and \(p>N\). Lifting at each \(\lambda_i\) in
the field containing \(\pi\) gives
\[
 x_i=\lambda_i-
       \frac{b_{k_*}\lambda_i^{k_*}}
            {N c_N\lambda_i^{N-1}}\,\pi^\nu
       +O(\pi^{\nu+1}).                                       \tag{3}
\]
Indeed all other interior terms have exponent at least
\(\nu+1\), the binomial derivative is a unit, and the quadratic
part of the Taylor correction has valuation at least
\(2\nu\ge\nu+1\) in powers of \(\pi\).
Divide (3) by \(\lambda_i\).
The relative correction is a unit multiple of
\(\lambda_i^{k_*}\), because \(\lambda_i^N=-c_0/c_N\).
These relative corrections have distinct reductions, since
\(k_*=1\) or \(N-1\) is coprime to \(N\).
For every two distinct upper roots, their quotient \(q\) therefore
satisfies
\[
          v_p(q-\omega)=\frac{\nu}{N},\qquad
          \omega=\lambda_i/\lambda_j\in\mu_N\setminus\{1\}.
 \tag{4}
\]
The scale \(\pi^{v-u}\) cancels in this quotient.
By (1), the number \(\nu/N\) has denominator exactly \(N\)
in lowest terms.

## A lower quotient cannot have this valuation

The lower roots are all units. Each root \(\alpha\) generates a
finite extension of \(\mathbb Q_p\) of degree at most \(n<p\).
Its ramification index \(e_\alpha\le n\) is prime to \(p\).
After a sufficiently large finite unramified extension \(K\)
of \(\mathbb Q_p\), one can write
\(\alpha\in K(p^{1/e_\alpha})\).
For this assertion, include the maximal unramified subfield of
\(\mathbb Q_p(\alpha)\) in \(K\). If \(\varpi\) is a uniformizer
of \(\mathbb Q_p(\alpha)\), then \(p/\varpi^{e_\alpha}\) is a unit.
Enlarge \(K\) to contain an \(e_\alpha\)-th root of its residue
and the needed roots of unity. Hensel lifting produces an
\(e_\alpha\)-th root of \(p\) in the compositum. The degree over
\(K\) is at most \(e_\alpha\) and is at least \(e_\alpha\) by
Eisenstein for \(x^{e_\alpha}-p\), so this generates the
compositum.

Use a common \(K\) for two lower roots \(\alpha,\beta\), also
containing \(\mu_N\). In the power basis of
\(p^{1/\operatorname{lcm}(e_\alpha,e_\beta)}\), the nonzero
coordinates of \(\alpha-\omega\beta\) have fractional
valuations in
\[
                 (1/e_\alpha)\mathbb Z\ \cup\
                 (1/e_\beta)\mathbb Z .
\]
Indeed \(\alpha,\beta\) use only the respective subfield
powers, multiplying by \(\omega\in K\) does not add powers,
and distinct basis powers have distinct fractional
valuations. A nonzero difference has the least valuation of
its nonzero coordinates. Since \(\beta\) is a unit,
a nonzero \(\alpha/\beta-\omega\) has valuation in the
same union. Neither index \(e_\alpha,e_\beta\le n\) can
supply a valuation with reduced denominator \(N\).
Thus (4) rules out every equality of an upper quotient of
distinct roots with a lower quotient.

For a separate local variant one can use a prime dividing some
cohort slopes. If \(p>N\), \(\tau\) is \(p\)-integral, and
all nonmarker factors at \(\tau\) are units, take
\(v_p(T-\tau)=e>0\). Each marker has weight \(e+v_p(s_i)\).
The proof applies with \(u,v\) the sums of these weights
in the two orientations, provided \(\gcd(N,|u-v|)=1\).
This can separate a pair even when the two unweighted
cohort counts agree.

## One progression and a finite exceptional set

Choose distinct primes satisfying the above requirements for
\(N=3,\ldots,J-1\) and combine the congruences modulo their
squares. The Chinese remainder theorem produces a full
progression. The lower/upper quotient separation at each
of its places rules out equality of arguments of two nonreal
roots, since such equality would give \(z/\bar z=w/\bar w\).
The even moment is strictly positive on the real line by its
product exponential integral, so the pair at terminal two
and common real arguments in any pair are excluded.
The absolute angular theorem covers \(N\ge J\).

For the finiteness assertion, take the degree-\(m\) polynomial
\(F=F_m(T,z)\), with endpoint coefficients \(f_0,f_m\), and form the
[ordered quotient resultant](MULTIVARIATE_SINGLE_MOVABLE_ROW_PADIC_PROGRESSION.md#the-one-variable-exceptional-bound)
\[
       D_m(T,y)=
       \frac{(-1)^m\operatorname{Res}_z(F(z),F(yz))}
            {f_0 f_m(y-1)^m}.
\]
The divisions are exact in the universal coefficient ring.
The power of \(y-1\) removes the diagonal terms of the root
product. The factors \(f_0,f_m\) divide the homogeneous Sylvester
resultant because setting either endpoint coefficient to zero
makes both polynomials vanish at the corresponding homogeneous
endpoint. In root notation the resulting polynomial is
\[
        D_m(T,y)=(f_0f_m)^{m-1}
                  \prod_{i\ne j}(y-z_j/z_i).
 \tag{5}
\]
Its degree in \(y\) is \(r_m=m(m-1)\) at every nonnegative
integer \(T\) because the endpoints are nonzero there.

For a degree bound in \(T\), let \(G_+,G_-\) be the sums of
the gaps of the increasing and decreasing moving rows, so
\(G=G_++G_-\). The coefficient at \(k\) has degree exactly
\(mG_-+k(G_+-G_-)\) in \(T\). Thus
\[
       T^{-mG_-} F_m\big(T,T^{-(G_+-G_-)}z\big)
\]
has a coefficientwise limit as \(T\to+\infty\) with both
endpoint coefficients positive. The roots of this rescaled
polynomial stay in a fixed compact annulus for large \(T\),
by the monic root bound and reversal. Its root quotients
are unchanged by rescaling and are therefore bounded.
In (5), the degree of \((f_0f_m)^{m-1}\) is \(Gm(m-1)\).
Every coefficient of the quotient product is bounded on the
positive ray, so \(\deg_T D_m\le G r_m\).

For a parameter in the congruence class all upper roots
are simple by lifting, and no ordered quotient of distinct upper
roots equals a lower quotient. An upper off-diagonal quotient
also cannot equal \(1\), even if the lower polynomial has
multiplicities. Thus the degree-\(r_n,r_{n+1}\) polynomials
have no common quotient there. It follows that
\(R_n(T)=\operatorname{Res}_y(D_n,D_{n+1})\) is nonzero as a
polynomial for each \(2\le n\le J-2\).

At an exact angular failure, the shared quotient \(w\ne1\) has
modulus one. Both \(w,w^{-1}\) occur, and if \(w=-1\) the two
ordered conjugate pairs give multiplicity at least two in
each quotient polynomial. Thus the Sylvester matrix of \(R_n\)
has corank at least two there, so the parameter is a zero of
multiplicity at least two. Constant row operations can make
two rows vanish at that parameter, and each such row contributes
one factor to the determinant. The Sylvester degree bound
gives \(\deg R_n\le2G r_n r_{n+1}\), proving (2).
In particular a computable root bound on the product of the
nonzero \(R_n\) bounds every exceptional integer \(T\).

Finally, ties in \(\tau_N\) can occur only on the hyperplanes
\((Na_i+1)s_j=(Na_j+1)s_i\) between positive moving directions.
Outside this finite collection, each upper cohort consists of
one row. There is also only one row in every cohort for the
one-row ray. For a common translation, the attaining cohort
is the minimum-baseline cohort independently of \(N\).
These observations give the examples in the theorem.

## Bounded exceptional fibers in the original exponent grid

Take any prescribed signed gaps and allow \(r\ge1\) smaller
exponents to vary, with all other rows fixed. One can construct
a full Cartesian product of integer congruence classes on
which every profile has all-degree angular separation. Choose
one varying base \(A_*\). At each terminal \(3\le N<J\), choose a
fresh prime \(p_N\) larger than \(N^2\) times the largest varying
gap and \(N\) times the largest fixed exponent. Impose
\[
               N A_*+1\equiv p_N\pmod{p_N^2}
       \qquad\hbox{and}\qquad
               A_i\equiv0\pmod{p_N}\quad(i\ne *)
\]
on the varying bases. All unselected rising factors reduce to
their nonzero indices at both orders, while the fixed factors
are units by the size condition. Thus the proof with one
marker applies at every terminal. One can refine any compatible
prescribed coordinate congruences by choosing the new primes
coprime to the prescribed moduli and using CRT.

Let \(d_*\) be the smallest gap among the varying rows. In the box
\(\{0,\ldots,H\}^r\), at most
\[
       d_*\!\sum_{n=2}^{J-2} n^2(n^2-1)\,(H+1)^{r-1}
\]
profiles have a pair of consecutive root arguments in common.
Fix all coordinates except the one with gap \(d_*\) and apply
(2) to that one-row ray. The bound is independent of the fixed
coordinates, dimension outside the box, and sizes of the other
gaps. It also applies to a box of translated integer intervals.
Thus the original profiles have all-degree angular separation
outside a set of density zero at any fixed signed gaps and
with at least one variable smaller exponent.
