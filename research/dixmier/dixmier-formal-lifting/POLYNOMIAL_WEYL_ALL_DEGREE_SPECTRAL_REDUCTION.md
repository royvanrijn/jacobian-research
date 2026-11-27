# A one-pole spectral reduction for the rank-one Dixmier conjecture

## Centered rectangles and the finite spectral pole

**Theorem PolynomialWeylSingleRegularPoleReduction.** If the rank-one
Dixmier conjecture fails over \(\mathbb C\), there is a polynomial
pair \(P,Q\) with \([Q,P]=1\) and support rectangles with corners
\((au,av),(bu,bv)\), where \(0<u<v\) and \(a,b\ge2\) are coprime.
In its differential presentation \(X=x,\ Y=D_x\), the leading
coefficients are proportional to \(x^{au},x^{bu}\).
Writing \(k=v-u\), the maximum of \(j-i\) over the support of
each operator is \(ak,bk\), respectively. For every spectral value,
\(P-\lambda\) and \(Q-\lambda\) have exactly one possible finite
singularity, at zero, and are regular singular there. Conversely,
a scalar-commutator pair with such rectangles and corners is proper.

Take a standard minimal counterexample using
[Guccione–Guccione–Valqui, Corollary 6.11](sources/guccione-guccione-valqui-1111.6100v3.pdf).
Proposition 6.3 gives nondividing Bernstein degrees. At weight
\((1,1)\) the corner monomials Poisson commute, so they are
proportional powers of a common monomial with exponent pair \((u,v)\).
Both coordinates are positive by subrectangularity. The negative
diagonal difference on the maximal \(X\)-face in Definition 6.10
implies \(u<v\).

The differential orders in \(Y=D_x\) are \(av,bv\), and their
common leading factor has degree \(u\). Since \(u<v\), the
[residue theorem](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md)
excludes two or more distinct roots. Translate its unique root
to zero. This preserves the rectangles and the corner terms.
Thus the maximal-\(Y\) faces are the single corners.

The following argument, used at degree forty-two in the
[regular finite pole proof](POLYNOMIAL_WEYL_FORTY_TWO_REGULAR_FINITE_POLE.md),
applies with arbitrary \(u,v,a,b\). Choose primitive weights
\((-p,q)\) with \(0<p<q\). For the first operator a maximizing
point \((i,j)\) satisfies
\[
       j-(p/q)i\ \ge\ av-(p/q)au,\qquad j-i\ge a(v-u)>0.
\]
The same holds with \(b\) for the mate. By Guccione–Guccione–Valqui,
Theorem 4.1, its leading symbol \(R\) has a polynomial bracket
mate \(H\) of weight \(q-p\) with \(\{R,H\}=cR\), \(c\ne0\).
Write \(T=X^qY^p\), \(R=X^iY^j L(T)\) with \(L(0)\ne0\),
and \(H=XYG(T)\). The form of \(H\) follows from integrality
since \(q\ge2\). Up to the sign of \(c\), the equation is
\[
   (i-j)LG+(q-p)TL'G+(pi-qj)TLG'=cL.
\]
If \(\ell=\deg L\) and \(m=\deg G>0\), its leading coefficient
at degree \(\ell+m\) contains the factor
\[
       (i+q\ell)-(j+p\ell)+(pi-qj)m .
\]
Both summands are negative by the maximizing-point inequality
and the positive weight. Thus \(G\) is constant. Now the
coefficient comparison forces \(L\) constant, since \(q-p>0\).
Every weight \((-s,1)\) for rational \(0<s<1\) therefore has
a unique maximizer. The single corner at \(s=0\) continues to
maximize up to \(s=1\): any earlier change would tie two
integral support points at rational \(s\). This proves the
claimed support bounds \(j-i\le ak,bk\). After multiplying
\(P-\lambda\) by \(x^{ak}\) it is a polynomial in \(x\) and
the Euler operator \(E=xD_x\), with constant nonzero leading
coefficient in \(E\). The analogous description holds for \(Q\).
This is the regular singular Euler form at zero.

For the converse it is enough to use the first Bernstein
monomial \(X^{m}Y^{n}\), with \(m=au>0\), \(n=av>m\).
The leading terms of the iterates of \(\operatorname{ad}_P\)
on \(Y\) have exponents \(((m-1)j,1+(n-1)j)\) at step \(j\).
The successive bracket factors are proportional to
\(-m+(m-n)j\) and never vanish. An image of a Weyl generator
in a generating pair has locally nilpotent adjoint action.
Thus the displayed pair is proper.

## A uniform lower bound on the diagonal gap

**Corollary PolynomialWeylDiagonalGapSix.** In the standard minimal
presentation one has \(v-u\ge6\). More generally, in a
polynomial scalar-commutator pair of nondividing differential orders
of gcd \(g\), a geometric root of the common leading factor has
multiplicity at most \(g-1\) or at least \(g+6\).

Consider the strict crossing from the
[translated-root argument](POLYNOMIAL_WEYL_TRANSLATED_ROOT_CROSSINGS.md#the-finite-list-for-a-fixed-differential-order-gcd),
or the standard minimal crossing from Proposition 7.2.
Put \(e=v_0-u_0>0\) at its normalized upper endpoint
\(C_0=(u_0,v_0)\), and write \(C_1=(r,t)\) for the other.
The diagonal deficiency strictly decreases as the incoming
endpoints move toward the crossing. Thus \(e\) is at most the
initial deficiency, itself at most \(v-u\) in the standard
presentation, or \(m-g\) in the translated-root presentation
at a root of multiplicity \(m>g\).

If the bracket-mate endpoint is \((f_1,f_2)=(k_0,k_0+j)\),
then \(u_0=k_0 e/j\), \(1\le j<e\), and the normal direction
before dividing by the gcd is \((k_0+j-1,1-k_0)\).
The unscaled weight at \(C_0\) is exactly \(e\).
The weight at \(C_1\) consequently gives
\[
                    e=(k_0-1)(r-t)+jr .
\]
In particular \(k_0+j-1\le e\), \(r\le e\) and \(u_0\le e^2\).
This gives a finite list in any fixed gap irrespective of the
ambient degrees. The
[integer check](verify_weyl_diagonal_gap.py) compares these
equations with an independent bounded equal-weight scan and
applies Proposition 7.3(6). For \(e\le5\) the sole remaining
rows are \((6,10),(8,12),(10,15)\).

The first row is excluded for arbitrary coprime quotients by
the [local scalar successor calculation](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_THIRTY_THREE.md#the-necessary-scalar-successor).
For \((8,12)\), the [sum-twenty proof](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_FORTY_ONE.md#sum-twenty)
has commuting-edge tests depending only on the common factor.
At all its endpoints, of lattice denominator \(l\), normalized
upper coordinate \(C_x\) and second coordinate \(\mu\),
the deficit \(\beta=l(\mu-C_x)\) is at most \(\mu\).
Thus its scalar argument also applies to arbitrary \(a,b\):
the reduced denominator is at least
\(((a+b)\mu-1)/\beta>a+b-1\) and coprime to \(\mu\).
It cannot divide either of \(a,b\), which is necessary for a
nonzero bracket at the origin. For \((10,15)\), the
[sum-twenty-five proof](POLYNOMIAL_WEYL_CONDITIONAL_THROUGH_FIFTY_ONE_FROM_F7.md#degree-fifty-and-the-sum-twenty-five-crossing)
has shifted state \(C=(6,7)\) with \(l=2,\beta=2\).
The same scalar argument applies. Its commuting multiple-root
rows all have integer multiplier \(M/\mu\), which is excluded
by the rational-primitive obstruction. Single-root shifts retain
the endpoint and lattice while the allowed slopes strictly
increase in a finite set. This excludes all three cases.
Exact multiplicity \(g\) is already excluded by the residue equation.

## The spectral monodromy in arbitrary orders

**Theorem PolynomialWeylSinglePoleSpectralMonodromy.** In the
centered rectangles above, let \(M_P(\lambda),M_Q(\lambda)\)
denote the monodromy around zero for the two spectral equations.
For each operator its monodromy conjugacy class is constant in
\(\lambda\). At each eigenvalue, every Jordan block size occurs
with multiplicity divisible by \(a\) for \(P\), and by \(b\) for
\(Q\). The monodromy fixed spaces consist of entire solutions
and have dimensions exactly \(ak,bk\), respectively.

Here are the changes needed to use the
[degree-forty-two spectral proof](POLYNOMIAL_WEYL_FORTY_TWO_REGULAR_SPECTRAL_MONODROMY.md)
in arbitrary orders. The leading negative Euler components are
\[
              P_0=x^{-ak}f(E),\qquad Q_0=x^{-bk}g(E).
\]
They commute because the commutator is scalar. The two operators
\(P_0^b,Q_0^a\) commute with \(P_0\) and have the same shift.
Their coefficient ratio in \(E\) is rational of period \(ak\)
by the commutation identity
\(F(E)x^{-j}=x^{-j}F(E-j)\). Such a rational function is constant
in characteristic zero. Thus \(Q_0^a=cP_0^b\) with \(c\ne0\).
Polynomiality gives \((E)_{ak}\mid f\) and \((E)_{bk}\mid g\),
where \((E)_j\) is the falling factorial.

On \(P-\lambda\) jets at an ordinary point, \(\partial_\lambda+Q\)
is flat with the differential equation. Indeed differentiation of
\((P-\lambda)y=0\) gives \((P-\lambda)\partial_\lambda y=y\),
and \((P-\lambda)Qy=-y\). Reducing the jets by the scalar
equation gives matrices rational in \(x\) with no poles except zero
and polynomial in \(\lambda\). This connection intertwines
analytic continuation around zero at all parameters. The same
argument for \(Q-\lambda\) uses \(\partial_\lambda-P\).

For \(P\), write \(\lambda=s^{ak}\), \(x=z/s\), and scale
\(P_s=s^{-ak}P(z/s,sD_z)\), \(Q_s=s^{-bk}Q(z/s,sD_z)\).
Their limits at \(s^{-1}=0\) are \(P_0,Q_0\) in \(z\).
The top coefficient of \(P_s-1\) is unchanged,
\(\alpha z^{au}\), so the equation on a fixed \(z\)-loop around
zero has analytic dependence on \(s^{-1}\), including at zero.
The spectral connection in jets at \(z=1\) is of the form
\[
        s^{(a+b)k-1}\bigl(B_0+O(s^{-1})\bigr),\qquad
        B_0=-ak[Q_0]_{P_0-1}.
\]
This is the chain rule applied to spectral transport, with the
additional term \(-zD_z/s\) of lower order. The constant matrix
\(B_0\) is semisimple because its appropriate power is a
nonzero scalar by \(Q_0^a=cP_0^b\).

Formal block diagonalization by the distinct leading scalars
uses only division by differences of those scalars, since the
coefficient of a derivative enters one step later than the
leading commutator in the formal recursion. A horizontal
endomorphism preserves each block: an off-diagonal entry with
nonzero Laurent head would have a leading right-hand term of
degree \((a+b)k-1\) higher than its head, whereas differentiation
lowers degree. Apply this to the scaled monodromy. Under
\(s\mapsto\zeta s\) for a primitive \(ak\)-th root \(\zeta\),
the connection leading scalars are multiplied by
\(\zeta^{(a+b)k}\), a primitive \(a\)-th root of unity.
The identification of jets is transport from \(z=1\) to
\(\zeta^{-1}\), followed by derivative rescaling.
It is analytic and invertible at \(s^{-1}=0\) and intertwines
monodromy and both connections. The same head comparison
makes it permute all the blocks in one cycle. They therefore
have identical Jordan partitions. Ranks of monodromy tests
computed over the formal field agree with generic analytic
ranks and hence with every finite spectral value by flatness.
For \(Q\) the scaling \(\lambda=s^{bk}\) gives a cycle of
length \(b\), using \(P_0\) in place of \(Q_0\).

For \(P-\lambda\), the output coefficient at \(x^{n-ak}\)
of a formal Taylor series uses diagonal \(f(n)\), plus finitely
many coefficients at lower indices. Above some finite \(n\)
the diagonal never vanishes and all solutions of the finite
initial system extend by recursion. If the initial unknowns
are indexed \(0,\ldots,N\), there are at most \(N-ak+1\)
equations because a polynomial differential operator has
no negative output coefficients on Taylor series. Hence
there are at least \(ak\) formal solutions. The coefficient
ratios in the tail are bounded since the diagonal has degree
\(av\) and the other coefficient polynomials have degree
at most \(av\). The finite-lag recurrence gives an exponential
bound, so these solutions converge. Zero is the only possible
finite singularity, and continuation makes these solutions entire.

In the homogeneous limit the equation is
\(f(zD_z)y=z^{ak}y\). A monodromy-fixed local solution is
meromorphic. The Euler companion equation has bounded
coefficients for \(zD_z\) near zero. Integration along
rays bounds the solution by a power of \(|z|^{-1}\),
uniformly over the circle for a single-valued solution.
Its Laurent series therefore has finite principal part
and satisfies \(f(n)y_n=y_{n-ak}\).
Each residue modulo \(ak\) contains an integer root of \(f\)
by \((E)_{ak}\mid f\). In each residue a nonzero solution
must start at the largest such root, since any earlier start
would run into a zero diagonal with nonzero preceding value.
The ensuing series converges by the recurrence. Thus the
limit fixed dimension is exactly \(ak\). Analytic
specialization at \(s^{-1}=0\) can increase the fixed dimension,
so the generic finite fixed dimension is at most \(ak\).
Together with the entire solutions this gives equality at
every parameter. Interchanging the operators gives the other count.

**Corollary PolynomialWeylSinglePoleIndicialPartitions.** In each
orbit of indicial roots modulo \(k\mathbb Z\), the polynomials
\(f,g\) have respective root counts \(a h,b h\) for an integer
\(h\ge0\), counting multiplicities. For a class modulo \(\mathbb Z\),
let \(r\) be the decreasing list of its positive \(h\)'s.
At the eigenvalue of that class, the normalized actual Jordan
partitions for \(P,Q\), with block multiplicities divided by
\(a,b\), both dominate \(r\).

Represent a \(k\mathbb Z\)-orbit by \(\rho+k\mathbb Z\) and
record multiplicities of \(f,g\) as Laurent polynomials \(F,G\)
in an indeterminate \(z\), using the integer index in this orbit.
The commutation of \(P_0,Q_0\) gives
\[
          (1-z^b)F=(1-z^a)G .
\]
The polynomials \(1+\cdots+z^{a-1}\) and \(1+\cdots+z^{b-1}\)
are coprime. Thus \(F=(1+\cdots+z^{a-1})H\) and
\(G=(1+\cdots+z^{b-1})H\) for an integral Laurent polynomial
\(H\), whose coefficients need not be positive. Put \(h=H(1)\).
Each residue modulo \(a\) in the orbit has \(h\) roots of \(f\),
and each residue modulo \(b\) has \(h\) roots of \(g\). In
particular \(h\ge0\).

For the homogeneous model of \(P\) in coordinate \(w=x^{ak}\),
an indicial class in \(w\) contributes one Jordan block of
size equal to its root count. To see this, divide by a chosen
complex exponent of \(w\). An eigenvector of monodromy then
has a Laurent series with finite principal part by regular
singular growth. Its one-step recurrence permits only the
largest indicial root in its class as a starting exponent.
The model algebraic multiplicity is the count of exponentiated
indicial roots, as follows by shrinking the loop for the
Euler companion system to its residue. Pullback to \(x\)
takes the \(ak\)-th power of monodromy and preserves block sizes
at nonzero eigenvalues. Thus each \(h>0\) contributes \(a\)
blocks of size \(h\) for \(P_0-1\) and \(b\) for \(Q_0-1\).
The model and actual characteristic polynomials coincide by the
analytic specialization in the preceding theorem. The ranks
of powers of \(M-\chi I\) for a fixed eigenvalue \(\chi\)
can only fall in the limit. Subtracting the common contribution
of the other eigenspaces and dividing by \(a\) or \(b\) gives
the asserted partition dominance.

## An algebraic model for the entire solution connection

**Corollary PolynomialWeylSinglePolePolynomialConnection.** The
entire solution spaces in the centered chart of \(P-\lambda\)
are the fibres of a polynomial connection of rank \(ak\).
It is the dual of \(\mathbb C[z]\) as a free
\(\mathbb C[\lambda]\)-module with \(\lambda\) acting by
\(P^\vee\) and connection operator \(-Q^\vee\), where
\(X^\vee=D_z,\ Y^\vee=z\) and \(\vee\) reverses products.
The \(Q\)-connection similarly has rank \(bk\), on the same
polynomial space with the roles exchanged and connection operator \(P^\vee\).

Pair a formal Taylor series \(y\) with \(z^n\) by \(y^{(n)}(0)\).
This pairing transposes \(X,Y\) to \(D_z,z\), respectively.
For an operator with rectangle \((au,av)\), the maximal
degree increase on polynomials in the adjoint presentation is
\(ak\). The coefficient of that increase on \(z^n\) is a polynomial
in \(n\) of degree \(au\): the corner term contributes the
unique term of that degree. Thus \(P^\vee-\lambda\) raises
degree by exactly \(ak\) on every polynomial of sufficiently
large degree, uniformly in \(\lambda\).

If there were a polynomial eigenvector for \(P^\vee\), successive
applications of \(Q^\vee\) would produce an infinite independent
family of generalized eigenvectors at the same value, by the
scalar commutator. But any generalized eigenvector has bounded
degree by the preceding degree increase. Hence there is no
polynomial eigenvector. This makes \(\mathbb C[z]\) torsion-free
over \(\mathbb C[\lambda]\), since any polynomial in \(P^\vee\)
factors into linear terms. It is finitely generated by degree
reduction with the high-degree monomials \(P^\vee z^n\).
It is therefore free. To compute its rank, take polynomial
spaces of degrees at most \(N\) and \(N+ak\), for \(N\) large.
The image of \(P^\vee-\lambda\) intersects the latter in
exactly its image on the former, and every higher polynomial
reduces to degree at most \(N+ak\). Its cokernel has dimension
\(ak\), giving the rank.

An \(\mathbb C[\lambda]\)-linear functional \(\phi\) on this
module gives the formal solution with derivatives
\(y^{(n)}(0,\lambda)=\phi(z^n)\), and every polynomial-parameter
formal solution comes from such a functional. This follows
by transposing the relation \(P-\lambda\).
The dual connection acts by
\(\phi\mapsto\partial_\lambda\phi+\phi\circ Q^\vee\),
which is exactly \(\partial_\lambda+Q\) on the solutions.
Specializing the free module and using convergence at the
regular pole identifies its fibres with the entire solutions.

## The obstruction required of a general proof

**Corollary PolynomialWeylSinglePoleSpectralOrbitBound.** In the
scaled \(P\)-spectral equation, the monodromy has a proper horizontal
fixed subspace of dimension \(ak<av\). If a formal refinement of
the spectral blocks is preserved by monodromy and permuted in
one orbit of length \(d\) by the deck identification, then
\[
                     d\mid a\gcd(u,v).
\]
For \(Q\) the corresponding bound is \(d\mid b\gcd(u,v)\).

Horizontality follows from the commutator equation for monodromy.
The deck identification intertwines the monodromy and carries
each block's fixed space to the next. Thus the common block
rank and common fixed rank give \(d\mid av\) and \(d\mid ak\).
The result follows from \(\gcd(v,k)=\gcd(u,v)\).
In particular, a full formal splitting into distinct
rank-one irregular types in a single deck orbit would exclude
the centered pair. A proof of such a splitting, or of an orbit
length incompatible with the displayed divisibility, would give
an all-degree route. The standard minimal reduction supplies the
centered hypotheses for a hypothetical counterexample.

The [fixed and moving spectral comparison](POLYNOMIAL_WEYL_FIXED_AND_MOVING_SPECTRAL_SPLITTINGS.md)
determines the irregular types after all lower terms are included.
In the moving presentation their transitive orbit always satisfies
the displayed divisor, with repeated rank-one factors. At a fixed
ordinary spatial point the full rank-one splitting is transitive,
and its use for DC1 now requires a formal comparison with the
entire-solution connection.

Any finite homogeneous differential operator \(T=x^{-j}r(E)\)
commuting with \(P_0\) has \(k\deg r=vj\). In fact the highest
differential commutator coefficient against the leading term
of \(P_0\) is proportional to
\(av(\deg r-j)-au\deg r\). Thus \(k/\gcd(u,v)\mid j\).
In the scaled spectral system a deck character given by such
a shift, \(\zeta^j\) for \(\zeta^{ak}=1\), already has
order dividing \(a\gcd(u,v)\). This applies to arbitrarily
many commuting homogeneous heads. In the
[degree-forty-two second splitting](POLYNOMIAL_WEYL_FORTY_TWO_SECOND_SPECTRAL_SPLITTING.md),
the refined orbit sizes six and twenty-one attain exactly this
divisor, since \((u,v)=(6,15)\).

The lowest Euler components alone permit the monodromy
divisibilities and the nontrivial fixed subspace. For example,
let \(\alpha\notin\mathbb Z\) and take
\[
        R=x^{-k}(E)_k(E-\alpha)^u,\qquad P_0=R^a,\quad Q_0=R^b .
\]
These are polynomial operators with the required rectangles.
For the model \(P_0-1\), the integer indicial roots are
\(0,\ldots,ak-1\), and the remaining roots are
\(\alpha+ik\), \(0\le i<a\), each repeated \(u\) times.
After \(w=z^{ak}\), each indicial class modulo \(ak\) contributes
one Jordan block. It has one eigenvector by the same one-step
recurrence after factoring out its complex exponent, and its
algebraic multiplicity is the number of exponentiated indicial
roots in the class. Thus the noninteger monodromy consists of
\(a\) blocks of size \(u\), with \(ak\) integer fixed solutions.
For \(Q_0\) the counts use \(b\).
This is commuting homogeneous data. A general exclusion needs
the constraints of the actual lower polynomial coefficients and
the scalar commutator on the subsequent spectral splitting.
