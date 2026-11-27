# Fixed and moving spectral splittings for polynomial Weyl pairs

Use the [centered Weyl reduction](POLYNOMIAL_WEYL_ALL_DEGREE_SPECTRAL_REDUCTION.md).
Over \(\mathbb C\), write \([Q,P]=1\), let \(v>u>0\), let \(k=v-u\),
and take the two rectangles \((au,av),(bu,bv)\) with \(\gcd(a,b)=1\).
The coefficient of \(D_x^{av}\) in \(P\) is \(\alpha x^{au}\), and
the coefficient of \(D_x^{bv}\) in \(Q\) is \(\beta x^{bu}\).
In particular the sole possible finite spatial pole is \(x=0\).
Write \(d=\gcd(u,v)=\gcd(k,v)\) and \(\ell=v/d\). We give the proof
for \(P-\lambda\) with \(\partial_\lambda+Q\). Interchanging \(P,Q\)
and using \(\partial_\lambda-P\) proves the claims for \(Q-\lambda\).

## Fixed ordinary spatial point

**Theorem PolynomialWeylFixedSpectralTransitivity.** Set \(n=av\)
and \(Z^n=\lambda\). At all but finitely many ordinary fixed points
\(x_*\in\mathbb C^\times\), the spectral connection splits over
\(\mathbb C((Z^{-1}))\) into \(n\) rank-one factors with pairwise
distinct irregular exponentials. The spectral deck group acts
transitively on the factors, so the connection at \(x_*\) is
irreducible over \(\mathbb C((\lambda^{-1}))\). On the spatial
cover \(x=w^v\), its exponential polynomials can be written
\[
  \Phi_i(w,Z)=F(\omega_i Z)+
       \frac{v}{k}\alpha^{-1/n}\omega_i w^kZ,\qquad
  \omega_i^n=1,\qquad F(T)\in T^\ell\mathbb C[T^\ell].
\]
Here a choice of \(\alpha^{-1/n}\) fixes the indexing.

The flatness of the two equations follows from \([Q,P]=1\), since
\((P-\lambda)\partial_\lambda y=y\) and
\((P-\lambda)Qy=-y\) on \(P-\lambda\) solutions. Use the jets
\(Z^{-j}D_x^jy\), \(0\le j<n\). Their spatial connection for
\(\partial_x\) has matrix \(ZA_0(x)+O(1)\). Its leading matrix is
the companion of \(T^n-(\alpha x^{au})^{-1}\), with simple
eigenvalues \(c_iw^{-u}\) for \(c_i=\alpha^{-1/n}\omega_i\).
Work over \(\mathbb C(w)((Z^{-1}))\) with
\(\partial_x=(vw^{v-1})^{-1}\partial_w\).
After using an eigenbasis for \(A_0\), a gauge with constant
term the identity diagonalizes the whole spatial connection.
At each successive \(Z^{-1}\) order the off-diagonal equations
divide by differences of the distinct leading eigenvalues.
Differentiation of the gauge preserves \(Z\)-degree and hence
enters a later step than commutation with \(ZA_0\). This is
a recursion over \(\mathbb C(w)\) with no poles at finite nonzero
\(w\), since the eigenvalue differences are nonzero monomials
and the original coefficient denominators are powers of \(w\).

Write the resulting spatial diagonal entries as
\(A_i=c_iw^{-u}Z+O(1)\), and write \(B\) for the transformed
spectral matrix for \(\partial_Z\). Flatness gives
\(\partial_x B_{ij}=(A_i-A_j)B_{ij}\) for \(i\ne j\).
A nonzero Laurent series cannot satisfy this equality, since
the right side would have \(Z\)-degree one higher than the
left side. Therefore \(B\) is diagonal and
\(\partial_x B_i=\partial_Z A_i\).
Its positive-degree coefficients in \(Z\) are constant
in \(w\), and its constant coefficient is
\((v/k)c_iw^k+C_i\) for \(C_i\in\mathbb C\).
Integrating the nonnegative \(Z\)-powers gives \(\Phi_i\), up to
additive constants which we take to be zero.

Use the unscaled jets to compare the gauges under deck
transformations. In these jets the projections to the spatial
lines are uniquely distinguished by the leading spatial
eigenvalues; a Laurent change of basis cannot intertwine two
different leading eigenvalues, by the same degree comparison.
A Laurent scalar change on a line has logarithmic \(Z\)
derivative \(O(Z^{-1})\). Thus the spectral deck
\(Z\mapsto e^{2\pi i/n}Z\) rotates the \(n\) labels in order and
pulls back their exponential polynomials accordingly.
All \(w\)-independent terms therefore have the form \(F(\omega_i Z)\)
with \(F(0)=0\). The spatial deck \(w\mapsto e^{2\pi i/v}w\)
also permutes the spatial lines, with the eigenvalue labels
rotated by a root of unity of order \(\ell\). On the spectral
exponentials it leaves the \(w\)-independent terms unchanged.
Consequently \(F\) is invariant under rotation of its argument
by this root, proving the stated form.

In particular \(\ell>1\), so at every nonzero \(w\) the linear
coefficients are distinct. Two distinct
exponential polynomials cannot be related by a Laurent scalar
gauge. The ramified differential module is a direct sum of
nonisomorphic simple modules, so any submodule is a sum of a
subset. Deck transitivity leaves no proper nonzero subset
that descends to the unramified field.

This calculation uses the same flat companion recursion as
[Beals and Sattinger, Lemma 9.3 and equation (9.9)](https://arxiv.org/abs/solv-int/9801010),
where the monic normalized string has spatial phase \(xzJ\).
Here the monomial leading coefficient gives \(w^kZ\), and the
spatial deck constrains the \(w\)-independent terms.

## Moving spatial point

**Theorem PolynomialWeylMovingSpectralMultiplicity.** Take
\(\lambda=s^{ak}\), \(x=z/s\) and evaluate the scaled connection
at \(z=1\). Over \(\mathbb C((s^{-1}))\) it is a sum of isotypic
rank-one blocks with distinct exponential polynomials
\(e_1(s),\ldots,e_r(s)\). The spectral deck identification
permutes these \(r\) types in one orbit, with
\[
                     a\mid r\mid ad,\qquad
   \text{multiplicity of every type}=\frac{av}{r}.
\]
The residues of all the rank-one factors agree in a gauge over
\(\mathbb C((s^{-1}))\). For \(Q\), using \(\lambda=s^{bk}\), the
corresponding formulas are \(b\mid r\mid bd\) and multiplicity \(bv/r\).
In particular each multiplicity is a multiple of \(\ell>1\).
When \(d=1\), there are exactly \(a\) types of rank \(v\) for
\(P\) and \(b\) types of rank \(v\) for \(Q\).

Set \(\tau=s^{-1}\), \(\Theta=zD_z=xD_x\),
\(P_s=s^{-ak}P(z/s,sD_z)\) and \(Q_s=s^{-bk}Q(z/s,sD_z)\).
The equation \(P_s-1\) gives in Euler jets a matrix \(A\)
for \(\Theta\), and the spectral chain rule gives a matrix \(B\)
for \(\partial_\tau\) at fixed \(z\). They are flat, with
\[
 A\in\operatorname{Mat}_n(\mathbb C(z)[[\tau]]),\qquad
 B\in\tau^{-(a+b)k-1}\operatorname{Mat}_n(\mathbb C(z)[[\tau]]).
\]
Indeed \(P_s,Q_s\) have only nonnegative \(\tau\)-powers by
the rectangles, and the top coefficient in the spatial
equation is unchanged. In Euler jets the coefficients of \(A\)
are holomorphic in \(z\) away from zero and analytic near
\(\tau=0\). The two contributions to the moving derivative
are the spectral derivative at fixed \(x\), multiplied by
\(\partial_\tau\lambda=-ak\tau^{-ak-1}\), and
\(\tau^{-1}\Theta\) from \(x=z\tau\).

First work over \(K=\mathbb C((h^{-1}))\) with \(h^v=z\).
At \(\tau=0\) the spatial equation is
\(f(\Theta)y=z^{ak}y\), where \(\deg f=n\) and the leading
coefficient of \(f\) is \(\alpha\). Shearing the Euler jets by
\(h^{-jk}\) gives a leading matrix at \(h=\infty\) with
simple eigenvalues \(c_i h^k\). Diagonalize at \(\tau=0\)
by successive orders of \(h^{-1}\). As in the fixed-point
calculation, \(\Theta\) preserves \(h\)-degree, so differentiating
a gauge term enters later than commutation with this leading
matrix. Continue the diagonalization order-by-order in \(\tau\)
with gauges over \(K[[\tau]]\). The off-diagonal operator
at each order is
\(\Theta-(A_{i,0}-A_{j,0})\) on \(K\), where \(A_{i,0}\) is
the \(i\)-th diagonal entry at \(\tau=0\). It is invertible:
the operator \((A_{i,0}-A_{j,0})^{-1}\Theta\)
decreases \(h\)-degree by at least \(k\), so its
geometric series converges formally. This produces a spatial
gauge \(G\in\mathrm{GL}_n(K[[\tau]])\), including the initial
jet shear and constant term gauge.

In this gauge flatness again makes \(B\) diagonal. The first
nonzero \(\tau\)-coefficient of an off-diagonal entry would
satisfy the homogeneous equation for the invertible operator
above. For each diagonal entry, \(\Theta B_i=\partial_\tau A_i\)
is regular in \(\tau\). Hence the entire polar part \(D_i(\tau)\)
of \(B_i\), including its residue, has coefficients in
\(\ker(\Theta:K\to K)=\mathbb C\).

Let \(\xi=e^{2\pi i/(ak)}\). The simultaneous transformation
\(\tau\mapsto\xi^{-1}\tau,\ z\mapsto\xi z\) preserves \(x\)
and \(\lambda\), and preserves the original jet presentation
since its jets are powers of \(xD_x\). Lift it using
\(h\mapsto e^{2\pi i/(akv)}h\). The leading spatial eigenvalues
are rotated by \(e^{2\pi i/(av)}\), so it permutes all \(n\)
labels in one cycle. Indeed the change between the two spatial
diagonal gauges is monomial on the matched labels. An entry
joining unmatched labels would satisfy a homogeneous spatial
equation whose first \(\tau\)-coefficient has invertible
operator \(\Theta\) minus a function of positive \(h\)-degree.
The nonzero entries of the change are \(\tau\)-adic units,
whose spectral logarithmic derivatives are regular in \(\tau\).
Thus the operation permutes the \(D_i(\tau)\,d\tau\) by
pullback, and in particular leaves their residues unchanged.
Since the labels form one cycle, all residues agree. Since
the operation on polar coefficients has order dividing \(ak\),
the number of distinct polar parts divides both \(av\) and
\(ak\). Equal residues mean polar parts differ precisely when
the exponential polynomials differ, by integrating the terms
of degree at most \(-2\). This proves the count and multiplicity
over the spatial formal field. The
[leading Euler splitting](POLYNOMIAL_WEYL_ALL_DEGREE_SPECTRAL_REDUCTION.md#the-spectral-monodromy-in-arbitrary-orders)
has \(a\) distinct leading scalars permuted in one deck orbit.
Projection of the full polar parts to their leading coefficient
is deck-equivariant and onto that orbit, so \(a\mid r\).

Here is the descent needed to read this formal calculation at
a finite ordinary \(z\). Put \(D=\operatorname{diag}(D_i)\).
The finite coefficient equations
\(BG-\partial_\tau G-GD=O(1)\) are linear in the necessary jets
of \(G\) over \(\mathbb C(z)\). They have a solution over \(K\)
with invertible constant matrix. Their solution space over
\(\mathbb C(z)\) extends to the one over \(K\) by row reduction.
The determinant of the constant matrix is therefore a nonzero
polynomial function on the solution space. Since \(\mathbb C(z)\)
is infinite, some solution over \(\mathbb C(z)\) has nonzero
determinant. Specialize these finite jets at an ordinary \(z_*\)
avoiding their poles and determinant zeros. The spectral
matrix there becomes \(D+O(1)\). Spatial transport from \(z_*\)
to \(1\) along a path in \(\mathbb C^\times\) is analytic and
invertible in \(\tau\) near zero, by the regular dependence
of \(A\), and is a horizontal identification in the spectral
direction by flatness. The same applies to the spatial leg
of the spectral deck identification.

Finally, a matrix with diagonal polar part \(D\) and regular
remainder separates formally into blocks for the distinct
exponential polynomials. To see this, successively use the
highest negative power at which polar coefficients differ to
separate the corresponding groups. At a pole of order \(m\ge2\),
division by differences of leading polar coefficients eliminates
an off-block term of order \(j\ge0\) with a gauge term of order
\(j+m\), while its derivative enters at strictly higher order.
Within each final block the whole polar part is scalar, including
the residue. After removing the scalar part the connection is
regular at \(\tau=0\), so a power series with invertible constant
term removes the remaining matrix by ordinary recursion.
Thus each repeated block is a sum of copies of one rank-one
module, and the displayed counts give the spectral types at
\(z=1\). The [degree-forty-two second splitting](POLYNOMIAL_WEYL_FORTY_TWO_SECOND_SPECTRAL_SPLITTING.md)
attains the maximum on both sides: \(v=15,d=3\) give
multiplicity five and \(r=6,21\) for \(a=2,b=7\).

## The same phase polynomial controls both presentations

**Theorem PolynomialWeylFixedMovingPhaseIdentity.** Use the
polynomial \(F\) and phases of the fixed-point theorem.
The moving exponential polynomials are exactly the distinct
polynomials
\[
                        F(\omega_i s^{k/v}),\qquad \omega_i^{av}=1,
\]
each counted with its number of occurrences. These are ordinary
polynomials in \(s\) because \(\ell=v/\gcd(k,v)\) divides every
power in \(F\). The analogous assertion holds for \(Q\).

To compare the two gauges, keep \(n=av\) and set
\(H=x^{k/v}Z\). Before choosing a spatial or spectral basepoint,
write the scalar spatial equation as
\[
      \left(f(\Theta)+\sum_{j>0}x^j f_j(\Theta)-H^n\right)y=0.
\]
Here \(\Theta\) on coefficients at fixed \(Z\) is
\(x\partial_x+(k/v)H\partial_H\). In the jets
\(H^{-j}\Theta^j y\), the spatial matrix has leading term
\(H A_{\rm sp}\) for a constant matrix with the distinct
eigenvalues \(c_i=\alpha^{-1/n}\omega_i\). Its other terms
are nonpositive in \(H\) and polynomial in \(x\).
Thus the spatial diagonal recursion produces a gauge
over \(\mathbb C[x][[H^{-1}]]\), with constant term an invertible
constant matrix. The recursive derivatives preserve polynomial
dependence on \(x\) and powers of \(H\). The diagonal spatial
coefficients are \(c_iH+O(1)\).

The matrix for \(Z\partial_Z\) at fixed \(x\) is a finite
Laurent polynomial in \(x,H\) before this gauge. Indeed jet
reduction of \(-n\lambda Q\) uses the displayed monic-in-order
spatial relation, with constant top coefficient, and the
derivative of the shear is diagonal. After gauge the matrix is
diagonal by flatness and the leading \(H\)-degree. Its diagonal
coefficients are Laurent series in \(H^{-1}\) with a uniform
lower bound on the powers of \(x\) that occur. In the flatness
equation
\[
                         \Theta B_i=(Z\partial_Z)A_i
\]
the right side has no negative \(x\)-powers. A monomial \(x^pH^q\)
on the left with \(p<0\) can therefore occur in \(B_i\) only if
\(p+(k/v)q=0\). Every such term has \(q>0\). It is precisely a
term of \(Z\partial_Z F(\omega_i Z)\) from the fixed-point
phase formula. The remaining positive \(H\)-power in \(B_i\)
is the derivative of the spatial term, \((v/k)c_iH\).

Now set \(x=z\tau,\ H=z^{k/v}\). The gauge is admissible as a
\(\tau\)-adic unit with coefficients in the spatial Laurent
field \(\mathbb C((h^{-1}))\), \(h^v=z\). Indeed the coefficient
of \(\tau^j\) is \(z^j\) times a Laurent series in \(h^{-k}\).
The same substitution works for the spectral matrix because
its \(x\)-degrees were bounded below. The moving operator
\(\tau\partial_\tau\) is
\(\Theta-(k/v)Z\partial_Z\). The two terms proportional to
\(c_iH\) cancel. The negative powers of \(\tau\) that remain
are exactly the \(\tau\partial_\tau\) derivatives of
\(F(\omega_i\tau^{-k/v})\). Finite polar descent and ordinary
spatial transport as in the moving theorem then give the stated
exponential polynomials at \(z=1\).

## Indicial classes bound algebraic spatial refinements

**Theorem PolynomialWeylIndicialClassRefinement.** The multisets
of roots of \(f\) and \(g\) in \(\mathbb C/k\mathbb Z\) are
\(a\mathcal R\) and \(b\mathcal R\) for a common multiset
\(\mathcal R\) of size \(v\). Let \(c\) be the order of its
stabilizer under the integer translations \(\mathbb Z/k\mathbb Z\).
Then \(c\mid d=\gcd(u,v)\). The rational spatial differential
modules at \(\tau=0\) have endomorphism algebras
\(\mathbb C^{ac}\) and \(\mathbb C^{bc}\).
Their primitive projectors have rank \(v/c\), and the rotation
of the spatial variable in the spectral deck acts transitively
on them. In particular the actual moving type counts satisfy
\[
                     a\mid r_P\mid ac,\qquad
                     b\mid r_Q\mid bc .
\]
This bound applies to every refinement with spectral-deck
covariant projectors that specialize to rational spatial
projectors at \(\tau=0\).

Work first with \(m=ak\) and the equation \(f(\Theta)y=z^m y\).
It descends via \(t=z^m\) to the generalized hypergeometric
equation \(f(m\theta)y=t y\), \(\theta=tD_t\).
This module is irreducible over \(\mathbb C(t)\). Indeed at
\(t=\infty\) the leading matrix in the sheared Euler jets has
\(n=av\) distinct eigenvalues on the degree-\(n\) ramification,
with distinct rank-one phases and a full deck cycle. The
fixed-point diagonal argument, now for the one spatial
connection, gives formal irreducibility.

A rank-one twist by \(t^\eta\) shifts the indicial roots
of the \(t\)-equation by \(\eta\). The twisted module is
isomorphic to the original precisely when their multisets
of roots agree modulo integers. Necessity follows from
regular-singular local monodromy at \(t=0\), whose eigenvalues
are the exponentials of the indicial roots. For sufficiency
one can move one root by one at a time. If \(L\) is
\((\theta-\rho)A(\theta)-c_0t\) with \(c_0\ne0\), and \(L_+\)
replaces \(\rho\) by \(\rho+1\), then
\[
                L_+(\theta-\rho)=(\theta-\rho-1)L .
\]
The operator \(\theta-\rho\) is invertible on solution modules,
with inverse \(c_0^{-1}t^{-1}A(\theta)\) there. Thus matching the
roots modulo integers suffices even for repeated or resonant
roots. This calculation is the specialization of the
[hypergeometric contiguity and Kummer methods of Katz, Chapter 3](https://web.math.princeton.edu/~nmk/wholebookcorrms.pdf)
to the operator with constant second polynomial.

Over \(\mathbb C(z)\) the endomorphism algebra of the pulled-back
module is graded by the characters of \(z\mapsto e^{2\pi i/m}z\).
Indeed in the \(\theta\)-jet presentation
\(\mathbb C(z)=\bigoplus_{j=0}^{m-1}z^j\mathbb C(t)\), and the
horizontal matrix equation in the \(j\)-th component is the
morphism equation to or from the \(t^{j/m}\)-twist.
Each such component is zero unless translation by \(j\)
preserves the original roots of \(f\) modulo \(m\), in which
case it is one-dimensional and generated by an invertible
endomorphism. Here irreducibility and complex constants give
the scalar endomorphisms on the base, or one can see this from
the full cycle of distinct formal types. The nonzero characters
form a cyclic subgroup of order \(R\), say. If \(J\) generates
their corresponding components by powers, then \(J^R\) is a
nonzero constant scalar. After constant normalization
\(J^R=1\). The algebra is \(\mathbb C[J]\cong\mathbb C^R\),
with primitive idempotents cyclically permuted and hence of
equal rank. The pulled-back module is semisimple: the sum
of the deck conjugates of any simple submodule descends and
must be the full base module by irreducibility. Its primitive
idempotents therefore give the simple spatial summands.

To compute \(R\), the identities for the two lowest Euler
operators can be read on their root divisors. From
\(P_0^b=\mathrm{const}\,Q_0^a\), reduction modulo \(k\)
gives \(b\bar f=a\bar g\), where bars denote root multisets.
Coprimality yields the common \(\mathcal R\) of size \(v\).
From \([P_0,Q_0]=0\), reduction of the root-divisor equation
modulo \(ak\) shows that the multiset for \(f\) is invariant
under translation by \(bk\), and hence by \(k\). Thus the roots
of \(f\) modulo \(ak\) are precisely the \(a\) uniform lifts
of \(\mathcal R\). Similarly the roots of \(g\) modulo \(bk\)
are its \(b\) uniform lifts. The two stabilizer orders are
\(ac,bc\). Integer translations act freely on each orbit of
\(\mathbb C/k\mathbb Z\), so \(c\mid v\) as well as \(c\mid k\).
This also gives the asserted common rank \(v/c\).

To pass from this spatial calculation to the moving phases,
return to finite polar descent before specializing \(z\).
The gauge polynomial solving the finite equations lives over
\(\mathbb C(z)\) and has invertible constant matrix. Formal
spectral block diagonalization for the distinct polar parts
then takes place over \(\mathbb C(z)[[\tau]]\). In this gauge
flatness forces the spatial connection to preserve those
blocks. An off-block entry of the spatial matrix would have a
first nonzero \(\tau\)-coefficient. Its equation under
\(\partial_\tau\) has on the other side a nonzero scalar pole
of order at least two times that coefficient, which cannot
balance the derivative. Thus the spectral projectors are
rational in \(z\) at \(\tau=0\) and give a coarsening of the
primitive projectors of the limiting module. The spectral
deck is the spatial rotation on the limiting module.
An equivariant partition of a single cyclic orbit into \(r_P\)
parts implies \(r_P\mid ac\). The lower divisibility and
the \(Q\) formula follow as before.

For the commuting model \(R_0=x^{-k}(E)_k(E-\eta)^u\) with
\(\eta\notin\mathbb Z\), \(\mathcal R\) has one copy of every
integer class and \(u\) copies of the single class \(\eta\).
Its stabilizer is trivial, even when \(d>1\).
For the degree-forty-two splitting with \((u,v)=(6,15)\),
the known phase counts force \(c=3\), since \(r_P=6\) and
\(c\mid3\). The new test distinguishes these two indicial
situations before any higher spectral coefficient is calculated.

## An integral constraint on the indicial means

**Corollary PolynomialWeylSpectralIndicialMeanIntegrality.**
Let \(\mu_P,\mu_Q\) be the means of the roots, counted with
multiplicity, of the indicial polynomials \(f,g\) of the two
lowest Euler components. Then the common residues in the
moving ramified connections are \(\mu_P,\mu_Q\), respectively,
and
\[
 \mu_*=\mu_P-\frac{(a-1)k}{2}
       =\mu_Q-\frac{(b-1)k}{2},
 \qquad k\mu_*\in\mathbb Z .
\]

First, at fixed \(x\), the matrix \(H(x,\lambda)\) for
spectral transport in Euler jets is polynomial in \(\lambda\)
and Laurent polynomial in \(x\), by reduction of \(Q\)
modulo the spatial equation. Hence in moving jets the trace is
\[
 \operatorname{tr}B=-ak\tau^{-ak-1}
    \operatorname{tr}H(z\tau,\tau^{-ak})+
    \tau^{-1}\operatorname{tr}A .
\]
A contribution from the first term to the \(\tau^{-1}\)
coefficient can arise only from a monomial \(x^m\lambda^j\)
with \(m=ak(j+1)>0\). It is therefore a polynomial in strictly
positive powers of \(z\). The second term contributes
\(\operatorname{tr}A(z,0)=n\mu_P\), from the companion
of \(f(\Theta)-z^{ak}\). The spatial formal gauge is a unit
in \(\tau\), so it does not change the trace residue.
That residue is constant in \(z\) by the moving theorem.
The polynomial in strictly positive powers of \(z\) must vanish,
and the common residue is consequently \(\mu_P\).

Next, the entire-solution polynomial connection embeds in the
moving jet connection over \(\mathbb C((\tau))\). To show this,
let \(t\) denote the polynomial variable in the
[dual model](POLYNOMIAL_WEYL_ALL_DEGREE_SPECTRAL_REDUCTION.md#an-algebraic-model-for-the-entire-solution-connection).
In a fixed free \(\mathbb C[\lambda]\)-basis of \(\mathbb C[t]\),
the coordinates of \(t^j\) have \(\lambda\)-degree at most
\((j+C)/(ak)\) for a constant \(C\) independent of \(j\).
Indeed \(P^\vee t^{j-ak}\) has degree \(j\) with nonzero
leading coefficient for large \(j\). Reduce \(t^j\) using
\(\lambda t^{j-ak}=P^\vee t^{j-ak}\) and lower powers of \(t\);
the claimed bound follows by induction from finitely many
initial coordinate vectors. Thus for each dual basis column
\(y\) there is a fixed integer \(C_0\ge C\) for which every
coefficient of the Taylor series in \(z\) of
\(\tau^{C_0}y(z\tau,\tau^{-ak})\) is a polynomial in \(\tau\).

These series converge locally jointly in \(z,\tau\) at \((0,0)\).
In fact the scaled equation is
\[
  \left(f(\Theta)+\sum_{j>0}(\tau z)^j f_j(\Theta)
                  -z^{ak}\right)y=0
\]
with a finite sum and \(\deg f_j\leq n\). Above the last
nonnegative integer root of \(f\), the Taylor recurrence
divides by \(f(m)\), uses finitely many preceding coefficients,
and has uniformly bounded coefficients for \(\tau\) in a fixed
compact disk, since \(\deg f=n\). The finite initial coefficients
are bounded there. Induction gives a uniform exponential
bound in the Taylor index, establishing local convergence.
On \(\mathbb C^\times\) ordinary spatial transport extends
this joint analytic evaluation to \(z=1\). Dividing out
\(\tau^{C_0}\) proves meromorphic evaluation there. It is
horizontal by the chain rule, and a maximal minor is nonzero
as a Laurent series because evaluation has full column rank for
each \(\tau\ne 0\) in the disk.

A differential submodule of a direct sum of simple modules
is semisimple. Every rank-one member of the moving direct sum
has residue \(\mu_P\) modulo integers, since a Laurent scalar
gauge changes the residue by its integer valuation.
It follows that the embedded rank-\(ak\) determinant has
residue \(ak\mu_P\) modulo integers. Its original connection
matrix is polynomial in \(\lambda\), whose trace under
\(\lambda=\tau^{-ak}\) has zero residue. Hence \(ak\mu_P\)
is an integer. The \(Q\) calculation gives \(bk\mu_Q\in\mathbb Z\).

Finally the leading negative Euler components obey
\(Q_0^a=cP_0^b\) with \(c\ne0\). The mean of the roots in
the coefficient of \(P_0^b\) is \(\mu_P+(b-1)ak/2\),
because \(F(E)x^{-j}=x^{-j}F(E-j)\). The mean for \(Q_0^a\)
is \(\mu_Q+(a-1)bk/2\). Equating the means gives the
displayed common \(\mu_*\). Since \(a(a-1)\) and \(b(b-1)\)
are even, the two integrality results give
\(ak\mu_*,bk\mu_*\in\mathbb Z\). Now use \(\gcd(a,b)=1\).

For the commuting Euler model
\(R=x^{-k}(E)_k(E-\eta)^u,\ P_0=R^a,\ Q_0=R^b\),
one has \(\mu_*=(k(k-1)/2+u\eta)/v\).
For \(u=6,v=15,k=9\), the choice \(\eta=1/2\) gives
\(k\mu_*=117/5\) and is excluded as the leading negative
Euler data of an actual Weyl pair. The choice \(\eta=2/3\)
gives \(k\mu_*=24\) and passes this condition. Both choices
meet the [homogeneous monodromy rank and block conditions](POLYNOMIAL_WEYL_ALL_DEGREE_SPECTRAL_REDUCTION.md#the-obstruction-required-of-a-general-proof).
Thus the mean constraint adds information beyond those
conditions while leaving some homogeneous models.

## Exact indicial means and the mixed Weyl coefficient

**Theorem PolynomialWeylExactIndicialMean.** For an actual centered Weyl
pair the means of the two lowest Euler polynomials satisfy
\[
             \mu_P=\frac{ak-1}{2},\qquad
             \mu_Q=\frac{bk-1}{2}.
\]
In particular \(\mu_*=(k-1)/2\). If the coefficient of
\(x^{au-1}D_x^{av-1}\) in \(P\) is \(\gamma\), then
\(\gamma=\alpha(au)(av)/2\). The corresponding coefficient of \(Q\)
is \(\beta(bu)(bv)/2\). Equivalently the mixed subprincipal
coefficients in symmetric Weyl ordering vanish. One translation of
\(D_x\) then makes the full symmetric subprincipal symbols of both
operators zero.

Write \(m=ak\) and \(n=av\). For each class \(r=0,\ldots,m-1\)
choose the largest nonnegative integer root \(\nu_r\) of \(f\)
congruent to \(r\) modulo \(m\). There is at least one in each
class because \((E)_m\mid f\). In the limiting spatial equation
\(f(\Theta)y=z^m y\) the derivatives of orders \(\nu_r\) at zero
give coordinates on the entire solutions. In each class the
Taylor recurrence \(f(j)y_j=y_{j-m}\) forces coefficients before
the last root to vanish. The coefficient at the last root is
free and all later coefficients are determined. The local Taylor
solutions converge by the same recurrence and extend to entire
solutions since there are no other finite singularities.

Consider the free polynomial module \(\mathbb C[t]\) with
\(\lambda=P^\vee\). Its submodule \(M'\) generated by the
\(t^{\nu_r}\) over \(\mathbb C[\lambda]\) is free of rank \(m\).
Indeed the coefficient of \(t^{q+m}\) in \(P^\vee t^q\) is
\(q!f(q+m)/(q+m)!\) by the derivative pairing, and all other
terms have smaller degree. It follows that the elements
\((P^\vee)^j t^{\nu_r}\), \(j\ge0\), have distinct degrees
\(\nu_r+jm\). Using these elements for leading-term division,
the quotient \(\mathbb C[t]/M'\) has basis the omitted degrees.
Its dimension is
\[
                     L=\sum_{r=0}^{m-1}\frac{\nu_r-r}{m}.
\]
Thus the coordinate determinant \(\Delta(\lambda)\) of the
\(t^{\nu_r}\) in any free polynomial basis has degree \(L\),
by the Smith form over \(\mathbb C[\lambda]\).

Take the dual basis and evaluate its entire sections at
\(x=z\tau,\ \lambda=\tau^{-m}\) as above. In the resulting \(m\)
columns the matrix of selected derivatives at \(z=0\) has
determinant
\[
            \tau^{\sum_r\nu_r}\Delta(\tau^{-m}).
\]
Its \(\tau\)-valuation is \(\sum_r r=m(m-1)/2\). This is also
the valuation of the exterior product of the columns in
ordinary Euler jets at \(z=1\). In detail, normalize the columns
by the inverse of the selected-derivative matrix. They are
meromorphic in \(\tau\) with coefficients holomorphic in \(z\)
near zero by the earlier joint convergence. If a normalized
column had a pole in \(\tau\), its leading coefficient would be
a solution of the limiting spatial equation whose selected
derivatives all vanish. The coordinate property just proved
excludes this. The normalized columns are therefore holomorphic
and specialize to an entire basis at \(\tau=0\). Transport to
\(z=1\) preserves this independence, so their exterior product
has valuation zero. Undoing the normalization gives the stated
valuation.

Use now the moving spectral gauge at \(z=1\), which is a power
series unit in \(\tau\). In this gauge the rank-one factors
have a common residue \(\mu_P\). The exterior product \(W\) of
the unnormalized evaluation columns satisfies
\[
 \partial_\tau W
      =\bigl(B_{\wedge^m}-\operatorname{tr}H_{\rm ent}\bigr)W ,
\]
where \(H_{\rm ent}\) is the polynomial source connection
pulled back from \(\lambda\) to \(\tau\). Its trace has zero
\(\tau\)-residue. Each nonzero component of \(W\) in the exterior
basis is a Laurent scalar whose logarithmic derivative
therefore has residue \(m\mu_P\). The minimum valuation of the
components is unchanged by a power series unit and is the
valuation computed above. Consequently
\(m\mu_P=m(m-1)/2\), giving the \(P\) formula. The same argument
applies to \(Q\).

Finally the top two terms of \(f(E)\) come from
\(\alpha(E)_n+\gamma(E)_{n-1}\), so
\(\mu_P=(n-1)/2-\gamma/(\alpha n)\). This proves the claimed
coefficient. The same subtraction is the Weyl-symbol correction
\(\gamma\mapsto\gamma-\alpha(au)(av)/2\) at the mixed corner.

For the full next differential coefficient, the rectangle and
the diagonal-gap bound allow only the powers \(x^{au-1}\) and
\(x^{au}\), and similarly for \(Q\). Write \(p,q\) for the
leading differential coefficients, \(n,N\) for the orders, and
write the next coefficients as \((n/2)p'+C_Pp\) and
\((N/2)q'+C_Qq\). In the order \(n+N-2\) coefficient of
\([Q,P]\), the terms containing \(C_P,C_Q\) are
\(C_Ppq'-C_Qqp'\), using \(Nqp'=npq'\) from the previous
order. All other terms have spatial degree at most
\(au+bu-2\). Since the commutator is scalar, it follows that
\(C_Pbu=C_Qau\), or \(C_P/n=C_Q/N\).
Translation \(D_x\mapsto D_x-C_P/n\) removes both constants.

In the commuting model \(R=x^{-k}(E)_k(E-\eta)^u\) discussed
above, exact equality requires \(\eta=(k-1)/2\). For \(k=9\)
the earlier integrality test accepts \(\eta=2/3\), whereas exact
equality excludes it, and indeed excludes every noninteger
\(\eta\) in this single-parameter model. The [conditional
subprincipal comparison in the Zheglov audit](ZHEGLOV_DC1_PRIMARY_SOURCE_AUDIT.md)
reaches the same Weyl coefficient under its later bridging
hypotheses. The determinant proof here uses the centered Weyl
pair and polynomial entire-solution module directly.

The same commuting powers permit the exact mean when their primitive
Euler polynomial has several extra roots. For \((u,v,k)=(6,15,9)\)
take the integer roots \(0,\ldots,8\) and six extra roots. Three
copies each of \(7/2,9/2\) have sum \(24\) and trivial stabilizer
modulo \(9\). The list \(1/2,7/2,13/2,3/2,9/2,15/2\) also
has sum \(24\), and its stabilizer has order three. Both
commuting models have \(\mu_*=4\), so the exact mean and
the indicial refinement still leave homogeneous leading data.

## The comparison required for DC1

At a fixed ordinary point all \(av\) types are distinct,
but for \(x=z/s\) the fixed-point term that separates them
becomes independent of \(s\): on a common cover,
\[
                          x^{k/v}Z=z^{k/v}.
\]
The moving theorem uses a separate spatial formal calculation
and finite polar descent. Substitution into an infinite
fixed-point Laurent gauge could mix its orders because its
coefficients have poles at \(x=0\).

The [entire-solution polynomial connection](POLYNOMIAL_WEYL_ALL_DEGREE_SPECTRAL_REDUCTION.md#an-algebraic-model-for-the-entire-solution-connection)
has rank \(ak\), strictly between zero and \(av\). Evaluation
into the full jets at a fixed ordinary \(x_*\) is horizontal
and injective on every finite spectral fibre by uniqueness of
the spatial initial-value problem. If it induced a nonzero
horizontal morphism of unramified formal differential modules
at \(\lambda=\infty\), the image would be a proper nonzero
submodule of the fixed-point formal connection, contradicting
its irreducibility. A rational horizontal comparison would
in particular give such a morphism. Finite-parameter
analyticity of evaluation supplies no formal Laurent comparison.
Producing one under the centered Weyl identities would exclude
the hypothetical counterexample and prove DC1.


The obstruction persists after any finite ramification. On a common
ramified cover the entire-solution connection embeds in the
moving direct sum and therefore has rank-one factors with
exponentials polynomial in \(s=\lambda^{1/(ak)}\).
No such polynomial has a term in \(\lambda^{1/(av)}\), since
\(0<ak<av\). At a nonzero fixed \(x\) every full-jet phase
has a nonzero term at precisely this exponent. No two of these
types, one from each presentation, can be isomorphic on a
further finite cover. Thus there is no nonzero formal horizontal
morphism between the entire-solution and fixed-jet connections
even after ramification.


### A growth criterion for the fixed-point comparison

**Corollary PolynomialWeylFixedEvaluationGrowth.** For any fixed ordinary
\(x_*\), the entire-solution evaluation matrix for \(P\), in polynomial
coordinates on the source and Euler jets on the target, satisfies a bound
\(C(1+|\lambda|)^N\exp(C|\lambda|^{1/(av)})\).
If a nonzero polynomial section of the source has a polynomially bounded
evaluation on even one ray in \(\lambda\), the centered pair is excluded
when \(av>2\). The \(Q\) criterion uses \(bv\) in place of \(av\).

In the moving chart the analytic evaluation proved above has a uniform
polynomial bound in \(|\lambda|\) at \(|z|=1\), because its pole order in
\(\tau\) is bounded and the circle in \(z\) is compact. For large
\(|\lambda|\), start at the point \(x\) of modulus
\(|\lambda|^{-1/(ak)}\) on the ray of \(x_*\) and transport to \(x_*\)
along that ray. Put \(r=|x|\) and \(R=|\lambda|^{1/n}r^{k/v}\) with
\(n=av\), so \(R\ge1\) along the segment. In the jets
\(R^{-j}\Theta^j y\), the equation for \(r\,d/dr\) has coefficients
bounded by \(C(1+R)\). Indeed the upper companion diagonal is \(R\),
the bottom term from \(\lambda x^{ak}\) has modulus
\(R^n/R^{n-1}=R\), and the other bottom coefficients are bounded:
they come from the polynomials \(f_j\) on a bounded \(x\)-segment
with shear powers \(R^{j-(n-1)}\), \(j\le n-1\).
The derivative of the shear adds only bounded diagonal entries
since \(r(\log R)'=k/v\). Integrating the matrix norm against
\(dr/r\) gives \(O(\log|\lambda|+|\lambda|^{1/n})\).
The initial evaluation and the final change of jets cost only
polynomial factors. This proves the bound, also for every polynomial
section of the entire-solution connection.

Suppose the evaluation of such a section is polynomially bounded
on a ray. Each entry is entire in \(\lambda\). To see that it is a
polynomial when \(n>2\), divide by \(\lambda^M\) for a sufficiently
large integer \(M\) on the slit exterior of the unit disk, using
the given ray as the slit. Choose \(1/n<\rho<1/2\) and take the
branch of the rotated \(\lambda^\rho\) with arguments
between \(-\pi\rho\) and \(\pi\rho\). Multiplication by
\(\exp(-\epsilon\lambda^\rho)\) in that branch suppresses the
outer circular boundary uniformly as its radius grows.
The two edges of the slit and the unit circle have bounded values
uniformly for \(0<\epsilon\le1\). The maximum principle and then
\(\epsilon\downarrow0\) give a global polynomial bound.
Cauchy's estimates give polynomial entries.

Covariant derivatives of this polynomial evaluation are again
polynomial, since the target spectral matrix at \(x_*\) is
polynomial in \(\lambda\). The source section and its successive
covariant derivatives generate a nonzero differential
submodule over \(\mathbb C(\lambda)\) of rank at most \(ak\).
Evaluation is horizontal, so its restriction to this submodule
is a nonzero rational differential morphism by the polynomial
images of the generators. This contradicts the fixed-point
phase obstruction. Thus a one-ray polynomial bound for one
nonzero polynomial section would prove DC1 if it could be
established under the centered Weyl hypotheses for one side
of every hypothetical counterexample. The diagonal gap bound
ensures \(av,bv>2\) throughout that reduction.

In the moving presentation a spectral decomposition by
distinct exponential types always has repeated isotypic blocks.
The decomposition alone cannot supply a preferred rank-one
refinement respected by the finite spatial monodromy and the
spectral deck. A refinement with one orbit on all \(av\) lines
would contradict the previously proved fixed-rank divisibility.
