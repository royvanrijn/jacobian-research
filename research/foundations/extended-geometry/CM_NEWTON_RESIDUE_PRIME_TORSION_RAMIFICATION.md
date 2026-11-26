# Ramification of residue-prime torsion in the split real-cubic case

Let \(A/K\) be an abelian threefold with full geometric rational
endomorphism algebra a totally real cubic field \(F\), acting over \(K\).
Work at a good absolute degree-one place \(v\) of characteristic
\(p>3\), split completely in \(F\) and prime to the endomorphism-order
index and a polarization degree. Exclude the finitely many ramified
rational primes of \(K\), so \(K_v=\mathbb Q_p\). Write the relative
Frobenius polynomial as \(X^2-a_v X+p\). The integral coefficient
idempotents split the local \(p\)-divisible group into three factors
of height two and dimension one. Write \(u_i\in\mathbb F_p^\times\)
for the reduction of \(a_v\) at an ordinary coefficient place. This is
the Frobenius eigenvalue on its étale group of order \(p\), since the
other root of the relative polynomial reduces to zero.
[Split local Newton data](CM_NEWTON_FROBENIUS_RAMIFICATION_AND_DENOMINATORS.md#inertia-and-the-frobenius-weyl-group).

## The full torsion field

**Theorem (Residue-prime ramification profile).** Let
\(J_v=K_v(A[p])\). The tame inertia order is \(m=p-1\) at an
ordinary place and \(m=p^2-1\) if any factor is nonordinary. The tame
inertia order over \(K_v(\mu_p)\) is therefore \(1\) or \(p+1\),
respectively. The wild inertia order \(w=p^r\) has \(r\) equal to the
number of distinct \(u_i^2\) among the ordinary factors whose local
mod-\(p\) extension of étale by connected points is nonsplit.
If \(w>1\), the sole positive upper ramification break over \(K_v\)
is \(1/(p-1)\). For all \(w\), the normalized different exponent is
\[
 \frac{v_{J_v}(\mathfrak D_{J_v/K_v})}{e(J_v/K_v)}
   =1-\frac1{mw}+\frac{w-1}{w(p-1)}
   =1+\frac1{p-1}-\frac1w\left(\frac1m+\frac1{p-1}\right).
\]
Here \(v_{J_v}\) takes value one on uniformizers of \(J_v\).

For an ordinary factor the connected and étale characters on generic
points are \(\bar\epsilon\alpha_i^{-1}\) and \(\alpha_i\), respectively,
where \(\alpha_i\) is unramified with arithmetic Frobenius value \(u_i\),
and \(\bar\epsilon\) is the mod-\(p\) cyclotomic character. Thus its
tame characters are \(\bar\epsilon\) and \(1\). For a nonordinary
factor the full group of order \(p^2\) comes from a one-dimensional
formal group of height two. Serre's formal group calculation over
absolute ramification index one identifies its \(p\)-torsion with a
one-dimensional \(\mathbb F_{p^2}\)-space on which inertia acts by
the fundamental character of order \(p^2-1\). Its norm is
\(\bar\epsilon\). The tame inertia quotient is procyclic, so its
order in the product is the least common multiple of these orders.
The cyclotomic assertion follows by taking the kernel of that norm
on the finite tame quotient.
[Serre, §1.9, Proposition 9, and §1.11, Propositions 11–12](https://www.college-de-france.fr/media/jean-pierre-serre/UPL5874918517843398173_Serre_proprie_te_s_galoisiennes_des_courbes_elliptiques.pdf).

For the wild image, take \(U/K_v\) unramified of degree dividing
\(p-1\) that kills all of the ordinary unramified characters.
Over \(O_U\) the connected and étale group schemes of order \(p\)
are \(\mu_p\) and the constant \(\mathbb Z/p\), respectively.
The inverse image of a generator of the étale quotient is a
\(\mu_p\)-torsor. Its class in \(O_U^\times/O_U^{\times p}\)
specifies the generic extension class by Kummer theory, because an
extension of vector spaces killed by \(p\) splits as soon as a lift
of the generator is fixed. The unit quotient is
\((1+p O_U)/(1+p O_U)^p \simeq O_U/pO_U\), using the \(p\)-adic
logarithm, or the first unit digit. Restriction to \(U\) is injective
on the original extension classes since \([U:K_v]\) is prime to \(p\).
A finite-flat ordinary extension of two rank-one group schemes is
being used here. For the elliptic case the same Kummer description
appears in [David–Weston, Lemma 3.2](https://arxiv.org/pdf/math/0701882).

Before restriction the coefficient character of the extension
class was \(\bar\epsilon\alpha_i^{-2}\). Its restricted Kummer class is
therefore a \(u_i^2\)-eigenvector for the natural Frobenius on
\(O_U^\times/O_U^{\times p}\). By a finite-field normal basis the
Frobenius on \(O_U/pO_U\) has each \([U:K_v]\)-th root of unity in
\(\mathbb F_p\) with multiplicity one. The span of the nonzero
classes has dimension \(r\) as asserted. After adjoining \(\mu_p\),
Kummer restriction is still injective, because this base change
also has degree prime to \(p\). The wild extension has group
\((\mathbb Z/p)^r\), as the different calculation below shows that
all of its nontrivial degree-\(p\) subextensions are ramified.
The nonordinary factors add only tame ramification.

For the break, represent any nonzero unit class by \(q=1+pb\)
with \(b\in O_U^\times\). The polynomial \((1+T)^p-q\) is Eisenstein
over \(U\). Its root field \(U_1\) has degree \(p\), and the
derivative at the uniformizer root has valuation \(p\) in \(U_1\).
Thus its different over \(U\) has exponent \(p\). Put
\(C=U(\mu_p)\) and \(V=U_1 C\). The degree-\(p-1\) extension
\(V/U_1\) is totally tamely ramified because \(v_{U_1}(\zeta_p-1)
=p/(p-1)\), where here the valuation is extended to an algebraic
closure. The tower through \(U_1\) gives
\[
 v_V(\mathfrak D_{V/U})=(p-1)p+(p-2).
\]
The tower through \(C\) then gives
\(v_V(\mathfrak D_{V/C})=2(p-1)\). The extension \(V/C\) is
totally ramified cyclic of degree \(p\), so the different formula
gives the upper and lower break \(1\). This applies to every cyclic
quotient of the elementary abelian Kummer extension over \(C\).
Herbrand's quotient theorem now gives a single upper break \(1\)
for that entire extension when \(r>0\), and it is totally ramified.
Lower numbering is compatible with subgroups. On viewing the
ordinary torsion extension over \(U\), its positive lower groups
agree with those over \(C\), so the corresponding upper break is
\(1/(p-1)\). Adding the nonordinary tame factors preserves the
positive upper filtration by the quotient theorem, since the
additional kernel has order prime to \(p\). Unramified base change
preserves the inertia groups and filtration. In the full torsion
extension the positive lower break is consequently \(m/(p-1)\)
if \(w>1\). The zeroth group contributes \(mw-1\) to the different
exponent and the remaining groups contribute
\((m/(p-1))(w-1)\). The same formula covers \(w=1\).
[Serre, Local Fields, Chapter III §4 and Chapter IV §§1 and 3](https://link.springer.com/book/10.1007/978-1-4757-5673-9).

## Coincidences of the ordinary residue eigenvalues

**Proposition (Low-weight residual coincidences).** Fix a nonzero
\(P(S,T)\in\mathbb Q[S,T]\) of weighted degree at most two,
where \(S\) has weight one and \(T\) weight two. Among the split
absolute degree-one places with exactly one nonordinary factor, the
places with
\[
                     P(u_1+u_2,u_1u_2)=0\quad\text{in }\mathbb F_p
\]
have relative density zero in the full split degree-one class.
Here \(u_1,u_2\) are the two ordinary residue eigenvalues.
In particular, the equality \(u_1^2=u_2^2\) has relative density
zero. Outside this exceptional locus the wild rank at a one-failure
place equals the number of nonsplit ordinary mod-\(p\) extensions.

Write \(s_1=\sum_i\sigma_i(a_v)\) and
\(s_2=\sum_{i<j}\sigma_i(a_v)\sigma_j(a_v)\), which are rational
integers. The zero residue trace at the failed factor means
\(s_1\equiv u_1+u_2\) and \(s_2\equiv u_1u_2\bmod p\).
Clear the fixed denominators of \(P\). The Weil bounds give
\(|s_1|\le6\sqrt p\) and \(|s_2|\le12p\), so divisibility forces an
exact equation \(P(s_1,s_2)=h p\) for an integer \(h\) in a finite
range depending only on \(P\). For a sufficiently large split
auxiliary prime, Lombardo's theorem gives the full rank-two
common-determinant image. On its connected algebraic group the
three traces are independent coordinates at fixed nonzero
determinant, as one sees by taking three companion matrices.
Thus for each \(h\), the equation
\(P(s_1,s_2)=h\det\) gives a proper algebraic subset. Each such
subset has measure zero on the compact open auxiliary image, and
finite-level Chebotarev, by approximating the zero locus in that
image, gives density zero for its Frobenius primes. The split class
has positive density, so the finite union also has relative
density zero in that class.
[Lombardo, Theorem 1.4](https://arxiv.org/pdf/1411.1703).

Apply the proposition to \(P=S\) and \(P=S^2-4T\) for opposite and
equal ordinary eigenvalues, respectively. A fixed nonzero rational
residue \(c\) for either ordinary eigenvalue likewise gives
\(P=T-cS+c^2\). For example a \(p\)-division point defined over
a local extension of degree at most two forces \(u_i=\pm1\) at
some ordinary factor. At a nonordinary factor every nontrivial
inertia scalar acts without a fixed nonzero vector, and restricting
to an extension of degree at most two keeps nontrivial inertia
since \(p^2-1>2\). In an ordinary factor the connected character
has inertia order \(p-1>2\), so a local point must have nonzero
étale image. Its Frobenius order divides the local residue degree,
which is at most two. Thus local \(p\)-torsion over any extension
of degree at most two cuts out a relative density-zero part of the
one-failure locus.

## Resolution needed from global discriminants

The exact normalized relative discriminant exponent of \(K(A[p])/K\)
at \(v\) equals the normalized local different in the theorem,
since for a Galois extension the residue degrees and number of places
cancel after dividing by the degree of the extension. For fixed
wild rank \(r\), the value with a nonordinary factor exceeds the
all-ordinary value by
\[
                              \frac{p}{p^r(p^2-1)}.
\]
The two values lie strictly below the all-ordinary value for
wild rank \(r+1\) when \(p>3\). Thus the exact coefficient at \(v\)
determines both the wild rank and whether any factor is nonordinary.
It gives a criterion at the residue prime itself, including cases
with full global image on \(A[p]\).

For a compositum of fields \(K(A[\ell])\) at distinct rational
primes \(\ell\), the other fields are unramified at a good place
of characteristic \(p\) by the good-reduction criterion on
prime-to-\(p\) torsion. They therefore leave this coefficient at
\(v\) unchanged. Its contribution to the logarithm of the absolute
root discriminant has weight \(\log p/[K:\mathbb Q]\).
At a fixed wild rank, changing the Newton type changes the
coefficient by \(O(p^{-1-r})\). Summed over \(p\le X\), all
these differences contribute at most \(O_K(\log X)\) to the
log root discriminant. The part where \(r\ge1\) is bounded
even as \(X\to\infty\), by summability of \((\log p)/p^2\).
The normalized discriminant is an exact local criterion, but its
leading aggregate asymptotic cannot separate the two cases.
[Serre–Tate, Good reduction of abelian varieties, Theorem 1](https://annals.math.princeton.edu/1968/88-3/p05).

## A degree-\(2(p+1)\) field for each factor

**Proposition (Projective line and determinant discriminants).**
Suppose one factor has full projective mod-\(p\) image
\(G=\operatorname{PGL}_2(\mathbb F_p)\), as it does for all
sufficiently large split \(p\) by Lombardo. In its projective
division field let \(K_1\) be fixed by a line stabilizer \(B\)
and \(K_2\) by \(B\cap\operatorname{PSL}_2(\mathbb F_p)\).
Then \([K_1:K]=p+1\) and \([K_2:K]=2(p+1)\).
Their relative discriminant exponents at \(v\) are
\[
\begin{array}{c|cc}
 &K_1/K&K_2/K\\ \hline
 \text{ordinary with trivial wild image}&p-2&2p-2\\
 \text{ordinary with nontrivial wild image}&p&2p+2\\
 \text{nonordinary}&p&2p
\end{array}
\]
In particular \(K_2\) determines the local type of this factor.

Write \(\Pi\) for the line permutation representation and \(\chi\)
for the quadratic determinant character of \(G\). The index-two
subgroup \(B\cap\ker\chi\) of \(B\) gives the permutation
representation \(\Pi\oplus(\Pi\otimes\chi)\) by induction in stages.
Thus the two relative discriminant exponents are the Artin conductors
of \(\Pi\) and \(\Pi\oplus(\Pi\otimes\chi)\), respectively.
For a local inertia subgroup, an orbit contributes a one-dimensional
invariant space in the twisted permutation representation precisely
when \(\chi\) is trivial on its point stabilizer. This is Frobenius
reciprocity on the permutation module of the orbit.
[Serre, Local Fields, Chapter VI §3, conductor-discriminant formula](https://link.springer.com/book/10.1007/978-1-4757-5673-9).

For tame ordinary inertia, the projective image is the split torus
represented by \(\operatorname{diag}(t,1)\). Its orbits on the line
space are the two eigenlines and their complement. The latter is a
regular orbit, and the character \(\chi\) is nontrivial on the torus
fixing each eigenline. The invariant dimensions are therefore
three for \(\Pi\) and one for \(\Pi\otimes\chi\), giving conductors
\(p-2,p\). For wild ordinary inertia the projective inertia group
is the affine group with translation subgroup of order \(p\).
Its orbits are the fixed line and its complement. The stabilizer
of the fixed line is the affine group itself, and the stabilizer
of a line in the other orbit is conjugate to the full split torus.
Thus both stabilizers have nontrivial \(\chi\). The inertia
invariant dimensions are two and zero. Wild inertia has the same
two orbits and trivial \(\chi\), so the wild invariant dimension
is two for both representations. The Swan term is
\((p-1)/(p-1)=1\) in either case by the single positive upper
break in the theorem. This gives conductors \(p,p+2\).

In the nonordinary case the projective image of inertia is the
nonsplit torus \(\mathbb F_{p^2}^\times/\mathbb F_p^\times\), acting
regularly on the set of \(\mathbb F_p\)-lines in \(\mathbb F_{p^2}\).
Its determinant modulo squares is the norm modulo squares, which
is surjective because the norm is surjective and takes a scalar
to its square. Both representations have one inertia invariant
from the regular orbit, and wild inertia is trivial. Their
conductors are both \(p\). Adding the two conductors in each case
gives the second column. The degree-\(p+1\) field has equal
exponents in the wild ordinary and nonordinary cases. The exact
coefficient of the larger field distinguishes them, but the
coefficient divided by its degree differs from \(p/(p+1)\) by
at most \(1/(p+1)\). The calculation supplies no estimate for the
frequency of any of these types as \(p\) varies.

## Small residue degrees in the full torsion field

**Proposition (Residue degree and a growing exclusion).** The residue
degree \(h_v\) of \(J_v/K_v\) is the least common multiple of the
orders of the ordinary eigenvalues \(u_i\), together with \(2\) if
there is any nonordinary factor. Among the split degree-one places
with exactly one nonordinary factor, the places with
\(h_v\le p^\theta\) have relative density zero for any fixed
\(0\le\theta<1/6\). After discarding the five exact trace loci
\(s_1\in\{-2,-1,0,1,2\}\), their count up to \(X\) is
\(O_{K,\theta}(X^{1/2+3\theta})\).

The unramified étale character of each ordinary factor gives the
corresponding divisor \(\operatorname{ord}(u_i)\) of \(h_v\).
For a nonordinary factor, Serre's semilinear description identifies
its decomposition image with the normalizer of the nonsplit
Cartan at residue degree one. Indeed inertia is the full Cartan,
and arithmetic Frobenius acts nontrivially on \(\mathbb F_{p^2}\)
by the \(p\)-power map. Its quotient by inertia has order two.
This gives the divisor \(2\) of \(h_v\) when present.
[Serre, §1.9, Proposition 9(i), and §1.11, Proposition 12(d)](https://www.college-de-france.fr/media/jean-pierre-serre/UPL5874918517843398173_Serre_proprie_te_s_galoisiennes_des_courbes_elliptiques.pdf).

Now let \(U/K_v\) be unramified of degree the stated least common
multiple, using degree one for the empty list. Over \(C=U(\mu_p)\)
the joint ordinary torsion field is the unit Kummer extension in
the theorem, so is totally ramified of \(p\)-power degree. If there
are several nonordinary factors, they define the same local
division field. To see this, choose coordinates on their
\(\mathbb F_{p^2}\)-lines so inertia acts by the same fundamental
character, conjugating the field coordinates where needed.
A fixed arithmetic Frobenius lift acts as \(z\mapsto c_i z^p\)
on each line. Its determinant is \(-N_{\mathbb F_{p^2}/\mathbb F_p}(c_i)\),
because \(z\mapsto z^p\) has eigenvalues \(1,-1\).
The common cyclotomic determinant makes these norms equal. The
ratios \(c_i/c_j\) are therefore of the form \(b/b^p\) in
\(\mathbb F_{p^2}^\times\). Conjugating by multiplication by \(b\)
identifies the Frobenius actions and commutes with inertia. Hence
the residual local representations are isomorphic and have the
same kernel.

When present, this one nonordinary division field after base
change to \(U\) is totally tame, since the base degree is even
and the residue quotient before base change had order two.
Its tame degree over \(C\) is \(p+1\), prime to the degree of the
ordinary Kummer field there. The ramification index of the
compositum is a multiple of both indices and hence equals the
full degree. The compositum over \(U\) is therefore totally
ramified, which gives the opposite divisibility for \(h_v\)
and proves the formula.

For the growing exclusion use \(u_1,u_2\) for the two ordinary
eigenvalues. On \(X/2<p\le X\), put \(d=h_v\le X^\theta\) and
introduce the monic integer polynomial
\[
 Q_d(Z)=\prod_{\zeta^d=\eta^d=1}(Z-\zeta-\eta).
\]
Its coefficients are Galois-invariant algebraic integers and its
degree is \(d^2\). Since \(d\mid p-1\), the two eigenvalues lift
to \(d\)-th roots of unity over \(\mathbb Z_p\). Their sum
reduces to \(s_1\bmod p\), so \(p\mid Q_d(s_1)\).
If \(Q_d(s_1)=0\), the rational integer \(s_1\) is a sum of
two complex roots of unity and has absolute value at most two.
The five exact trace loci have density zero by the same algebraic
Chebotarev argument as above, with equations \(s_1=j\). If
\(Q_d(s_1)\ne0\), then
\[
                   |Q_d(s_1)|\le (6\sqrt X+2)^{d^2}.
\]
At most \(O(d^2)\) distinct primes in \((X/2,X]\) divide such
an integer, since their product is bounded by its absolute value.
Summing over \(O(\sqrt X)\) possible integers \(s_1\) and
\(d\le X^\theta\) gives \(O(X^{1/2+3\theta})\) places in the
interval, using at most \([K:\mathbb Q]\) places above each
rational prime. Summing the dyadic intervals proves the count.
Since \(1/2+3\theta<1\), it is \(o(X/\log X)\).
