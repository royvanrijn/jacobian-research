# Low-degree totally real lifts from quadratic norm-one traces

Every abc triple has a totally real lift of degree linear in its
logarithmic radical, with bounded normalized cusp radical. Its
log root discriminant recovers the rational logarithmic radical
after removing the contribution of at most one explicit auxiliary
prime. The lift comes from a single norm-one trace over a real
multiquadratic field of bounded degree. Its solvable Galois closure
is a compositum of trace splitting fields. The possible auxiliary
prime is the least quadratic nonresidue at a single prime on a
summand. For a fixed odd level, a correction prime removes the
class number restriction. Comparing two such lifts and their
compositum cancels their auxiliary discriminant contributions.

## Splitting a prime support using at most eight quadratic fields

Let \(S\) be any finite set of odd primes. There are at most eight
imaginary quadratic fields such that each \(p\in S\) splits in
at least one of them. Their discriminants are supported on
\(S\cup\{2\}\) except possibly for one additional prime \(q\).
The exception arises only if there is a prime \(p_0\in S\)
with \(p_0\equiv3\bmod4\) at which every prime in
\((S\setminus\{p_0\})\cup\{2\}\) is a quadratic residue.
In that case \(q\) may be taken to be the least positive
quadratic nonresidue modulo \(p_0\), which is prime. Put
\(q=1\) if no additional prime is needed.

Use \(\mathbb Q(i)\) for \(p\equiv1\bmod4\). For the set
\(S_-\) of primes \(p\equiv3\bmod4\), let \(A_{pp'}\) be
the indicator that \(p'\) is a quadratic nonresidue modulo
\(p\) for \(p\ne p'\), and put \(A_{pp}=0\). Work over
\(\mathbb F_2\). By quadratic reciprocity,
\(A_{pp'}+A_{p'p}=1\) for distinct indices. Define
\[
              B=A+\operatorname{diag}(A\mathbf1).
\]
For a vector \(x\in\mathbb F_2^{S_-}\), one has
\((Bx)_p=\sum_{p'}A_{pp'}(x_{p'}+x_p)\). If this equals
one, the product \(m\) of the primes on the side of the
cut opposite \(p\) is a nonresidue modulo \(p\). Thus
\(p\) splits in \(\mathbb Q(\sqrt{-m})\). Each vector \(x\)
requires at most two imaginary quadratic fields, from the
products on its two sides.

There is a set of three such vectors that covers every
nonzero row of \(B\). For a subset \(U\) of \(u\) row
indices, the symmetrization of its principal submatrix
is \(I_u+J_u\), where \(J_u\) is all ones. This has rank
\(u\) for even \(u\) and \(u-1\) for odd \(u\). As
\(\operatorname{rank}(B_U+B_U^\top)\le2\operatorname{rank}(B_U)\),
the full rows indexed by \(U\) have rank at least
\(\lfloor u/2\rfloor\). On the set of nonzero rows every
subset has size at most three times its rank, including
singletons. Edmonds' matroid partition theorem therefore
partitions these rows into three independent sets.
On each set solve the linear system \((Bx)_p=1\).
The solutions give the required cuts and at most six
quadratic fields.
[Edmonds, Theorem 1](https://nvlpubs.nist.gov/nistpubs/jres/69B/jresv69Bn1-2p67_A1b.pdf).

A zero row has no nonresidue among the other members of
\(S_-\). There is at most one such row, since the two
opposite entries for distinct primes sum to one. If it
has a nonresidue \(q'\) in \(S\cup\{2\}\), cover it with
\(\mathbb Q(\sqrt{-q'})\). Otherwise it is the exceptional
\(p_0\), and \(\mathbb Q(\sqrt{-q})\) covers it. No odd prime
outside the listed supports is needed in the discriminants.
Counting the first field, the six cut fields and the
possible last field proves the bound eight.

For the odd support of a primitive positive \(a+b=c\), an
exceptional \(p_0\) lies on \(a\) or \(b\).
Indeed if \(p_0\mid c\), then \(a\equiv-b\bmod p_0\)
and \(a,b\) are units there. Their quadratic symbols are
opposite, whereas if every prime factor in \(ab\) is a
residue they are both one.

## A totally real trace field with prime-order inertia

Let \(J\) be one of the imaginary quadratic fields, and
assign it a nonempty set \(S_J\subset S\) of primes that
split there. Pick one prime ideal above each of these
primes. Let \(\mathfrak a\) be their product, let \(H_J\)
be the class number, choose \(\beta\) with
\((\beta)=\mathfrak a^{H_J}\), and put \(w=\beta/\bar\beta\).
Then \(w\bar w=1\). At each \(p\in S_J\) its two valuations
are \(H_J\) and \(-H_J\), and it is a unit at all other
finite primes of \(J\).

Let \(\ell\) be an odd prime coprime to
\(H_J\operatorname{disc}(J)\prod_{p\in S_J}p\). For
\(z^\ell=w\) and a primitive \(\ell\)-th root of unity
\(\zeta\), take the traces
\[
                  x_i=\zeta^i z+\zeta^{-i}z^{-1}.
\]
These are the roots of \(P_\ell(x)-(w+w^{-1})\), where
\(P_\ell(z+z^{-1})=z^\ell+z^{-\ell}\). They are distinct
and real: all \(\zeta^i z\) are on the unit circle, and
equality between distinct traces would require
\(w=w^{-1}\), incompatible with its nonzero valuations.
At each \(p\in S_J\), a local valuation of \(x_i\) is
\(-H_J/\ell\), since one summand strictly dominates.
Thus each trace field has degree \(\ell\) and is totally
ramified at every assigned \(p\).

The normal closure of a trace field is totally real with
Galois group \(\operatorname{AGL}_1(\mathbb F_\ell)\).
In fact \(W=J(\zeta,z)\) is Galois over \(\mathbb Q\).
The intersection \(J\cap\mathbb Q(\zeta)\) is trivial
since \(\ell\nmid\operatorname{disc}(J)\), and the local
valuations give \([W:J(\zeta)]=\ell\). If an automorphism
of \(W\) sends \(\zeta\mapsto\zeta^a\) and
\(z\mapsto\zeta^b z^s\) for \(s=\pm1\), its action on
the traces is \(i\mapsto s(ai+b)\). All affine actions
occur. The order-two kernel is generated by complex
conjugation. The trace splitting field \(M_J\) contains
\(\mathbb Q(\zeta+\zeta^{-1})\) and has degree \(2\ell\)
over this common real cyclotomic field.

At an odd \(p\ne\ell\) outside \(S_J\) the roots of the
unit \(w\) are unramified over the corresponding local
field of \(J(\zeta)\). Thus inertia is trivial away from
\(\operatorname{disc}(J)\), and has order two at odd
primes of \(\operatorname{disc}(J)\). In the latter case
it fixes \(\zeta\), exchanges the two embeddings of \(J\),
and acts as a reflection on the traces. Such an action
has one fixed point and \((\ell-1)/2\) transpositions.
Consequently its normalized logarithmic discriminant
contribution to the trace field is
\((1-1/\ell)(\log p)/2\). At each \(p\in S_J\), inertia
acts by translations of full order \(\ell\). These
statements use the tame different and the orbit count
on the roots.
[Chebyshev traces and local ramification](ABC_TOTALLY_REAL_RAMIFICATION_LIFTS.md#the-two-endpoint-marks-for-fixed-real-covers).

The affine group agrees with Girstmair's De Moivre
classification. His irreducibility criterion uses the
same quadratic norm quotient. Here the prescribed
valuations also determine inertia at the assigned primes.
[Girstmair, Theorems 1 and 2](https://arxiv.org/abs/2007.14183).

## One auxiliary prime and a polynomial degree bound

Write \(R=\operatorname{rad}(abc)\), \(A=\log R\), and
take \(S\) to be the odd primes in \(R\). If \(S\) is
empty use \(\mathbb Q\); the bounds below then hold with
\(A\le\log2\). Otherwise cover \(S\) as above, assign each
prime to one splitting field, and discard unused fields.
If \(q\ne1\), assign its field only \(p_0\). The fields
chosen before adding \(q\) already cover the other primes.
Every auxiliary quadratic discriminant has size at most
\(4R\), since \(q\le p_0\) if an extra prime is needed.
The imaginary quadratic class numbers are at most
\(O(\sqrt{|D|}(1+\log|D|))\) for discriminant \(D\).
For example, the Minkowski bound gives an integral ideal
of norm \(O(\sqrt{|D|})\) in each class, and the count of
ideals of a given norm is at most the divisor function.
The product of \(R\), the quadratic discriminants and
these class numbers therefore has log size \(O(A+2)\).
The prime number theorem supplies a common prime
\(\ell\asymp A+2\) avoiding that product, by taking a
dyadic interval with a sufficiently large absolute
constant multiplying \(A+2\).

For each assigned field take one trace field \(L_J\) and
its trace splitting field \(M_J\). Let \(L\) and \(M\) be
their respective composita. They are totally real,
\([L:\mathbb Q]\le\ell^8\), and \(M/\mathbb Q\) is solvable
Galois with degree at most
\((\ell-1)(2\ell)^8/2\), by their common real cyclotomic
subfield. Write \(d(F)=\log\operatorname{rd}(F)\), and
let \(r_F(t)\) denote the normalized logarithmic cusp
radical of \(t=a/c\) over \(F\). Then
\[
 \begin{aligned}
 r_L(t),\ r_M(t)&=O(1),\\
 d(L),\ d(M)&=A+\tfrac12\log q+O(\log(A+2)).
 \end{aligned}                                                   \tag{NT1}
\]
At every odd prime of \(R\) the local ramification
indices in either compositum are all divisible by
\(\ell\), since that prime is totally ramified in an
assigned \(L_J\). The corresponding radical weight is
at most \(A/\ell\). The \(2\)-part is at most \(\log2\).
The same indices are tame: at every prime away from
\(2\ell\) all the individual closures are tame, and
their compositum is tame. Thus the log discriminant
weight over the support is between
\((1-1/\ell)\log(R/\gcd(R,2))\) and
\(\log(R/\gcd(R,2))\).

Outside this support and \(2\ell\), the only possible
ramification is at the auxiliary \(q\ne1\). It is
quadratic in \(M\), and its normalized contribution to
\(d(L)\) is between \((1-1/\ell)\log q/2\) and
\(\log q/2\), by inclusion of its own trace field and
by monotonicity of the local normalized discriminant
in a tower. The difference between these bounds is
\(O(1)\). Finally the local different bound of (PD2)
gives \(O(\log\ell)\) at \(2\ell\), since the order of
the Galois compositum is polynomial in \(\ell\).
[Normalized local-different bound](ABC_PROJECTIVE_DIVISION_FIELD_DISCRIMINANTS.md#projective-division-fields-at-every-level).

There is also a version unramified outside \(2\ell R\).
When \(q\ne1\), omit its assigned trace field and retain
those covering all \(p\in S\setminus\{p_0\}\).
The only odd support prime where ramification can fall
below \(\ell\) is \(p_0\), and its index is one or two.
Write \(M_0\) for this solvable Galois compositum. Then
\[
 d(M_0)+r_{M_0}(t)=A+O(\log(A+2)),\qquad
 r_{M_0}(t)=\frac{\log p_0}{e_{p_0}(M_0/\mathbb Q)}+O(1).
                                                               \tag{NT2}
\]
Here there is no \(p_0\) term when \(q=1\) and \(M_0=M\).
At each odd support prime, the tame discriminant and
remaining radical sum to its original weight. Only \(2\ell\)
can contribute an additional \(O(\log(A+2))\) term.
The alternative (NT2) confines the potentially large
normalized odd-prime conductor to one rational prime.

For example, \((a,b,c)=(1,7,8)\) has \(p_0=7\),
\(q=3\), and \(J=\mathbb Q(\sqrt{-3})\). The generator
\(\beta=2+\sqrt{-3}\) of a prime over \(7\) gives
\(w+w^{-1}=2/7\). Already the order-five trace field
has degree five, all real embeddings, full inertia at
\(7\), and reflection inertia at \(3\) in its closure.

## A single trace gives an optimal linear degree bound

Let \(J_1,\ldots,J_s\) be the distinct quadratic fields
retained above, with a nonempty disjoint assigned support
for each. Set \(B=J_1\cdots J_s\), and let \(F=B^+\) be its
maximal real subfield. In fact \([F:\mathbb Q]\le32\).
For any of the three cuts, the two associated negative
squareclasses multiply to the same squareclass
\(\prod_{p\in S_-}p\). The fields from the cuts therefore
span a squareclass space of dimension at most four.
The field \(\mathbb Q(i)\) and the possible last field add
at most two dimensions in total.

Use the same common \(\ell\), and write \(w_i=z_i^\ell\)
for the norm-one element of \(J_i\) and its chosen
root. Define
\[
 w=\prod_{i=1}^s w_i,\qquad z=\prod_{i=1}^s z_i,\qquad
 K=F(z+z^{-1}).
\]
The product has norm one from \(B\) to \(F\).
Every embedding of \(B\) commutes with complex
conjugation, so all conjugates of \(w\) lie on the unit
circle and all conjugates of the trace are real.
The trace satisfies \(P_\ell(X)=w+w^{-1}\) over \(F\).
In fact \(w+w^{-1}\) generates \(F\) over \(\mathbb Q\).
If an automorphism of \(B\) preserves this coefficient,
it sends \(w\) to \(w\) or \(w^{-1}\). Comparing at one
assigned prime for each \(J_i\) forces it to fix every
\(w_i\) or to invert every \(w_i\), respectively.
These actions are the identity and complex conjugation
on \(B\). Hence the trace generates \(K\) over \(\mathbb Q\).

At every prime over an odd \(p\mid R\), precisely one
factor \(w_i\) has nonzero valuation in \(B\). Its
valuation is \(\pm H_{J_i}e(B/J_i)\), prime to
\(\ell\). With the valuation on \(B\) normalized to
have value group \(\mathbb Z\), the trace thus has
valuation \(-|v(w)|/\ell\). Normalizing on \(F\)
changes this by a factor of one or two. Consequently
\(K/F\) has degree \(\ell\) and is totally ramified
at every prime of \(F\) above the odd support of \(R\).
For each such rational \(p\) the ramification index
in \(K/\mathbb Q\) is \(\ell e_p(F/\mathbb Q)\).
In particular \([K:\mathbb Q]\le32\ell=O(A+2)\).

To see that \(K\subset M\) for the same Galois compositum
as in (NT1), put \(U=B(\zeta)\) and \(W=U(z_1,\ldots,z_s)\).
The field \(W\) is the compositum of the individual Galois
fields \(J_i(\zeta,z_i)\). An element of
\(\operatorname{Gal}(W/M)\) restricts on each such field
to the identity or to complex conjugation. Since the
restriction to the common \(\zeta\) determines which,
\(\operatorname{Gal}(W/M)\) consists precisely of the
identity and complex conjugation. In particular it
fixes \(F\) and the trace \(z+z^{-1}\).

If \(q\ne1\) and \(s\ge2\), its inertia in \(B\) is
nontrivial only on the one \(J_i=\mathbb Q(\sqrt{-q})\).
It therefore survives in \(F\), since complex conjugation
is nontrivial on every \(J_i\). In this case the normalized
log discriminant contribution at \(q\) is exactly
\(\tfrac12\log q\) in both \(K\) and \(M\). If \(s=1\),
\(K\) is the single trace field and its contribution is
\((1-1/\ell)\tfrac12\log q\).
Using \(K\subset M\) to bound the contribution at \(2\ell\)
and tame inertia to compute the odd support contribution
gives
\[
 [K:\mathbb Q]\le32\ell,\qquad r_K(t)=O(1),\qquad
 d(K)=A+\tfrac12\log q+O(\log(A+2)).                    \tag{NT3}
\]
The degree order for bounded normalized radical is optimal
over all number fields. For any number field \(K'\) of
degree \(n\), the contribution to \(r_{K'}(t)\) at each
rational \(p\mid R\) is at least \((\log p)/n\), since all
places above it lie on the same cusp and the sum of their
residue degrees is at least one. Thus \(r_{K'}(t)\ge A/n\).

When \(q\ne1\), form \(K_0\) by the same product after
omitting the factor and field assigned to \(p_0\). If none
remains use \(\mathbb Q\). Then \(K_0\subset M_0\), so it is
unramified outside \(2\ell R\) and its indices over \(p_0\)
are one or two. At the other odd support primes the trace
still gives full inertia of order divisible by \(\ell\).
The tame identity between discriminant and residual radical
therefore gives
\[
 d(K_0)+r_{K_0}(t)=A+O(\log(A+2)).
                                                               \tag{NT4}
\]
Here the contribution to \(r_{K_0}(t)\) at \(p_0\) is between
\(\tfrac12\log p_0\) and \(\log p_0\), and the total at all
other primes is bounded. When \(q=1\) use \(K_0=K\).

For \((a,b,c)=(1,343,344)\) the exceptional prime is
\(7\) and \(q=3\). Assign \(7\) to \(\mathbb Q(\sqrt{-3})\)
and \(43\) to \(\mathbb Q(\sqrt{-7})\), using generators
\(2+\sqrt{-3}\) and \(6+\sqrt{-7}\) of the respective
prime ideals. With \(\ell=5\), the trace over
\(F=\mathbb Q(\sqrt{21})\) satisfies
\[
 \left(P_5(X)-\frac{58}{301}\right)^2
        -21\left(\frac{96}{301}\right)^2=0,\qquad
 P_5(X)=X^5-5X^3+5X.
\]
It has degree ten over \(\mathbb Q\), with ramification
indices ten at \(7\), five at \(43\), and two at \(3\).

### The exact cost of taking the Galois closure

Let \(C=\operatorname{Gal}(B/\mathbb Q)\) with complex
conjugation \(c_B\), and let \(\chi_i:C\to\{1,-1\}\)
be the character of \(J_i\). The \(\ell\)-Kummer classes
of the \(w_i\) in \(U^\times/U^{\times\ell}\) are
independent. For each \(i\), a valuation over any of its
assigned primes annihilates all the other classes and
has nonzero value modulo \(\ell\) on \(w_i\).
Also \(B\cap\mathbb Q(\zeta)=\mathbb Q\) by ramification
at \(\ell\). It follows that the full group of \(M\) is
\[
 \operatorname{Gal}(M/\mathbb Q)
 =\mathbb F_\ell^s\rtimes
 \bigl((C\times\mathbb F_\ell^\times)/\langle(c_B,-1)\rangle\bigr).
                                                               \tag{NT5}
\]
In \(W\), a pair \((g,a)\) in the complementary group
lifts by \(z_i\mapsto z_i^{\chi_i(g)}\) and
\(\zeta\mapsto\zeta^a\). These lifts exist by Kummer
independence and multiply as prescribed.
The action on translations in coordinate \(i\) is by
\(a\chi_i(g)\). Its kernel before the quotient is
\(\langle(c_B,-1)\rangle\), since the \(\chi_i\) separate
the elements of \(C\). In particular the action in
(NT5) is faithful.

The translation subgroup fixes \(K\) along the single
hyperplane \(\sum_i v_i=0\). Its conjugates intersect
trivially, because the distinct characters \(\chi_i\)
are linearly independent over \(\mathbb F_\ell\).
If a normal subgroup intersects the translations trivially,
its commutator with them is trivial. Faithfulness of
the action then puts the subgroup inside the translations,
so it is trivial. Applying this to the largest normal
subgroup fixing \(K\) shows that \(M\) is its normal closure.
More generally, any Galois subfield of \(M/\mathbb Q\)
that keeps the \(\ell\)-inertia at one assigned prime per
\(J_i\) must be all of \(M\). At such a prime the
\(\ell\)-inertia is the corresponding coordinate axis.
Any invariant subspace of the translations is a sum of
coordinate axes, by the projections obtained from the
distinct \(\chi_i\). Thus the normal subgroup defining
such a subfield intersects the translations trivially.

### Modular base change and a single Hecke twist

The normal closure of \(K/F\) is the real subfield
\(N\) of \(B(\zeta,z)\). The proof for one quadratic
norm-one trace applies over \(F\): the local valuations
give the degree \(\ell\) Kummer extension, and \(B/F\)
is quadratic and disjoint from \(F(\zeta)/F\).
The group of \(N/F\) is \(\operatorname{AGL}_1(\mathbb F_\ell)\).
Its translation subgroup \(T\) fixes the field
\(U^+=(B(\zeta))^+\), abelian over \(\mathbb Q\).
Write \(\rho\) for the permutation action on the
\(\ell\) conjugates of the trace over \(F\). For any
nontrivial character \(\psi\) of \(T\),
\[
 \rho=\mathbf 1\oplus
 \operatorname{Ind}_{T}^{\operatorname{AGL}_1(\mathbb F_\ell)}\psi.
\]
Indeed \(T\) acts regularly on the roots, and the
multiplicative group acts transitively on the nontrivial
characters in its regular representation. By induction
of Euler factors, an elliptic curve \(E/\mathbb Q\)
therefore satisfies
\[
 L(E/K,s)=L(E/F,s)\,L(E_{U^+}\otimes\psi,s).
                                                               \tag{NT6}
\]
Here \(\psi\) is the Hecke character of the cyclic
extension \(N/U^+\). The same identity applies to the
adjoint Euler product. This uses all local factors,
including at ramified primes.
[Artin formalism, Dokchitser–Evans–Wiersema, §§2.5–2.6](https://arxiv.org/abs/1905.04282).

The character \(\psi\) has conductor exponent one at
each prime of \(U^+\) above the odd support of \(R\).
At each such prime the degree-\(\ell\) relative inertia
is full and tame. Away from \(2\ell R\) the character
is unramified. At \(q\), the inertia in \(M\) has order
two, while \(N/U^+\) has odd prime degree.
Thus the nonabelian Euler product in (NT6) needs only
an abelian base change and a single cyclic twist.
Artin conductor induction also gives the exact identity
of relative discriminant and conductor ideals
\[
 \mathfrak d_{K/F}
 =\mathfrak d_{U^+/F}\,
   \mathrm N_{U^+/F}\bigl(\mathfrak f(\psi)\bigr).
                                                               \tag{NT7}
\]
Indeed the conductor of \(\rho\) is the relative field
discriminant, and the trivial summand has unit conductor.
The twist accounts for the odd support, and the auxiliary
\(q\) occurs only in the discriminant of the abelian base.

The rational Frey curve is modular over \(K\) and \(K_0\).
For nonempty odd support it has multiplicative reduction
at an odd prime and hence has no complex multiplication.
Its Tate representation stays absolutely irreducible
over finite extensions by the open image theorem.
Start with rational modularity and base change up a
cyclic tower to the common solvable Galois closure.
Descend along a cyclic tower to \(K\) or \(K_0\).
At each descent the upstairs representation is invariant
because it matches the restriction of the Tate
representation. Cyclic descent gives a lower form that
matches after restriction. Absolute irreducibility and
Schur's lemma show that the lower Galois representations
differ by a character of the cyclic quotient. Twisting
corrects this difference before going down the next step.
An empty odd support uses the rational field.
[Rational modularity, Theorem A](https://math.stanford.edu/~conrad/papers/tswfinal.pdf).
[Serre, Théorème 2](https://www.college-de-france.fr/media/jean-pierre-serre/UPL5874918517843398173_Serre_proprie_te_s_galoisiennes_des_courbes_elliptiques.pdf).
[Langlands, §§2 and 11](https://publications.ias.edu/sites/default/files/base-change-for-gl2_rpl_2.pdf).
[Solvable descent for irreducible restrictions, Dieulefait, p. 1024](https://annals.math.princeton.edu/wp-content/uploads/annals-v176-n2-p07-p.pdf).

## The remaining cost in the pointwise height problem

There is no auxiliary cost for every support without
an exceptional \(p_0\). Even on the remaining inputs
Burgess' least-nonresidue bound gives, for any
\(\epsilon>0\),
\[
 d(M)\le A+
    \left(\frac{1}{8\sqrt{\mathrm e}}+\epsilon\right)\log p_0
    +O_\epsilon(\log(A+2)).
\]
The same estimate holds for \(L\) and \(K\). Under GRH for Dirichlet
\(L\)-functions, the bound
\(q\ll(\log p_0)^2\) makes the auxiliary cost \(O(\log(A+2))\).
The construction and (NT1) are unconditional.
[Burgess, Theorem 2](https://www.cs.umd.edu/users/gasarch/TOPICS/res/burgess.pdf).
[Conditional quadratic nonresidue bound](https://doi.org/10.1090/mcom/4154).

Unconditionally, define \(D_{\rm tr}(t)=d(K)-\frac12\log q\)
for the linear-degree trace lift. By (NT3),
\(D_{\rm tr}(t)=\log\operatorname{rad}(abc)+O(\log(h(t)+2))\).
Thus abc is equivalent to the sharp pointwise bound
\(h(t)\le(1+\epsilon)D_{\rm tr}(t)+O_\epsilon(1)\) on these
lifts. The same equivalence over \(M\) follows from (NT1).
Any choices in the construction have the same error
bounds. The primitive Frey curve is multiplicative
precisely over the odd support of \(R\), so its normalized
odd-prime log conductor is bounded over \(K\) and \(M\).
The one possibly subtracted prime \(q\) lies outside
\(2R\), so the original elliptic curve has good reduction
there. This growing-level height estimate requires
degree uniformity and a sharp discriminant coefficient.
The fixed-level construction below permits constants
depending on each level. Pasten's Theorem 18.10 has a
parametrization degree term and a field-degree-dependent
coefficient on the absolute discriminant and conductor,
with hypotheses of Shimura level one. Theorem 18.13
also needs the modular-degree conjecture and allows the
coefficient to depend on the degree. Automorphic base
change and (NT6) provide the form and its Euler factors.
A uniform comparison of the integral periods or
congruence invariants with the Faltings height is still
needed for this route.
[Pasten, Theorems 18.10 and 18.13](https://arxiv.org/abs/1705.09251).
The complete log root discriminant yields the same
unconditional reformulation when \(q=1\). Under the indicated
GRH bound it does so on all the constructed lifts. For the
construction supported on \(2\ell R\), (NT4) gives the
unconditional equivalent bound
\[
 h(t)\le(1+\epsilon)\bigl(d(K_0)+r_{K_0}(t)\bigr)+O_\epsilon(1).
\]
One can also use \(M_0\) by (NT2).
Here only one rational prime can carry an unbounded
normalized cusp weight. Both equivalences absorb the
\(O(\log(h(t)+2))\) error by applying the proposed bound
with a smaller positive \(\epsilon\).

## Fixed odd levels with a single ideal-class correction prime

Fix any odd integer \(N\ge3\), independently of the triple.
Put \(S=\{p:p\mid R,\ p\nmid 2N\}\). There is a totally
real lift \(K_N\) of degree at most \(32N\), with solvable
normal closure of degree bounded in terms of \(N\) alone.
Away from \(2N\), it ramifies only on \(S\) and at most two
primes of good reduction for the original Frey curve.
One is the previous cover prime \(q\). The other, \(P\),
corrects every ideal class simultaneously and satisfies
\(\log P=O(A+\log(2N))\) with an absolute constant.
One or both primes may be absent.

Use the same imaginary quadratic cover on \(S\).
If \(S\) is empty, take \(K_N=\mathbb Q\). Otherwise let
\(J_1,\ldots,J_s\) be its distinct used fields, \(B\) their
compositum and \(F=B^+\). Choose one prime ideal over each
assigned rational prime in its \(J_i\), and let
\(\mathfrak a_i\) be their product. Write
\[
 A_i=\operatorname{Cl}(J_i)/N\operatorname{Cl}(J_i).
\]
It suffices to correct the class of \(\mathfrak a_i\) in
\(A_i\). A further factor with ideal exponent divisible by
\(N\) can then make the ideal principal. Its primes do not
introduce tame Kummer ramification.

### Simultaneous classes at one split prime

Let \(H_i/J_i\) be the unramified class field with group
\(A_i\). Since \(B/J_i\) has degree a power of two and
\(A_i\) has odd order, \(BH_i/B\) also has group \(A_i\).
The extensions \(BH_i/B\) are unramified. Their compositum
\(H/B\) has group equal to the whole product \(\prod_i A_i\).
Indeed its group is a subgroup of the product, projects
onto each factor, and is stable under
\(C=\operatorname{Gal}(B/\mathbb Q)\). The action on the
\(i\)-th factor is by the quadratic character \(\chi_i\) of
\(J_i\), since complex conjugation inverts ideal classes.
The characters are distinct and \(|C|\) is invertible
modulo \(N\). On a module killed by \(N\), the averaging
operator
\[
 \frac1{|C|}\sum_{g\in C}\chi_i(g)g
\]
projects onto the \(i\)-th factor. Applying it to the
subgroup proves surjectivity onto the product.

Choose the target element \((-[\mathfrak a_i])_i\) in
this product. If it is zero, take \(P=1\) and omit all
correction factors. Otherwise choose a degree-one prime
\(\mathfrak P\) of \(B\) with this Artin symbol in \(H/B\)
and lying over \(P\nmid 2NRq\). Since \(B/\mathbb Q\)
is Galois and unramified at \(P\), this rational prime
splits completely in \(B\). The prime
\(\mathfrak P_i=\mathfrak P\cap J_i\) has class
\(-[\mathfrak a_i]\) in \(A_i\), by compatibility of Artin
symbols with the norm and residue degree one.

For a quantitative choice, apply Thorner and Zaman's
Theorem 1.1 over \(B\) to the abelian extension
\(H B(\zeta_m)/B\), where
\[
 m=256\!\prod_{\substack{p\mid NRq\\p\ {\rm odd}}}p^2 .
\]
Choose any automorphism restricting to the target in
\(H/B\). At every excluded odd prime, the cyclotomic
ramification has an odd factor \(p\), which cannot be
absorbed by \(B\). At \(2\) its ramification index is
at least \(128>[B:\mathbb Q]\). Thus all excluded primes
ramify in the enlarged extension over \(B\).
The unramified class fields contribute unit conductor.
The cyclotomic characters pulled back by norm have
conductors dividing \(m\mathcal O_B\). The norm of this
modulus is at most \(m^{64}\). Also \(B\) has degree at
most \(64\) and \(\log|\operatorname{disc}(B)|=O(A+1)\),
as its odd inertia is at most quadratic and supported on
\(Rq\). In Theorem 1.1 the parameter \(\mathcal Q\) is the
largest character conductor norm, and the guaranteed
prime ideal has rational prime norm. These facts give
\(\log P=O(A+\log(2N))\), with no degree factor coming from
the unramified class field.
[Thorner–Zaman, Theorem 1.1](https://msp.org/ant/2017/11-5/ant-v11-n5-p04-s.pdf).

For example, \(\mathbb Q(\sqrt{-23})\) and
\(\mathbb Q(\sqrt{-31})\) have class number three.
The split primes over \(13\) and \(5\), respectively,
are nonprincipal. One rational prime \(P=41\) corrects
both classes. The generators
\(21+2\sqrt{-23}\) and \(9+2\sqrt{-31}\) have norms
\(13\cdot41\) and \(5\cdot41\), respectively, and give
the chosen orientations of the prime ideals. At order
\(N=3\) the signed correction weights \(1,2\) give full
relative inertia at one prime over \(41\) and trivial
relative inertia at the other. Both assigned primes
still have full relative inertia.

### Ramification and the fixed-degree height obligation

Take \(\epsilon_i=1\) if \([\mathfrak a_i]\ne0\) in
\(A_i\), and zero otherwise. Choose \(\mathfrak b_i\)
such that
\[
 (\beta_i)=\mathfrak a_i\mathfrak P_i^{\epsilon_i}
                           \mathfrak b_i^N
\]
is principal, where a correction factor of exponent
zero is omitted. Define \(w_i=\beta_i/\bar\beta_i\),
\(w=\prod_i w_i^{2^{i-1}}\), and \(z^N=w\).
Then \(K_N=F(z+z^{-1})\) is totally real. At every
place of \(B\) above any assigned \(p\in S\), the
valuation of \(w\) is prime to \(N\). Precisely one
component contributes a nonzero residue modulo \(N\),
equal up to sign to \(2^{i-1}e(B/J_i)\). These
integers are powers of two. At such a place the
trace therefore has a valuation with denominator \(N\),
which proves \([K_N:F]=N\) and total relative
ramification there. At any place of \(B\) away from
\(S\) and \(P\), the valuation of \(w\) is divisible by
\(N\), and the same holds for every \(w_i\).
After adjoining \(\zeta_N\), an \(N\)-th root of a unit
is unramified at residue characteristics prime to \(N\).
Consequently there is no further tame ramification
outside that of \(B\).

For the normal closure bound, adjoin \(\zeta_N\) and
an \(N\)-th root of each \(w_i\) to \(B\). This field is
Galois over \(\mathbb Q\), since conjugating \(J_i\)
inverts \(w_i\). Its group is an extension of an
abelian group by an abelian group, and has order at
most \(64\varphi(N)N^8\). Complex conjugation is central:
it inverts each chosen root and \(\zeta_N\), and
commutes with every automorphism on \(B\). Thus the
maximal real subfield is Galois, contains \(K_N\),
and supplies the asserted bound. The normalized log
different above \(2N\) is \(O(\log(2N))\) by (PD2).
The non-CM modularity argument above applies to this
solvable closure as well.
[Normalized local-different bound](ABC_PROJECTIVE_DIVISION_FIELD_DISCRIMINANTS.md#projective-division-fields-at-every-level).

If \(N\) is prime and \(N>255\), the discriminant cost
at \(P\ne1\) is exactly \((1-1/N)\log P\).
Indeed at each prime over \(P\) the valuation of \(w\)
is congruent to a nonempty signed sum of distinct
powers of two among \(1,\ldots,128\). The sum is
nonzero and has absolute value at most \(255\).
The trace is then totally ramified of index \(N\)
above every prime of \(F\) over \(P\). For general
odd \(N\) its local cost at \(P\) can be removed
by restricting the discriminant to the original
support. Define
\[
 d_S(K_N)=\frac{1}{[K_N:\mathbb Q]}
       \sum_{p\in S} v_p(|\operatorname{disc}(K_N)|)\log p .
\]
Thus the omitted rational primes divide \(2N\) or
have good reduction for the original Frey curve.
For all odd \(N\ge3\),
\[
 (1-1/N)(A-\log(2N))\le d_S(K_N)\le A,\qquad
 d_S(K_N)+r_{K_N}(t)=A+O(\log(2N)).
                                                               \tag{NT8}
\]
At every \(p\in S\), this is the tame identity:
discriminant and residual radical contribute
\((1-1/e_p)\log p\) and \((1/e_p)\log p\), with
\(e_p=N e_p(F/\mathbb Q)\). At the omitted primes
of the original support the total radical weight
is at most \(\log(2N)\). In particular
\(r_{K_N}(t)\le A/N+\log(2N)\).
The complete log root discriminant is \(O(A+\log(2N))\)
by the bound on \(P\), the bound on \(q\), and the
wild different estimate.

Full abc is equivalent to a pointwise bound
\[
 h(t)\le b_N d_S(K_N)+C_N
                                                               \tag{NT9}
\]
on the constructed fields at an unbounded sequence of
fixed odd levels \(N\), with \(b_N\to1\).
The additive constant may depend arbitrarily on the fixed
level. In one direction use \(d_S(K_N)\le A\) and take
one sufficiently large \(N\) for the chosen exponent.
Conversely abc and the lower bound in (NT8) give (NT9)
with \(b_N=(1+1/N)/(1-1/N)\) for every odd \(N\ge3\).
This applies to any choices meeting the construction.
The pointwise comparison must use the discriminant
supported at primes of multiplicative reduction, or
control the auxiliary good-reduction contributions
separately. A bound for full root discriminant with
the same leading coefficient does not remove the
repair prime. The available estimate for its logarithm
is proportional to \(A+\log(2N)\).

### Cancelling the auxiliary discriminants by a second lift

Construct two lifts at the same fixed odd level \(N\)
with disjoint auxiliary ramification away from \(2N\).
First choose \(q_1,P_1\) as above. If the cover has an
exceptional \(p_0\), replace its field in the second
cover by \(\mathbb Q(\sqrt{-q_2})\), choosing \(q_2\)
a quadratic nonresidue modulo \(p_0\) and outside
\(2NRq_1P_1\). One can impose a nonzero nonresidue at
\(p_0\) and the residue \(1\) at all other primes in
this excluded product. The cyclotomic case of
Thorner and Zaman's Theorem 1.1 gives such a \(q_2\)
with \(\log q_2=O(A+\log(2N))\). In the absence of an
exception take \(q_2=1\).
For the second correction \(P_2\), use the same
Artin-symbol prescription with \(q_1,P_1\) also
excluded from the cyclotomic modulus. The norm
and discriminant bounds remain \(O(A+\log(2N))\)
on the logarithmic scale. A correction prime may
again be omitted if its target class vector is zero.

Choose the totally real solvable Galois fields
\(M_1,M_2\) containing the two trace fields by the
normal closure construction above, and put
\(M_{12}=M_1M_2\). All three degrees are bounded
solely in terms of \(N\). Set
\[
 D_N^{\cap}(t)=d(M_1)+d(M_2)-d(M_{12}).
\]
At an odd prime away from \(N\), tame inertia in
the compositum has index equal to the least common
multiple of the two inertia indices. This follows
by taking both quotients of the same procyclic tame
inertia group. If the prime ramifies only in one
field, its three contributions cancel exactly.
If \(p\in S\), its two indices are \(N\) or \(2N\):
each is divisible by \(N\) by the trace construction,
has factor at most two from \(B\), and has tame
Kummer inertia of exponent dividing \(N\).
One index therefore divides the other.
The combined contribution at \(p\) is
\((1-1/\min(e_1,e_2))\log p\).
Only \(2N\) can contribute an error to this tame formula,
and its combined normalized different is \(O(\log(2N))\).
Consequently
\[
 (1-1/N)A-O(\log(2N))\le D_N^{\cap}(t)
                         \le A+O(\log(2N)).             \tag{NT10}
\]

The difference uses ordinary root discriminants of
three totally real fields and treats every rational
prime by the same rule. Its tame cancellation at
auxiliary primes does not require a bound on the size
of each prime. Full abc is equivalent to
\(h(t)\le b_N D_N^{\cap}(t)+C_N\) on these pairs at
an unbounded sequence of fixed odd levels with
\(b_N\to1\). This is also the difference between
the log root discriminant of \(M_1\) and the
normalized relative log discriminant of
\(M_{12}/M_2\). A comparison that bounds the two
absolute discriminants separately at a common
leading coefficient does not give the required
relative saving.
The same class-correction construction also gives a
conductor comparison for mixed Hecke twists over one
abelian field, recovering the full radical at a single
fixed odd level.
[Mixed twist conductors](ABC_PAIRED_TWIST_CONDUCTORS.md).
