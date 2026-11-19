# Alexander-module bounds for protected braid closures

The ordinary Alexander module cannot certify the full displayed strand
number for a closure with protected blocks. Each block supplies a
relation that removes one generator over the integral Laurent ring.
There is a second independent relation per block at all field
specializations except \(t=59\) in characteristic \(83\).

## The Alexander lower bound

Let \(R=\mathbb Z[t,t^{-1}]\) and \(G=\pi_1(S^3\setminus K)\).
The Alexander module \(A_K\) is the abelianization of the kernel of the
exponent-sum map \(G\to\mathbb Z\), with \(t\) acting by conjugation by a
positive meridian. If \(m\) meridians \(x_1,\ldots,x_m\) generate \(G\),
the corresponding exponent-sum kernel in the free group maps onto
this kernel. Schreier rewriting with representatives \(x_1^j\) gives
the free kernel generators \(x_1^j(x_i x_1^{-1})x_1^{-j}\), where
\(2\le i\le m\) and \(j\in\mathbb Z\). After abelianization these give
\(m-1\) module generators. For any field \(F\) and \(\alpha\ne0\), with
\(F_\alpha\) the \(R\)-module defined by \(t\mapsto\alpha\), this yields
\[
 \operatorname{mr}(K)\ge 1+g(A_K)\ge 1+\dim_F(A_K\otimes_R F_\alpha).
\]
The related elementary-ideal formulation appears in Proposition 4 of
Boden and Gaudreau, *Bridge numbers for virtual and welded knots*.
[Boden and Gaudreau, §2](https://math.mcmaster.ca/~boden/preprints/hub29.pdf)

For a braid word \(\beta\) on \(b\) strands, let \(B(\beta)\) be the
matrix with rows given by the abelianized Fox derivatives of the
bottom meridian words in terms of the top meridians. For the signed
Hurwitz convention \(T(a,c)=(c,c^{-1}ac)\), a positive crossing changes
the two rows \((v,w)\) to \((w,t^{-1}v+(1-t^{-1})w)\).
A negative crossing changes them to \(((1-t)v+tw,v)\).
These formulas follow from the Fox product rule since both input words
have exponent sum one. Start with \(I_b\) and apply the updates in word
order. The closure relators identify every bottom meridian with the
corresponding top meridian, so their Fox rows form \(D=B(\beta)-I_b\).

The boundary out of \(R^b\) in the infinite cyclic cover is
\((t-1)\sum_i v_i\). Its kernel is the free \(R\)-module
\(S=\{v:\sum_i v_i=0\}\), and the rows of \(D\) present \(A_K\) as a
quotient of \(S\). Deleting the last coordinate identifies \(S\) with
\(R^{b-1}\). This presentation specializes by tensoring, including at
\(\alpha=1\). As the rows of \(D\) sum to zero coordinatewise within
each row, deleting that coordinate preserves row rank. Thus
\[
 \dim_F(A_K\otimes_R F_\alpha)=b-1-\operatorname{rank}_F D(\alpha).
\]
The right nullity of the full \(b\times b\) Fox matrix is one larger
than the module dimension.

## The protected relations

Use the signed block
\[
 U=(1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3).
\]
The prefix \(P_k\) puts copies of \(U\) on \(k\) disjoint consecutive
sets of four strands. The tail \(\tau\) fixes the two outer free
meridian words of each block. For a knot closure \(K=\widehat{P_k\tau}\),
its two outer closure relations in each block therefore have the Fox
rows of \(B(U)-I_4\) in the four columns of that block, with zeros
elsewhere. This applies to any such tail.
[Protected braid construction, §§1–3](MULTIPLE_PROTECTED_OBSTRUCTIONS.md)

**Theorem (integral generator bound and specialization ranks).**
For every such closure, with \(k\ge1\) and \(b\ge4k\),
\[
 g(A_K)\le b-1-k.
\]
For every field \(F\) and \(\alpha\in F^\times\), the two outer rows
have rank \(2\) in each block unless \(\operatorname{char}F=83\)
and \(\alpha=59\) in the prime subfield. At that specialization they
have rank \(1\). Consequently the specialized Alexander-module
dimension is at most \(b-1-2k\) in the first case and \(b-1-k\) in
the second case.

**Proof.** Let \(D_U=B(U)-I_4\). The crossing recurrence gives the first
and fourth entries of its first row as \(t^{-4} f(t)\) and
\(-t^{-1}(t-1)^2h(t)\), where
\[
 f=t^6-3t^5+4t^4-7t^3+7t^2-4t+1,\qquad
 h=t^3-t^2+2t-1.
\]
These entries generate the unit ideal in \(R\). In fact, \(f(1)=-1\)
shows that \(t-1\) is a unit modulo \(f\). The identities
\[
 \begin{aligned}
 f&=(t^3-2t^2-2)h+3t^2-1,\\
 3h&=(t-1)(3t^2-1)+7t-4,\\
 3(7t+4)(7t-4)-49(3t^2-1)&=1
 \end{aligned}
\]
show that \(f,h\) generate the unit ideal already in \(\mathbb Z[t]\).
It follows that \(f,(t-1)^2h\) generate the unit ideal there as well.

In the sum-zero submodule on a block, use the basis
\(e_1-e_3,e_2-e_3,e_4-e_3\) in local indices. The protected first row
has coordinates \((a,c,d)\) with \(ua+vd=1\) for some \(u,v\in R\).
It extends to a basis by adjoining \(e_2-e_3\) and
\(-v(e_1-e_3)+u(e_4-e_3)\), whose determinant with the first row is one.
The block submodules together form a direct summand of \(S\): summing
coordinates within each block and keeping each extra coordinate gives
a split map to the sum-zero module on the \(b-3k\) resulting parts.
A section picks one coordinate in each part for its sum. We can
therefore quotient \(S\) by the \(k\) first protected rows and obtain
a free module of rank \(b-1-k\). Its further quotient \(A_K\) is
generated by at most that many elements.

For the local rank statement, take the six \(2\times2\) minors of the
two outer rows of \(D_U\), clearing Laurent monomials in each. The
numerators for column pairs \((1,3)\) and \((2,3)\) have resultant
\(-83\). Thus a simultaneous zero in a field can only occur in
characteristic \(83\). Over \(\mathbb F_{83}[t]\) the monic gcd of all
six numerators is \(t+24\), so the unique common zero in any extension
field is \(59\). The unit ideal calculation guarantees rank at least
one there. These are exact polynomial calculations obtained from the
signed word in the [checker](verify_protected_alexander.py). The same
checker independently differentiates the literal free words of \(U\)
to verify the recurrence. Disjoint column supports make the block
ranks add, and the formula for the module dimension gives the two
bounds. \(\square\)

In particular \(g(A_K)+1\le b-k\), the same numerical value as the
[two-sided augmentation upper bound](TWO_SIDED_PROTECTED_OBSTRUCTIONS.md).
The Alexander-generator lower bound and that augmentation upper bound
alone cannot certify a positive gap for this architecture. For
\(\alpha\ne1\), the dimension bound \(b-1-k\) is also the bound from
[the rank-one cohomology argument](RANK_ONE_COHOMOLOGY_COVECTOR_BARRIER.md)
with one-dimensional coefficients. The integral generator bound and
the stronger rank-two calculation use the explicit protected Fox rows.

## Search on the two-block closures

Use the explicit mobile-band tail \(R_{2,b}^{\,1+(b-4)\ell}\) of the
protected construction. The search covers all \(8\le b\le12\),
\(0\le\ell\le15\), all primes \(p\le101\), and every
\(\alpha\in\mathbb F_p^\times\). It checks each closure permutation is a
single cycle. It then computes \(B(R_{2,b})^{1+(b-4)\ell}B(P_2)-I_b\)
over the prime field and evaluates the quotient dimension by modular
elimination.

**Proposition (bounded prime-field search).** Across the 80 braid words and
90,800 evaluations,
the maximum Alexander-module dimension for each of the 80 words
is \(2\). Every word attains it over \(\mathbb F_{83}\) at \(t=59\).
For example, the eight-strand word with \(\ell=0\) also attains it
over \(\mathbb F_5\) at \(t=2\). The resulting meridional lower bound
is \(3\). The actual meridional ranks and the prime hyperbolic
gap-two realization remain open.

The executable checker compares the two Fox computations on the complete
eight-strand word with \(\ell=0\) and tests the trefoil at \(t=1,-1\)
in characteristic \(3\) as a dimension control. The finite search uses
prime fields only. The proof of the obstruction applies to every field
and every tail with the stated free-word protection.

From the repository root, using Python 3.12.3 and SymPy 1.14.0, run

```sh
.venv/bin/python research/meridional-rank/representations/verify_protected_alexander.py \
  --write-receipt research/meridional-rank/representations/protected_alexander_verification.json
```

[Exact verification receipt](protected_alexander_verification.json)
