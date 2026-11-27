# Composita of dihedral fields with twisted quadratic resolvents

Fix an odd prime \(\ell\) and a positive squarefree integer \(t>1\).
Consider two dihedral Galois closures \(N_1,N_2\) over \(\mathbf Q\)
with groups \(D_\ell=C_\ell\rtimes C_2\) of order \(2\ell\)
and quadratic resolvents
\[
            K_1=\mathbf Q(\sqrt d),\qquad
            K_2=\mathbf Q(\sqrt{td})
\]
for signed squarefree \(d\) with both fields quadratic.
Let \(E_i\) be a degree-\(\ell\) subfield of \(N_i\).
The quadratic resolvents are distinct, so \(N_1,N_2\) are
linearly disjoint: the only nontrivial proper quotient of
\(D_\ell\) is \(C_2\), and either nontrivial common quotient
would identify the quadratic resolvents. Thus \(E=E_1E_2\)
has degree \(\ell^2\) and normal closure \(N=N_1N_2\) of degree
\(4\ell^2\). Its permutation action is the product of the
two degree-\(\ell\) actions.

**Theorem TD1.** Suppose that over either prescribed radicand
sign the mixed class-group moment satisfies
\[
 \sum_{\substack{|d|\le X\\d\ {\rm squarefree}}}
   |\operatorname{Cl}(K_1)[\ell]|\,|\operatorname{Cl}(K_2)[\ell]|
                   \ll_{t,\ell,\epsilon}X^{\Theta+\epsilon},
 \qquad \Theta\ge1 .                                             \tag{1}
\]
Then the number of degree-\(\ell^2\) composita \(E\) of
absolute discriminant at most \(Y\) in this family is
\(O_{t,\ell,\epsilon}(Y^{2\Theta/(\ell^2-1)+\epsilon})\).
The corresponding normal closures \(N\) of absolute
discriminant at most \(Y\) are
\(O_{t,\ell,\epsilon}(Y^{\Theta/(2\ell^2)+\epsilon})\).
These counts include arbitrary finite relative ramification
in \(N_i/K_i\).

## Shared and private relative ramification

Let the fixed bad rational primes be those dividing \(2\ell t\).
At a prime outside this set dividing \(d\), both quadratic
fields ramify. In each \(N_i/\mathbf Q\) the cyclic tame inertia
surjects onto the quadratic quotient and hence is generated
by a reflection. A cyclic inertia group in \(D_\ell\)
that contains a reflection has order two. The inertia in
\(N\) therefore acts by a pair of reflections. On the
\(\ell^2\) geometric embeddings of \(E\) this fixes exactly
one point, giving discriminant exponent \((\ell^2-1)/2\).
On the regular action of \(N\) the exponent is \(2\ell^2\).

At a prime outside the fixed set and not dividing \(d\),
both quadratic fields are unramified. Inertia for \(N\)
is a cyclic subgroup of \(C_\ell\times C_\ell\), so if
nontrivial it has order \(\ell\). Its generator has a
nontrivial rotation on at least one factor, hence has
no fixed points on the product action. The exponent is
\(\ell^2-\ell\), whether one or both of the two relative
extensions ramify there. In the regular action it is
\(4\ell(\ell-1)\). These exponents follow from the tame
degree-minus-orbits formula for permutation Artin conductors.
[Jones–Roberts, §§3B–3C](https://msp.org/ant/2014/8-3/ant-v8-n3-p04-p.pdf)

There is also a necessary compatibility at a prime ramifying
relatively in both factors. Frobenius conjugates a tame
inertia generator by the \(p\)-th power. In \(D_\ell\),
the conjugation on the rotation group is multiplication
by \(+1\) or \(-1\) according to the quadratic Frobenius.
Thus \(p\bmod\ell\) must equal that sign for each ramifying
factor. In particular both factors can ramify relatively at
such a prime only if \(p\) splits in the fixed twist field
\(\mathbf Q(\sqrt t)\). The upper count below can allow
every squarefree relative support, without using this
additional restriction.

Write \(b\) for the squarefree product of the rational
primes outside the fixed set where at least one relative
extension ramifies. The resulting bounds are
\[
 \begin{split}
 |\operatorname{Disc}(E)| &\gg_{t,\ell}
       |d|^{(\ell^2-1)/2} b^{\ell^2-\ell},\\
 |\operatorname{Disc}(N)| &\gg_{t,\ell}
       |d|^{2\ell^2} b^{4\ell(\ell-1)}.                            \tag{2}
 \end{split}
\]
Since \(d\) is squarefree, omitting the fixed primes from
its product changes only the constant.

## Uniform counting of relative conductors

For fixed \(d,b\), bound the two cyclic extensions by
characters of ray class groups over \(K_i\), allowing
all primes above \(b\) and the fixed primes. Away from
the fixed primes the conductor exponent at a ramified
prime is one since the extension is tame. At the
fixed primes the conductor exponent is bounded uniformly,
because the characters have order \(\ell\) and a
sufficiently deep local unit group consists of \(\ell\)-th
powers by Hensel's lemma in local fields of degree at
most two. The real places contribute no odd-order kernel.
The ray class exact sequence then bounds the total
number of ordered pairs of extensions by
\[
          \ll_{\ell,t} \ell^{4\omega(b)}
                 |\operatorname{Cl}(K_1)[\ell]|\,
                 |\operatorname{Cl}(K_2)[\ell]|.
 \tag{3}
\]
Indeed there are at most two primes in each \(K_i\) at
each rational \(p\mid b\), with cyclic residue unit group
at each place. Allowing all the ray characters bounds
the subset with the required inversion action. The
passage from normal extensions to their degree-\(\ell\)
subfields costs only constants depending on \(\ell\).

Use (1) and (2) at fixed \(b\). For degree \(\ell^2\)
the exponent in the \(b\)-sum is
\((\ell^2-\ell)2\Theta/(\ell^2-1)=2\ell\Theta/(\ell+1)>1\).
For regular degree \(4\ell^2\) it is
\(4\ell(\ell-1)\Theta/(2\ell^2)=2(\ell-1)\Theta/\ell>1\).
The factor in (3) costs at most an arbitrarily small
power of \(b\), so both sums converge. This proves
the two bounds with all the allowed relative conductors.

## Upper and lower counts

**Corollary TD2.** For every fixed squarefree \(t>1\), the
restricted families in TD1 at \(\ell=3,5\) have the following
upper exponents, each up to \(\epsilon\), for either sign.

| \(\ell\) | Degree \(\ell^2\) | Regular degree \(4\ell^2\) |
| --- | --- | --- |
| 3 | \(17/56\) | \(17/252\) |
| 5 | \(1/8\) | \(3/100\) |

On the imaginary side the corresponding lower exponents,
without an \(\epsilon\) loss, are \(1/12,1/54\) for
\(\ell=3\) and \(1/60,1/250\) for \(\ell=5\).

For the upper bounds, Cauchy–Schwarz and the involution
\(d\mapsto\operatorname{sf}(td)\) reduce (1) to the
quadratic second moments at \(\ell=3,5\). The involution
changes the size window by at most the fixed factor \(t\).
The [signed quadratic \(3\)-torsion moments](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md)
and the [higher odd prime-path moments](ABC_COHERENT_PRIME_PATH_HIGHER_ODD_TORSION.md)
give \(\Theta=17/14\) and \(\Theta=3/2\), respectively.

For the lower bounds use the
[simultaneous small-norm class construction](QUADRATIC_TWIST_SHARED_SMALL_TORSION_NORMS.md).
It gives \(\gg_{t,\ell}X^{1/\ell}\) negative squarefree
radicands up to \(X\) with nontrivial \(\ell\)-torsion
in both fields. Select a nonzero character from each
ordinary class group to obtain cyclic unramified
\(N_i/K_i\), dihedral over \(\mathbf Q\) by inversion.
For these composita \(b=1\), and both discriminants are
bounded above by fixed multiples of the corresponding
powers of \(|d|\) in (2). At the fixed primes one uses
the bounded local discriminant exponents in fixed degree.
The normal closure has only boundedly many quadratic
subfields, so it accounts for only boundedly many
radicands of this construction. The normal closure is
determined by \(E\), so the same bounded multiplicity
holds for degree-\(\ell^2\) composita. Substitution of
the two discriminant powers proves the lower exponents.

In particular there are \(\gg Y^{1/10}\) pairs of
dihedral quintic fields with negative quadratic resolvents
of radicands \(d,5d\) prime to 5 in the first component,
\(|\operatorname{Disc}(E_1)|\le Y\), and
\[
                     |\operatorname{Disc}(E_2)|
                     =25|\operatorname{Disc}(E_1)|.
\]
Use the prime-to-five coefficient-one lower construction
in the [shared norm rigidity result](ABC_FIFTH_SHARED_NORM_RIGIDITY.md).
For unramified \(N_i/K_i\), the permutation character
discriminant relation gives
\(|\operatorname{Disc}(E_i)|=|\operatorname{Disc}(K_i)|^2\),
and the two quadratic discriminants have ratio 5.
The same quadratic resolvent cannot occur for distinct
first radicands, so the lower count survives.
[The individual dihedral discriminant relation](GENERALIZED_DIHEDRAL_COUNTS_FROM_QUADRATIC_TORSION.md)

These are counts for a prescribed twist between the
quadratic resolvents. The broader compositum method of
Wang and the counting problem for individual dihedral quintics
address different field families.
[Wang, Theorem 1.1 and §2](https://wangjiuya.github.io/research/malle.pdf)
[Individual dihedral quintics](https://msp.org/involve/2012/5-1/involve-v5-n1-p10-s.pdf)

## Exact cubic control

The irreducible cubics \(x^3-x^2+x+2\) and \(x^3-x^2-5\)
have discriminants \(-139\) and \(-695\), both fundamental.
Their normal closures are unramified over the quadratic
resolvents by the individual dihedral discriminant relation.
On the nine embeddings of the compositum, the shared prime
139 gives the pair of reflections with one fixed point.
At 5 only the second factor has a reflection, giving three
fixed points. Thus the degree-nine absolute discriminant
is \(139^4\cdot 5^3=46662630125\).
The same value follows by computing the field discriminant
of the polynomial compositum in PARI/GP.
