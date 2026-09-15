# Two structural gates for a rational rank-19 marking

These gates exclude 119 further rows of the retained foundry queue. They
do not construct a new MW17 surface. The first applies uniformly to a
family of transcendental lattices; the second uses the proved real-period
criterion and a fixed finite congruence audit.

## 1. The scaled split family

For N a positive integer, let

\[
T_N=2(U\oplus\langle2N\rangle),\qquad
G_N=\begin{pmatrix}0&0&2\\0&4N&0\\2&0&0\end{pmatrix}.
\]

The full canonical marked curve is **X_0(4N) over Q**. In particular,
for **N at least five**, no K3 over Q of geometric Picard rank 19 has
this actual T and its full geometric NS rationally marked.

This generalizes the [determinant-800 proof](DET800_FULL_MARKING_OBSTRUCTION_2026-09-14.md),
including its arithmetic descent, rather than generalizing just its
norm-one calculation.

### Integral and local groups

The isometries of T_N and T_N/2 are the same. The even Clifford order of
the primitive form ef+Nh² is the split Eichler order R^0(N). Its local
normalizer consists of units and, at each prime dividing N, the coset
interchanging the two endpoint maximal orders. Consequently its rational
positive normalizer has the exact-divisor Atkin–Lehner cosets. This is the
same order-normalizer argument used in the determinant-800 proof; it
applies to every prime power in N, not just to squarefree N.

For g=[[a,Nb],[c,d]], delta=ad-Nbc, the action on T is

\[
A(g)=\delta^{-1}\begin{pmatrix}
a^2&-2Nab&-Nb^2\\
-ac&ad+Nbc&bd\\
-Nc^2&2Ncd&d^2
\end{pmatrix}.
\]

The dual-generator denominators are (2,4N,2). At 2, let delta be a unit.
Stability means that (A-I)G_N^-1 is 2-integral. Its entries involving
the h-dual vector give ab, cd, and bc even. Since ad-Nbc is odd, bc even
forces a,d odd, and the first two conditions force b,c even. Conversely
these parity conditions make every entry 2-integral. Thus the stable
unit subgroup is exactly b=c=0 modulo two, for every N.

There are no stable determinant-minus-one orthogonal elements at 2.
For -A to be stable, its action modulo two must again be identity. If N
is odd, the trace-zero adjoint representation of GL2(F2) is faithful, so
g must be scalar modulo two. If N is even, the off-diagonal h-coordinates
force b,c even directly, since a,d are already odd. In both cases a,d
are odd and b,c are even. Then

\[
A_{hh}=1+2Nbc/\delta\equiv1\pmod{4N\mathbf Z_2},
\]

so -A cannot act trivially on the cyclic 2-primary h component, whose
order is at least four.

If N is even, the remaining 2-adic normalizer coset exchanges e and f
modulo two, whereas every Eichler-unit action fixes both e and f there
(it may add h-coordinates). The exchange cannot be removed by a unit or
by changing the common sign. Hence that coset has no stable element.
If N is odd the local order at 2 is maximal, so there is no such coset.

At an odd prime p dividing N, units act trivially on A_(T,p); the endpoint
exchange acts as -1 on its cyclic p-primary discriminant group. Thus no
Atkin–Lehner coset survives in SO at an odd bad prime. At the other
primes the level is maximal.

### Arithmetic descent and exclusion

The absence of any stable determinant-minus-one element at 2 kills the
single compatible global orientation character of a full rational NS
marking. One cannot restore an odd-prime exchange by independently
choosing a minus sign at that prime. The stable group therefore consists
of units everywhere, with the above parity condition at 2. All determinant
units occur: diagonal matrices diag(u,1), for odd u, already verify this
at 2. At odd primes it is the ordinary Eichler-unit determinant image.

In Gamma_0(N) coordinates the resulting level is b even and c divisible
by 2N. Rational conjugation by diag(1/2,1) changes it into the usual
Gamma_0(4N) level, including its determinant units. Hence the canonical
curve is X_0(4N), not an unspecified form of its complex curve. The
canonical arithmetic period-map and stable-group inputs are the
Dolgachev/Rizov results cited in the determinant-800 proof.

The rational cyclic-isogeny classification permits only degrees
1 through 19,21,25,27,37,43,67,163. Its only positive multiples of four
are 4,8,12,16. See
[Banwait–Najman–Padurariu, Theorem 1.1 and Table 1.1](https://arxiv.org/html/2206.08891v3).
Therefore X_0(4N)(Q) has no noncuspidal point for N at least five.
The small values N=1,2,3,4 are not excluded by this theorem; this note
does not claim an MW17 fibration for any of them.

### Exact catalogue transports

The certificate finds 84 retained rows with literal content two, even
primitive half-lattice, and a retained primitive isotropic vector of
divisibility one. For each, it constructs an integral vector w with
(e,w)=1, sets f=w-(w²/2)e, and takes the primitive generator h of the
orthogonal complement of e,f. It checks the resulting basis has
determinant plus or minus one and transports the literal Gram matrix
to G_N. Their N values lie between 35 and 303.

All 84 are excluded. Determinant 800 was already excluded, so this
removes **83 additional queue rows**, including determinants 560 and 656.
The calculation uses actual integral isometries, not determinant labels
or similarity classes. It does not claim that the retained isotropic
vectors find every split half-lattice in the catalogue.

## 2. Congruence obstruction to the required square-two vector

The [rank-19 real-period theorem](REAL_RANK19_PERIOD_OBSTRUCTION_2026-09-14.md)
requires a rationally anisotropic T to represent 2 integrally. The audit
tests the fixed moduli 8,16,3,5,7,11,13 on the 218 anisotropic rows in the
721-row frontier before this note. For each modulus it enumerates all
three-coordinate residues until a witness exists or that finite space is
exhausted. A missing residue is an exact integral-representation obstruction.
Passing every modulus makes no claim of an integral representation.

Thirty-six rows fail: 17 at modulus 8, five at 16, thirteen at 3, and one
at 7. For each excluded row the checker independently recomputes the
rational diagonal form and its Clifford Hilbert symbol at the retained
anisotropy prime, obtaining -1. Thus both necessary premises are checked:
rational anisotropy and failure to represent 2 integrally.

These 36 rows are disjoint from the split-family exclusions. The combined
frontier has **602 unresolved rows**, no surviving member of the original
21-row low-genus shortlist, and no new positive marking witness.

## Replay and limits

```sh
sage -python research/elkies-k3/scripts/certify_marking_structural_gates.py --check
```

The [packet](../artifacts/generated-results/elkies-k3-marking-structural-gates-v1.json)
retains every integral transport, local obstruction, input hash, and the
ordered surviving frontier. The 13 local-level controls include N=1
through 8 and N=10,12,35,41,50; the uniform claim uses the written parity
and normalizer argument, not an extrapolation from those controls.
The tests are exact and cheap, with no rootless-frame, equation, or
specialization campaign. No independent implementation, formal
verification, external review, or novelty claim is made.
