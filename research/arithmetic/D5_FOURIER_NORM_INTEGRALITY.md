# Integral generators in the dihedral quintic norm equation

Larson and Rolen [Theorem 1.1 and §§4–5](https://msp.org/involve/2012/5-1/involve-v5-n1-p10-s.pdf)
attach a quadratic norm triple to a small trace-zero integral generator of
a dihedral quintic field and reconstruct at most one such field from the
triple. Their reconstruction bounds the number of fields by a count of
integral norm triples. We use the same variables to characterize when the
reconstructed roots themselves are integral.

Let \(z=\zeta_5\), \(F=\mathbf Q(\sqrt5)\), and
\(\sqrt5=1+2(z+z^{-1})\). Use \(A'\) and \(B'\) for conjugation in \(F\),
and write \(\operatorname{Tr}\) and \(N\) for the trace and norm of \(F\).
The norm equation is
\[
  N(B^2-4A^2 A')=5C^2,\qquad A,B\in F,\quad C\in\mathbf Q.
  \tag{1}
\]

## Rational reconstruction

**Theorem DF1.** For each solution of (1) with \(A C\ne0\), the Fourier
equations have one dihedral orbit of ordered trace-zero root vectors.
The reconstructed monic quintic has rational coefficients. If irreducible
its Galois group is \(C_5\) or \(D_5\).

To prove this, let \(V_j=\sum_{i=0}^4 z^{ij}x_i\) for \(1\le j\le4\) and
put
\[
\begin{split}
 A&=V_2V_3,\qquad A'=V_1V_4,\\
 (M_1,M_2,M_3,M_4)&=(V_1 V_2^2,\ V_3^2V_4,\ V_2 V_4^2,\ V_1^2V_3),\\
 B&=M_1+M_2,\qquad B'=M_3+M_4,\qquad
 \sqrt5 C=(M_1-M_2)(M_3-M_4).
\end{split}
\tag{2}
\]
Given (1), choose \(\delta_1\) with
\(\delta_1^2=B^2-4A^2 A'\) and set \(\delta_2=\sqrt5 C/\delta_1\);
then \(\delta_2^2=(B')^2-4(A')^2 A\). The pairs of \(M_i\) are
\((B\pm\delta_1)/2\) and \((B'\pm\delta_2)/2\), in that order. Their
products are \(A^2A'\) and \((A')^2 A\), so all four are nonzero.
Take
\[
 V_1^5=M_1M_4^2/A^2,\qquad
 V_3=M_4/V_1^2,\qquad V_4=A'/V_1,\qquad V_2=A/V_3 .
 \tag{3}
\]
Substitution recovers (2). The five choices of \(V_1\) differ by the
rotation \(V_j\mapsto z^j V_j\). Reversing both signs of the \(\delta_i\)
gives the reflection \(V_j\mapsto V_{5-j}\). Inverse Fourier transform
\(x_i=5^{-1}\sum_{j=1}^4 z^{-ij}V_j\) gives the claimed orbit.

For rational descent, replacing \(z\) by \(z^2\) at a fixed ordered
\(x\)-vector sends \((A,B,C)\) to \((A',B',C)\): it sends the two
differences in (2) to the second difference and the negative of the
first, respectively. Replacing \(z\) by \(z^4\) leaves the triple
unchanged. Therefore for any absolute Galois automorphism \(\sigma\),
the \(x\)-vector \(\sigma(x)\), evaluated with \(z\), has the same
invariants as \(x\). The sole dihedral orbit is stable even when some
entries coincide. If the polynomial is irreducible, the permutation
action on its distinct roots is transitive and lies in \(D_5\).
For a genuine \(D_5\) quintic generator, \(A\ne 0\) by Larson–Rolen
Lemma 4.3. Also \(C\ne0\): otherwise both differences are zero,
and (3) gives at most five ordered vectors, whereas a \(D_5\) generator
has ten distinct ordered vectors in its orbit.

## The exact integral conditions

**Theorem DF2.** Assume \(A,B\in\mathcal O_F\) and \(C\in\mathbf Z\)
in Theorem DF1. Define \(D,E\in F\) by
\[
  D=\frac{BB'-\sqrt5 C}{2A},\qquad
  E=\frac{B'D-B(A')^2}{A}.
  \tag{4}
\]
The Fourier coordinates are all algebraic integers if and only if
\(E\in\mathcal O_F\), which also implies \(D\in\mathcal O_F\).
With integral Fourier coordinates, the four conditions
\[
\begin{split}
 \operatorname{Tr}(A)&\equiv0\pmod 5,\\
 \operatorname{Tr}(B)&\equiv0\pmod {25},\\
 \operatorname{Tr}(A^2)-N(A)-\operatorname{Tr}(D)
      &\equiv0\pmod {125},\\
 \operatorname{Tr}(E)-5\operatorname{Tr}(AB)+5\operatorname{Tr}(A'B)
      &\equiv0\pmod {3125}
\end{split}
\tag{5}
\]
are equivalent to integrality of all five roots. In particular (4)
gives two nested divisibilities in \(\mathcal O_F\), and (5) is a test
in the rational integers.

For the first statement, the identities
\[
 AD=M_1M_4+M_2M_3,\qquad
 D=V_1^3V_2+V_3V_4^3,\qquad E=V_1^5+V_4^5
 \tag{6}
\]
follow by expanding (2) and (4). If \(E\) is integral, \(V_1^5,V_4^5\)
are roots of the monic polynomial with sum \(E\) and product \((A')^5\).
The conjugate polynomial gives the same assertion for \(V_2^5,V_3^5\).
For this last identity one can replace \(z\) by \(z^2\) in (6), which
conjugates \(D,E\) from (4). Thus all \(V_j\) are integral. Equation
(6) proves the necessity and the assertion about \(D\).

For (5), with \(P_k=\sum_i x_i^k\), the inverse Fourier formula gives
\(5^{k-1}P_k\) as the sum of the ordered degree \(k\) monomials with
index sum zero modulo \(5\). For \(k=2,3,4,5\) these numerators are
respectively
\[
  2\operatorname{Tr}(A),\quad 3\operatorname{Tr}(B),\quad
  4\operatorname{Tr}(D)+6\operatorname{Tr}(A^2)+24 N(A),\quad
  \operatorname{Tr}(E)+20\operatorname{Tr}(AB)+30\operatorname{Tr}(A'B).
\]
Newton's identities give the four nonleading coefficients of the
trace-zero polynomial, in descending degree:
\[
 -\frac{\operatorname{Tr}(A)}5,\quad
 -\frac{\operatorname{Tr}(B)}{25},\quad
 \frac{\operatorname{Tr}(A^2)-N(A)-\operatorname{Tr}(D)}{125},\quad
 -\frac{\operatorname{Tr}(E)-5\operatorname{Tr}(AB)+5\operatorname{Tr}(A'B)}{3125}.
 \tag{7}
\]
A monic rational polynomial has integral roots exactly when its
coefficients are integers.

There is also a local version of (5). In any field containing the
\(V_j\) and \(z\), put \(\pi=1-z\). At every place above \(5\) the
inverse Fourier transform is integral exactly when
\(\sum_{j=1}^4 j^k V_j\in\pi^{4-k}\mathcal O_{\mathrm{loc}}\) for
\(k=0,1,2,3\), where \(\mathcal O_{\mathrm{loc}}\) is the local
valuation ring. Expand \((1-\pi)^{-ij}\) modulo \(\pi^4\). As a
polynomial of degree at most three in \(i\), its four values at
\(i=0,1,2,3\) detect every coefficient modulo \(\pi^4\), since
their Vandermonde determinant \(12\) is a unit. The degree \(k\)
coefficient of the sum over \(j\) is \(\pi^k\sum j^kV_j\) times a
unit, by the leading term of each binomial polynomial. Also
\((5)=(\pi^4)\) locally. At every other finite place the Fourier
transform and its inverse preserve integrality.

## The two divisibilities and the prime support

Outside primes dividing \(AA'\) the \(V_j\) are integral: each \(M_i\)
is integral, being a root of \(T^2-BT+A^2A'\) or its conjugate,
and the fifth power reconstruction (3) and its conjugates and
reflections have unit denominators.

At a prime \(\mathfrak p\nmid 10\) of \(F\) at which
\(v_{\mathfrak p}(A)=a>0\) and \(A',B,B'\) are units,
the \(V_j\) are integral locally if and only if
\[
 \sqrt5 C\equiv BB'\pmod{\mathfrak p^a}.
 \tag{8}
\]
Necessity follows from (4). For sufficiency put \(S=\sqrt5 C\).
The norm equation gives
\[
 (BB'-S)(BB'+S)=
 4AA'\big(A(B')^2+A'B^2-4A^2(A')^2\big).
\]
If (8) holds, \(BB'+S\) is a unit and
\(D\equiv (A')^2B/B'\pmod{\mathfrak p^a}\). Thus
\(B'D-B(A')^2\) is divisible by \(A\) locally. The conjugate formula
for \(E'\) has unit denominator there, since both \(2\) and \(A'\)
are units. Both fifth power equations now give integrality of the
\(V_j\) locally. This test chooses one of the two unit square-root
signs of \(5C^2\) modulo \(\mathfrak p^a\). The second divisibility
adds no local restriction under these unit hypotheses.

## Counterfamilies and exact controls

**Theorem DF3.** Fix a \(D_5\) quintic field. There are infinitely
many solutions of (1) with \(A,B\in\mathcal O_F\), nonzero
\(A,C\), and \(C\in\mathbf Z\) whose reconstructed polynomial is
irreducible and generates this field, in each of the following two
families: \(D\notin\mathcal O_F\); and \(D\in\mathcal O_F\) but
\(E\notin\mathcal O_F\).

Let \(p\ne 5\) split completely in the composite of the normal closure
of the fixed field with \(\mathbf Q(z)\). These primes are infinite by
Chebotarev. At a chosen prime of the composite above \(p\), take any
\(D_5\) ordering of the five embeddings. The five residues of an
integral element \(\gamma\) modulo \(p^{k+2}\) can be prescribed
independently by the Chinese remainder theorem. Fourier transform
over \(\mathbf Z/p^{k+2}\mathbf Z\) is invertible, since \(p\ne5\)
and \(z\) reduces to a primitive fifth root. Choose \(\gamma\)
whose four nontrivial Fourier modes have valuations
\((0,k+1,k+1,k+1)\), with \(k=2\) or \(3\). For instance choose their
residues to be \(1,p^{k+1},p^{k+1},p^{k+1}\).
Let \(\beta=5\gamma-\operatorname{Tr}_{K/\mathbf Q}(\gamma)\).
This ensures trace zero and scales the four modes by a unit at \(p\).
The element \(\alpha=\beta/p\) has Fourier valuations
\((-1,k,k,k)\) at the chosen prime, and permutations thereof at
every prime above \(p\). Indeed the Galois group acts on the
indices by the dihedral permutations and the cyclotomic
permutations, with factors that are fifth roots of unity.
The modes are integral away from \(p\).
The element \(\alpha\) is nonzero by the prescribed modes, so
trace zero makes it nonrational and hence a generator of the
prime-degree field.

At \(p\), every monomial of \(A,A'\) and of the \(M_i\) is
integral, as (2) has no monomial with a power greater than \(2\)
of a single mode. Equation (2) then gives \(A,B\in\mathcal O_F\)
and \(C\sqrt5\in\mathcal O_F\). The latter gives \(C\in\mathbf Z\)
because \(C\) is rational and the valuation of \(\sqrt5\) at the
prime above \(5\) in \(F\) is \(1\). At the chosen prime both
\(M_1-M_2\) and \(M_3-M_4\) are nonzero by unequal valuations.
For \(k=2\) the two terms of \(D\) in (6) have valuations
\(-1\) and \(8\). For \(k=3\) every summand of \(D\) is integral
at every place, while the two terms of \(E\) have valuations
\(-5\) and \(15\) at the chosen place. In both cases the root vector
fails integrality only above \(p\), so distinct primes give distinct
triples by Theorem DF1. This argument allows arbitrarily large
denominators while the field is fixed.

An exact small control uses \(u=-1-z-z^2+z^3\), with norm \(11\) in
\(\mathbf Q(z)\). Put \(V_j=\sigma_j(v)\), where
\(\sigma_j(z)=z^j\). The choices \(v=11^k/u^{k+1}\) for \(k=2,3\)
give the same valuation pattern above \(11\) as the two
counterfamilies. All of \(A,B,C\) from (2) are integral and
\(AC\ne0\). The coordinate denominators of \(D,E\) in the integral
basis \(1,z,z^2,z^3\) are respectively \((11,11^5)\) and
\((1,11^5)\). The choice \(v=u\) makes all Fourier modes integral
and keeps \(AC\ne0\), but \(x_0=\operatorname{Tr}_{\mathbf Q(z)/\mathbf Q}(u)/5
=-3/5\), so (5) adds a genuine condition. These three controls
use rational reconstructed roots. The fixed-field construction
above works, for example, with the field generated by a root of
\(T^5-2T^4+2T^3-T^2+1\), which has Galois group \(D_5\) and field
discriminant \(2209\).

## Relation to counting

For Larson–Rolen's archimedean generator box, bounding every
\(|x_i|\) by a constant times \(Y^{1/8}\) bounds both embeddings of
\(D,E\) by constants times \(Y^{1/2},Y^{5/8}\), respectively.
Equations (4)–(5) can therefore be imposed together with the
original box, and (7) permits an irreducibility and discriminant
test on the survivors. The counterfamilies show a strict
difference between integral norm triples and integral generators.
They do not estimate the proportion of inadmissible triples in
that box. A power saving would require controlling the signed
divisibility (8) or the nonunit sectors of (4) as \(A,B\) vary.
