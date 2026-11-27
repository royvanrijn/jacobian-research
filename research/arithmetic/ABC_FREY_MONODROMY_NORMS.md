# Primitive Frey monodromy norms at primes in the level

Let \(a+b=c\) be primitive and positive, \(R=\operatorname{rad}(abc)\),
\(A=\log R\), \(h=\log c\), and \(E\) its Frey curve as in the
[degree comparison](ABC_CUBIC_ADJOINT_SHIMURA_PROJECTORS.md). Let
\(S_{\rm odd}\) be the odd primes of \(R\), \(k=|S_{\rm odd}|\), and
\(n_r=2v_r(abc)\) for \(r\in S_{\rm odd}\). Write \(N=N_E\). For
\(k\ge2\) take an odd admissible discriminant \(D\) and \(M=N/D\)
with at least two odd primes. Put \(J=J_0^D(M)\) and let
\(q:J\to E_0\) be the optimal elliptic quotient in the isogeny class
of \(E\). Write \(\delta=\delta_{D,M}\), so \(qq^\vee=[\delta]\)
with the canonical principal polarizations.

## Exact integral projector on the torus

Fix any \(p\in S_{\rm odd}\). The Jacobian has semistable reduction
at \(p\), by Čerednik–Drinfeld for \(p\mid D\) and by the
two-component reduction at \(p\parallel M\). The quotient \(E_0\)
has multiplicative reduction. Work on geometric character and
component groups, and let \(X=X_p(J)\). Let \(v\) be a primitive
generator of the saturation of \(q^*X_p(E_0)\) in \(X\). Define
\[
 H_{D,M,p}=\langle v,v\rangle,\qquad
 g_p=\gcd\{ \langle x,v\rangle:x\in X\},\qquad
 n_p^0=v_p(\Delta_{E_0}),
\]
using the monodromy pairing. Write \(i_p,j_p\) for the orders of the
image and cokernel of \(\Phi_p(J)\to\Phi_p(E_0)\). If \(P\) denotes the
rational projector \(q^\vee q/\delta\), let \(d_X(P)\) be the least
positive integer clearing its denominator on \(X\). Then
\[
 \boxed{\quad g_p=i_p,\qquad
 q^*X_p(E_0)=j_p\mathbb Z v,\qquad
 H_{D,M,p}=\frac{\delta n_p^0}{j_p^2}
          =\frac{\delta i_p^2}{n_p^0},\qquad
 d_X(P)=\frac{H_{D,M,p}}{i_p}=\frac{\delta}{j_p}.\quad}
 \tag{MN1}
\]

Indeed for \(L=\mathbb Z v\), the saturated character line,
the image of the map \(X\to\operatorname{Hom}(L,\mathbb Z)\)
defined by the pairing is \(g_p\mathbb Z\) in the dual basis.
Its cokernel has order \(g_p\), and the image of \(L\) has
index \(H_{D,M,p}/g_p\) inside the image of \(X\).
Conrad–Stein Proposition 6.5 identifies this cokernel with the
image on components, and their Theorem 6.1 gives
\(n_p^0/\delta=g_p^2/H_{D,M,p}\). Here their \(m_A\)
is \(\delta\) because the induced polarization on \(E_0\) is
multiplication by \(\delta\). If a generator of \(X_p(E_0)\)
pulls back to \(s v\) with \(s>0\), then the index of the
pairing image of \(\mathbb Z s v\) in that of \(X\), still
evaluating on \(L\), is \(sH_{D,M,p}/g_p\). Applying
Theorem 6.1 with this sublattice and \(L\) shows
\(s=n_p^0/g_p=j_p\), using their Lemma 6.3 to identify
the cokernel for \(q^*X_p(E_0)\) with \(\Phi_p(E_0)\). Functoriality
makes the projector self-adjoint for monodromy. Its action is
\(x\mapsto v\langle x,v\rangle/H_{D,M,p}\). Since \(v\) is
primitive, its denominator is \(H_{D,M,p}/g_p\).
[Conrad–Stein, Lemma 6.3, Theorem 6.1 and Proposition 6.5](https://wstein.org/papers/compgrp/compgrp.pdf).

At an odd prime \(p\mid D\), Pasten's Theorem 6.17 bounds \(j_p\)
by an absolute constant. At \(p\parallel M\), Pasten's Lemma
6.14 bounds \(i_p\) by an absolute constant, with
\(S=\{2\}\). Both factors are supported at coefficient primes
at most \(163\), and \(n_p^0/n_p\) has height at most \(\log 163\)
with the same support by his Lemma 6.8. In particular,
\[
 H_{D,M,p}\asymp
 \begin{cases}
   \delta_{D,M} n_p,&p\mid D,\\
   \delta_{D,M}/n_p,&p\parallel M.
 \end{cases}
 \tag{MN2}
\]
The quotient of either side by its indicated comparison is a
positive rational number with uniformly bounded height and the
specified prime support. At a level prime with large \(n_p\),
the lattice projector can therefore have a substantially smaller
denominator than the endomorphism projector.
[Pasten, Lemmas 6.8 and 6.14 and Theorem 6.17](https://arxiv.org/abs/1705.09251).

## Two vertices suffice for the norm bound

On a geometric semistable dual graph, \(X\) is the integral cycle
lattice, with pairing given by the sum of edge products weighted
by node thicknesses. An edge of thickness \(w\) can equivalently
be subdivided into a chain of \(w\) unit edges, with the same cycle
coefficient on the entire chain. In the case \(p\parallel M\),
the graph before these subdivisions has two vertices and one edge
per supersingular point of the Shimura curve of level \(M/p\)
at \(p\). For \(D=1\) this is the modular reduction graph.
Orient all the edges from one vertex to the other. Writing \(v_e\)
for the cycle coordinates and \(w_e\) for the thicknesses, we get
\[
       \sum_e v_e=0,\qquad \gcd_e v_e=1,\qquad
       H_{D,M,p}=\sum_e w_e v_e^2.
 \tag{MN3}
\]
The cycle lattice is saturated in the edge lattice, so its
primitive vector has content one in edge coordinates. Taking
differences of two edges as the cycle generators also gives
\[
   i_p=\gcd_{e,e'}(w_e v_e-w_{e'}v_{e'}).
\]
In particular for \(\ell>163\) the weighted coordinates are
nonconstant modulo \(\ell\). Divisibility of the norm at such
a prime records isotropy of a nonzero cycle with nonzero pairing
functional modulo \(\ell\). Rational Jacquet–Langlands and strong multiplicity one identify the line
by the good Hecke eigenvalues of \(E\). The multiplicity of this
full-level newform is one on the toric character space. Thus the
input to the norm is a rational Hecke line in the degree-zero
supersingular module, with its primitive integral normalization.
[Boutot–Carayol on uniformization and Illusie, equations (2.1.1), (2.1.3), and Theorem 2.7](https://www.numdam.org/issue/AST_1991__196-197__1_0.pdf).
[Conrad–Stein, §7.2 on the classical modular reduction](https://wstein.org/papers/compgrp/compgrp.pdf).
[Longo–Rotger–Vigni, §§2.2 and 3.5 on the supersingular module at \(p\) exactly dividing the level](https://arxiv.org/abs/0910.3391).

Every \(n_p\) is \(O(h+1)\), since \(p^{n_p/2}\mid abc<c^3\).
Equations (MN1) and (SP9) in the
[degree comparison](ABC_CUBIC_ADJOINT_SHIMURA_PROJECTORS.md#uniform-shimura-transfer)
give, uniformly for every permitted \(D,M,p\),
\[
 \log H_{D,M,p}=h+A+O\!\left(\frac{h}{\log(h+2)}\right).
 \tag{MN4}
\]
In fact \(1\le j_p\le n_p^0\ll h+1\), which suffices for this
logarithmic conclusion. For any one choice of \(D,M,p\) per triple,
or the minimum over any nonempty collection, the uniform estimate
\(H_{D,M,p}\ll_\varepsilon R^{2+\varepsilon}\) for every
\(\varepsilon>0\) is equivalent to abc. In the reverse implication,
bound the error in (MN4) by \(\eta h + O_\eta(1)\) and choose
\(\eta\) and the exponent slack small enough for the desired abc
exponent. For one explicit selection keep the two smallest odd
primes in \(M\) when \(k\) is even and the three smallest when \(k\)
is odd, and take \(p\) to be the smallest. This uses only graphs of
the form (MN3), including \(D=1\) when \(k=2,3\). The cases with
at most one odd prime follow from the
[elementary prime-power reduction](ABC_CUBIC_ADJOINT_SHIMURA_PROJECTORS.md#a-small-level-projector-target).

## Gcd across the two-vertex graphs

Let \(m=2\) when \(k\) is even and \(m=3\) when \(k\) is odd.
Move over every \(m\)-element subset \(U\subset S_{\rm odd}\),
set \(M=N_2\prod_{r\in U}r\) with \(N_2\) the \(2\)-part of \(N\),
and take every \(p\in U\). Define
\[
 W=\frac{\delta_{1,N}}{\prod_{r\in S_{\rm odd}} n_r},\qquad
 G=\gcd_{U,\,p\in U} H_{N/M,M,p},\qquad
 T=\gcd_{\substack{V\subset S_{\rm odd}\\|V|=m-1}}
                  \prod_{r\in V} n_r.
\]
The exact comparison for each choice, from Pasten's global degree
transfer and (MN2), is
\[
 H_{N/M,M,p}=W\left(\prod_{r\in U\setminus\{p\}}n_r\right)c_{U,p}
 \tag{MN5}
\]
where \(c_{U,p}>0\) is rational of logarithmic height \(O(k+1)\),
supported on coefficient primes at most \(163\). Thus
\[
 G=W T\,\xi,\qquad
 h_{\mathbb Q}(\xi)=O(k+1),\qquad
 \operatorname{supp}(\xi)\subset\{\ell:\ell\text{ prime},\ \ell\le163\}.
 \tag{MN6}
\]
Indeed every \((m-1)\)-element subset occurs as \(U\setminus\{p\}\).
For each individual graph, dividing (MN5) by its spectator product
recovers every valuation of \(W\) at coefficient primes above \(163\).
At each coefficient prime \(\ell\), pass to the minimum of the
valuations in (MN5). The difference from the minimum without the
\(c_{U,p}\) terms is at most the maximum of
\(|v_\ell(c_{U,p})|\). Summing over the fixed support proves
the height bound, even when \(W\) has a denominator.
[Pasten, Theorem 6.1(b.2) and Corollary 6.2](https://arxiv.org/abs/1705.09251).

For even \(k\), \(T=\gcd_r n_r\) divides an absolute power of two.
Pasten's Lemma 6.12 rules out every odd prime factor. Lemma 6.10
with \(S=\{2\}\) and \(L=8\) rules out \(8\mid T\) outside a
finite set of Frey curves, whose contribution to \(T\) is
bounded. Consequently \(v_\ell(G)=v_\ell(W)\) for every
\(\ell>163\) when \(k\) is even. For odd \(k\), write
\(e_{\ell,1}\le e_{\ell,2}\le\cdots\le e_{\ell,k}\) for the
ordered list of valuations \(v_\ell(n_r)\). Then
\(v_\ell(T)=e_{\ell,1}+e_{\ell,2}\). For \(\ell>163\) the
first term vanishes and
\[
 (k-1)e_{\ell,2}\le\sum_{r\in S_{\rm odd}}v_\ell(n_r).
\]
This last inequality holds at every prime, and the product
contributed by the \(e_{\ell,1}\) is uniformly bounded by the
same argument as in the even case. It follows that
\(\log T\ll 1+A/(k-1)\) unconditionally using Pasten's
product-of-valuations bound (with a fixed exponent slack).
In both parity cases one also has \(\log T\ll\log(h+2)\)
from the definition by products of at most two \(n_r\)'s.
[Pasten, Lemmas 6.10 and 6.12 and Theorem 16.7](https://arxiv.org/abs/1705.09251).

In particular \(G=W\exp(O(k+\log(h+2)))\) and
\[
 \log G=h+A+O\!\left(\frac{h}{\log(h+2)}\right).
 \tag{MN7}
\]
Here one uses the classical degree estimate and the elementary
divisor bound on \(\sum_r\log n_r\) as in (SP9).
Thus a sharp bound \(G\ll_\varepsilon R^{2+\varepsilon}\)
is itself equivalent to abc. This gcd uses \(m\binom{k}{m}\)
norms. It requires a bound on shared prime powers of the norms
and permits a different choice of \(U,p\) at each coefficient
prime. On families with \(k\to\infty\), (MN6) has total
multiplicative error \(R^{o(1)}\) relative to \(W\), since
\(k=o(A)\) and \(\log T=o(A)\).

There is also a single small-pivot option with no gcd. By
Theorem 16.7, \(\sum_{r\in S_{\rm odd}}\log(r n_r)\ll A\).
Choose the \(m\) primes with least \(r n_r\), then choose any
\(p\) in the chosen set. This gives
\[
 M\prod_{r\in U}n_r\le R^{O(1/k)},\qquad
 p n_p\le R^{O(1/k)}.
\]
The bounded \(2\)-part is absorbed since \(A/k\ge\log 3\).
Equation (MN5) gives \(H_{N/M,M,p}=W\exp(O(k+A/k))\)
for this prescription. These bounds concern \(M\) and the
removed local factors, while the size of the primitive
supersingular cycle still requires an arithmetic estimate.


## A modular example

For \((a,b,c)=(1,80,81)\), the Frey curve has conductor \(240\)
and optimal quotient \(240\mathrm{d}1\). Here \(D=1\). At \(p=3\)
and \(p=5\) the supersingular modules have 24 and 32 coordinates,
respectively. Intersecting the good Hecke eigenspaces gives
primitive vectors of squared norm \(16\) in both modules, with
pairing gcd \(1\). Each optimal discriminant exponent is \(1\).
An independent calculation of the degree within the same Sage
session, via the modular-symbol lattice, also gives \(16\).
The following is an exact check for this example.
[Cremona's elliptic curve data and optimality convention](https://johncremona.github.io/ecdata/).

```python
from sage.all import BrandtModule, EllipticCurve, QQ, VectorSpace, gcd

frey = EllipticCurve([0, 79, 0, -80, 0])
optimal = EllipticCurve([0, 1, 0, 0, -12])
assert frey.conductor() == optimal.conductor() == 240
degree = frey.modular_symbol_space(sign=0).abelian_variety().modular_degree()
assert degree == 16
for p in (3, 5):
    brandt = BrandtModule(p, 240 // p)
    space = VectorSpace(QQ, brandt.dimension())
    for ell in (7, 11, 13, 17, 19, 23):
        operator = brandt.hecke_matrix(ell)
        space = space.intersection((operator - frey.ap(ell)).left_kernel())
        if space.dimension() == 1:
            break
    assert space.dimension() == 1
    v = space.basis()[0]
    v *= v.denominator()
    v /= gcd(list(v))
    weighted = [weight * entry for weight, entry in zip(brandt.monodromy_weights(), v)]
    assert sum(v) == 0 and gcd(list(v)) == 1
    image = gcd([entry - weighted[0] for entry in weighted])
    norm = sum(entry * product for entry, product in zip(v, weighted))
    assert (norm, image, optimal.discriminant().valuation(p)) == (16, 1, 1)
    print(p, ell, brandt.dimension(), norm, image)
```

## A limitation of spectral bounds on a cycle lattice

Even on the cycle lattice of a connected bipartite graph with
unit weights and maximum degree four, a bounded, sparse integral
self-adjoint operator can have a simple eigenvalue with a fixed
spectral gap and an exponentially large primitive eigenvector
norm. For \(t\ge2\), replace each edge in a path of length \(t\)
by two parallel oriented edges. Its cycle basis consists of the
differences of the two parallel edges at each step, with pairing
matrix \(2I_t\). Put \(B:\mathbb Z^t\to\mathbb Z^{t-1}\) with
\((Bx)_j=x_{j+1}-2x_j\). The endomorphism \(B^{\mathsf T}B\)
has simple kernel generated by \((1,2,\ldots,2^{t-1})\).
Its other eigenvalues are those of \(BB^{\mathsf T}\), which
lie in \([1,9]\) by the diagonal \(5\) and off-diagonals \(-2\).
The primitive squared norm is \(2(4^t-1)/3\) and the denominator
of its orthogonal projector is \((4^t-1)/3\).
This construction concerns abstract cycle lattices and supplies
no arithmetic example. A bound on the shared prime powers in
(MN6), or one on a primitive cycle in (MN3), must use additional
arithmetic information beyond these spectral properties.

