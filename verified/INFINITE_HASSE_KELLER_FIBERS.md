# Infinitely many Hasse failures in one fixed Keller map

## 1. The theorem

There is one polynomial map

\[
 F:\mathbb A^3_{\mathbb Q}\longrightarrow\mathbb A^3_{\mathbb Q},
 \qquad \det DF=-2,
\]

and an infinite set of rational targets \(y_\ell\), indexed by the rational
primes

\[
 \ell\equiv1\pmod {27},
\]

such that every complete fiber \(F^{-1}(y_\ell)\) is a reduced scheme of
degree five satisfying

\[
 F^{-1}(y_\ell)(\mathbb Q_v)\ne\varnothing
 \quad\text{for every place }v\text{ of }\mathbb Q,
 \qquad
 F^{-1}(y_\ell)(\mathbb Q)=\varnothing.
\]

Thus one fixed Keller map has infinitely many everywhere locally soluble
rational target fibers with no rational point.

This converts the classical Berend--Bilu example from one isolated fiber
into an arithmetic progression of fibers of a single map.  The key is to
use the quadratic-tilt Keller pencil rather than the linear-tilt pencil.

## 2. The fixed map

Put

\[
 G(S)=S^5-\frac32S^4+\frac32S^3-\frac54S^2+\frac9{16}S.
\]

Its relevant coefficients are

\[
 g_1=\frac9{16},\qquad g_2=-\frac54,\qquad
 g_3=\frac32,\qquad g_4=-\frac32,\qquad g_5=1.
\]

For source coordinates \((x,y,z)\), set

\[
 t=1+xy,\qquad
 q=t^2z+\frac38y^2(1+3t),
\]

and define

\[
\begin{aligned}
 F_1={}&tq,\\
 F_2={}&y+8xq-\frac{40}{9}tq
          -\frac{32}{3}t^2x^2q^4
          +\frac{80}{9}t^2x^3q^5,\\
 F_3={}&x(5-3t)-\frac83x^3z
          +\frac{16}{3}(xq)^4-\frac{16}{3}(xq)^5.
\end{aligned} \tag{2.1}
\]

This is formula (26) of the root-engineered quadratic gauge specialized to
\(G\).  Its determinant identity gives

\[
 \boxed{\det DF=-2}. \tag{2.2}
\]

For a target denoted \((P,B,C)\), the inverse polynomial on the chart
\(P\ne0\) is

\[
 E_{P,B,C}(S)
 =G_P(S)-\frac9{32}(BS^2+C), \tag{2.3}
\]

where \(G_1=G\).  On the source incidence,

\[
 \frac{\partial E_{P,B,C}}{\partial S}=\frac9{16}D,
 \qquad D=\frac1t. \tag{2.4}
\]

Consequently every simple root of (2.3) reconstructs one source point.
When \(P=1\), every source point in the fiber satisfies \(tq=1\), so there
is no omitted boundary chart.

## 3. The infinite target line

Write \(X=S-\tfrac12\).  The elementary identity

\[
\begin{aligned}
 X^3(X^2+X+1)
 &=(S-\tfrac12)^3(S^2+\tfrac34)\\
 &=G(S)-\frac3{32}
\end{aligned} \tag{3.1}
\]

shows why the quadratic pencil is the correct one.  For a rational number
\(a\), take

\[
 y_a=(P,B,C)
 =\left(1,\frac{32a}{9},\frac{8a+1}{3}\right). \tag{3.2}
\]

Substitution in (2.3) gives

\[
\boxed{
 E_{y_a}(S)
 =\left((S-\tfrac12)^3-a\right)
  \left(S^2+\frac34\right).
} \tag{3.3}
\]

After the translation \(X=S-\tfrac12\), this is exactly

\[
 (X^3-a)(X^2+X+1). \tag{3.4}
\]

For a prime \(\ell\), the two factors in (3.4) are irreducible over
\(\mathbb Q\): the cubic has no rational root and the quadratic has
discriminant \(-3\).  They are coprime unless \(\ell=1\), and both are
separable.  Hence (3.3) is squarefree of degree five and has no rational
root.  Equations (2.3)--(2.4) therefore identify the complete fiber as

\[
 F^{-1}(y_\ell)
 \simeq
 \operatorname{Spec}\mathbb Q[S]/(E_{y_\ell}), \tag{3.5}
\]

and prove \(F^{-1}(y_\ell)(\mathbb Q)=\varnothing\).

## 4. Local solubility

Let \(\ell\equiv1\pmod {27}\) be prime.  We prove that (3.4) has a root in
every completion.

* Over \(\mathbb R\), the cubic \(X^3-\ell\) has a root.
* Over \(\mathbb Q_2\), \(X^3-\ell\) has a root by Hensel's lemma: modulo
  \(2\), \(X=1\) is a root and the derivative \(3X^2\) is a unit.
* Over \(\mathbb Q_3\), apply the strong form of Hensel's lemma at \(X=1\).
  Since \(v_3(1-\ell)\ge3>2v_3(3)\), the cubic has a root.
* Let \(p\ne2,3,\ell\).  If \(p\equiv1\pmod3\), then
  \(X^2+X+1\) has two simple roots modulo \(p\), hence a root in
  \(\mathbb Q_p\).  If \(p\equiv2\pmod3\), cubing is a bijection of
  \(\mathbb F_p^\times\), so \(X^3-\ell\) has a simple root modulo \(p\)
  and hence in \(\mathbb Q_p\).
* At \(p=\ell\), the congruence \(\ell\equiv1\pmod3\) makes
  \(X^2+X+1\) split with distinct roots modulo \(\ell\), so it has a root
  in \(\mathbb Q_\ell\).

This covers every place.  Dirichlet's theorem supplies infinitely many
primes \(\ell\equiv1\pmod {27}\), and the targets (3.2) are pairwise
distinct.  The theorem follows.

## 5. Height asymptotics

Use the standard multiplicative height on \(\mathbb A^3(\mathbb Q)\),
obtained from the embedding

\[
 (P,B,C)\longmapsto[1:P:B:C]\in\mathbb P^3.
\]

For \(\ell\equiv1\pmod {27}\), the target (3.2) has primitive integral
coordinates

\[
 [9:9:32\ell:24\ell+3].
\]

Consequently

\[
 H(y_\ell)=32\ell.
\]

If \(N_{\rm HP}(B)\) denotes the number of targets in this constructed
family with height at most \(B\), then

\[
\begin{aligned}
N_{\rm HP}(B)
 &=\pi(B/32;27,1)\\
 &\sim \frac1{\varphi(27)}\frac{B/32}{\log(B/32)}
 \sim\boxed{\frac{B}{576\log B}}. \tag{5.1}
\end{aligned}
\]

Thus the construction gives not only infinitude but an exact
prime-progression counting problem and an unconditional height asymptotic.
It is a lower bound for the total number of Hasse-failing targets of the
fixed map, not an asymptotic for that larger set.

## 6. An exceptional finite-field line

The same target line has deterministic local behavior at good primes.  In
the \(X\)-coordinate its inverse polynomial is

\[
 (X^3-a)(X^2+X+1).
\]

Let \(p>3\).  If \(p\equiv1\pmod3\), the quadratic factor splits into two
distinct linear factors over \(\mathbb F_p\).  If \(p\equiv2\pmod3\), cubing
is a bijection on \(\mathbb F_p\), so the cubic has exactly one root for
every \(a\); it is simple when \(a\ne0\).  Hence, away from the single
degenerate parameter when necessary, every point of this target line lies
in the finite-field image.

This does not contradict the ambient \(S_5\) fixed-point law.  The line is
contained in the reducibility locus, so its specialized monodromy is much
smaller than the generic symmetric monodromy.  It gives a concrete
arithmetic-statistical problem: compare the random-permutation distribution
on the full target with the deterministic local solubility on exceptional
thin subvarieties.

## 7. What this changes

The fixed-map problem does not require the genus-one common-discriminant
slice of the linear pencil.  That slice remains useful for understanding
Hasse failures inside the weighted family, and its Jacobian has positive
rank, but ramified decomposition groups still have to be controlled there.

The quadratic-tilt pencil absorbs the entire one-parameter intersective
family because translating \(X^2+X+1\) removes its linear term.  More
generally, any family

\[
 (X^3-a)Q(X)
\]

with a fixed quadratic \(Q\) can be tested for occurrence in one
quadratic-tilt Keller pencil after centering \(Q\).  This turns the search
for fixed-map Hasse families into a compatibility problem between:

1. an intersective polynomial family affine-linear in its parameter; and
2. the horizontal coordinate of a polynomializable marked-line Keller
   incidence.

The exact symbolic replay is

```bash
.venv/bin/python scripts/verify_infinite_hasse_keller_fibers.py
```

## 8. Arithmetic references

* D. Berend and Y. Bilu,
  [*Polynomials with roots modulo every integer*](https://doi.org/10.1090/S0002-9939-96-03210-8),
  Proc. Amer. Math. Soc. **124** (1996), 1663--1671.
* J. Sonn,
  [*Polynomials with roots in \(\mathbb Q_p\) for all \(p\)*](https://arxiv.org/abs/math/0612528),
  Proc. Amer. Math. Soc. **136** (2008), 1955--1960.

