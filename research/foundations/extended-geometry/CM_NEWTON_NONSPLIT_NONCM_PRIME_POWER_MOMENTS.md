# Prime-power moments in connected non-CM covers of unbounded order

An order-\(N\) cyclic cover can have a nonzero Frobenius return
on the coherent line used by a sphere moment, even when a different
primitive character orbit has zero stable rank at the same prime.
In the family below the two conclusions come from the even and odd
characters of a refinement to order \(2N\). A semistable degeneration
at a second prime forces the primitive order-\(N\) factor to contain
a non-CM subvariety of growing dimension.
At \(p=5,c=3\) the exact primitive \(5\)-rank comes from a digit
condition and is asymptotic to half the dimension.

## The covers and the moment

Fix a rational prime \(p>3\), with \(p\equiv2\pmod3\), and an odd
integer \(r\ge3\). Put
\[
                 m=p-1,\qquad N=\frac{p^r-1}{p-1}.
\]
Then \(N\) is odd, \(3\nmid N\), and \(p\nmid N\).
Take \(c\in\mathbb Q\) with \(v_3(c)=1\) and with \(c,1-c\)
units at \(p\). In particular \(c=3\) works for every indicated \(p\).
Write \(C\) and \(D\) for the smooth projective normalizations of
\[
 C:\quad W^N=z^{N-2}(z^2-1)(z^2-c),\qquad
 D:\quad V^N=u^{N-2}(u-1)^2(u-c)^2.
\]
For either cover, write \(J^{[N]}\) for the Jacobian factor cut out by
the primitive order-\(N\) idempotent of the indicated deck group.

**Theorem.** The connected cover \(C\) has a primitive factor
of dimension \(2\varphi(N)\), with a factor \(J^{[N]}(D)\) of
dimension \(\varphi(N)\) that is totally toric at \(3\).
Thus \(J^{[N]}(C)\) is not potentially CM. At \(p\), it is
nonordinary, and the \(r\)-fold coherent Frobenius is invertible
on its character-one line.

On the sphere \(2xy+t^2=1\), use the literal Laurent chart in
\(x,t\), and let \(z=x/t\). Every phase \(f\) with highest literal
weight \(j+i=N\) and top angular polynomial
\[
                      z^{-2}(z^2-1)(z^2-c)
\]
has least sphere degree \(N\) and is an absolute nonpower in the
sphere function field. Suppose the lower coefficients belong to a
number field and are integral at a place above \(p\). For every
integer \(k\ge1\), put \(q=p^{kr}\) and \(m_k=(q-1)/N\). Then
\(q\mathcal A(f^{m_k})\) is integral at that place and reduces to one.
In particular \(\mathcal A(f^{p-1})\ne0\), and each fixed phase in
the theorem has infinitely many nonzero moments certified at this
same prime.

## The refinement and the toric factor

The cover \(C\) is connected because the exponent at each root of
\(z^2-1\) is one. Its function field has the equivalent presentation
\[
                     W^{2N}=u^{N-2}(u-1)^2(u-c)^2,\qquad u=z^2 .
\]
Indeed \(z^{N-2}=W^N/((u-1)(u-c))\) recovers \(z\) together
with \(z^2=u\), since \(N\) is odd. The refined cyclic cover of
degree \(2N\) is connected since \(\gcd(N-2,2N)=1\).
The original generator is the square of a refined generator.
With \(V=W^2\), the map to \(D\) is the quotient by the
order-two subgroup of the refined deck group. Thus its cohomology
pulls back to the even characters of the refined group. This is
an isogeny factor of the Jacobian, since the norm composed with
pullback is multiplication by two.

To study \(D\) at \(3\), blow up the collision of \(u=0\) and
\(u=c\), which separates the two marked sections because \(v_3(c)=1\).
Near the resulting node use \(us=c\). The equation of the cover
there becomes
\[
                       (V/u)^N=(u-1)^2(1-s)^2 .
\]
The right side is a unit and \(3\nmid N\), so the normalization
is étale over the node, with \(N\) geometric points there.
On the main component the reduced equation is
\((V/u)^N=(u-1)^2\). On the exceptional component, putting
\(U=u/c\), the reduced equation is
\((V/c)^N=U^{N-2}(U-1)^2\).
Both components are connected and genus zero. Each is a cyclic
cover with two fully ramified points, since \(\gcd(N,2)=1\).
Away from the node, the distinct tame branch sections have smooth
normalization. This gives semistable reduction, after an unramified
extension if a splitting field for the node points is needed.

The dual graph has two vertices joined by \(N\) edges, on which
the deck group acts regularly. Its graph homology, the torus
character lattice of the Jacobian, is the augmentation representation.
The primitive part has rank \(\varphi(N)\), and the entire
Jacobian of \(D\) has totally toric reduction. Toric rank is
preserved under finite extension within semistable reduction and
under isogeny. A positive-dimensional CM abelian variety acquires
good reduction after finite extension, so the primitive factor of
\(D\) cannot be potentially CM. Its pullback also excludes
potential CM for the larger primitive factor of \(C\).
[SGA 7 I, Exposé IX, §12.3.7 and §12.4](https://grothendiecksga.com/read/sga7/en/9-12.html)
· [Serre–Tate, Theorem 6](https://annals.math.princeton.edu/1968/88-3/p05)

## The failing character orbit

The refined branch exponents at \(0,1,c,\infty\) are
\(N-2,2,2,N-2\). For an original primitive character
\(i\in(\mathbb Z/N)^\times\), put \(x=i/N\) using the representative
of \(i\) in \([1,N-1]\). The two refined labels are determined by
their parity \(e\in\{0,1\}\). Their coherent dimensions are
\[
      h_x(e)=-1+2\{e/2-x\}+2x,\qquad
      h_x(0)=1,\qquad h_x(1)=2\mathbf 1_{x>1/2}.
\]
For example, an eigenfunction of refined exponent
\(\ell\in[1,2N-1]\) over the affine line is generated by
\(W^\ell/\prod_{a=0,1,c}(u-a)^{\lfloor \ell e_a/(2N)\rfloor}\),
where \(e_a\) are the finite branch exponents. At infinity
the transition exponent is the ceiling of the sum of the
finite fractional parts. The cohomology of the resulting
line bundle on \(\mathbb P^1\) has dimension that ceiling
minus one, or the sum of all four fractional parts minus one.
The fractions at \(0,\infty\) are \(\{e/2-x\}\), and those
at \(1,c\) are \(x\), which gives the display. The even
primitive part has dimension \(\varphi(N)\).
For the full original primitive part, the dimensions
of \(i\) and \(-i\) sum to four, giving \(2\varphi(N)\).
The character calculation applies at \(p\), where all four
branch points of the refined cover remain distinct and
the cover is tame.

Absolute Frobenius sends \((i,e)\) to \((pi\bmod N,e)\), by
\(p\)-semilinearity on the deck eigenvalues. For \(e=1\), start
at the unit \(i=2\), where \(h_{2/N}(1)=0\). Write \(X=p^{r-1}\).
The representative of \(p^{r-1}i\bmod N\) is \(2X-N\), and
\[
                \frac{2X-N}{N}=\frac{(p-2)X+1}{pX-1}
                        \in(1/2,1).
\]
That component has dimension two. Every vector there eventually
maps through the zero space on its orbit under Frobenius, so the
stable rank on this orbit is zero. On an ordinary abelian variety
absolute Frobenius on coherent first cohomology is bijective.
The primitive factor in the theorem is therefore nonordinary.
[Cyclic Hasse–Witt decomposition and formula, Lin–Mantovan–Singhal, §2.4](https://arxiv.org/html/2406.09632v1#S2.SS4)

## The surviving line and the isolated moment

The original character-one line of \(C\) has Čech generator \(W/z\).
Indeed all finite branch exponents of the original cover are less
than \(N\), and their sum \(N+2\) gives transition exponent two at
infinity for this character. At \(p\), the \(r\)-fold Frobenius
preserves the line because \(p^r\equiv1\bmod N\). If
\(G(z)=z^{N-2}(z^2-1)(z^2-c)\), it sends the representative to
\(W G(z)^m/z^{p^r}\). The coefficient of \(W/z\) is consequently
\[
 \begin{aligned}
 [z^{p^r-1}]G(z)^m
   &=[u^m](u-1)^m(u-c)^m\\
   &=\sum_{a=0}^{p-1} c^a
     =(c-1)^{p-1}=1\qquad\bmod p .
 \end{aligned}
\]
The second equality in the residue field uses
\((u-b)^{p-1}=(u^p-b^p)/(u-b)=\sum_{a=0}^{p-1}b^{p-1-a}u^a\).
The last equality uses \(c\in\mathbb Q\) and \(c-1\) a unit.
In particular the selected line is invertible under the return.

The same return can be computed one step at a time on \(D\).
At step \(a\) for \(0\le a<r\), the character is
\((N+p^a)/2\), with finite branch residues
\((N-p^a,p^a,p^a)\) modulo \(N\). Raising the normalized
affine eigensheaf generator to the \(p\)-th power gives the
coefficient at \(u^{p-1}\) of a polynomial with branch
exponents \(\lfloor p(N-p^a)/N\rfloor\) at \(0\) and
\(\lfloor p^{a+1}/N\rfloor\) at each of \(1,c\). The exponents are
\((p-1,0,0)\) through step \(r-2\), and \((0,p-1,p-1)\)
at step \(r-1\). Thus each individual step is invertible.

A regular lift of the prescribed face is
\[
                4c\,y^2t^{N-2}-(1+c)t^N+x^2t^{N-2}.
\]
The relation \(2xy+t^2=1\) makes the maximal \(j+i\) in each
nonzero angle equal to the degree of its highest normal-basis term.
This proves least degree \(N\) for every compatible lower face.
In the function field with coordinates \(z,t\), a power exponent
would divide the pole order \(N\) at \(t=\infty\), and it would
also divide the pole order two at \(z=0\) of the leading
coefficient there. Since \(N\) is odd, the phase is an absolute
nonpower.

At moment order \(m=p-1\), the angular constant term has height
at most \(Nm=p^r-1\). Integration of a monomial of even height
\(j\) divides by \(j+1\); odd heights integrate to zero. The only
denominator in this range divisible by \(p^r\) is the top one,
and its coefficient comes from products of face terms alone.
Thus multiplication by \(p^r\) and reduction leaves the same
coefficient as the Frobenius return. Since the \(r\)-fold return
fixes the basis class \(W/z\), so does every \(k\)-fold iterate of
it. At order \(m_k\), the highest height is \(q-1\), isolating the
coefficient of the \(kr\)-fold return in the same way. This proves
all the moment congruences. In particular the connected conductors
and the least sphere degrees both grow for a fixed certifying
prime, without a density assertion about the primes of any one
fixed cover.
[Prime-power interval isolation](CM_NEWTON_INFINITE_MOMENT_FAMILIES.md#isolating-the-coefficient-with-a-prime-power)

## Exact primitive rank at five

**Proposition.** In the theorem take \(p=5,c=3\). Let \(S_r\) be the
units \(1\le i<N\) whose \(r\)-digit base-five block for \(4i\)
has every digit at least two, counting leading zeros. The \(5\)-rank
is \(\varphi(N)+2\lvert S_r\rvert\) in dimension \(2\varphi(N)\).
The ratio tends to one half along odd \(r\to\infty\).

Use the finite eigensheaf generators of the refined \(2N\)-cover
and their coherent basis classes \(w_\ell/u^j\). Set
\(x=i/N\), \(a=\lfloor5x\rfloor\), \(x'=5x-a\), and
\(Q=(u-1)(u-3)=u^2+u+3\) over \(\mathbb F_5\). The Frobenius step
multiplies the target generator by \(u^gQ^a\), with
\(g=\lfloor5\{e/2-x\}\rfloor\) for parity \(e\). The three exponents
in that multiplier are the floors of five times the finite fractional
parts at the source. Thus the entry from \(j\) to \(j'\) is
\([u^{5j-j'-g}]Q^a\). The \(p\)-rank is the dimension of the stable
image of these maps, summed over their orbits.
[Bouw, §2 and Lemma 5.1](https://doi.org/10.1023/A:1017513122376).

For even parity \(g=4-a\), and all line maps are units given by
\([u^a]Q^a=1,1,2,4,1\) for \(a=0,1,2,3,4\), respectively. The
even rank is \(\varphi(N)\). An odd orbit meeting \(x<1/2\) has
zero stable rank. When \(x,x'>1/2\), both spaces have dimension two
and \(g=6-a\). For \(a=2,3,4\) their Frobenius matrices are
\[
  \begin{pmatrix}4&0\\0&1\end{pmatrix},\qquad
  \begin{pmatrix}2&1\\2&3\end{pmatrix},\qquad
  \begin{pmatrix}2&4\\3&3\end{pmatrix},
\]
all with determinant four modulo five. An orbit with all \(x>1/2\)
therefore contributes two per label. Since \(x=4i/(5^r-1)\),
multiplication by five modulo one rotates its repeating base-five
block. Every rotation exceeds \(1/2=0.\overline{2}_5\) precisely
when all digits are at least two, since the all-two block is excluded
by odd \(N\). Lastly \(\lvert S_r\rvert\le3^r\) and
\(\varphi(N)\ge N^{3/4}\), using \(1-1/\ell\ge\ell^{-1/4}\)
for every prime \(\ell\mid N\), all of which are at least five.
This makes the odd contribution \(o(\varphi(N))\).

For \(N=31\), the digit blocks come from the three permutations of
each of \((2,2,4)\) and \((2,3,3)\), giving \(5\)-rank \(42\) in
dimension \(60\). The fourth moment remains nonzero with normalized
residue one for the degree-31 phase \(t^{29}(x^2+12y^2-4t^2)\)
plus any strict lower-degree perturbation integral at \(5\).

The [exact arithmetic replay](verify_cm_newton_nonsplit_noncm_prime_power_moments.py)
computes the individual Frobenius coefficients on \(D\), the dimensions
on the indicated odd orbit, and the full rational sphere moment with
strict lower terms added. It gives normalized residue one at
\((p,N)=(5,31),(5,781),(11,133),(17,307)\).
It checks the \(5\)-rank \(42\) at \(N=31\) by coherent orbits and
the digit count, and independently expands the top coefficient at
order \(504=(5^6-1)/31\), giving residue one.
