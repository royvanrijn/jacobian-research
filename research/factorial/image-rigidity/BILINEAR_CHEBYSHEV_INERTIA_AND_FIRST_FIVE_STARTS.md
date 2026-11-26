# Compressed Chebyshev inertia and the first five pure-cross three-windows

Let \(\Lambda(X^iY^j)=i!j!\). This note proves that, over every
characteristic-zero field, a phase
\[
                         f=BXY+dX+eY,\qquad Bde\ne0,
\]
cannot have \(\Lambda(f^4)=\Lambda(f^5)=\Lambda(f^6)=0\).
It does not prove the unrestricted later-start pure-cross theorem.

The mathematical input is the native resonance reduction, Bessel
real-double exclusion and compressed-matrix identity in
[Bilinear real-conjugate ratio three-window](BILINEAR_REAL_CONJUGATE_RATIO_THREE_WINDOW.md),
Parts A--B. Their exact scopes are recalled below. The matrix-inertia
argument and the one finite scalar certificate used here are proved
locally.

## 1. Exact matrix inertia

For \(N\ge2\), let \(x_1<\cdots<x_N\) be the roots of \(T_N\), choose
\(c=x_h\), and omit index \(h\). Put \(n=N-1\) and
\[
 G_n(c,b)=\sum_{j=0}^n U_n^{(j)}(c)b^j,\qquad
 K_{ii}=\frac{x_i}{2(1-x_i^2)}+\frac1{x_i-c},\quad
 K_{ij}=\frac1{x_i-x_j}\quad(i\ne j).
\]
The imported exact interpolation/cofactor identity is
\[
 \frac{G_n(c,b)}{U_n(c)}=\det(I-bK).
\]
For the real diagonal \(J=\operatorname{diag}(x_i-c)\), direct entry
calculation gives
\[
 JK+K^tJ=H,\qquad
 H={\bf1}{\bf1}^t+
       \operatorname{diag}\left(\frac{1-cx_i}{1-x_i^2}\right)>0. \tag{1}
\]
This excludes imaginary-axis eigenvalues: an eigenvector \(v\ne0\)
with eigenvalue \(\lambda\) would give
\(v^*Hv=2\Re\lambda\,v^*Jv\).

**Lemma 1.** The algebraic multiplicities of the eigenvalues of \(K\)
in the right and left half-planes are respectively \(N-h\) and \(h-1\).

Indeed let \(E_+\) and \(E_-\) be the sums of generalized eigenspaces
with respectively positive and negative real parts. There are no
imaginary-axis eigenvalues, so these spaces have total dimension
\(N-1\). For nonzero \(v\in E_+\), set \(w(t)=e^{-tK}v\).
Jordan normal form shows \(w(t)\to0\) exponentially up to a fixed
polynomial factor. From (1),
\[
 \frac{d}{dt}(w(t)^*Jw(t))=-w(t)^*Hw(t),\qquad
                  v^*Jv=\int_0^\infty w(t)^*Hw(t)\,dt>0.
\]
Thus \(J\) is positive definite on \(E_+\). Similarly, for
\(v\in E_-\), use \(w(t)=e^{tK}v\) to obtain
\[
                         -v^*Jv=\int_0^\infty w(t)^*Hw(t)\,dt>0.
\]
The positive and negative indices of \(J\) are \(N-h\) and \(h-1\).
A positive subspace has dimension at most its positive index, since
projection onto the positive coordinate space is injective; the
negative case is identical. Both dimension inequalities must be
equalities because their sums equal \(N-1\). This proves the lemma,
including algebraic multiplicity through generalized eigenspaces.

Inversion \(\lambda\mapsto b=1/\lambda\) preserves the sign of the real
part. It also preserves nonzero root multiplicity in the determinant
identity. Consequently \(G_n(c,\cdot)\) has exactly \(N-h\) roots in
the right half-plane and \(h-1\) in the left, counted algebraically.

A repeated nonreal root brings its distinct conjugate repeated root
in the same half-plane, using at least four algebraic dimensions
there. Hence, if
\[
                              N-h\le3,\qquad h-1\le3,          \tag{2}
\]
there are no repeated nonreal roots. The imported Bessel argument
excludes every repeated nonzero real root at a \(T_N\) resonance.
Also \(G_n(c,0)=U_n(c)\ne0\). Thus (2) implies squarefreeness of the
entire polynomial \(G_n(c,\cdot)\). Accretivity alone was not used to
assert simplicity; the dimension count and real-double exclusion
are both necessary.

## 2. The primitive quartic factor for \(T_6\)

A length-three zero block starting at power four has middle index
\(n=5\). The exact native reduction requires
\[
                         T_6(c)=0,\qquad G_5(c,b)=G_{5,b}(c,b)=0.
\]
Factor
\[
 T_6(c)=(2c^2-1)(16c^4-16c^2+1).                              \tag{3}
\]
The two central \(T_6\) nodes are \(c=\pm\cos(5\pi/12)\); their
ordered indices are \(h=3,4\). Both satisfy (2), since their
half-plane counts are \(3,2\) or \(2,3\). Thus \(G_5\) is squarefree
at either central node.

The quartic \(16c^4-16c^2+1\) is irreducible over \(\mathbb Q\).
To see this let \(u=2c\); its equation is \(u^4-4u^2+1=0\).
Inside \(\mathbb Q(u)\),
\[
                       (u^2-2)^2=3,\qquad (u^3-3u)^2=2.
\]
Hence \(\mathbb Q(u)\) contains square roots of both 3 and 2.
These generate a degree-four extension: a hypothetical
\(\sqrt2=a+b\sqrt3\), \(a,b\in\mathbb Q\), forces \(ab=0\) and then
a rational square equal to 2 or \(2/3\), both impossible.
It follows that \(u\) has degree four, proving irreducibility.

Every root of this quartic is therefore conjugate to a central node.
A repeated root \(b\) of \(G_5(c,\cdot)\) is algebraic over
\(\mathbb Q(c)\). An embedding sending \(c\) to the central node
extends to its algebraic coefficient field and preserves both
\(G_5=0\) and \(G_{5,b}=0\). This contradicts squarefreeness at that
node. Thus the quartic factor in (3) is completely excluded, not only
the two central real nodes.

## 3. The remaining quadratic factor: one literal finite certificate

Suppose \(c^2=1/2\), and put \(b=cz\). Direct differentiation of
\(U_5(c)=32c^5-32c^3+6c\) gives
\[
 G_5(c,b)=-2c-2b+128cb^2+768b^3+3840cb^4+3840b^5
                         =2c F(z),
\]
\[
                  F(z)=480z^5+480z^4+192z^3+32z^2-z-1.        \tag{4}
\]
Over \(\mathbb F_7[z]\), the literal identity
\[
 (6+5z+5z^2+6z^3)F(z)
 +(3z+3z^2+z^3+3z^4)F'(z)=1                                 \tag{5}
\]
is checked by ordinary multiplication and reduction of the integer
coefficients. Both degrees are preserved modulo 7: their leading
coefficients are \(480\not\equiv0\) and \(2400\not\equiv0\).
Thus their integer Sylvester resultant is nonzero modulo 7, and in
particular is a nonzero integer. Therefore \(F,F'\) are coprime over
\(\mathbb Q\) and \(F\) is squarefree over every characteristic-zero
field. Equation (4), with \(c\ne0\), proves the needed squarefreeness
of \(G_5(c,\cdot)\).

This is an exact fixed polynomial certificate in a uniform
coefficient argument, not a sample of moment values and not an
extrapolation to other indices.

## 4. Characteristic-zero phase conclusion and limits

For coefficients in an arbitrary characteristic-zero field, extend
to an algebraic closure and choose \(s^2=de\). The native reduction
uses \(c=(d+e)/(2s)\), \(b=B/(2s)\). A hypothetical block at powers
\(4,5,6\) would make \(c\) a root of (3), hence algebraic over
\(\mathbb Q\), and \(b\) a root of the nonzero \(G_5(c,\cdot)\), hence
algebraic as well. Embed this finite algebraic field in \(\mathbb C\).
Sections 2 and 3 exclude every possible resonant \(c\), giving a
contradiction.

Together with the proof for starts 1,2,3 in
[Bilinear real-conjugate ratio three-window](BILINEAR_REAL_CONJUGATE_RATIO_THREE_WINDOW.md),
Section C.2, this establishes the first four positive starting
positions for every pure-cross three-term phase.

## 5. The first five starts and the remaining terminal orders

**Corollary.** Over every characteristic-zero field, for
\(f=BXY+dX+eY\) with \(Bde\ne0\), and every integer \(1\le q\le5\),
at least one of
\[
                  \Lambda(f^q),\quad\Lambda(f^{q+1}),
                  \quad\Lambda(f^{q+2})
\]
is nonzero.

Starts 1,2,3 are the exact result just cited. Start 4 is proved in
Sections 1--4. Start 5 is the \(p=7\) instance of the
[prime-terminal theorem](BILINEAR_PRIME_TERMINAL_THREE_WINDOWS.md),
Section 6. This deduction combines proved windows; it does not
extrapolate from a finite scan.

More generally, the prime-terminal theorem excludes every start
\(q=p-2\) with \(p\) an odd prime. Thus any remaining counterexample
must start at \(q\ge6\), and its terminal \(q+2\) cannot be an odd
prime. The previously proved restriction that its squared cross
ratio is algebraic with no real conjugate is unchanged. No
every-start conclusion for unrestricted ratios, and no full binary
quadratic SFC, SFC(2), or SFC(3) closure is claimed.
