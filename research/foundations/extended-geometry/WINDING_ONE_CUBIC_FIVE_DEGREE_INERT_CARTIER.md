# Inert-prime constraints for the degree-five cubic boundary

Use the sphere functional \(\Lambda\) and the regular winding-one phase
\[
 f=\frac{v}{x}R(x,vt),\qquad
 v=x-r,\qquad R(x,T)=U_0+U_1 T+U_2 T^2+UT^3,\qquad U=e(x-s)
\]
of least sphere degree seven. Here \(e,r,s\ne0\), \(r\ne s\), the field
has characteristic zero, and we extend it to name the roots if needed.
This form of \(f\) follows from the
[logarithmic boundary theorem](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md)
in the eventual-purity case with \(E=v^4U\).

## Proposition (Inert-prime fixed deficits)

Suppose that \(\Lambda(f^m)=0\) for all sufficiently large \(m\).
Take scalars \(a,b\) from the
[marker dependence](WINDING_ONE_CUBIC_FIVE_DEGREE_MARKER_DEPENDENCE.md)
and put \(A=ax+bx^2\). Thus \(\Lambda(v(vt-A)f^m)=0\) for every \(m\ge0\).
Define
\[
 \alpha=U_2/U,\qquad \beta=U_1/U,\qquad \gamma=U_0/U,\qquad
 H_{j,d}=[z^d](1+\alpha z+\beta z^2+\gamma z^3)^{j-2/3},
\]
with \(H_{j,-1}=0\), and let
\[
 \Psi_j=(vU/x)^j(H_{j,3j}-A H_{j,3j-1}),\qquad j\ge0.
\]
Write
\[
 \mathcal L_2(h)=h'-\frac{2h}{3}\left(\frac2x+\frac1v+\frac1{x-s}\right).
\]
For every \(j\ge0\) there are unique \(c_j,d_j\) and \(h_j\in K(x)\)
with
\[
 x\Psi_j=c_j x+d_j\,x/v+\mathcal L_2(h_j).
 \tag{1}
\]
Then \(d_j=0\). In particular \(x\Psi_1\), with
\[
 \Psi_1=\frac{vU}{3x}\left(q-(A+\alpha/3)p_0\right),\qquad
 p_0=\beta-\alpha^2/3,\qquad q=\gamma-\alpha\beta/3+2\alpha^3/27,
\]
has representative a scalar multiple of \(x\).

### Prime selection and the rational reduction

Over a number field fix \(j\) and use sufficiently large good primes
of residue characteristic \(p=3k+2\). Set \(m=k+j\), so that the
height of \(v(vt-A)f^m\) is at most \(p-1+3j\). Its only height
with interval denominator divisible by \(p\) is \(p-1\). This
height is even. Thus \(p\Lambda(v(vt-A)f^m)=0\) reduces to
\[
 [x^{k+j}]v^{p+k+j}
       \left([T^{p-2}]R^m-A[T^{p-1}]R^m\right)=0.
\]
Here \(v(vt-A)=v(T-A)\), and changing the coefficient variable
from \(t\) to \(T=vt\) adds \(p-1\) to the power of \(v\).
Replace \(v^p\) by \(x^p-r^p\). Since \(k+j<p\) and the other
factors are polynomial, the \(x^p\) term contributes zero. The
two deficits from the top of \(R^m\) are \(3j\) and \(3j-1\),
and \(m\equiv j-2/3\pmod p\). The finite binomial expansions
therefore give
\[
 0=[x^k](vU)^k\Psi_j=[x^{p-1}] F_0^k x\Psi_j,\qquad F_0=x^2vU.
 \tag{2}
\]
In particular, at \(j=0\) this gives \([x^k](vU)^k=0\).
The displayed expression for \(\Psi_1\) follows from
\(H_{1,2}=p_0/3\) and \(H_{1,3}=q/3-\alpha p_0/9\).

The poles of \(x\Psi_j\) lie among \(0,r,s,\infty\).
At \(0,r,s\), the respective residues of the multiplier in
\(\mathcal L_2\) are \(4/3,2/3,2/3\). As with the
[split-prime reduction](WINDING_ONE_CUBIC_FIVE_DEGREE_BOUNDARY_MOMENTS.md#reduction-by-rational-derivatives),
remove finite poles of order at least two by subtracting images of
the appropriate principal parts. The leading term at infinity of
\(\mathcal L_2(x^n)\) is \((n-8/3)x^{n-1}\), so remove polynomial
terms of positive degree. Then use affine \(h\) with prescribed
values at \(0,s\) to remove the two simple poles there. The
resulting representatives \(1,1/v\) are independent by the same
leading-term and affine-value tests. The identity
\[
 \mathcal L_2(x(x-s))=s/3-2x/3-2(r-s)x/(3v)
\]
replaces that basis by \(x,x/v\), whose change-of-basis determinant
is \(rs/2\). Also
\(\mathcal L_2(h)=0\) has only the zero rational solution, by the
pole and infinity tests. This proves (1).

Modulo \(p\) for \(k\equiv-2/3\), \(F_0^k\mathcal L_2(h)=(F_0^k h)'\).
Its coefficient at \(x^{p-1}\) is zero. The \(c_j\) term in (1)
is killed by the \(j=0\) instance of (2). Write \(\mu=r/s\).
After scaling \(x=rX\) and extracting nonzero constants, the
\(j=0\) coefficient and the \(d_j\) coefficient are proportional,
respectively, to
\[
 H_p(\mu)=(-1)^k\sum_{i=0}^k\binom ki^2\mu^i,\qquad
 \frac{\mu}{k}H_p'(\mu).
\]
These formulas use the binomial comparisons of the
[split-prime calculation](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_MARKERS.md#the-unmarked-and-marked-coefficients),
which hold also at \(p=3k+2\). The coefficient recurrence proves
the same second-order equation for \(H_p\). At \(\mu\ne0,1\),
that equation and \(k<p\) exclude a common zero of \(H_p,H_p'\).
Thus (2) gives \(d_j=0\) at all sufficiently large good primes.
This proves the identity over a number field. Over any
characteristic-zero field, a nonzero \(d_j\) would persist under
specialization of a finitely generated coefficient algebra to a
number field, together with the eventual-purity and marker identities.

## Theorem (Complex multiplication restriction)

For an eventually pure phase as above let \(Y\) be the smooth
projective curve with function field \(W^3=F_0\). It is bielliptic,
and the elliptic quotient described below has complex multiplication
by an order in \(\mathbb Q(\sqrt{-3})\).

### Ordinary and supersingular progressions

The four branch points \(0,r,s,\infty\) are fully ramified, so
Riemann–Hurwitz gives genus two. The independent holomorphic
forms are \(\omega=dx/W\) and \(\omega'=x\,dx/W^2\).
At \(0,r,s,\infty\), the respective orders of \(W\) are
\(2,1,1,-4\), and those of \(dx\) are \(2,2,2,-4\). These
give regularity of both forms.

The Möbius involution \(\tau(x)=s(x-r)/(x-s)\) swaps \(0,r\)
and \(s,\infty\). It lifts to
\[
 \iota(x,W)=(\tau(x),\delta\,x(x-r)/((x-s)W)),
      \qquad \delta^3=e^2s^3(s-r)^2.
\]
Indeed \(F_0(\tau(x))F_0(x)=(\delta x(x-r)/(x-s))^3\) and
\(x(x-r)/(x-s)=x\tau(x)/s\) is \(\tau\)-invariant. The lift
has order two and conjugates the order-three deck transformation
to its inverse. Consequently it exchanges the lines spanned by
\(\omega,\omega'\). The invariant part has dimension one,
giving an elliptic quotient, after adjoining a quotient point
if needed. Choices of \(\delta\) give conjugate involutions by
deck transformations. The geometric quotient therefore depends
only on \(r/s\), since scaling \(x=rX\) and rescaling \(W\) put
the branch polynomial in the form \(X^2(X-1)(X-1/\mu)\).

For \(p=3k+2\), write \(\omega'=(W^{-1})^pF_0^k x\,dx\).
The only exponent congruent to \(p-1\) modulo \(p\) in the
polynomial row lies at \(p-1\), as its exponents range from
\(2k+1\) to \(4k+1<2p-1\). Its coefficient is zero by (2) at
\(j=0\). Thus Cartier kills \(\omega'\) and, by \(\iota\),
also \(\omega\). The pullback of the elliptic differential is
injective and commutes with Cartier, so the elliptic quotient is
supersingular.

At \(p=3k+1\), \(\omega=(W^{-1})^pF_0^k\,dx\) has only
the possible Cartier coefficient at \(p-1\). It is a nonzero
multiple of \(H_p(\mu)\), which is nonzero because unmarked
purity forces \(H_p'(\mu)=0\) and the two have no common zero.
The same nonvanishing holds on \(\omega'\) by \(\iota\).
Cartier preserves the two deck eigenlines at split \(p\), so it
is invertible on holomorphic differentials and the elliptic
quotient is ordinary. These statements hold at every place
above all sufficiently large rational primes in the respective
classes. They can be checked after any finite extension of
the residue field.

### Identification of the CM field

Work first over a number field, and enlarge it to a finite
Galois extension \(L/\mathbb Q\) containing a field of definition
for the elliptic quotient \(E_1\) and \(\zeta_3\). Suppose \(E_1\)
has no complex multiplication. By Serre's open-image theorem its
\(\ell\)-adic image over \(L\) is Zariski dense in \(\mathrm{GL}_2\).
[Serre, introduction and §4.2](https://www.college-de-france.fr/media/jean-pierre-serre/UPL5874918517843398173_Serre_proprie_te_s_galoisiennes_des_courbes_elliptiques.pdf)

Induce \(V_\ell(E_1)\) from \(G_L\) to \(G_{\mathbb Q}\), with
representation \(\rho\). Let \(H\) be the Zariski closure of
\(\rho(G_L)\) over \(\overline{\mathbb Q}_\ell\). It acts blockwise,
and its projection onto every two-dimensional block is
\(\mathrm{GL}_2\) by the same density for each conjugate elliptic
curve. Its identity component is reductive, since its unipotent
radical projects to the trivial normal unipotent subgroup in every
block. The semisimple group \(S=(H^\circ)^{\mathrm{der}}\) projects
onto \(\mathrm{SL}_2\) in every block. Each simple Lie algebra ideal
of \(S\) projects either trivially or isomorphically to a given block's
\(\mathfrak{sl}_2\). Exactly one projects nontrivially for each
block, by surjectivity and commutation of distinct ideals. Thus
the connected simple normal factors have disjoint nonempty block
supports. In particular every element of \(H\) preserves each such
factor under conjugation, since conjugation in \(H\) is blockwise.

Choose \(\sigma\in\operatorname{Gal}(L/\mathbb Q)\) acting
nontrivially on \(\zeta_3\), of order \(d\), and write \(g\)
for the matrix of a lift to \(G_{\mathbb Q}\). Conjugation by
\(g\) permutes the simple factors of \(S\). Let \(S_0\) be the
one projecting onto the original \(E_1\) block and \(c\) its
orbit length. Then \(c\mid d\) because \(g^d\in H\). The
automorphism \(\operatorname{Ad}(g)^c\) of \(S_0\) is inner,
as \(S_0\) has type \(A_1\). Replace \(g\) by \(g'=g h\)
with \(h\in S_0\) chosen to cancel this action. Indeed the other
conjugates of \(h\) met on the first \(c-1\) steps lie in factors
commuting with \(S_0\). For \(t\) in a maximal torus of \(S_0\),
project \((g't)^d\) to the original block. All factors
\(\operatorname{Ad}(g')^i(t)\) outside \(S_0\) project trivially,
and each of the \(d/c\) factors in \(S_0\) equals \(t\).
Also \((g')^d\) centralizes \(S_0\), so it projects to a scalar.
It follows that the squared trace divided by determinant on this
block varies as \((z^{d/c}+z^{-d/c})^2\), where \(z,z^{-1}\)
are the two eigenvalues of the projected \(t\).

For a supersingular elliptic curve over a finite field of size
\(q\), the squared Frobenius trace divided by \(q\) belongs to
\(\{0,1,2,3,4\}\). Indeed the normalized Frobenius eigenvalues
are roots of unity, so the ratio is a rational algebraic integer
in \([0,4]\). It follows from the preceding torus variation that
the regular polynomial
\[
 \prod_{u=0}^4\left(\operatorname{tr}(B)^2-u\det B\right),
 \qquad B=(G^d)|_{V_\ell(E_1)},
\]
does not vanish identically for \(G\in gH\). By density it is
nonzero for some matrix of an actual lift of \(\sigma\), and by
continuity it is nonzero on a neighborhood specified by a finite
quotient of \(G_{\mathbb Q}\). Chebotarev gives rational
Frobenius elements in this neighborhood, choosing representatives
of their conjugacy classes. Their residue characteristic is
\(p\equiv2\pmod3\), and their \(d\)-th powers are Frobenius
elements of places of \(L\) with residue size \(p^d\). The
nonvanishing contradicts supersingularity at those places.
This proves complex multiplication, including when inert rational
primes have no degree-one places in \(L\).

By Deuring's reduction theorem a CM elliptic curve is ordinary
at sufficiently large rational primes split in its CM field and
supersingular at inert primes. The two progressions above identify
that field as \(\mathbb Q(\sqrt{-3})\), because distinct quadratic
fields give different splitting characters by Chebotarev in
their biquadratic compositum.
[Deuring, introduction and §4](https://download.uni-mainz.de/mathematik/Algebraische%20Geometrie/Lehre/WS23.Padische.1941.Deuring.pdf)

In a general characteristic-zero coefficient field the marker and
ratio analysis makes \(\mu\) algebraic. A number-field
specialization over \(\mathbb Q(\mu)\) retains the required
nonzero data and all pure moment identities. It retains the
geometric elliptic quotient because that quotient depends only on
\(\mu\). The restriction therefore holds over the original field
as well. It leaves CM orders in \(\mathbb Q(\sqrt{-3})\) with
their conductors to be determined.
