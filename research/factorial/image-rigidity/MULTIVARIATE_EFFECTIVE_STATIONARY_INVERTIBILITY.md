# Quantitative invertibility of stationary fibers on interior slabs

Use the notation of the
[effective stationary-fiber localization](MULTIVARIATE_EFFECTIVE_STATIONARY_FIBER_LOCALIZATION.md),
with \(z=2q-1\), \(W=\sum_i w_i\), and
\[
 X(z)=S(z)+\log(1-z)-\log(1+z),\qquad
 S(z)=\sum_i d_i\log(1+r_i z).
\]
Here \(X'(z)\) denotes the derivative in \(z\). Let
\[
 K=\{w:|\Re w|\le M,\ -\pi+\epsilon\le\Im w\le-\epsilon\},
 \qquad 0<\epsilon\le\pi/2.
\]

## Direct bounds when critical values stay away

Fix a positive curvature band and \(0<r\le\epsilon/2\).
Apply the effective localization to the rectangle with real bound
\(M+r\) and strip margin \(\epsilon-r\), giving \(R,h\).
If the disk \(|w-w_0|<r\) contains no critical value, with
\(w_0\in K\), each inverse branch over that disk is analytic and
bounded by \(R\). Indeed the proper stationary map is an unramified
finite covering over this disk, and the disk is simply connected.
Cauchy's estimate gives
\[
                  |X'(z)|\ge r/R
                  \qquad \bigl(X(z)=w_0\bigr).             \tag{1}
\]
The [one-sheet count](MULTIVARIATE_UNIFORM_SADDLE_MAP_COMPACTNESS.md)
and [shared supercritical geometry](MULTIVARIATE_SADDLE_MAP_TWO_SHEETS.md)
allow \(r=\epsilon/2\) throughout those families.
For the two-sheet subcritical family, the same estimate applies
whenever the disk avoids the unique fold value.

Two distinct preimages at such a field also have separation at least
\[
               \frac{2r}{R A_2},\qquad
               A_2=\frac{W_+ R^2/2+2}{h^2}.                \tag{2}
\]
Indeed the segment between them stays inside \(|z|\le R,\ \Im z\ge h\).
For real \(a\), \(|1+a z|\ge \Im z/|z|\), by minimizing over \(a\).
This bounds the row terms of \(X''\) by \(W_+ R^2/(2h^2)\).
The two base logarithms contribute at most \(2/h^2\).
Taylor's formula along the segment, using equal endpoint values
and (1), gives (2).

## A finite interval certificate for a shared slab

In the two-sheet case a single contour can test the whole slab for
critical values, including profiles with a fold. Choose
\(0<\eta\le\epsilon/2\). Use the effective localization for the
enlarged target rectangle \(K^+\) with real bound \(M+\eta\)
and strip margin \(\epsilon-\eta\), giving \(R,h\).
Take \(\Gamma\) to be the positive boundary of
\[
          |\Re z|\le R+1,\qquad h/2\le\Im z\le R+1 .
\]
Every boundary point is outside the localization for \(K^+\).
Consequently \(|X(z)-w|\ge\eta\) for every \(z\in\Gamma,w\in K\),
uniformly in the profiles. Define
\[
 S_k(w)=\frac1{2\pi i}\int_\Gamma
               z^k\,\frac{X'(z)}{X(z)-w}\,dz,\qquad
 {\cal D}(w)=2S_2(w)-S_1(w)^2.
\]
The residue theorem gives \(S_k=z_1^k+z_2^k\) with multiplicity.
Thus \({\cal D}=(z_1-z_2)^2\) is analytic in the field and vanishes
precisely at a critical value. This formula also applies to the
retained and diffuse closure of the profile family.

For explicit interval tests put \(B=2(R+1)\),
\(\Lambda=\operatorname{length}\Gamma\), and
\(A=(W_+B+4)/h\). This bounds \(|X'|\) on the contour by
\(|1+r_i z|\ge\Im z/|z|\) and the base logarithmic derivatives.
If \(Y\) approximates a profile map \(X\) on \(\Gamma\) within
\(\delta_0<\eta/2\) in value and \(\delta_1\) in derivative, its
contour denominator is at least \(\eta/2\).
The trace integral changes by at most
\[
  \frac{\Lambda B^k}{2\pi}
      \left( \frac{2\delta_1}{\eta}
                    +\frac{2A\delta_0}{\eta^2}\right).
                                                               \tag{3}
\]
This follows by writing the difference of the two trace integrands
with the derivative error in the first term and the denominator
error in the second term. Each trace also has field derivative bounded by
\[
                M_k=\frac{\Lambda B^k A}{2\pi\eta^2}.
                                                               \tag{4}
\]
For a target cell one may therefore integrate only at its center,
then add the bound from (4) times its radius.
The integral on each segment of \(\Gamma\) has a continuous
analytic integrand in a neighborhood of the segment, so complex
interval quadrature or Riemann sums with derivative enclosures
can certify it.

The explicit finite profile nets in the
[localization note](MULTIVARIATE_EFFECTIVE_STATIONARY_FIBER_LOCALIZATION.md#finite-analytic-approximation-of-the-profile-family)
bound \(\delta_0,\delta_1\) in (3). In particular a cutoff on
\(|r_i|\), a mesh for the remaining diffuse mass, and a target mesh
give entirely finite interval tests. A certified \(d_0>0\) lower
bound for \(|{\cal D}|\) on such cells proves
\[
 |X'(z)|\ge
 \left(\frac{M_1}{2}+\frac{2M_2+4R M_1}{4\sqrt{d_0}}\right)^{-1}.
                                                               \tag{5}
\]
Indeed \(z=(S_1\pm\sqrt{\cal D})/2\) locally and
\(|S_1|\le2R\). If the closed profile and field set excludes
critical values, compactness and continuity give a positive
lower bound for \(|{\cal D}|\), so arbitrarily fine finite meshes
on that set eventually certify (5).
For the full curvature band one can enumerate retained lists and
mesh each diffuse interval within the stated curvature bounds.
Those lists with allowed diffuse mass are themselves in the
retained and diffuse closure, keeping the certificate on the
same compact family.

An unrestricted slab can contain genuine folds. For example the
diffuse profile of mass one has \(S(z)=z/2\) and a critical point
at \(z=i\sqrt3\) with target \(i(\sqrt3/2-2\pi/3)\).
The interval test requires a critical-value-free slab or a
specified subregion away from the folds. The explicit radius
in the localization note can be large as the strip margin shrinks,
so these estimates give finite certification bounds without a
practical complexity claim.

## The fold and the interior derivatives at vanishing curvature

Suppose \(0<W\le1/4096\) and put \(u=\sqrt W\).
Every row is shared, and \(|r_i|\le W/2\) by nonzero integer gaps.
Also \(\max\rho(q)\le W/[4(1-W/2)]<1\) for the real stationary
curvature function in the two-sheet classification.
Its unique upper critical point \(z_c\) has field \(w_c=X(z_c)\).
Then
\[
 \left|z_c-\frac{2i}{\sqrt W}\right|<4,\qquad
 \bigl|w_c+i\pi-2i\sqrt W\bigr|<8W.                        \tag{6}
\]

For the proof put \(v=u z/2\). The critical-point equation
\((1-z^2)S'(z)-2=0\) becomes
\[
 -2(v^2+1)+u^2/2+(u^2/2-2v^2)E(v)=0,\qquad
 E(v)=\frac2W S'(2v/u)-1 .
\]
On and inside \(|v-i|=2u\), one has \(|v|\le33/32\) and
\(\max_i|r_i z|\le u|v|<1/60\). Since
\(S'(z)=\tfrac12\sum_i w_i/(1+r_i z)\), it follows that
\(|E(v)|<2u\). The perturbation of \(-2(v^2+1)\) is
less than \(7u\) there. On the circle the main term has
modulus at least \(4u(2-2u)>7u\), so Rouché gives one
critical point with \(|v-i|<2u\). The classification identifies
it as \(z_c\), and the first bound in (6) follows.

At this point
\(|S(z)-Wz/2|\le W^2|z|^2/4<2u^2\),
using \(|\log(1+s)-s|\le |s|^2\) for \(|s|\le1/2\).
Also \(|Wz/2-iu|<2u^2\). On the upper half-plane where
\(|z|>1\), the base term equals
\[
  -i\pi+\log(1-1/z)-\log(1+1/z).
\]
This follows by evaluating on the positive imaginary ray and
then by analytic continuation. Here \(|v|\ge31/32\)
and hence \(|1/z|<u\). The two small logarithms differ from
\(-2/z\) by less than \(3u^3<u^2\), by the odd power series.
Moreover \(|-2/z-iu|<3u^2\) from \(|v-i|<2u\).
Adding these errors gives the second bound in (6).

In particular, if \(u\le\epsilon/8\), all the disks of radius
\(\epsilon/2\) about \(K\) avoid the fold, since its height
above \(-\pi\) is at most \(2u+8u^2<\epsilon/2\).
The inverse bound can also be made polynomial in \(W\) in this
regime. For any interior rectangle with real bound \(M'\) and
margin \(\epsilon'\), a preimage with \(|z|\ge8/\epsilon'\)
satisfies
\[
  S(z)=X(z)+i\pi-\log(1-1/z)+\log(1+1/z),
  \qquad
  |\log(1-1/z)-\log(1+1/z)|\le4/|z|.
\]
Thus \(S\) at this preimage has imaginary part between
\(\epsilon'/2\) and \(\pi-\epsilon'/2\), and real part of modulus
at most \(M'+1\). In the coordinate \(v=Wz\), the same sum
\(S\) has slopes \(r_i/W\) and masses \(w_i/W\).
The slopes are bounded by \(1/2\), the total mass is one, and
the gaps remain integers. The radius proof in the effective
localization note used only these slope, mass, imaginary-part,
and real-part bounds for \(S\), after removing the base
logarithms. It supplies an explicit
\(C=C(M',\epsilon')\) bounding \(|v|\).
Together with the preimages of modulus less than \(8/\epsilon'\)
this puts all fibers in \(|z|\le\max(C,8/\epsilon')/W\).
Apply (1) with \(M'=M+\epsilon/2\), \(\epsilon'=\epsilon/2\)
and \(r=\epsilon/2\). Uniformly in the small-curvature slab,
\[
    |X'(z)|\ge
    \frac{\epsilon\,W}{2\max(C(M+\epsilon/2,\epsilon/2),16/\epsilon)}.
                                                               \tag{7}
\]
In particular the explicit bounds control stationary invertibility
up to vanishing curvature on each fixed interior slab. They do
not yet quantify root separation at the critical field or the
activation of the actual moment contours.
