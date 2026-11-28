# Compact Gaussian score contours and branches at infinity

The [exterior algebraic-arc obstruction](THREE_CELL_CIRCULAR_INTERFACE_OBSTRUCTION.md#algebraic-carriers-and-the-conic-consequence)
uses a branch of the complexified score equation on which one
isotropic coordinate has a pole and the other stays bounded. Compact
indicator sources with analytic closed score contours need not have
such a branch. The following example has the same two positive source
weights as an exterior contrast of a three-cell partition. A separate
bound below controls balls contained in actual bounded winner regions.

## An indicator source with no meromorphic branch at infinity

Fix \(0<\rho<1\), \(a>R>0\), and let
\[
 f=2\,1_{B(ae_1,R)}+1_{B(-ae_1,R)},\qquad \sigma^2=1-\rho^2.
\]
For a complex observation \((z_1,z_2)\) use
\(u=z_1+iz_2\) and \(v=z_1-iz_2\). With
\(\eta=\rho/(2\sigma^2)\), put \(\zeta(y)=\eta(y_1-iy_2)\) and
\[
 L(u,v)=\int_{\mathbb R^2} f(y)e^{u\zeta(y)+v\overline{\zeta(y)}}
                         e^{-|y|^2/(2\sigma^2)}\,dy .
\]
The complexified score equation \(T_\rho f=c>0\) is
\(L(u,v)=2\pi\sigma^2 c\,e^{\alpha uv}\) for
\(\alpha=\rho^2/(2\sigma^2)\).

Put \(b=\eta a\), \(\ell=a/\rho\), and \(c_0=e^{-a^2/(2\sigma^2)}\).
Changing variables by the disk centers and integrating in angle gives
\[
\begin{split}
 L(u,v)=c_0\big[&2e^{b(u+v)}H((u-\ell)(v-\ell))\\
                &+e^{-b(u+v)}H((u+\ell)(v+\ell))\big],\\
 H(w)&=2\pi\int_0^R s e^{-s^2/(2\sigma^2)}
            \sum_{n\ge0}\frac{(\eta^2s^2w)^n}{(n!)^2}\,ds .
\end{split}                                                   \tag{CB1}
\]
For the first disk the Gaussian weight supplies the shifts
\(u\mapsto u-\ell\) and \(v\mapsto v-\ell\), with the opposite shifts
for the second. The angular integration keeps exactly the terms with
equal powers of the two plane coordinates in the exponential series.
This gives the series in (CB1). Its absolute sum is at most
\(e^{2\eta R\sqrt{|w|}}\), by comparing its diagonal terms with the
square of the series for \(e^{\eta R\sqrt{|w|}}\). In particular \(H\)
is entire, nonzero and satisfies \(|H(w)|\le C e^{C\sqrt{|w|}}\).

Suppose that a germ satisfying the score equation has both isotropic
coordinates meromorphic at its endpoint, with at least one pole.
They cannot both have poles. If the pole orders are \(d,e>0\), the
leading term of \(\alpha uv\) has order \(d+e\) and is positive real
to leading order along a suitable parameter ray. The bound
\(|L(u,v)|\le C\exp(C(|u|+|v|))\) from compact support contradicts the
equation there. If exactly one coordinate has a pole, interchange
\(u,v\) if needed and choose a parameter \(z\) near infinity with
\(u=z^m\), \(m\ge1\). This absorbs the analytic root of the leading
unit of the pole. Now \(v(z)=v_0+O(z^{-1})\) is holomorphic near
infinity.

After substitution in (CB1), move the right side of the score
equation to the left. The three constant slopes of \(z^m\) in the
resulting exponential sum are \(b,-b,\alpha v_0\). The coefficients
multiplying the exponentials are holomorphic for large \(|z|\).
Each satisfies \(|h(z)|\le C\exp(C|z|^k)\) for a common \(k<m\).
For the disk terms this follows from the bound on \(H\) with argument
\(O(|z|^m)\), and for the right side it follows from
\(z^m(v(z)-v_0)=O(|z|^{m-1})\).

Distinct constant slopes with coefficients of this growth cannot
cancel identically. Combine equal slopes and discard zero
coefficients. An exposed vertex of the remaining finite set of
slopes dominates throughout a fixed interval of ray angles, since
\(e^{im\arg z}\) sweeps every direction. On that interval the sum
identity bounds its coefficient by \(\exp(-c_1|z|^m)\) for large
\(|z|\), whereas on the remaining angles it has the original upper
bound. For a nonzero coefficient, the circular mean of
\(\log|h(re^{i\theta})|\) is finite and convex in \(\log r\) on the
exterior annulus. Convexity follows by rotational averaging of the
subharmonic function \(\log|h|\), since a radial subharmonic function
is convex in the logarithmic radius. Zeros on a circle have integrable
logarithmic singularities. By fixing two radii, convexity gives an
affine lower bound in \(\log r\) for all larger radii. The negative
upper bound of order \(-r^m\) from the fixed angular interval
contradicts that lower bound.

If \(\alpha v_0\notin\{b,-b\}\), the right side's slope cannot have
zero coefficient, a contradiction. If \(\alpha v_0=b\), the
negative disk slope stands alone and its coefficient would have to
vanish identically. Since \(\alpha\ell=b\), we have \(v_0=\ell\), and
the argument \((u+\ell)(v+\ell)\) is an unbounded holomorphic function
of \(z\). Its composition with the nonzero entire function \(H\)
cannot vanish identically by the open mapping and identity theorems.
When \(\alpha v_0=-b\), the same argument applies to the positive
disk slope and \((u-\ell)(v-\ell)\). Thus no such meromorphic germ
exists at any positive score level.

Nevertheless sufficiently small positive levels are closed analytic
contours. The real score is strictly positive, tends uniformly to
zero at infinity, and its outward radial derivative is strictly
negative on \(|x|>(a+R)/\rho\). In fact differentiating the conditional
Gaussian density in direction \(x/|x|\) multiplies it by
\(\rho(y\cdot x/|x|-\rho|x|)/\sigma^2<0\) throughout the source
support. Choose a level smaller than the positive minimum on the
closed disk of radius \((a+R)/\rho\). Each ray meets it once and the
implicit function theorem gives a regular analytic boundary there.

The two disks can both have Gaussian mass \(1/3\), so their
complement also has mass \(1/3\). Fix
\(0<\varepsilon<\Phi^{-1}(2/3)\), where \(\Phi\) is the standard
normal distribution function, and use \(a=R+\varepsilon\). The
positive disk is described by
\((y_1-\varepsilon)^2+y_2^2<2R(y_1-\varepsilon)\). Its mass
is continuous in \(R\), tends to zero as \(R\downarrow0\), and
increases to the mass of the half-space \(y_1>\varepsilon\) as
\(R\to\infty\). This limiting mass exceeds \(1/3\), giving the
desired radius by continuity. Reflection gives the other mass.

The corresponding scores can also have two genuine transverse
triples on the same unbounded exterior winner component. Write
\(U_1,U_2\) for the smoothed disk indicators and \(U_3=1-U_1-U_2\).
At \((0,y_0)\), with \(y_0>R/\rho\) sufficiently large, reflection
gives \(U_1=U_2=t\in(0,1/3)\). Set \(\beta_1=\beta_2=0\) and
\(\beta_3=1-3t\). All adjusted scores tie at \((0,\pm y_0)\).
Along the vertical axis the horizontal derivative of \(U_1-U_2\)
is strictly positive, since the positive disk lies in \(y_1>0\)
and the negative disk is its mirror. The vertical derivative of
\((U_1+U_2)/2-U_3\) is nonzero at the two points because the
conditional mean is above or below the full source support in that
coordinate. Since \(U_1-U_2\) vanishes identically on the axis,
the two score differences have independent gradients and give all
three local winners. Label 3 strictly wins on the vertical rays
beyond these points, where \(U_1=U_2<t\), and connects the two rays
outside a large disk, where it also wins strictly. On the local
\(1/3\) arcs the exterior score equation is \(T_\rho f=3t>0\),
so the branch obstruction applies. The nearby strict winners of
labels 1 and 2 have \(|x_2|>R\) and differ from the original disk
indicators. Thus the example has balanced compact sources and the
local score geometry of two exterior triples. Global agreement of
the sources with the winner rule and stability remain independent
inputs to a branch theorem for maximizing partitions.

## A ball comparison for bounded winner regions

Let \(A\) be a three-label partition in the planar Gaussian space
whose indicators coincide almost everywhere with the strict winner
rule for \(T_\rho1_{A_i}-\beta_i\). Assume each label has positive
mass. For label \(i\) define
\(\theta_i=\max_{j\ne i}(\beta_i-\beta_j)\). Each price difference
lies in \((-1,1)\) because the two labels are strict winners on sets
of positive mass and their smoothed score difference has magnitude
strictly below one. Let \(R_i>0\) be the unique solution of
\[
 \mathbb P\bigl(|\rho R_i e_1+\sigma Z|<R_i\bigr)
                    =\tfrac12(1+\theta_i),\qquad Z\sim\gamma_2.   \tag{CB2}
\]
The probability is strictly increasing and continuous from zero to
one as the radius increases, since the event is membership in
\((R_i/\sigma)B(-\rho e_1,1)\) and this ball contains the origin
in its interior before dilation.

If \(A_i\) contains a ball \(B(p,R)\) up to null sets with \(R>R_i\),
choose \(1<\lambda<1/\rho\) sufficiently close to one that
\(\mathbb P(|\rho\lambda R e_1+\sigma Z|<R)>(1+\theta_i)/2\).
For \(x\in B(p/\rho,\lambda R)\), the conditional Gaussian mean
is at most \(\rho\lambda R\) from \(p\). Its probability of the
original ball is at least the just-displayed probability. Indeed
the Gaussian probability of a centered ball decreases with the
norm of the mean. After rotating the mean onto the first axis,
this follows by conditioning on the other coordinate and using the
monotonicity for symmetric intervals. Hence
\(T_\rho1_{A_i}(x)>(1+\theta_i)/2\), so its unadjusted contrast
against any other label exceeds \(\theta_i\) and \(i\) strictly wins.

The same \(\lambda\) works for every larger ball radius, by the
nested dilations of \(B(-\rho\lambda e_1,1)\). Induction gives balls
\(B(p\rho^{-n},\lambda^nR)\) in the strict winner region. Therefore
a bounded label has inradius at most \(R_i\). This controls the
size of pure balls in a compact competitor without assuming any
regularity of the interface.

The same comparison applies to a union of two labels. Let
\(S=A_1\cup A_2\), \(V_i=T_\rho1_{A_i}-\beta_i\), and suppose there
is an active tie \(V_1=V_2>V_3\). For actual planar maximizers
the [complete adjacency and interface regularity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)
provide such points. With
\[
 p_S=\frac{2+\beta_1+\beta_2-2\beta_3}{3},\qquad
 s(x)=T_\rho1_S(x),
\]
we have
\[
 \max(V_1(x),V_2(x))-V_3(x)
   \ge \tfrac32(s(x)-p_S).
\tag{CB3}
\]
Here \(V_3=1-s-\beta_3\), and the maximum is at least the
mean of the two inside scores. At the active tie this maximum
equals the mean, so \(s>p_S\) there. Thus \(p_S<1\), while
\(\beta_i-\beta_3>-1\) for \(i=1,2\) gives \(p_S>0\).
Wherever \(s>p_S\), the winner belongs to the union almost
everywhere. The ball iteration with threshold \(p_S\) in place of
\((1+\theta_i)/2\) therefore gives inradius at most the unique
solution of (CB2) with right side \(p_S\) whenever \(S\) is bounded.
