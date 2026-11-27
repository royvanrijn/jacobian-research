# Excluding the maximal-z quintic Veronese leading part in dimension three

## The degree-four source and the quintic scope

Wang and Yan, [*Quartic Keller Maps and Tameness in Three Variables*](https://borowiecki.dev/pdf/2610.36495), Theorems 1.1 and 3.1, prove invertibility for a polynomial map over \(\mathbb C\) in three variables with nonzero constant Jacobian determinant and ordinary total degree at most four. The ordinary degree is the maximum over the three coordinates. The characteristic-zero invertibility consequence follows by embedding the finitely generated coefficient field in \(\mathbb C\) and using the unique formal inverse at the origin after translation. The authors' separate Theorem 4.1 concerns tameness over \(\mathbb C\) in that same degree range. Their Proposition 3.11 handles independent quartic leading parts and Proposition 3.18 handles dependent ones. Neither result has degree-five hypotheses.

Their Lemma 3.6 excludes the quartic leading part \(g(x^2,xy,y^2)\) with \(\deg g=\deg_z g=2\). Its adjugate and highest-\(z\) coefficient strategy gives the following restriction in degree five. The exponent changes in the highest-\(z\) calculation and the surviving \(z^2\) term in the homogeneous cubic need to be retained.

## The maximal \(z\)-degree case

**Proposition QuinticVeroneseResidual.** Let \(F\) over a characteristic-zero field have three coordinates, ordinary total degree five and nonzero constant Jacobian determinant. Write \(H_j\) for its homogeneous part of degree \(j\). Assume
\[
H_5=gN,\qquad N=(x^2,xy,y^2),\qquad \deg_z g=\deg g=3.
\]
Rescale the target to make \([z^3]g=1\), and translate \(x,y\) in the source as below. After those operations \([z^3]H_4=0\), \([z^3]H_3=0\) and \([z^2]H_2=0\). There are constants \(a,b\) such that
\[
[z^2]H_3=DN(a,b),\qquad [z]H_1=N(a,b),\qquad (a,b)\ne(0,0).
\]
Here \(DN(a,b)=aN_x+bN_y\), with the derivatives evaluated at \((x,y)\). In particular, the case \([z^2]H_3=0\) after the indicated normalization is excluded.

Put \(\nu=(y^2,-2xy,x^2)\) and \(X=(x,y,z)^{\mathsf T}\). Using \(N_x\times N_y=2\nu\), \(N_y\times N=-x\nu\) and \(N\times N_x=-y\nu\) in the three rows of the adjugate gives
\[
\operatorname{adj}(JH_5)=g(5g e_3-g_z X)\nu^{\mathsf T}.
\]
Euler's identity gives the last row in the form \(g(5g-zg_z)\nu^{\mathsf T}\). The contribution using two columns of \(JH_5\) and one of \(JH_j\) is the trace of this adjugate times \(JH_j\). For \(S=\nu\cdot H_j\), it equals \(g(5g S_z-jg_z S)\). If \(\deg_z S=e\) and \(5e-3j\ne0\), its top coefficient is \((5e-3j)[z^e]S\), at \(z\)-degree \(e+5\).

Write \(\Delta_d\) for the degree \(d\) part of \(\det JF\). A contribution from homogeneous indices \(i,j,\ell\) has \(i+j+\ell=d+3\) and \(z\)-degree at most \(\deg_z H_i+\deg_z H_j+\deg_z H_\ell-1\), because one of the three columns differentiates with respect to \(z\). Every \(\Delta_d\) for \(d>0\) vanishes. At \(d=11\) the only triple is \((5,5,4)\). Since \(5e-12\) never vanishes for integral \(e\), \(\nu\cdot H_4=0\).

For any column \(H\) with \(\nu\cdot H=0\), reduction modulo \(x\) and \(y\) shows that \(H_1\) is divisible by \(x\) and \(H_3\) by \(y\), where the subscripts here index coordinates. The middle coordinate then gives \(H=DN(U,V)\) for polynomials \(U,V\). Thus \([z^3]H_4=DN(U_3,V_3)\) for constants \(U_3,V_3\). Translate \((x,y)\) by \((-U_3,-V_3)\). This shifts \([z^3]H_4\) by \(DN(-U_3,-V_3)\) from \(H_5\), since \([z^3]H_5=N\). The resulting \(H_4\) has \(\deg_z H_4\le2\).

At \(d=10\) the triple \((5,4,4)\) has \(z\)-degree at most six. For the other triple \((5,5,3)\), the formula with \(j=3\) forces \(\deg_z(\nu\cdot H_3)\le1\). No nonzero constant column is annihilated by \(\nu\) as a polynomial identity, so \([z^3]H_3=0\). The syzygy for the linear column \([z^2]H_3\) gives \(DN(a,b)\). At \(d=9\) the triples \((5,4,3)\) and \((4,4,4)\) again have \(z\)-degree at most six. Applying the formula to \((5,5,2)\) excludes a nonzero \(z^2\) coefficient in \(\nu\cdot H_2\) and hence in \(H_2\).

At \(d=8\) the triples \((5,4,2)\) and \((4,4,3)\) have \(z\)-degree at most five. The coefficient at \(z^6\) in the \((5,3,3)\) contribution depends only on the already prescribed \(z^3\) and \(z^2\) terms. Compare it to
\[
zN(zx+a,zy+b)=z^3N(x,y)+z^2 DN(a,b)+zN(a,b).
\]
This polynomial map has zero determinant since its three derivative columns lie in the span of the two derivatives of \(N\) at \((zx+a,zy+b)\), by the degree-two Euler identity for \(N\). The contribution of \((5,3,3)\) at \(z^6\) therefore cancels that of \((5,5,1)\) with \([z]H_1=N(a,b)\). Subtracting that vanishing coefficient from the actual coefficient gives \(2\nu\cdot([z]H_1-N(a,b))=0\). The difference is a constant column, so it vanishes. Finally \(F_z(0,0,z)=N(a,b)\), since the identified higher \(z\) powers are zero on the axis. The Keller condition forces this column to be nonzero.

## Laurent exclusion of the remaining coefficients

**Theorem QuinticVeroneseMaximalZExclusion.** Over a characteristic-zero field there is no three-variable Keller map with ordinary degree five and leading part \(H_5=g(x^2,xy,y^2)\), where \(g\) is homogeneous cubic and \(\deg_z g=3\).

Use the normalization and coefficient identities of the proposition. On \(z\ne0\), make the Laurent change of source coordinates
\[
 (x,y,z)=((u_1-a)/z,(u_2-b)/z,z).
\]
Its Jacobian determinant is \(z^{-2}\). The transformed map has determinant \(\kappa z^{-2}\) for some \(\kappa\in k^\times\), and has the form
\[
 zN(u_1,u_2)+M_0(u_1,u_2)+z^{-1}M_1(u_1,u_2)+\cdots .
\]
The only terms of positive \(z\)-degree before collecting are \(z^3N(x,y)\), \(z^2DN(a,b)\), and \(zN(a,b)\). Together they become \(zN(u_1,u_2)\). All other terms have \(z\)-degree at most zero in the new coordinates by the bounds on the \(z\)-degrees of the \(H_j\). All Laurent coefficients are polynomial in \(u=(u_1,u_2)\).

Write \(E=u_1\partial_{u_1}+u_2\partial_{u_2}\) and \(\nu=(u_2^2,-2u_1u_2,u_1^2)\) in these coordinates. If the next coefficient after \(zN(u)\) is \(z^{-r}M(u)\), \(r\ge0\), the determinant coefficient at \(z^{1-r}\) is
\[
 -\nu\cdot(E+2r)M.
\]
Indeed, taking the \(z\)-derivative on \(z^{-r}M\) contributes \(-2r\nu\cdot M\). Taking the \(u_1\)- or \(u_2\)-derivative on \(z^{-r}M\) contributes \(-u_1\nu\cdot M_{u_1}-u_2\nu\cdot M_{u_2}\) by the two cross products in the proposition. Two uses of \(M\) have \(z\)-degree at most \(-2r<1-r\), and a lower Laurent term contributes only below \(1-r\).

For \(r=0,1,2\) the coefficient at \(z^{1-r}\) must vanish. Decompose \(M\) by total degree in \(u\). Each homogeneous piece \(M_j\) satisfies \((j+2r)\nu\cdot M_j=0\), because the corresponding homogeneous products have distinct degrees. If \(r=0\), remove the constant piece by translating the target. Every other piece satisfies \(\nu\cdot M_j=0\). The syzygy in the proposition gives \(M=DN(u)(U(u),V(u))\) for polynomials \(U,V\). Precompose by
\[
 (u,z)\longmapsto\left(u-z^{-r-1}(U(u),V(u)),\,z\right).
\]
Expanding the quadratic \(N\) shows that this cancels \(z^{-r}M\), leaving \(zN(u)\) and terms of \(z\)-degree strictly less than \(-r\). This precomposition has determinant \(1+O(z^{-r-1})\) as a Laurent polynomial in \(z^{-1}\). Therefore the determinant still has leading Laurent term \(\kappa z^{-2}\). These three steps are substitutions of polynomials in \(u\) with Laurent coefficients; global invertibility of the substitutions is unnecessary.

For \(r=3\), take \(M\) to be the coefficient at \(z^{-3}\) after the three cancellations, allowing \(M=0\). The same determinant expansion gives
\[
 \kappa=-\nu\cdot(E+6)M .
\]
The right side vanishes at \(u=0\), contradicting \(\kappa\ne0\).

The [quadratic-in-\(z\) proof](QUINTIC_VERONESE_QUADRATIC_Z_EXCLUSION.md)
excludes \(\deg_z g=2\), and the
[weighted square exclusion](QUINTIC_VERONESE_SQUARE_EXCLUSION.md)
completes \(\deg_z g=1\). The
[binary factor exclusion](QUINTIC_VERONESE_BINARY_RANK_EXCLUSIONS.md)
settles \(\deg_z g=0\), completing the Veronese leading class.
The [quintic leading-part reduction](QUINTIC_KELLER_LEADING_PART_RESIDUAL.md) gives the last independent normal form and the dependent branch. The [completion exclusion](QUINTIC_KELLER_TERNARY_PENCIL_TRANSVERSE_PARAMETER.md#independent-degree-five-leading-parts) eliminates that independent form.

## Collision construction targets

The current construction bounds are [dimension twelve and support
thirty-six](../cubic-reduction/PROOF12.md) for integral cubic collisions,
[dimension eighteen](../subtree-shuffling/CUBIC_INDEX_SEVENTEEN_SHUFFLE_FAILURE.md)
for a rational homogeneous cubic collision, and
[dimension eighty-nine](../cubic-reduction/POWER89.md) for an integral
cube-linear map with a rational collision. Support counts nonzero
nonlinear output-monomial pairs in the twelve-variable example. Each
optimization must preserve the full constant Jacobian identity and an
explicit collision in its own class. All four figures are construction
upper bounds. The maximal-\(z\) exclusion concerns a dimension-three
quintic class and leaves these construction bounds as separate targets.
