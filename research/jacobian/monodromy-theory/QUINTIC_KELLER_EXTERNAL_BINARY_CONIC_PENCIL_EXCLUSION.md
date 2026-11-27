# The external linear factor of a binary quintic conic pencil

**Theorem QuinticKellerExternalBinaryPencil.** Let \(p,q\) be coprime homogeneous quadrics in two independent linear forms in three variables, and suppose the linear form \(\ell\) is independent of those two forms. There is no degree-five polynomial Keller map with leading part
\[
                            (\ell p^2,\ \ell pq,\ \ell q^2).
\]
The result holds over \(\mathbb C\) and hence over characteristic-zero fields by embedding the finitely generated coefficient field in \(\mathbb C\).

## Weighted normalization

The projective pencil of coprime binary quadrics contains two independent squares. The square locus in the projective plane of binary quadrics is a nonsingular conic. A pencil line meets it twice, and it is tangent at a square \(u^2\) only if both generators are divisible by \(u\). Coprimality therefore gives two distinct squares. A pencil change induces a constant target change on the three Veronese components. Choose binary coordinates \(x,y\) and third coordinate \(z=\ell\) to get
\[
                        H_5=z N(x^2,y^2),\qquad N(u,v)=(u^2,uv,v^2).
\]

Write \(p_0=x^2,\ q_0=y^2,\ \nu=(q_0^2,-2p_0q_0,p_0^2)\). For a homogeneous column \(U\) of degree \(j\), put \(S=\nu\cdot U\). The determinant contribution using two columns of \(JH_5\) and one of \(JU\) is
\[
             \frac{z\{p_0,q_0\}}2(5z S_z-jS),
 \tag{1}
\]
where \(\{\cdot,\cdot\}\) is the binary Jacobian. Indeed, if \(D=N(p_0,q_0)\), then \(D_x\times D_y=2\{p_0,q_0\}\nu\) and \(4D=xD_x+yD_y\). The rows of the adjugate of \(J(zD)\) are therefore the three multiples of \((D_x\times D_y)^{\mathsf T}\) with coefficients \(-zx/4,-zy/4,z^2\). Taking the trace against \(JU\) and using Euler's identity for \(U\) gives (1). At ordinary degree eleven, only \(U=H_4\) contributes. The integer \(5k-4\) is nonzero at every \(z\)-power \(k\), so \(\nu\cdot H_4=0\). The coprimality syzygy gives \(H_4=DN(p_0,q_0)(U_2,V_2)\) for ordinary homogeneous quadratics \(U_2,V_2\). In particular \(\deg_z H_4\le2\).

For source weights \((1,1,2)\), the highest term of \(F(tx,ty,t^2)\) is \(t^6D_*\), with
\[
                 D_*=N(x^2,y^2)+[z^2]H_4+[z^3]H_3 .
\]
The three binary terms have degrees four, two and zero, respectively, so \(D_*\) is even under \((x,y)\mapsto(-x,-y)\). The determinant in \(x,y,t\) is \(2\kappa t^3\), where \(\kappa=\det JF\ne0\). Its \(t^{17}\) coefficient gives \(\det((D_*)_x,(D_*)_y,D_*)=0\). Homogenize \(D_*\) to degree four and apply Wang and Yan's [Lemma 2.4, printed pages 6–7](https://arxiv.org/pdf/2610.36495v1). Its three entries are independent and have no common homogeneous factor, since their degree-four parts are \(x^4,x^2y^2,y^4\), with gcd one. The factorization is therefore in two linear forms with binary quartic entries, or in two quadratic forms with binary quadratic entries.

In the linear-argument case the two linear forms at infinity are independent. Choosing their basis and comparing the degree-four terms gives
\(D_*=N((x+a)^2,(y+b)^2)\). Evenness of the first and third coordinates forces \(a=b=0\).

In the quadratic-argument case, the three independent binary quadratics are a constant invertible linear combination of \(N\). At infinity the two arguments have nonconstant ratio and the resulting projective image is the standard Veronese conic. Hence the constant target transformation preserves that conic. Its projective action is induced by a change of the binary arguments because the conic parametrization by \(\mathbb P^1\) is an isomorphism. Absorb the scalar by rescaling those arguments over \(\mathbb C\). Comparing the degree-four terms now gives \(N(x^2,y^2)\) as their leading Veronese vector, so their degree-two parts are \(x^2,y^2\) up to a common sign. The square of each argument is even, and its nonzero even leading term precludes an odd argument. Thus the arguments have no linear terms.

The two factorization cases together give \(D_*=N(x^2+\alpha,y^2+\beta)\). Equivalently, in the original variables the leading part in weights \((1,1,2)\) is
\[
           \mathcal H=(z p^2,z pq,z q^2),\qquad
           p=x^2+\alpha z,\quad q=y^2+\beta z .
 \tag{2}
\]

## A nonzero shift

Assume \((\alpha,\beta)\ne(0,0)\). Let \(R\ne0\) be weighted homogeneous of positive weight \(d\le11\) with \(\det J(\mathcal H_1,\mathcal H_2,R)=0\). In the chart \((tx,ty,t^2)\), put
\[
       P=x^2+\alpha,\quad Q=y^2+\beta,\quad r=Q/P,\quad w=t^3P,\quad R_*=R(x,y,1).
\]
The three functions in the determinant become \(w^2,rw^2,t^d R_*\). Thus
\[
             dw\wedge dr\wedge d(t^d R_*)=0,\qquad
             d(R_*^3/P^d)\wedge dr=0 .
\]
The second equality follows by taking the degree-zero ratio \((t^dR_*)^3/w^d\) and using \(w_t\ne0\). The generic projective conic \(Q-\lambda P=0\) over \(\mathbb C(\lambda)\) is smooth, since the ternary quadratic matrix in \(x,y,s\) is diagonal with entries \(-\lambda,1,\beta-\lambda\alpha\). It has a rational point from any common zero of \(P,Q\), so its function field is rational over \(\mathbb C(\lambda)\). Thus \(\mathbb C(r)\) is relatively algebraically closed in \(\mathbb C(x,y)\), and \(R_*^3/P^d\in\mathbb C(r)^\times\).

At either point at infinity of the generic conic after extending the constants, its intersection with \(s=0\) is simple and \(x\ne0\). The rational function \(P=(x^2+\alpha s^2)/s^2\) has a pole of order two there. Valuations give \(3\mid 2d\), so \(3\mid d\). Write \(d=3k\). Applying the same differential argument to \(t^dR_*/w^k\) gives \(R_*/P^k\in\mathbb C(r)\). The change \((x,y)\mapsto(-x,-y)\) fixes \(r,P\) and acts on \(R_*\) by \((-1)^d\), so \(d\) is even. Consequently \(d=6\), and restoring the weights gives
\[
                          R=z p^2 a(q/p),\qquad a\in\mathbb C(T).
\]
The two polynomials in (2) are coprime and neither \(p\) nor a member \(q-cp\) divides \(z\). Therefore \(a\) has no finite pole and is a polynomial of degree at most two. For a finite pole, polynomiality would force \(q-cp\mid z\), and at degree greater than two it would force \(p\mid z\). Thus \(R\) belongs to the span of \(\mathcal H\).

Now \(G=F_2^2-F_1F_3\) has weight at most eleven after leading cancellation and satisfies \(\det J(F_1,F_2,G)=-\kappa F_1\), just as in the
[ordinary conic-pencil exclusion](QUINTIC_KELLER_CONIC_PENCIL_FIRST_INTEGRALS.md#the-veronese-equation).
If \(G\) has positive leading weight \(d\), the part of weight \(d+8\) of this determinant is \(\det J(\mathcal H_1,\mathcal H_2,G_d)\). It vanishes because the right side has weight six. The first-integral conclusion shows \(d=6\) and \(G_d\) lies in the span of \(\mathcal H\). Subtract the corresponding linear combination of \(F\) from \(G\). Its weight is at most five and its determinant with \(F_1,F_2\) still has weight six on the right. It cannot have a positive-weight leading term by the same argument, and a constant term alone gives zero determinant. This contradiction excludes every nonzero shift.

## Zero shifts

If \(\alpha=\beta=0\), then \([z^2]H_4=[z^3]H_3=0\) and \(F=C+zB+z^2A\), with \(A\) affine and \(B_4=N(x^2,y^2)\). The three entries of \(B_4\) are independent, as are their transforms by any invertible constant target change. The following rank arguments use only this independence and the degree bounds, as in the
[scalar weighted exclusion](QUINTIC_VERONESE_SCALAR_WEIGHTED_EXCLUSION.md).

Rank three for \(A\) is impossible by the \(z^5\) determinant coefficient \(2\det(A_x,A_y,A)\), since three independent affine binary forms have nonzero such determinant. In rank one put \(A=(0,0,h)\). The highest binary part of the \(z^3\) coefficient is
\[
        (2-\deg(h)/4)h_{\rm top}\{B_{1,4},B_{2,4}\}.
\]
Indeed the coefficient is \(2h\{B_1,B_2\}+B_1\{B_2,h\}-B_2\{B_1,h\}\), and Euler's identity gives the displayed top part. It is nonzero by independence.

In rank two write \(A=(a,b,0)\) and \(B=(p_1,q_1,r_1)\). If \(a,b\) have independent linear parts, normalize them to \(x,y\). The \(z^4\) coefficient gives \((1-2E)r_1=0\), where \(E=x\partial_x+y\partial_y\). This forces \(r_1=0\), contrary to \(B_4\) independent. Otherwise normalize to \(A=(x,1,0)\). The \(z^4\) coefficient now forces \((r_1)_y=0\), and \(z^3\) gives
\[
        (C_3)_y=r_1'(p_1)_y+(r_1/2-xr_1')(q_1)_y .
\]
Here \(C_3\) is the third entry of \(C\); the formula follows by expanding the determinant with the three \(z^2\) coefficients \(x,1,0\). The degree-seven part forces \((q_{1,4})_y=0\) because \(r_{1,4}\) is a nonzero multiple of \(x^4\) and all other terms have degree at most six. Hence \(q_{1,4}\) and \(r_{1,4}\) are proportional, again contradicting independence. Thus \(A=0\).

The affine-in-\(z\) map has \(\det(B_x,B_y,B)=0\) from the \(z^2\) coefficient and \(B\) has no common affine zero. Homogenize \(B\) to degree four and use the same factorization as for \(D_*\). Two affine linear arguments either have a common zero or have dependent highest parts, so the linear-argument case is impossible. The quadratic case has \(B=N(P_2,Q_2)\) for affine quadratics with leading parts \(x^2,y^2\), by the same matching at infinity used for \(D_*\). These quadratics have a common projective zero by Bézout after homogenization, and none on the line at infinity. The resulting common affine zero contradicts the nonvanishing of the \(z\) column. This finishes the zero-shift case.
