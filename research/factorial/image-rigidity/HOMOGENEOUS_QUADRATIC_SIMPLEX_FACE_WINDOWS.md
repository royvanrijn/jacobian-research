# Uniform factorial windows for homogeneous quadratics in arbitrary dimension

For \(m\ge1\) set \(R=2^m-1\) and \(W_m=2R+R^2(R+1)/2\).
For every nonzero homogeneous quadratic \(F\) in \(m\) variables
over a characteristic-zero field, every window of \(W_m\) consecutive
positive orders contains a nonzero ordinary factorial moment.
The bound is independent of the rank, the coefficients, and the
starting order. It can exceed the actual support size.

The same window holds for a homogeneous \(F\) of any degree if its
restriction \(q\) to the simplex \(\sum_i X_i=1,\ X_i\ge0\) has the
following property. On the complex affine span of each nonempty face
\(P\), including the whole simplex, \(q_P\) is constant, admits a
constant vector \(v\) with \(\partial_v q_P=1\), or has the form
\(c+Q(x-x_0)\) with \(x_0\) in that affine span and \(Q\) homogeneous
of positive degree.

## A differential system on the faces

Use Euclidean measure on each face in affine coordinates, with the
full simplex parametrized by \(x_1,\ldots,x_{m-1}\) and
\(x_m=1-\sum_{i<m}x_i\). Dirichlet integration identifies
\[
 \int_\Delta q(x)^n\,dx=\frac{\Lambda(F^n)}{(n\deg F+m-1)!},\qquad
 G(t)=\int_\Delta\frac{dx}{1-tq(x)}=zJ_\Delta(z),\quad z=1/t ,
\]
where \(J_P(z)=\int_P(z-q_P)^{-1}\) is defined initially for large \(|z|\).
A vertex is given the point mass of mass one.

If \(q_P\equiv c_P\), then \(J'_P=-(z-c_P)^{-1}J_P\).
If \(q_P\) admits the translation vector \(v\), the divergence theorem
applied to \(v/(z-q_P)\) gives
\[
 J'_P=-\sum_{Q\text{ facet of }P} b_{P,Q}J_Q,\qquad
 b_{P,Q}=v\cdot n_{P,Q},
 \tag{1}
\]
where \(n_{P,Q}\) is the outward unit normal in the span of \(P\).
For \(q_P=c_P+Q_e(x-x_P)\) with \(Q_e\) homogeneous of degree \(e\),
use \(V=x-x_P\). Its normal component \(b_{P,Q}=V\cdot n_{P,Q}\)
is constant on each facet. Euler's identity gives
\[
 J'_P=\frac{k-e}{e(z-c_P)}J_P
      -\sum_{Q\text{ facet of }P}\frac{b_{P,Q}}{e(z-c_P)}J_Q,
 \qquad k=\dim P .
 \tag{2}
\]
Indeed the divergence is
\((k-e)/(z-q_P)+e(z-c_P)/(z-q_P)^2\), and
\(J'_P=-\int_P(z-q_P)^{-2}\).
These formulas give a rational system \(Y'=A(z)Y\) with one
coordinate for every nonempty face.
Integration by parts on the simplex also underlies one-loop
relations that express derivatives through lower-dimensional
quadratic-denominator integrals.
[Bern, Dixon and Kosower, equations (3.1)–(3.4)](https://arxiv.org/pdf/hep-ph/9306240)

Every affine quadratic satisfies the face condition. Write its
gradient in affine coordinates as \(Hx+\ell\) for \(H\) symmetric.
If \(\ell\in\operatorname{im} H\), translating to a solution of
\(Hx+\ell=0\) gives a constant plus a homogeneous quadratic.
Otherwise some \(v\in\ker H\) has \(\ell\cdot v\ne0\), since
\((\ker H)^\perp=\operatorname{im}H\) for the complex bilinear
pairing. Rescaling \(v\) gives the translation case.

## Decay of every full-face output

A constant-zero face can have coordinate \(C/z\) in the system.
If a nonconstant \(P\) has \(c_P=0\) in (2), each of its
constant-zero facets \(Q\) has \(b_{P,Q}=0\). Its complex affine
hyperplane is contained in the zero set of \(Q_e(x-x_P)\).
If that plane did not contain \(x_P\), its dilates about \(x_P\)
would contain the complement of a linear hyperplane translated
through \(x_P\). Homogeneity would give \(Q_e=0\) identically.
So \(x_P\) lies in the facet hyperplane and the normal component
is zero.

Let \(E\ge1\) be the largest degree in the chosen equations (2),
or \(1\) when there are none, and put \(\varepsilon=1/E\).
On a fixed small ray \(z=\rho e^{i\vartheta}\), every solution of
the whole system has, on faces with \(q_P\not\equiv0\),
\[
 J_P=O\!\left(\rho^{-1+\varepsilon}
                         (1+|\log\rho|)^{\dim P}\right).
 \tag{3}
\]
Prove this by induction through the faces. A nonzero constant face
has bounded coordinate. In (1), or (2) with \(c_P\ne0\), the
coefficients are bounded near zero and every forcing term is
\(O(\rho^{-1})\) by induction. Integrating with a bounded
integrating factor gives \(O(1+|\log\rho|)\).
For (2) with \(c_P=0\), the constant-zero facets have zero
coefficient. Write \(\alpha=k/e-1\ge-1+\varepsilon\). The integrating
factor is \(z^{-\alpha}\), and the forcing times the integrating
factor is \(O(\rho^{-2+\varepsilon-\alpha}(1+|\log\rho|)^{k-1})\).
If \(\alpha=-1+\varepsilon\), integration adds at most one logarithm.
Otherwise it contributes \(O(\rho^{-1+\varepsilon-\alpha}
(1+|\log\rho|)^{k-1})\). Restoring \(z^\alpha\) gives (3)
in both cases, including the homogeneous solution \(Cz^\alpha\).

A nonzero homogeneous phase has \(q\not\equiv0\), by scaling to
the affine hyperplane outside \(\sum X_i=0\).
Thus \(zJ_\Delta\to0\) for every solution continued to a small
ray. A full-face output cannot be a nonzero polynomial in \(1/z\).

## Coefficient-uniform width

Set \(b(z)=z e_\Delta^{\mathsf T}\), \(b_0=b\), and
\(b_{j+1}=b'_j+b_j A\). Let \(r\) be the first dependent index
of the rows \(b_0,\ldots,b_r\) over \(\mathbb C(z)\), so
\(1\le r\le R\). The resulting scalar equation of order \(r\)
has exactly the output space of \(Y'=AY\) near an ordinary
point. The first \(r\) rows have rank \(r\) there, so any scalar
initial jet through order \(r-1\) comes from a vector initial
value, and uniqueness for the scalar equation identifies the
solutions. If a polynomial in \(1/z\) were a solution, continue
its output realization along a path to a ray at zero avoiding
the finite set of system poles and use (3). The polynomial
would have to be zero.

Let \(D(z)\) be the product of distinct factors \(z-c_P\) used in
the face system, and let \(d=\deg D\le R\).
The matrix \(DA\) is polynomial of degree at most \(d\).
The rows \(L_j=D^j b_j\) are polynomial of degree at most \(1+jd\)
by the numerator \(DL'_j-jD'L_j+L_j(DA)\) at the next step.
In \(r\) columns where the first \(r\) rows are independent,
signed minors of \(L_0,\ldots,L_r\) supply a dependence
\(\sum_j c_j L_j=0\) with \(c_r\ne0\).
It holds on all columns because restriction to those columns
is injective on the span of the first \(r\) full rows and that
span contains \(L_r\). The scalar equation is therefore
\[
        \sum_{j=0}^r c_jD^j\,\partial_z^j G(1/z)=0.
\]
The degrees of its coefficients are at most
\(K=r+d\,r(r+1)/2\): the minor omits the row of degree at
most \(1+jd\), and the power \(D^j\) restores at most \(jd\).

With \(\theta=t\partial_t\), replacing \(z\) by \(1/t\) takes
\(\partial_z^j\) to
\((-1)^j t^j\theta(\theta+1)\cdots(\theta+j-1)\).
Multiplying by \(t^K\) gives a nonzero operator
\(\mathcal L=\sum_{i=0}^{K+r} t^i P_i(\theta)\) for \(G\)
with no nonzero polynomial solution. Suppose that a block of
\(W_m\) zero coefficients of \(G\) starts at \(n\ge1\), and let
\(B(t)\) be the Taylor truncation before degree \(n\).
Then \(\mathcal L(G-B)\) has order at least \(n+W_m\) and
\(\mathcal L B\) has degree at most \(n-1+K+r<n+W_m\).
The equality \(\mathcal L B=-\mathcal L(G-B)\) gives
\(\mathcal L B=0\), contrary to \(B(0)=G(0)=\operatorname{vol}(\Delta)\).
Here \(K+r\le2R+R^2(R+1)/2\).
Finally the field generated over \(\mathbb Q\) by the finitely
many coefficients of a phase embeds in \(\mathbb C\), which
preserves all factorial moment equalities. The same proof works
for the stated face conditions over any characteristic-zero field
by including their finitely many data in the field of definition.
