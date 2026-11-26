# Pure-cross three-windows for real-conjugate cross ratios

Let \(\Lambda(X^iY^j)=i!j!\). The following every-start result is proved
by an exact Bessel equation and a compressed Chebyshev differentiation
matrix.

**Theorem.** Let \(K\) have characteristic zero and
\(f=BXY+dX+eY\in K[X,Y]\), \(Bde\ne0\).
If \(r=B^2/(de)\) is transcendental over \(\mathbb Q\), or is algebraic
and has a real conjugate, then for every \(q\ge1\) at least one of
\[
                  \Lambda(f^q),\quad\Lambda(f^{q+1}),
                                     \quad\Lambda(f^{q+2})
\]
is nonzero. In particular the conclusion holds over \(\mathbb C\)
whenever \(r\in\mathbb R\setminus\{0\}\), and it holds for any algebraic
\(r\) of odd degree. For arbitrary \(B,d,e\), the windows starting at
\(q=1,2,3\) are also proved below.

A hypothetical counterexample at any later start must have algebraic
ratio \(r\) with no real conjugates. No exclusion of that remaining
case is claimed. The full arbitrary-start binary quadratic problem
and the full two-variable Strong Factorial Conjecture remain outside
the conclusion.

The proofs in Parts A and B are over \(\mathbb C\) with a real
normalized ratio. Part C gives the characteristic-zero transfer and
the three initial starts.

## A.1. The exact Chebyshev reduction

Choose \(s^2=de\), and write
\[
 c=\frac{d+e}{2s},\quad b=\frac B{2s},\quad
 P_j=\frac{\Lambda((f/s)^j)}{j!},\quad
 G_n(c,b)=\sum_{i=0}^n U_n^{(i)}(c)b^i.                       \tag{A.2}
\]
Changing the choice of \(s\) changes both \(c,b\) to their negatives. The formal native expansion is
\[
 \sum_{j\ge0}P_jT^j
 =\sum_{i\ge0}\frac{i!(2b)^iT^i}{(1-2cT+T^2)^{i+1}},
\]
so differentiation of the generating function of \(U_j\) gives \(P_j=G_j\).
The product rule and \(T'_{n+1}=(n+1)U_n\) give the polynomial identities
\[
 G_{n-1}=(c-nb)G_n+b^2\partial_bG_n-T_{n+1}(c),\qquad
 P_{n+1}=2(n+1)bP_n+P_{n-1}+2T_{n+1}(c).                    \tag{A.3}
\]
The same native reduction is also recorded in the local Chebyshev
transfer note; the formulas above are derived here and no physical
source pin is needed as an additional mathematical premise.

Consequently a three-zero block starting at a positive index, with middle index \(n\ge2\), is equivalent to
\[
                T_{n+1}(c)=0,\qquad G_n(c,b)=\partial_bG_n(c,b)=0,
                \qquad b\ne0.                              \tag{A.4}
\]
The roots of \(T_{n+1}\) are real and strictly between \(-1\) and \(1\).
At such a root set
\[
 A=1-c^2>0,\quad H=U_n(c),\quad C=AH.
\]
The elementary Chebyshev identities
\[
 T_{n+1}^2+(1-c^2)U_n^2=1,\qquad
 (1-c^2)U_n'-(cU_n)=-(n+1)T_{n+1}
\]
show \(H\ne0\), \(C\ne0\), and \(C\) is real.

## A.2. A universal polynomial differential equation

The Chebyshev equation and its \(j\)-th derivative are
\[
 (1-c^2)U_n''-3cU_n'+n(n+2)U_n=0,
\]
\[
 (1-c^2)U_n^{(j+2)}-(2j+3)cU_n^{(j+1)}
       +[n(n+2)-j(j+2)]U_n^{(j)}=0.
\]
Multiplying by powers of \(b\) and collecting coefficients, including the two initial coefficients, proves for \(W(b)=G_n(c,b)\)
\[
 -b^4W''-(3b^3+2cb^2)W'
   +[n(n+2)b^2-cb+1-c^2]W
  =(1-c^2)U_n(c)-(n+1)bT_{n+1}(c).                          \tag{A.5}
\]
At a resonance in (A.4) the right side is the nonzero constant \(C\).
In particular every repeated nonzero root of \(W\), if one exists, has multiplicity exactly two:
at \(W=W'=0\), equation (A.5) gives \(-b^4W''=C\ne0\).
An inhomogeneous equation of order two does not, by itself, exclude this double root.

For \(t\ne0\), put
\[
 R(t)=t^nW(1/t),\qquad V(t)=e^{-ct}R(t).
\]
Direct substitution in (A.5), at resonance, yields
\[
 tR''-(2ct+2n+1)R'
       +[(c^2-1)t+(2n+1)c]R=-Ct^{n+1},
\]
\[
                 V''-\frac{2n+1}{t}V'-V=-Ct^ne^{-ct}.       \tag{A.6}
\]
Writing \(V=t^{n+1}Z\) identifies the homogeneous equation with the
modified Bessel equation of order \(n+1\). This is an exact polynomial
change followed by an exponential gauge; it does not use a moment EGF
as a convergent germ at a singular parameter.

## A.3. No real normalized cross coefficient can be a repeated root

Let \(t>0\), \(\nu=n+1\), and
\[
 \psi(t)=t^{n+1}K_{n+1}(t),\quad p(t)=t^{-2n-1},\quad
 \mathcal W(t)=p(t)\big(V'(t)\psi(t)-V(t)\psi'(t)\big).
\]
The standard modified Bessel equation gives
\(\psi''-(2n+1)t^{-1}\psi'-\psi=0\). Equation (A.6) therefore gives
\[
                   \mathcal W'(t)=-C e^{-ct}K_{n+1}(t).     \tag{A.7}
\]
For real \(t>0\), \(K_{n+1}(t)>0\), as follows from its integral
representation. Its large-positive-argument asymptotic, and the same
asymptotic for its derivative, imply \(\mathcal W(t)\to0\) as
\(t\to+\infty\): \(R,R'\) are polynomials, and the exponential factor is
\(e^{-(1+c)t}\), with \(1+c>0\).

If a real \(b_0>0\) were a repeated root, \(t_0=1/b_0>0\) would satisfy
\(V(t_0)=V'(t_0)=0\), hence \(\mathcal W(t_0)=0\). Integrating (A.7) gives
\[
                0=C\int_{t_0}^{\infty}e^{-ct}K_{n+1}(t)\,dt,
\]
which is impossible. For negative \(b\), use the exact parity
\(G_n(-c,-b)=(-1)^nG_n(c,b)\), preserving the resonance. Thus
\[
 T_{n+1}(c)=0,\ b\in\mathbb R\setminus\{0\}
       \quad\Longrightarrow\quad (G_n,\partial_bG_n)(c,b)\ne(0,0). \tag{A.8}
\]

For a phase with \(B^2/(de)>0\), a hypothetical triple would force
\(c\in(-1,1)\), while \(b^2=B^2/(4de)>0\) forces \(b\) real and nonzero.
Equation (A.8) excludes it. Therefore for every integer \(q\ge1\), at
least one of
\(\Lambda(f^q),\Lambda(f^{q+1}),\Lambda(f^{q+2})\) is nonzero.
This conclusion is stated over \(\mathbb C\) with the specified real
ratio; it is not an assertion for arbitrary coefficients over every
characteristic-zero field. The separate opposite-axis case \(e=-d\)
is already closed for all complex \(B\) by the existing Robin-resolvent
argument.

The Bessel inputs used here were checked in the actual primary DLMF:
[modified Bessel equation 10.25.1](https://dlmf.nist.gov/10.25.E1),
[integral representation 10.32.9](https://dlmf.nist.gov/10.32.E9), and
[large-argument formulas 10.40.2 and 10.40.4](https://dlmf.nist.gov/10.40).
They apply for the fixed integer order and positive real argument used
in this proof.

## A.4. The precise nonreal obstruction

For \(\Re t>0\) the horizontal-ray integral
\[
 J_{\nu,c}(t)=\int_t^{\infty}e^{-cu}K_\nu(u)\,du
\]
is analytic, with the ray directed toward \(+\infty\). The same decay
calculation gives
\[
                       \mathcal W(t)=C J_{\nu,c}(t).         \tag{A.9}
\]
Thus a repeated root with \(\Re(1/b)>0\) must make both \(V(t)\) and
\(J_{\nu,c}(t)\) vanish. A zero of \(J\) alone is not sufficient.
Substitution of the Bessel integral representation and absolute
Fubini integration gives
\[
 J_{\nu,c}(t)=\int_0^\infty
     \frac{\cosh(\nu v)}{c+\cosh v}\,
                   e^{-t(c+\cosh v)}\,dv,\qquad \Re t>0.      \tag{A.10}
\]
The integrand has positive weight, but this does not establish complex
zero-freeness. For example the Laplace transform \(1+2e^{-t}\) has zeros
\(\log2+(2j+1)\pi i\) in the right half-plane. A half-plane zero-free
theorem is not imported or asserted. Parity moves any \(b\) with nonzero
real part to this half-plane; purely imaginary \(b\) is a separate
boundary case. The full pure-cross three-window remains open here.

## A.5. Two necessary restrictions on possible shortcuts

The resonance condition is essential even for real \(c\). Explicitly
\[
 G_2(c,b)=4c^2-1+8cb+8b^2
        =8(b+c/2)^2\quad\hbox{if }c^2=1/2,
\]
whereas \(T_3(c)=-c\ne0\).

Nor does real-rootedness of an arbitrary input polynomial prohibit a
nonreal multiple root of its factorial derivative transform. The
polynomial
\[
 p(x)=1+4x+\frac{31}{10}x^2+\frac{11}{15}x^3
                      +\frac{121}{2400}x^4
\]
has four distinct negative real roots: \(2400p\) has alternating signs
at \(0,-1,-2,-5,-10\). Nevertheless
\[
                \sum_{j=0}^4p^{(j)}(0)b^j
                           =(1+2b+\tfrac{11}{10}b^2)^2,
\]
and the quadratic factor has nonreal roots. This is a counterexample
to that broad transform lemma, not to any factorial-window conjecture.

## A.6. Classical derivative-transform context and scope

Csordas--Charalambides--Waleffe,
[A new property of a class of Jacobi polynomials](https://people.math.wisc.edu/~fwaleffe/submitted/Jacobi_CCW.pdf),
Theorem 1, proves left-half-plane stability for derivative transforms
of \(P_n^{(\alpha,\beta)}\) evaluated at the endpoint \(1\), for
\(-1<\alpha\le1\), \(\beta>-1\). Its weighted-integration and Gaussian
quadrature proof is a possible model for a further argument.
The endpoint hypothesis does not include the interior resonant nodes
used above. No endpoint theorem is applied as an interior theorem.

The retained Filaseta--Rocks 2016 theorem concerns the different family
\(\sum_j\binom nj(mj)!x^j\), with a density conclusion on irreducibility.
It does not supply the required all-index Chebyshev-resonant
squarefreeness. No novelty or literature-completeness claim is made.

## B.1. A compression identity for an arbitrary polynomial with simple roots

First let \(P\) be any degree-\(N\) polynomial with distinct roots
\(x_1,\ldots,x_N\), and choose \(c=x_h\). Put
\[
 S_i=\sum_{j\ne i}\frac1{x_i-x_j}=\frac{P''(x_i)}{2P'(x_i)},\qquad
 M_{ii}=S_i,\quad M_{ij}=\frac1{x_i-x_j}\quad(i\ne j).          \tag{B.2}
\]
All arguments may be made over any splitting field of characteristic
zero. On the polynomials of degree at most \(N-1\), use coordinates
\(q(x_i)/P'(x_i)\). The Lagrange interpolation identity
\[
 q(x)=\sum_i \frac{q(x_i)}{P'(x_i)}\frac{P(x)}{x-x_i}
\]
shows that differentiation has matrix \(M\). In particular \(M^N=0\)
and \(\det(I-bM)=1\).

For the finite polynomial resolvent
\(u(x)=\sum_{j=0}^N b^jP^{(j)}(x)\), one has \(u-bu'=P\).
Writing \(q=u-P\), \(a_i=q(x_i)/P'(x_i)\), and \({\bf1}\) for the
all-one column, interpolation gives
\[
                 (I-bM)a=b{\bf1}.
\]
Thus
\[
 \frac{u(c)}{bP'(c)}
  =e_h^t(I-bM)^{-1}{\bf1}
  =\sum_{j=0}^{N-1}\frac{P^{(j+1)}(c)}{P'(c)}b^j.             \tag{B.3}
\]
The quotient at \(b=0\) is understood as its polynomial continuation.

Decompose the coordinate space as
\(\langle{\bf1}\rangle\oplus\ker(e_h^t)\). The lower-right compression
of \(M\) in this decomposition is
\[
 A=(I-{\bf1}e_h^t)M\big|_{\ker(e_h^t)}.
\]
In the coordinates indexed by \(i\ne h\), \(A_{ij}=M_{ij}-M_{hj}\).
The cofactor formula for the upper-left entry of
\((I-bM)^{-1}\), together with \(\det(I-bM)=1\), proves
\[
                      e_h^t(I-bM)^{-1}{\bf1}=\det(I-bA).     \tag{B.4}
\]
This is a cofactor identity, so no inverse is being assumed at a
specialized \(b\): the equality is a polynomial identity.

Conjugate \(A\) by the diagonal matrix \(D_{ii}=c-x_i\), \(i\ne h\).
For \(i\ne j\),
\[
 M_{ij}-M_{hj}
       =\frac{c-x_i}{(x_i-x_j)(c-x_j)}.
\]
Consequently \(K=D^{-1}AD\) has
\[
 K_{ii}=S_i+\frac1{x_i-c},\qquad
 K_{ij}=\frac1{x_i-x_j}\quad(i\ne j,\quad i,j\ne h).           \tag{B.5}
\]

## B.2. Specialization to the Chebyshev nodes

Take \(P=T_N\). Its \(N\) roots \(x_i\) are distinct and belong to
\((-1,1)\). The elementary differential equation
\[
              (1-x^2)T_N''-xT_N'+N^2T_N=0
\]
gives
\[
                    S_i=\frac{x_i}{2(1-x_i^2)}.
\]
Also \(P'=NU_n\). Equations (B.3)--(B.5) therefore prove the exact formula
\[
 \boxed{\quad
 \frac{G_n(c,b)}{U_n(c)}=\det(I-bK),\qquad
 K_{ii}=\frac{x_i}{2(1-x_i^2)}+\frac1{x_i-c},\quad
 K_{ij}=\frac1{x_i-x_j}\ (i\ne j).
 \quad}                                                     \tag{B.6}
\]
Here the indices omit the selected root \(c\), and \(U_n(c)\ne0\)
because \(T_N\) has simple roots.

## B.3. An exact positive energy identity

Let
\[
 J=\operatorname{diag}(x_i-c),\qquad
 H={\bf1}{\bf1}^t+
       \operatorname{diag}\!\left(\frac{1-cx_i}{1-x_i^2}\right).
\]
Direct calculation gives
\[
                              JK+K^tJ=H.                    \tag{B.7}
\]
Indeed the off-diagonal entry is
\((x_i-c-(x_j-c))/(x_i-x_j)=1\), and the diagonal entry is
\(2(x_i-c)[x_i/(2(1-x_i^2))+1/(x_i-c)]\).
All diagonal summands in the second term defining \(H\) are strictly
positive, since \(|c|,|x_i|<1\). Thus \(H\) is positive definite on
complex vectors.

If \(Kv=\lambda v\), \(v\ne0\), then
\[
                 2\Re(\lambda)\,v^*Jv=v^*Hv>0.              \tag{B.8}
\]
No eigenvalue of \(K\) can lie on the imaginary axis. A nonzero root
\(b\) of \(G_n(c,\cdot)\) has \(1/b\) as an eigenvalue of \(K\) by (B.6).
This proves the claimed absence of purely imaginary \(b\).
It also proves, for example, that at the largest node \(c\) all roots
of \(G_n(c,\cdot)\) lie in the open left half-plane, because \(J\) is
then negative definite. At the smallest node they lie in the open
right half-plane.

Equation (B.7) is an accretivity identity in an indefinite metric. It
is not selfadjointness in a Pontryagin space and it does not establish
simplicity of nonreal eigenvalues. Positive definite symmetric part,
even with a definite metric, does not in general exclude Jordan
blocks. No such inference is made here.

## B.4. The exact coefficient sector closed for pure-cross phases

For a phase with nonzero real ratio, choose \(s^2=de\), and set
\(c=(d+e)/(2s)\), \(b=B/(2s)\). The native Chebyshev reduction says that
a three-zero block at positive indices would require
\[
                    T_{n+1}(c)=0,\qquad G_n=G_{n,b}=0,
                    \qquad n\ge2,\ b\ne0.
\]
In particular \(c\in(-1,1)\). If \(B^2/(de)<0\), then \(b\) is purely
imaginary and (B.6)--(B.8) exclude even a single zero of \(G_n\).
If \(B^2/(de)>0\), then \(b\) is real, and the
Bessel argument in Sections A.2--A.3, excludes its being a double root. Hence every phase with nonzero real ratio has at least one nonzero moment in each positive-start block
of length three.

The condition is unchanged by choosing the other square root \(s\).
The statement is over \(\mathbb C\) with an explicit reality condition,
not over arbitrary coefficient fields without a specified embedding.
All phases with a nonreal value of \(B^2/(de)\) remain outside this
proof. The independent opposite-axis theorem \(e=-d\), valid for all
complex cross coefficients, remains available.

## C. Algebraic real conjugates and exact initial starts

### C.1. Transfer to an arbitrary real conjugate

Let \(r=B^2/(de)\) be algebraic over \(\mathbb Q\), with a real conjugate
\(r_0\ne0\). Embed \(\mathbb Q(r)\) in \(\mathbb C\) sending \(r\) to
\(r_0\). This embedding extends to the finitely generated coefficient
field \(\mathbb Q(B,d,e)\): choose an algebraically independent image
of a transcendence basis over the image of \(\mathbb Q(r)\), then
extend the resulting algebraic embedding into the algebraically
closed field \(\mathbb C\). Every moment is an integer polynomial in
\(B,d,e\). A zero block is preserved by this embedding and would
contradict Part B, Section B.4. This proves the theorem over \(K\).
An irreducible polynomial of odd degree over \(\mathbb Q\) has a real
root, which proves the odd-degree specialization.

Conversely, a hypothetical triple at any start already forces \(r\)
to be algebraic, without an algebraicity hypothesis on the original
coefficients. Work in an algebraic closure of the coefficient field.
The reduction gives \(T_{n+1}(c)=0\), so \(c\) is algebraic over
\(\mathbb Q\). Since \(G_n(c,b)\) has degree \(n\) with nonzero leading
coefficient \(2^n n!\), its root \(b\) is algebraic over
\(\mathbb Q(c)\). Therefore \(r=4b^2\) is algebraic. If it had a real
conjugate, the preceding argument would give a contradiction. Any
counterexample ratio thus has no real conjugates and has even
minimal-polynomial degree.

### C.2. The first three starts without any coefficient restriction

Embed the finitely generated coefficient field in \(\mathbb C\).
At a hypothetical triple the middle index is \(n=q+1\), and \(c\)
is real. Part A excludes real repeated roots; Part B excludes purely
imaginary roots. A nonreal repeated root of the real polynomial
\(G_n(c,\cdot)\) comes with its distinct conjugate, also repeated.
For \(n=2,3\), degree alone rules this out. These are the windows
starting at \(q=1,2\).

For \(n=4\), such a pair would exhaust all four roots and force
\[
 G_4(c,b)=384(b^2+u b+v)^2,\qquad u,v\in\mathbb R.
\]
Native differentiation of \(U_4=16c^4-12c^2+1\) gives
\[
 G_4=384b^4+384cb^3+(192c^2-24)b^2
              +(64c^3-24c)b+16c^4-12c^2+1 .
\]
Matching the coefficients of \(b^3,b^2,b\) gives
\[
       u=c/2,\qquad v=c^2/8-1/32,\qquad
                         cv=c^3/6-c/16.
\]
If \(c\ne0\), these equations force \(c^2=3/4\), \(v=1/16\).
But the constant coefficient of the square would be \(384v^2=3/2\),
whereas \(U_4(c)=1\). If \(c=0\), the square would have constant
\(384/32^2=3/8\), again different from \(U_4(0)=1\).
Thus \(G_4(c,\cdot)\) has no repeated nonreal root for any real \(c\).
Part A handles a real repeated root at the resonant \(c\), proving
the window starting at \(q=3\).

These are exact separate low-degree arguments, not an extrapolation
from a finite scan to all starts.

## D. Exact supplemental controls and remaining issue

The complete uniform Bessel and compression proofs above are the
mathematical premises. Exact supplemental controls independently
reconstruct the Chebyshev polynomials, check 75 universal polynomial
identities and 50 reciprocal/gauge identities for indices 1--25,
and test 81 arbitrary rational-node compression instances and 3,894
energy entries. They reject five identity corruptions and four
diagonal corruptions. Finite replay is a transcription check only.

No general simplicity theorem is inferred from strict accretivity.
The remaining problem is a possible pair of nonreal conjugate
double roots of \(G_n(c,\cdot)\) at a real Chebyshev resonance, for
\(n\ge5\), whose squared-root parameter has no real conjugate.
The endpoint Jacobi stability result and positivity of the
incomplete-\(K\) integral do not settle that issue.
