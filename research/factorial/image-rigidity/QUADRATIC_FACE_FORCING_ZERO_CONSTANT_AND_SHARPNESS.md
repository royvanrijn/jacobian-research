# Zero-constant quadratic forcing and sharpness of the three-window

Let K be a characteristic-zero field, Lambda(X^d)=d!,
h_j(r,z)=Lambda((rX^2+X+z)^j), and T_j=h_j+jh_(j-1), j>=1.
The zero-constant boundary has a shorter forcing window than the universal
three-window.

**Theorem 1 (zero-constant two-window).** For every r in K and every k>=1,
T_k(r,0) and T_(k+1)(r,0) are not both zero.

The universal cubic syzygy, for i>=2 and D_i=2i-1+2z, is

 z^(i-1)Pi_i(z)=D_i T_(i+1)
  +2(i+1)[1-r(2i+1)D_i]T_i
  -i(i+1)(D_i+2)(1-4rz)T_(i-1),
 Pi_i(z)=2(z+i)^3+3(z+i)^2+i.                                  (1)

For a self-contained derivation, the factorial derivative identity gives

 h_n=2rn(2n-1)h_(n-1)+n(n-1)(1-4rz)h_(n-2)+z^(n-1)(z+n), n>=2.

Write T_i-i(1-4rz)T_(i-1)=2riD_i h_(i-1)+z^(i-1)(z+i).
Combine this equation at i and i+1 using h_i=T_i-i h_(i-1);
the h terms cancel and give(1), without a coefficient-dependent division.
At z=0 and ONLY for i>=2, it becomes the homogeneous backward relation

 T_(i+1)=2(i+1)[r(2i+1)-1/(2i-1)]T_i
          +[i(i+1)(2i+1)/(2i-1)]T_(i-1).                       (2)

The coefficient of T_(i-1) is a nonzero rational integer fraction. Thus if
a consecutive pair starting at k>=2 were zero, equation(2) at i=k would
force T_(k-1)=0. Repeating reaches T_1=T_2=0. A pair starting at k=1 is
already in that base case. Direct native evaluation gives

 T_1(r,0)=2r+2,
 T_2(r,0)=24r^2+16r+4.

The first equation forces r=-1, where T_2=12!=0. This proves Theorem1
without dividing by r or any coefficient-dependent pivot. The restriction
i>=2 matters: at i=1 the left side of(1) is Pi_1(0)=6, not zero, so no
homogeneous recurrence at i=1 is used.

**Corollary 2 (all-constant three-window).** The universal no-three-forcing
statement extends to every z, including z=0. For z!=0 it is proved in
Corollary5 of [the direct three-window proof](QUADRATIC_FACE_FORCING_THREE_WINDOW_AND_SYMMETRIC_FIVE_WINDOW.md), whose global proof combines the native
cubic valuation alternatives and the three local valuation branches;
at z=0 Theorem1 is stronger.

The exact direct source is
[the direct three-window proof](QUADRATIC_FACE_FORCING_THREE_WINDOW_AND_SYMMETRIC_FIVE_WINDOW.md),
SHA256 2572d6e869692ab8afe2c8b9783209fdf5170a7e24421194c571a8a8c990937a.

## 1. Symmetric zero-constant moments

Let f=a(x^2+y^2)+bxy+e(x+y), with e!=0. Put M_j=Lambda_2(f^j),
H_j=Lambda_1((aX^2+eX)^j). The already proved symmetric recurrence is

 M_(j+1)=(2a+b)(j+1)^2 M_j+j(j+1)e^2 M_(j-1)
                              +H_(j+1)+(j+1)eH_j.              (3)

It is Equation10 of [the symmetric quadratic proof](FULL_SYMMETRIC_BINARY_QUADRATIC_STRONG_FACTORIAL.md),
SHA256 c29933ff9bfee00d55703e049cde90c7ceb1cc577ada03a34dce4b9413530c00.
It is a polynomial identity and remains valid at c=0 and at every vanishing
quadratic pivot.

If M_n,M_(n+1),M_(n+2),M_(n+3) were zero for n>=1, apply(3) at j=n+1,n+2.
It gives the two adjacent face-forcing zeros at indices n+2,n+3. Dividing
by e to the respective powers identifies them with
T_(n+2)(a/e,0),T_(n+3)(a/e,0), contradicting Theorem1. Hence this displayed
symmetric zero-constant family has no four consecutive zero moments.

When b=0 and ae!=0, its actual four monomials x^2,y^2,x,y satisfy their
Strong Factorial window. When abe!=0, its five actual monomials satisfy
the stronger four-zero exclusion. No shorter support window is inferred
merely by setting additional coefficients to zero.

## 2. Exact forcing-window sharpness, without a moment-window claim

Take r=-2/5 and z=-6/5. Direct factorial evaluation gives

 h_0=1, h_1=-1, h_2=2, h_3=-1986/125,
 T_1=0, T_2=0, T_3=-1236/125!=0.                               (4)

Both r and z are nonzero, so the auxiliary quadratic phase has all three
terms. Thus the universal three-consecutive-forcing exclusion cannot be
replaced by a universal two-forcing exclusion, even in that dense face
class. This proves sharpness of the FORCING sequence window only.
It does not construct a symmetric binary phase with four zero moments,
or prove optimality of a factorial-moment SFC window.

All conclusions are algebraic over Q and hold in every characteristic-zero
field. No general nonsymmetric binary quadratic, full SFC(2), or SFC(3)
claim is made. The zero-constant two-window and the explicit forcing
sharpness example are proved locally above; the all-constant three-window
uses the linked nonzero-constant theorem, and the symmetric application
uses only the displayed pinned interior recurrence.
