# Fully shared terminal ten from outer mean .52 to .76

Let two distinct ordinary exponent columns have all unequal rows shared.
Write \(u\) for their normalized order-nine outer reflected mean.
If \(.52\le u\le .76\) and all reciprocals \(d/(9\ell+1)\) are at most
\(1/8\), the order-nine and order-ten factorial moment polynomials have
disjoint root arguments under independent positive rescalings.
In particular they are coprime over every characteristic-zero field.

## Exact retained enumeration

In each band let \(c\) be the cap of the
[shared light-row proof](MULTIVARIATE_TERMINAL_TEN_SHARED_LIGHT.md).
If all reciprocals are at most \(c\), that proof applies. Otherwise
retain every row with reciprocal greater than \(c\). The individual
squared order-nine outer mean \(u_{\ell,d}^2\) obeys
\[
 u_{\ell,d}^2\le(1+8c/(1+c))^{-d}
 \quad\text{when}\quad 9\ell+1<d/c .
\]
Each individual row and each product of retained rows has mean
at least \(u_-\), the band's lower endpoint. These inequalities
bound the gaps and bases in the candidate list, after restricting
reciprocals to at most \(1/8\). The checker computes each candidate's
exact squared mean and exhausts its multiplicities by multiplying
the squared means until the product falls below \(u_-^2\).
It identifies each signed list with its global column reversal.

The following lists remain. The pair notation denotes the two
exponents of a row, and the orientation relations are up to global
column reversal.

| Band | Light cap \(c\) | Nonempty retained lists |
| --- | ---: | --- |
| \([.52,.62]\) | \(99/1000\) | One \((1,2)\), two \((1,2)\) in the same orientation, two \((1,2)\) in opposite orientations, or one \((2,4)\) |
| \([.62,.70]\) | \(3/40\) | One \((1,2)\) |
| \([.70,.76]\) | \(11/200\) | One \((1,2)\) |

## Fixed core and residual replay

Let \(v_C\) be the retained order-nine outer mean and \(v=u/v_C\)
the residual mean. Its range is contained in
\([u_-/v_C,\min(1,u_+/v_C)]\). The
[exact checker](verify_multivariate_terminal_ten_retained_mid_bands.py)
brackets these endpoints outward on a rational grid by comparing
squares against the exact \(v_C^2\). It tests coverage of this
bracket by the frozen ray table of each retained list. All
residual reciprocals are at most \(c\). Within each residual
slab the checker implements the ordinary row cap calculation
from the light proof. For the singleton \((1,2)\), the first
effective residual caps in the three bands are \(2/37,1/19,1/40\),
respectively. All three other lists have effective cap at most
\(1/40\) throughout.

For each order \(j\), normalize the retained multiplier \(F_{j,k}\)
at coefficient \(k\) to one at both endpoints. Its \(j\)-th power
is a positive rational product of the retained row factorials.
Integer comparisons enclose \(F_{j,k}\) on a dyadic grid.
Multiply each coefficient vector of the
[matched moving-ray model](MULTIVARIATE_TERMINAL_TEN_SHARED_LIGHT.md#grouped-moment-and-moving-ray-comparison)
by the same fixed factor \(F_{j,k}\) in degree \(k\).
In the base and shared-displacement groups this preserves the
common affine parameter exactly. The retained core fixes the
column orientation, so the residual signed parameter ranges over
all of \([-1,1]\). The singleton uses the grouped cubic exponential
model. The other lists use the quadratic model.

The independent residual coefficient caps verify the annulus
\(1/12<|z|<12\), including the exact retained factors. At fixed
retained list the logarithmic homotopy scales the residual skew
and defects to zero, then moves the outer ratio to \(81/80\).
The resulting Gaussian-with-core references move to the lower
endpoint of their respective residual bands. The ray cosines
remain ordered and move affinely between consecutive knots,
with strict Bernstein projection signs along every verified
cell. At residual mean one, the actual and Gaussian-with-core
arrays both equal the exact retained core.

Each lower endpoint has independent exact rational Rouché disks.
They give alternating upper sector orders \(10,9,10,9,10,9,10,9,10\)
and one extra negative real root of order nine. The last
sector joined to its reflection has only this one nonic root,
which stays real. Strict positivity of the actual tenth moment
on the real axis follows from its Gamma integral. Hence the
tenth-order pair in that sector is nonreal, proving the claimed
disjointness of arguments.

## Reproduction

Run the retained checker from the repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_retained_mid_bands.py --native
~~~

The optional --native flag uses the same outward dyadic
C++17 projection routine as the light checker. It enumerates
four, one, and one retained lists in the three bands and replays
each list. In the first band the three non-singleton lists use
797, 782, and 694 radial cells in table order. The marked
singleton uses 4,257, 3,112, and 2,127 cells in its three bands.
