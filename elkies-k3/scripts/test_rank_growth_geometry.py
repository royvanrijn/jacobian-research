from __future__ import annotations

import math
import numpy as np

from rank_growth import cascade_metrics, extension_metrics, matrix_numerical_rank


# Construct an exact Euclidean toy model in R^21.
# e0..e16 are the old rank-17 core.  Each next vector contains a small new
# orthogonal increment, and later vectors couple strongly to the newest one.
e = np.eye(21)
rows = [e[i].copy() for i in range(17)]

q18 = 2.0 * e[0] + 0.20 * e[17]
q19 = -1.5 * e[1] + 0.95 * e[17] + 0.10 * e[18]
q20 = 0.75 * e[2] - 0.20 * e[17] + 0.90 * e[18] + 0.08 * e[19]
q21 = -0.50 * e[3] + 0.10 * e[17] - 0.15 * e[18] + 0.95 * e[19] + 0.06 * e[20]
rows.extend([q18, q19, q20, q21])
A = np.asarray(rows)
G = A @ A.T

assert matrix_numerical_rank(G) == 21

m18 = extension_metrics(G[:18, :18], 17, 17)
assert math.isclose(m18.orthogonal_height, 0.20**2, rel_tol=1e-10, abs_tol=1e-12)

m19 = extension_metrics(G[:19, :19], 18, 18)
c19 = cascade_metrics(G[:19, :19], 17, 18, 18)
assert m19.orthogonal_height > 0
assert math.isclose(c19.last_increment_corr, 1.0, rel_tol=1e-12)

m20 = extension_metrics(G[:20, :20], 19, 19)
c20 = cascade_metrics(G[:20, :20], 17, 19, 19)
assert m20.orthogonal_height > 0
assert c20.last_increment_corr > 0.95

m21 = extension_metrics(G, 20, 20)
c21 = cascade_metrics(G, 17, 20, 20)
assert m21.orthogonal_height > 0
assert c21.last_increment_corr > 0.95

print("PASS")
print(f"delta18={m18.orthogonal_height:.12g}")
print(f"corr19={c19.last_increment_corr:.12g}")
print(f"corr20={c20.last_increment_corr:.12g}")
print(f"corr21={c21.last_increment_corr:.12g}")
