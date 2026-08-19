from pathlib import Path
import argparse
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument(
    "gram",
    help="Gram matrix with first 17 vectors the known rank17 basis"
)

args = parser.parse_args()

S = np.loadtxt(args.gram, dtype=float)

assert S.shape[0] == S.shape[1]
assert S.shape[0] >= 18

H = S[:17,:17]

Hi = np.linalg.inv(H)

print("ambient rank =", len(S))
print("extra vectors =", len(S)-17)
print()

for j in range(17,len(S)):

    b = S[:17,j]

    norm = S[j,j]

    coeff = Hi @ b

    projection_norm = float(b @ coeff)

    delta = float(norm - projection_norm)

    nearest = np.rint(coeff)

    coord_error = float(
        np.linalg.norm(coeff-nearest)
    )

    print(
        f"EXT|j={j}"
        f"|norm={norm:.15g}"
        f"|projection={projection_norm:.15g}"
        f"|delta={delta:.15g}"
        f"|coord_err={coord_error:.8g}"
        f"|coeff={coeff.tolist()}"
    )
