from sage.all import *
from pathlib import Path
import re, sys

BASE = Path(__file__).resolve().parents[1]
AMBIENT = BASE / "niemeier-all/24A1/ambient.txt"
OUT = BASE / "seeds/24A1-target-genus-pneighbor-grams"

# Prefer newest target-genus frontier first.
DIRS = [
    BASE / "seeds/24A1-target-genus-frontier-g5",
    BASE / "seeds/24A1-target-genus-frontier-g4",
    BASE / "seeds/24A1-target-genus-frontier-g3",
    BASE / "seeds/24A1-target-genus-frontier-g2",
    BASE / "seeds/24A1-target-genus-r1-n1311",
]

wanted = {1, 6, 11, 16, 21, 23}


def ambient_gram():
    lines = AMBIENT.read_text().splitlines()
    i = lines.index("GRAM 24")

    rows = []
    for line in lines[i+1:]:
        x = re.findall(r'-?\d+', line)
        if len(x) == 24:
            rows.append([ZZ(v) for v in x])
        if len(rows) == 24:
            break

    return matrix(ZZ, rows)


def basis(path):
    text = path.read_text(errors="replace")

    m = re.search(
        r'basis_coords\s*=\s*\n(\[\[.*?\]\])',
        text, re.S
    )
    if not m:
        return None

    rows = []
    for line in m.group(1).splitlines():
        x = re.findall(r'-?\d+', line)
        if len(x) == 24:
            rows.append([ZZ(v) for v in x])

    return matrix(ZZ, rows) if len(rows) == 7 else None


def root_id(path):
    text = path.read_text(errors="replace")
    m = re.search(r'^rootmask = (\d+)$', text, re.M)
    if not m:
        return None

    x = int(m.group(1))

    if x == 0 or x & (x-1):
        return None

    return x.bit_length() - 1


G = ambient_gram()
OUT.mkdir(parents=True, exist_ok=True)

found = set()

for d in DIRS:
    if not d.exists():
        continue

    for p in sorted(d.glob("*.txt")):
        r = root_id(p)

        if r not in wanted or r in found:
            continue

        B = basis(p)
        if B is None:
            continue

        C = matrix(ZZ, (B*G).right_kernel().basis_matrix())
        H = C * G * C.T

        assert H.nrows() == 17
        assert H.det() == 948
        assert all(H[i,i] % 2 == 0 for i in range(17))

        out = OUT / f"root-{r}-gram.txt"

        out.write_text(
            "\n".join(
                " ".join(map(str,row))
                for row in H.rows()
            ) + "\n"
        )

        print(r, p, "->", out)
        found.add(r)

print("FOUND:", sorted(found))
print("MISSING:", sorted(wanted-found))
