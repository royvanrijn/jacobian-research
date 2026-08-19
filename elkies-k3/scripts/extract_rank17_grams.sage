from sage.all import *
from pathlib import Path
import re, sys

POOL = Path(sys.argv[1])
AMBIENT = Path(sys.argv[2])
OUT = Path(sys.argv[3])

def ambient_gram(path):
    lines = path.read_text().splitlines()
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
    m = re.search(r'basis_coords\s*=\s*\n(\[\[.*?\]\])', text, re.S)
    if not m:
        return None
    rows = []
    for line in m.group(1).splitlines():
        x = re.findall(r'-?\d+', line)
        if len(x) == 24:
            rows.append([ZZ(v) for v in x])
    return matrix(ZZ, rows) if len(rows) == 7 else None

G = ambient_gram(AMBIENT)
OUT.mkdir(parents=True, exist_ok=True)

for p in sorted(POOL.glob("root-*.txt")):
    B = basis(p)
    if B is None:
        print("PARSE FAIL", p)
        continue

    C = matrix(ZZ, (B * G).right_kernel().basis_matrix())
    H = C * G * C.T

    assert H.rank() == 17
    assert H.det() == 948
    assert all(H[i,i] % 2 == 0 for i in range(17))

    out = OUT / (p.stem + "-gram.txt")
    out.write_text(
        "\n".join(" ".join(map(str,row)) for row in H.rows()) + "\n"
    )

    print(p.name, "->", out)
