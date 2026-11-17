"""Check the local manuscript bundle and selected exact Gaussian identities.

Run with the repository's .venv Python. This is a finite algebra and artifact
check, not a verification of the formal-frame or analytic classification proofs.
No existing audit receipts are overwritten. Use --build to compile the current
source and compare every rendered page with the retained PDF. Add --refresh to
install the rebuilt PDF and update its artifact identities in verification.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from importlib.metadata import version
from pathlib import Path

import sympy as S

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_tree():
    """Follow the actual TeX input graph, not historical version directories."""
    pending, sources = [HERE / "manuscript.tex"], {}
    while pending:
        path = pending.pop()
        require(path.parent == HERE, f"nonlocal TeX input: {path}")
        if path.name in sources:
            continue
        require(path.is_file(), f"missing TeX input: {path}")
        body = re.sub(r"(?<!\\)%[^\n]*", "", path.read_text())
        sources[path.name] = body
        for name in re.findall(r"\\input\{([^}]+)\}", body):
            pending.append(HERE / (name if name.endswith(".tex") else name + ".tex"))
    all_text = "\n".join(sources.values())
    labels = re.findall(r"\\label\{([^}]+)\}", all_text)
    require(len(labels) == len(set(labels)), "duplicate TeX label")
    references = set(re.findall(r"\\(?:eqref|ref|pageref)\{([^}]+)\}", all_text))
    require(not references - set(labels), f"missing labels: {references - set(labels)}")
    bibliography = set(re.findall(r"\\bibitem\{([^}]+)\}", all_text))
    citations = set()
    for group in re.findall(r"\\cite(?:\[[^\]]*\])*\{([^}]+)\}", all_text):
        citations.update(key.strip() for key in group.split(","))
    require(not citations - bibliography, f"missing bibliography: {citations - bibliography}")
    return {name: digest(HERE / name) for name in sorted(sources)}


def wick(poly, variables, circular_pairs, real_variables=()):
    """Independent monomial contraction, without a generating-function oracle."""
    value = 0
    for powers, coefficient in S.Poly(S.expand(poly), *variables).terms():
        exponents = dict(zip(variables, powers))
        moment = coefficient
        for z, w in circular_pairs:
            if exponents[z] != exponents[w]:
                moment = 0
                break
            moment *= S.factorial(exponents[z])
        for x in real_variables:
            k = exponents[x]
            if k % 2:
                moment = 0
                break
            moment *= S.factorial2(k - 1) if k else 1
        value += moment
    return S.simplify(value)


def exact_checks():
    z, w, x, z2, w2, t = S.symbols("z w x z2 w2 t")
    p3 = (1 + z) * (w - (2 + z) * x**2 / 2)
    p4 = (1 + z2) * (w * (1 - z) + w2)
    require(S.Poly(p3, z, w, x).total_degree() == 4, "quartic degree")
    require(S.Poly(p4, z, w, z2, w2).total_degree() == 3, "cubic degree")
    for m in range(1, 7):
        for phase, marker, variables, pairs, real in (
            (p3, z, (z, w, x), ((z, w),), (x,)),
            (p4, z2, (z, w, z2, w2), ((z, w), (z2, w2)), ()),
        ):
            require(wick(phase**m, variables, pairs, real) == 0, f"pure witness m={m}")
            require(wick(marker * phase**m, variables, pairs, real) == S.factorial(m),
                    f"mixed witness m={m}")

    # Exact rational identities behind the all-order contractions in Section 5.
    g = t / (1 - t)
    require(S.cancel(1 + t * (1 + g) * (2 + g) - (1 - t)**-2) == 0,
            "quartic radicand identity")
    require(S.cancel(1 / ((1 - t) * (1 + t * (1 + g))) - 1) == 0,
            "four-variable prefactor identity")

    a, c, d, alpha, beta, gamma = S.symbols("a c d alpha beta gamma", nonzero=True)
    u = c + d * (a * t + beta)
    contracted_exponent = (alpha * (a * t + beta) - t * (a * t + beta) * u**2 / (2 * a)
                           + (gamma + t * u)**2 / 2)
    displayed_exponent = alpha * (a * t + beta) + gamma**2 / 2 + gamma * t * u - beta * t * u**2 / (2 * a)
    require(S.expand(contracted_exponent - displayed_exponent) == 0, "shear source identity")
    exceptional = w + z * x - z**3 / 2
    for q in range(1, 5):
        require(wick(w**q * exceptional**(3 * q), (z, w, x), ((z, w),), (x,))
                == S.factorial(3 * q) * (-S.Rational(1, 2))**q, f"sharpness q={q}")

    # The frame test must retain every lower-degree coefficient.
    coordinates = S.symbols("x0:3")
    direction = S.symbols("v0:3")
    cubic = quadratic = 0
    for degree in (2, 3):
        term = 0
        for i in range(degree + 1):
            for j in range(degree - i + 1):
                k = degree - i - j
                term += S.Symbol(f"c{degree}_{i}{j}{k}") * coordinates[0]**i * coordinates[1]**j * coordinates[2]**k
        if degree == 3:
            cubic = term
        else:
            quadratic = term
    directional = lambda f: sum(v * S.diff(f, y) for v, y in zip(direction, coordinates))
    substitution = dict(zip(coordinates, direction))
    rhs = 2 * sum(S.diff(cubic, y).subs(substitution) * y for y in coordinates) + 2 * quadratic.subs(substitution)
    require(S.expand(directional(directional(cubic + quadratic)) - rhs) == 0, "frame polarization identity")

    # Existing symbolic checker is read-only and does not mutate its old receipt.
    from check_formal_identities import check
    check()
    return {
        "pure_and_mixed_witness_powers": [1, 6],
        "sharpness_multiplier_degrees": [1, 4],
        "universal_identities": ["counterexample contraction factors", "shear source exponent",
                                 "general ternary frame polarization", "Airy and mixed-phase eliminations"],
    }


def build(install):
    # Keep auxiliary files out of the bundle; plain --build is read-only there.
    with tempfile.TemporaryDirectory(prefix="gaussian-boundary-") as directory:
        command = ["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error",
                   f"-outdir={directory}", "manuscript.tex"]
        process = subprocess.run(command, cwd=HERE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        log = Path(directory, "manuscript.log")
        log_text = log.read_text(errors="replace") if log.exists() else process.stdout
        require(process.returncode == 0, "LaTeX build failed:\n" + process.stdout[-8000:])
        bad = re.findall(r"(?m)^.*(?:Overfull|undefined|multiply defined|Rerun to get).*$", log_text)
        require(not bad, "LaTeX diagnostics require review:\n" + "\n".join(bad))
        rebuilt = Path(directory, "manuscript.pdf")
        if install:
            shutil.copyfile(rebuilt, HERE / "manuscript.pdf")
        else:
            import pymupdf
            with pymupdf.open(rebuilt) as new, pymupdf.open(HERE / "manuscript.pdf") as old:
                require(len(new) == len(old), "rebuilt page count differs; review and --refresh")
                for index, (new_page, old_page) in enumerate(zip(new, old), 1):
                    require(new_page.get_text() == old_page.get_text(),
                            f"rebuilt text differs on page {index}; review and --refresh")
                    new_image, old_image = new_page.get_pixmap(), old_page.get_pixmap()
                    require((new_image.width, new_image.height, new_image.samples)
                            == (old_image.width, old_image.height, old_image.samples),
                            f"rebuilt rendering differs on page {index}; review and --refresh")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true", help="compile with latexmk and compare rendered pages; add --refresh to install")
    parser.add_argument("--refresh", action="store_true", help="record current source/PDF hashes after checking")
    args = parser.parse_args()
    if args.refresh and not args.build:
        parser.error("--refresh requires --build to compile the recorded source/PDF pair")
    hashes = source_tree()
    checks = exact_checks()
    if args.build:
        build(install=args.refresh)
    import pymupdf
    pdf_path = HERE / "manuscript.pdf"
    with pymupdf.open(pdf_path) as pdf:
        pages = len(pdf)
        text = "\n".join(page.get_text() for page in pdf)
    require("??" not in text, "unresolved references in PDF")
    require("The dimension" in text, "unexpected PDF")
    record_path = HERE / "verification.json"
    record = json.loads(record_path.read_text())
    metadata = json.loads((HERE / "metadata.json").read_text())
    abstract = (HERE / "manuscript.tex").read_text().split("\\begin{abstract}", 1)[1].split("\\end{abstract}", 1)[0].strip()
    require(metadata["abstract"] == abstract, "metadata abstract differs from manuscript")
    if args.refresh:
        record.update(source_hashes=hashes, pdf_sha256=digest(pdf_path), pages=pages)
        record["build_environment"] = {
            "python": sys.version.split()[0],
            "sympy": version("sympy"),
            "pymupdf": version("pymupdf"),
            "latexmk": subprocess.check_output(["latexmk", "-v"], text=True).strip().splitlines()[0],
            "pdflatex": subprocess.check_output(["pdflatex", "--version"], text=True).splitlines()[0],
        }
        record["reproducible_check"] = {
            "command": "../../.venv/bin/python verify.py --build",
            "checker_sha256": digest(Path(__file__)),
            "support_hashes": {"check_formal_identities.py": digest(HERE / "check_formal_identities.py")},
            "scope": "Local source closure, TeX references, clean LaTeX build when requested, artifact identities, and selected exact algebraic identities.",
            "exact_checks": checks,
            "does_not_verify": ["general formal-frame constructions", "analytic estimates",
                                "full ternary-cubic classification", "independent external review"],
        }
        record["verification_note"] = "Source/PDF identities and selected algebra replayed with verify.py; mathematical review scope and remaining review are recorded separately."
        record_path.write_text(json.dumps(record, indent=2) + "\n")
    require(record["source_hashes"] == hashes, "source hashes differ; review changes, rebuild and use --refresh")
    require(record["pdf_sha256"] == digest(pdf_path), "PDF hash differs")
    require(record["pages"] == pages, "PDF page count differs")
    require(record["reproducible_check"]["checker_sha256"] == digest(Path(__file__)), "checker identity differs")
    for name, expected in record["reproducible_check"]["support_hashes"].items():
        require(digest(HERE / name) == expected, f"support checker identity differs: {name}")
    print(json.dumps({"status": "passed", "tex_sources": len(hashes), "pages": pages,
                      "build_replayed": args.build, "selected_exact_checks": checks,
                      "full_mathematical_proof_verified": False}, indent=2))


if __name__ == "__main__":
    main()
