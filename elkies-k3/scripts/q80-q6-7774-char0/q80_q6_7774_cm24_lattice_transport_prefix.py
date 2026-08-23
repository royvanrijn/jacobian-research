#!/usr/bin/env python3
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path


def find_repo() -> Path:
    env = os.environ.get("JACOBIAN_RESEARCH")
    candidates = []
    if env:
        candidates.append(Path(env).expanduser())

    cwd = Path.cwd().resolve()
    candidates.extend([cwd, *cwd.parents])

    home = Path.home()
    candidates.extend([
        home / "jacobian-research",
        home / "GitHub" / "jacobian-research",
        home / "github" / "jacobian-research",
        home / "Developer" / "jacobian-research",
        home / "Development" / "jacobian-research",
        home / "Projects" / "jacobian-research",
        home / "projects" / "jacobian-research",
        home / "code" / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "Documents" / "GitHub" / "jacobian-research",
    ])

    marker = Path("elkies-k3/scripts/classify_kumar_cm_frame_extensions.sage")
    seen = set()
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if (candidate / marker).exists():
            return candidate

    for depth1 in home.iterdir():
        if not depth1.is_dir() or depth1.name.startswith("."):
            continue
        if depth1.name in {"Library", "Movies", "Music", "Pictures"}:
            continue
        direct = depth1 / "jacobian-research"
        if (direct / marker).exists():
            return direct.resolve()
        try:
            for depth2 in depth1.iterdir():
                if not depth2.is_dir() or depth2.name.startswith("."):
                    continue
                if depth2.name == "jacobian-research" and (depth2 / marker).exists():
                    return depth2.resolve()
        except PermissionError:
            pass

    raise SystemExit(
        "Could not locate jacobian-research. "
        "Set JACOBIAN_RESEARCH=/path/to/jacobian-research and rerun."
    )


SAGE_CODE = '\nfrom __future__ import print_function\n\nimport csv\nimport sys\nfrom pathlib import Path\n\nfrom sage.all import (\n    QQ, ZZ, block_diagonal_matrix, gcd,\n    matrix, vector, xgcd\n)\n\nROOT = Path(sys.argv[1]).resolve()\nSCRIPTS = ROOT / "elkies-k3" / "scripts"\nDATA = ROOT / "elkies-k3" / "data" / "fibrations"\n\nCLASSIFY = SCRIPTS / "classify_kumar_cm_frame_extensions.sage"\nassert CLASSIFY.exists()\n\n_saved_argv = list(sys.argv)\n_saved_file = globals().get("__file__", None)\n\n# We only need the deterministic construction prefix through q80_cm24.\n# The repository file also contains later exploratory level-79 searches;\n# those are unrelated to this specialization replay and may carry their own\n# assertions.  Execute the Sage-preparsed prefix only.\nfrom sage.repl.preparse import preparse_file\n\nsource = CLASSIFY.read_text()\ncut_marker = (\n    "assert q80_cm24_profiles[3] == q80_cm24_optimal[1]"\n)\ncut = source.find(cut_marker)\nassert cut >= 0\ncut += len(cut_marker)\nprefix = source[:cut] + "\\n"\n\n# Sanity: the excluded tail contains the later exploratory machinery that\n# caused the previous run to abort.\ntail = source[cut:]\nassert "level79_weyl_orbits" in tail\nassert "q80_cm24, q80_embedding24" in prefix\nassert "q80_child" in prefix\n\nsys.argv = [\n    str(CLASSIFY),\n    "--search-max-q", "0",\n    "--orbit-counts-only",\n]\nglobals()["__file__"] = str(CLASSIFY)\ntry:\n    exec(preparse_file(prefix), globals())\nfinally:\n    sys.argv = _saved_argv\n    if _saved_file is None:\n        globals().pop("__file__", None)\n    else:\n        globals()["__file__"] = _saved_file\n\nassert \'q80_child\' in globals()\nassert \'q80_cm24\' in globals()\nassert \'q80_embedding24\' in globals()\nassert \'root_invariants\' in globals()\nassert \'root_components\' in globals()\n\nassert (ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_2.txt").exists()\nprint(\n    f"Q807774CMTRANSPORT|classify_source={CLASSIFY}|"\n    f"prefix_chars={len(prefix)}|tail_skipped_chars={len(tail)}|"\n    "status=PASS_REPOSITORY_PREFIX_LOAD",\n    flush=True,\n)\n\nHYPER = matrix(ZZ, [[0, 1], [1, 0]])\n\n\ndef load_matrix(path):\n    return matrix(\n        ZZ,\n        [\n            [ZZ(value) for value in line.split()]\n            for line in path.read_text().splitlines()\n            if line.strip() and not line.startswith("#")\n        ],\n    )\n\n\ndef bezout_vector(pairings):\n    current = ZZ(0)\n    coeffs = [ZZ(0)] * len(pairings)\n    for i, pairing in enumerate(pairings):\n        pairing = ZZ(pairing)\n        if pairing == 0:\n            continue\n        g, left, right = xgcd(current, pairing)\n        coeffs = [left*c for c in coeffs]\n        coeffs[i] += right\n        current = g\n    assert abs(current) == 1\n    if current == -1:\n        coeffs = [-c for c in coeffs]\n    return vector(ZZ, coeffs)\n\n\ndef neighbor(parent, qnorm, a, b, coordinates):\n    ns = block_diagonal_matrix(HYPER, -parent)\n    coordinates = vector(ZZ, coordinates)\n    fiber = vector(ZZ, [a, b] + list(coordinates))\n\n    assert ZZ(a)*ZZ(b) == ZZ(qnorm)\n    assert coordinates * parent * coordinates == 2*ZZ(qnorm)\n    assert fiber * ns * fiber == 0\n    assert gcd([abs(ZZ(x)) for x in ns*fiber]) == 1\n\n    mate = bezout_vector(list(ns*fiber))\n    mate -= ZZ(mate*ns*mate)//2 * fiber\n\n    assert fiber*ns*mate == 1\n    assert mate*ns*mate == 0\n\n    complement = matrix(\n        ZZ, [list(fiber*ns), list(mate*ns)]\n    ).right_kernel_matrix()\n\n    child = -(complement*ns*complement.transpose())\n    transport = matrix(\n        ZZ,\n        [list(fiber), list(mate)] + [list(row) for row in complement.rows()],\n    )\n    assert abs(transport.det()) == 1\n    return child, transport\n\n\ndef enhance_neighbor(transport, embedding, special_frame):\n    lifted = [\n        vector(\n            ZZ,\n            list(row[:2]) + list(vector(ZZ, row[2:]) * embedding),\n        )\n        for row in transport.rows()\n    ]\n\n    special_ns = block_diagonal_matrix(HYPER, -special_frame)\n    fiber, mate = lifted[:2]\n\n    assert fiber * special_ns * fiber == 0\n    assert mate * special_ns * mate == 0\n    assert fiber * special_ns * mate == 1\n\n    complement = matrix(\n        ZZ,\n        [list(fiber*special_ns), list(mate*special_ns)],\n    ).right_kernel_matrix()\n\n    basis = matrix(\n        ZZ,\n        [list(fiber), list(mate)] + [list(row) for row in complement.rows()],\n    )\n    assert abs(basis.det()) == 1\n\n    child = -(complement*special_ns*complement.transpose())\n    inverse_basis = basis.inverse()\n\n    embedding_rows = []\n    for row in lifted[2:]:\n        coordinates = row * inverse_basis\n        assert coordinates[0] == coordinates[1] == 0\n        embedding_rows.append(list(coordinates[2:]))\n\n    child_embedding = matrix(ZZ, embedding_rows)\n    return child, child_embedding\n\n\ndef ade_name(invariant):\n    rank, roots, det = map(int, invariant)\n    if rank == 0:\n        return "rootless"\n    if roots == rank*(rank+1) and det == rank+1:\n        return f"A{rank}"\n    if rank >= 4 and roots == 2*rank*(rank-1) and det == 4:\n        return f"D{rank}"\n    if (rank, roots, det) == (6, 72, 3):\n        return "E6"\n    if (rank, roots, det) == (7, 126, 2):\n        return "E7"\n    if (rank, roots, det) == (8, 240, 1):\n        return "E8"\n    raise ArithmeticError(f"unknown ADE invariant {invariant}")\n\n\ndef special_signature(frame):\n    components = root_components(frame)\n    names = [ade_name(component) for component in components]\n    names.sort()\n\n    root_rank = sum(int(component[0]) for component in components)\n    root_count = sum(int(component[1]) for component in components)\n    root_det = ZZ(1)\n    for component in components:\n        root_det *= ZZ(component[2])\n\n    assert frame.nrows() == 18\n    mw_rank = 18-root_rank\n\n    return {\n        "ade": "+".join(names) if names else "rootless",\n        "root_rank": root_rank,\n        "root_count": root_count,\n        "root_det": int(root_det),\n        "mw_rank": mw_rank,\n        "components": tuple(tuple(map(int, component)) for component in components),\n    }\n\n\ngeneric = load_matrix(DATA / "kumar_q80_e6_d5_a3_mw3_frame.txt")\nassert generic == q80_child\n\nspecial = q80_cm24\nembedding = q80_embedding24\nassert generic == embedding * special * embedding.transpose()\n\nprint(\n    f"Q807774CMTRANSPORT|start_generic_det={generic.det()}|"\n    f"start_special_det={special.det()}|"\n    f"embedding_shape={embedding.nrows()}x{embedding.ncols()}|"\n    "status=PASS_PINNED_CM24_EXTENSION",\n    flush=True,\n)\n\nwith (DATA / "kumar_q80_to_rootless_path.tsv").open() as handle:\n    canonical = list(csv.DictReader(handle, delimiter="\\t"))\nassert len(canonical) >= 2\n\nsteps = []\nfor index in (0, 1):\n    row = canonical[index]\n    steps.append({\n        "name": f"prefix_q4_{index+1}",\n        "q": ZZ(row["q"]),\n        "a": ZZ(row["a"]),\n        "b": ZZ(row["b"]),\n        "v": tuple(map(ZZ, row["v"].split(","))),\n    })\n\nsteps.extend([\n    {\n        "name": "escape",\n        "q": ZZ(6), "a": ZZ(2), "b": ZZ(3),\n        "v": (-5,-3,6,6,-8,-4,2,4,-1,8,-16,-1,0,3,5,-2,-2),\n    },\n    {\n        "name": "orbit424",\n        "q": ZZ(4), "a": ZZ(2), "b": ZZ(2),\n        "v": (32,48,-21,28,8,-52,-34,0,18,5,-23,43,9,-18,16,-6,-6),\n    },\n    {\n        "name": "orbit1222",\n        "q": ZZ(4), "a": ZZ(2), "b": ZZ(2),\n        "v": (10,53,-192,-114,29,-256,-170,-12,-14,74,-32,-14,-6,-26,-58,84,-28),\n    },\n    {\n        "name": "q6_7774",\n        "q": ZZ(6), "a": ZZ(2), "b": ZZ(3),\n        "v": (85,2699,1257,7718,3756,-41,3077,-4614,-6615,6032,2584,-1678,121,-736,-913,1,1165),\n    },\n])\n\norbit1222_signature = None\nq67774_signature = None\n\nfor step_index, step in enumerate(steps, 1):\n    child, transport = neighbor(\n        generic,\n        step["q"], step["a"], step["b"],\n        vector(ZZ, step["v"]),\n    )\n\n    special_child, child_embedding = enhance_neighbor(\n        transport, embedding, special\n    )\n\n    assert child == (\n        child_embedding\n        * special_child\n        * child_embedding.transpose()\n    )\n\n    sig = special_signature(special_child)\n\n    print(\n        f"Q807774CMTRANSPORT|step={step_index}|name={step[\'name\']}|"\n        f"special_ade={sig[\'ade\']}|"\n        f"root_data=({sig[\'root_rank\']},{sig[\'root_count\']},{sig[\'root_det\']})|"\n        f"MW={sig[\'mw_rank\']}|"\n        f"components={sig[\'components\']}|"\n        "status=PASS_SPECIAL_STEP",\n        flush=True,\n    )\n\n    if step["name"] == "orbit1222":\n        orbit1222_signature = sig\n    if step["name"] == "q6_7774":\n        q67774_signature = sig\n\n    generic = child\n    special = special_child\n    embedding = child_embedding\n\nassert orbit1222_signature is not None\nassert (\n    orbit1222_signature["root_rank"],\n    orbit1222_signature["root_count"],\n    orbit1222_signature["root_det"],\n    orbit1222_signature["mw_rank"],\n) == (15, 90, 392, 3)\n\nprint(\n    "Q807774CMTRANSPORT|orbit1222_regression="\n    "2A6+3A1/MW3|status=PASS_ORBIT1222_CM24_REGRESSION",\n    flush=True,\n)\n\nassert q67774_signature is not None\npredicted = (\n    q67774_signature["root_rank"],\n    q67774_signature["root_count"],\n    q67774_signature["root_det"],\n)\n\nequation_candidates = {\n    (16,190,84): ("pair1", (-5,4,0), "A13+A2+A1", 2),\n    (15,170,32): ("pair3", (-1,-4,-1), "D7+D7+A1", 3),\n    (14,130,55): ("pair4", (-1,-4,0), "A10+A4", 4),\n}\n\nprint(\n    f"Q807774CMTRANSPORT|predicted_root_data={predicted}|"\n    f"predicted_ade={q67774_signature[\'ade\']}|"\n    f"predicted_MW={q67774_signature[\'mw_rank\']}|"\n    "status=PASS_LATTICE_SPECIALIZATION_PREDICTION",\n    flush=True,\n)\n\nassert predicted in equation_candidates, (\n    "CM24 lattice prediction did not match any equation candidate",\n    predicted,\n    equation_candidates,\n)\n\ncandidate, mw_coordinates, expected_ade, expected_mw = equation_candidates[predicted]\nassert q67774_signature["mw_rank"] == expected_mw\n\nprint(\n    f"Q807774CMSELECT|candidate={candidate}|"\n    f"mw={mw_coordinates}|"\n    f"root_data={predicted}|"\n    f"special_ade={expected_ade}|MW={expected_mw}|"\n    "status=PASS_UNIQUE_Q6_7774_CM24_MARKING",\n    flush=True,\n)\n\nprint(\n    "Q807774CMTRANSPORT|next=pin_selected_q6_7774_CM24_child_and_"\n    "continue_q4_1938_equation_level|"\n    "status=PASS_Q6_7774_CM24_LATTICE_EQUATION_MATCH",\n    flush=True,\n)\n'


def main():
    repo = find_repo()
    sage = shutil.which("sage") or "/usr/local/bin/sage"
    if shutil.which("sage") is None and not Path(sage).exists():
        raise SystemExit("sage not found")

    print(f"repo={repo}", flush=True)
    print(f"sage={sage}", flush=True)

    with tempfile.TemporaryDirectory() as td:
        script = Path(td) / "q80_q6_7774_cm24_lattice_transport.sage"
        script.write_text(SAGE_CODE)
        subprocess.run([sage, str(script), str(repo)], check=True)


if __name__ == "__main__":
    main()
