#!/usr/bin/env sage -python
"""Enumerate section-translation orbits of lattice bisections in rootless MW17.

For the rootless fibration the Neron--Severi lattice is ``U + (-M)``.  With
``F=(1,0,0)`` and ``O=(-1,1,0)``, a degree-two (-2)-class has the form

    B_w=((w.M.w-2)/4, 2, w),       w.M.w == 2 (mod 4).

Translation by the section indexed by ``x in M`` sends ``w`` to ``w+2x``.
Thus its translation orbit is a class of ``M/2M``.  Moreover

    B_w . S_x = (w-2x).M.(w-2x)/4 - 5/2.

Since the fibration is rootless, the class is nonnegative on every section
exactly when its coset has no norm-2 or norm-6 representative.  This script
enumerates all such cosets, proves that each has minimum norm 10, and emits
one norm-10 representative per orbit.  On any K3 realization of this
rootless fibration, the section-nonnegative condition also forces every such
effective class to be an irreducible smooth rational bisection.  It does not
construct their equations, branch divisors, quadratic extensions, or
Mordell--Weil sections after base change.
"""

from sage.all import QQ, ZZ, QuadraticForm, matrix, pari, vector

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RANK17 = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
SHORT_COORDS = ROOT / "elkies-k3/data/lattice/short_vector_basis_coords.txt"
SHORT_GRAM = ROOT / "elkies-k3/data/lattice/short_vector_basis_gram.txt"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-rootless-bisection-orbits.json"
DEFAULT_ORBITS = ROOT / "artifacts/generated-results/elkies-k3-rootless-bisection-orbits.tsv"
ALTERNATE_FRAME = (
    ROOT / "artifacts/generated-results/"
    "q80-alternate-fifth-q6-rootless-transport.json"
)
ALTERNATE_OUTPUT = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-q80-alternate-rootless-bisection-orbits.json"
)
ALTERNATE_ORBITS = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-q80-alternate-rootless-bisection-orbits.tsv"
)


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mask(vector_value):
    result = 0
    for index, value in enumerate(vector_value):
        if ZZ(value) % 2:
            result |= 1 << index
    return result


def quadratic_form(gram):
    coefficients = []
    for row in range(gram.nrows()):
        for column in range(row, gram.ncols()):
            coefficients.append(
                gram[row, row] // 2 if row == column else gram[row, column]
            )
    return QuadraticForm(ZZ, gram.nrows(), coefficients)


def norm(value, gram):
    return ZZ(value * gram * value)


def entries(value):
    return " ".join(str(ZZ(entry)) for entry in value)


def display_path(path):
    """Return a repository-relative path for either relative or absolute CLI input."""
    return str(Path(path).resolve().relative_to(ROOT))


def integer_quadratic_norm(value, gram):
    """Evaluate a row-vector norm without creating a Sage vector at every leaf."""
    return sum(
        value[left]*gram[left, right]*value[right]
        for left in range(len(value))
        for right in range(len(value))
    )


def streaming_short_vectors(gram, *, bound, representative_key=None):
    """Stream every vector of integral norm at most ``bound``.

    PARI's ``qfminim`` materializes all short vectors.  The alternate q80
    rootless frame has nearly three million signed vectors through norm ten,
    which exceeds its vector-list stack despite being small enough to visit.
    This reduced-basis Fincke--Pohst traversal stores only the parity masks and
    canonical representatives needed by the bisection quotient.  Every leaf
    is rechecked with the integral Gram form.  A quarter-unit floating guard
    is harmless because the lattice is even; the signed shell counts are then
    compared with PARI's exact count-only qfminim output below.
    """
    from math import ceil, sqrt
    from collections import Counter

    dimension = gram.nrows()
    lower = matrix(QQ, dimension, dimension)
    diagonal = []
    for index in range(dimension):
        lower[index, index] = 1
        value = QQ(gram[index, index])-sum(
            lower[index, prior]**2*diagonal[prior]
            for prior in range(index)
        )
        assert value > 0
        diagonal.append(value)
        for row in range(index+1, dimension):
            lower[row, index] = (
                QQ(gram[row, index])-sum(
                    lower[row, prior]*lower[index, prior]*diagonal[prior]
                    for prior in range(index)
                )
            )/value

    lower_float = [
        [float(lower[row, column]) for column in range(dimension)]
        for row in range(dimension)
    ]
    diagonal_float = [float(value) for value in diagonal]
    coordinates = [0]*dimension
    signed_counts = Counter()
    masks_by_norm = {value: set() for value in range(2, bound+1, 2)}
    representatives = {}
    representative_keys = {}
    multiplicities = Counter()
    # The exact norm is integral and even.  Thus the 1/4 guard cannot include
    # a norm above ``bound``; it only makes floating branch pruning generous.
    floating_bound = float(bound)+0.25
    epsilon = 1.0e-10

    def visit(index, used, exact_partial):
        if index < 0:
            exact_value = exact_partial
            if not exact_value or exact_value > bound:
                return
            assert exact_value % 2 == 0
            signed_counts[exact_value] += 1
            orbit = mask(coordinates)
            masks_by_norm[exact_value].add(orbit)
            if exact_value == bound:
                multiplicities[orbit] += 1
                candidate = tuple(coordinates)
                candidate_key = (
                    candidate if representative_key is None
                    else representative_key(candidate)
                )
                prior_key = representative_keys.get(orbit)
                if prior_key is None or candidate_key < prior_key:
                    representatives[orbit] = candidate
                    representative_keys[orbit] = candidate_key
            return

        center = sum(
            lower_float[row][index]*coordinates[row]
            for row in range(index+1, dimension)
        )
        radius = sqrt(max(0.0, (floating_bound-used)/diagonal_float[index]))
        lower_bound = ceil(-center-radius-epsilon)
        upper_bound = int(-(-(-center+radius+epsilon)//1))
        cross = sum(
            gram[index, row]*coordinates[row]
            for row in range(index+1, dimension)
        )
        for entry in range(lower_bound, upper_bound+1):
            cost = diagonal_float[index]*(entry+center)**2
            if used+cost > floating_bound+epsilon:
                continue
            coordinates[index] = entry
            visit(
                index-1,
                used+cost,
                exact_partial+gram[index, index]*entry**2+2*entry*cross,
            )

    visit(dimension-1, 0.0, ZZ.zero())
    expected_signed_count = ZZ(pari(gram).qfminim(bound, 1)[0])
    assert sum(signed_counts.values()) == expected_signed_count
    assert all(value % 2 == 0 for value in signed_counts)
    assert all(count % 2 == 0 for count in multiplicities.values())
    return {
        "signed_counts": signed_counts,
        "masks_by_norm": masks_by_norm,
        "representatives": representatives,
        "unoriented_multiplicities": {
            orbit: count//2 for orbit, count in multiplicities.items()
        },
        "pari_signed_count": expected_signed_count,
    }


def alternate_rootless_main(arguments):
    """Enumerate the second rootless q80 lattice without materializing its shell."""
    frame_data = json.loads(arguments.frame_artifact.read_text())
    frame = matrix(ZZ, frame_data["rootless_frame"])
    assert frame.nrows() == frame.ncols() == 17
    assert frame.is_positive_definite() and frame.det() == 948
    assert pari(frame).qfminim(2)[0] == 0

    # ``LLL_gram`` acts on columns, while all lattice coordinates in this
    # script are rows.  Its transpose is therefore the row-coordinate change.
    change = frame.LLL_gram().transpose()
    gram = change*frame*change.transpose()
    assert abs(change.det()) == 1
    assert pari(gram).qfminim(2)[0] == 0
    streaming = streaming_short_vectors(gram, bound=10)
    masks_by_norm = streaming["masks_by_norm"]
    multiplicities = streaming["unoriented_multiplicities"]
    representatives = streaming["representatives"]
    dimension = gram.nrows()
    parity_vectors = {
        value: vector(ZZ, [(value >> index) & 1 for index in range(dimension)])
        for value in range(1 << dimension)
    }
    residue_two = {
        value for value, representative in parity_vectors.items()
        if norm(representative, gram) % 4 == 2
    }
    excluded = masks_by_norm[2] | masks_by_norm[6]
    candidates = residue_two-excluded
    assert not masks_by_norm[2]
    assert candidates <= set(representatives)
    assert all(norm(vector(ZZ, representatives[orbit]), gram) == 10 for orbit in candidates)

    arguments.orbits_output.parent.mkdir(parents=True, exist_ok=True)
    with arguments.orbits_output.open("w") as stream:
        stream.write(
            "orbit_mask\thex\tmin_norm\tminimal_unoriented_count\t"
            "short_basis_w\talternate_rank17_w\tlattice_bisection\n"
        )
        for orbit in sorted(candidates):
            short_vector = vector(ZZ, representatives[orbit])
            original_vector = short_vector*change
            assert norm(original_vector, frame) == 10
            stream.write(
                f"{orbit}\t0x{orbit:05x}\t10\t{multiplicities[orbit]}\t"
                f"{entries(short_vector)}\t{entries(original_vector)}\t"
                f"2 2 {entries(original_vector)}\n"
            )

    multiplicity_histogram = {
        str(count): frequency
        for count, frequency in sorted(Counter(multiplicities[orbit] for orbit in candidates).items())
    }
    payload = {
        "schema": "elkies-k3.rootless-bisection-orbits.v1",
        "status": "PASS_ALTERNATE_ROOTLESS_LATTICE_BISECTION_ORBITS",
        "scope": (
            "Exact alternate-rootless Neron--Severi lattice enumeration only; no "
            "characteristic-zero bisection equation, branch divisor, quadratic-extension "
            "hash, collision, or base-changed height matrix is asserted."
        ),
        "input": {
            "alternate_rootless_transport": display_path(arguments.frame_artifact),
            "alternate_rootless_transport_sha256": digest(arguments.frame_artifact),
            "rootless_frame_sha256": frame_data["rootless_frame_sha256"],
            "rank": 17,
            "determinant": 948,
            "row_short_basis_change": [list(map(int, row)) for row in change.rows()],
        },
        "streaming_enumeration": {
            "method": "LLL-reduced Fincke-Pohst traversal with exact leaf norms",
            "bound": 10,
            "floating_guard_above_integral_bound": "1/4",
            "pari_exact_signed_count_through_bound": int(streaming["pari_signed_count"]),
            "signed_shell_counts": {
                str(value): int(streaming["signed_counts"].get(value, 0))
                for value in range(2, 11, 2)
            },
            "parity_cosets_hit_by_shell": {
                str(value): len(masks_by_norm[value])
                for value in range(2, 11, 2)
            },
        },
        "translation_action": {
            "class": "B_w=((w.M.w-2)/4,2,w)",
            "section_translation": "w -> w + 2x",
            "quotient": "M/2M",
            "fibrewise_inversion": {
                "action_on_bisection_vector": "w -> -w",
                "action_on_translation_orbits": "identity",
                "proof": "-w-w=-2w belongs to 2M",
                "collision_search_consequence": (
                    "Fibrewise inversion creates no second translation orbit; "
                    "it cannot be counted as an independent equal-extension collision."
                ),
            },
            "total_cosets": 1 << dimension,
            "residue_two_mod_four_cosets": len(residue_two),
            "rootless_fibre_component_orbits": "trivial",
        },
        "section_nonnegative_filter": {
            "excluded_norms": [2, 6],
            "excluded_cosets": len(excluded),
            "surviving_translation_orbits": len(candidates),
            "minimum_norm_of_every_survivor": 10,
            "minimum_section_intersection": "0",
        },
        "minimal_representatives": {
            "norm": 10,
            "lattice_bisection": [2, 2, "w"],
            "bisection_square": -2,
            "fibre_degree": 2,
            "intersection_with_zero_section": 0,
            "unoriented_norm_ten_pairs_in_survivors": sum(multiplicities[orbit] for orbit in candidates),
            "multiplicity_histogram": multiplicity_histogram,
        },
        "orbits_tsv": display_path(arguments.orbits_output),
        "orbits_tsv_sha256": digest(arguments.orbits_output),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
    print(
        "ALTROOTLESSBISECT|translation_orbits={}|minimum_norm=10|"
        "minimal_pairs={}|fiber_components=trivial|status=PASS_ALTERNATE_ROOTLESS_LATTICE_BISECTION_ORBITS".format(
            len(candidates), sum(multiplicities[orbit] for orbit in candidates)
        ), flush=True,
    )
    print(f"ALTROOTLESSBISECT|output={arguments.output}|sha256={digest(arguments.output)}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--orbits-output", type=Path, default=DEFAULT_ORBITS)
    parser.add_argument(
        "--frame-artifact",
        type=Path,
        help="enumerate the alternate rootless frame exported by this transport artifact",
    )
    arguments = parser.parse_args()

    if arguments.frame_artifact:
        if arguments.output == DEFAULT_OUTPUT:
            arguments.output = ALTERNATE_OUTPUT
        if arguments.orbits_output == DEFAULT_ORBITS:
            arguments.orbits_output = ALTERNATE_ORBITS
        alternate_rootless_main(arguments)
        return

    # The complete norm-10 shell has more than one million unoriented vectors.
    # PARI's default Sage stack is insufficient for this exact enumeration.
    pari.allocatemem(4 * 1024**3)

    # The short basis is an integral change of the pinned rootless frame.  It
    # makes the norm-10 shell enumeration practical while representatives are
    # also exported back in pinned rank17_gram coordinates.
    pinned = load_matrix(RANK17)
    change = load_matrix(SHORT_COORDS)
    gram = load_matrix(SHORT_GRAM)
    assert pinned.nrows() == 17 and pinned.det() == 948
    assert pinned.is_positive_definite()
    assert abs(change.det()) == 1
    assert gram == change * pinned * change.transpose()
    assert all(pinned[index, index] % 2 == 0 for index in range(17))
    assert pari(gram).qfminim(2)[0] == 0, "the input must be rootless"

    form = quadratic_form(gram)
    dimension = gram.nrows()
    parity_vectors = {
        value: vector(ZZ, [(value >> index) & 1 for index in range(dimension)])
        for value in range(1 << dimension)
    }
    residue_two = {
        value for value, representative in parity_vectors.items()
        if norm(representative, gram) % 4 == 2
    }
    assert len(residue_two) == 65792

    # Sage returns one representative of every +/- pair in each shell.  Signs
    # have the same class modulo 2M, which is exactly the quotient used here.
    short = form.short_vector_list_up_to_length(6, True)
    roots = tuple(short[1])
    norm_six = tuple(short[3])
    norm_ten = tuple(short[5])
    assert all(norm(value, gram) == 2 for value in roots)
    assert all(norm(value, gram) == 6 for value in norm_six)
    assert all(norm(value, gram) == 10 for value in norm_ten)
    assert not roots

    excluded = {mask(value) for value in roots + norm_six}
    candidates = residue_two - excluded
    assert len(excluded) == 26672
    assert len(candidates) == 39120

    representatives = {}
    multiplicities = Counter()
    for value in norm_ten:
        orbit = mask(value)
        if orbit not in candidates:
            continue
        multiplicities[orbit] += 1
        representative = vector(ZZ, value)
        prior = representatives.get(orbit)
        if prior is None or tuple(representative) < tuple(prior):
            representatives[orbit] = representative

    assert set(representatives) == candidates
    assert sum(multiplicities.values()) == 806238
    assert all(norm(value, gram) == 10 for value in representatives.values())

    # A norm-10 representative has B=(2,2,w): B^2=-2, B.F=2, B.O=0.
    # Its minimum intersection with a section is 10/4-5/2=0.
    min_pairing = QQ(10) / 4 - QQ(5) / 2
    assert min_pairing == 0

    # Geometric consequence on any K3 realizing this rootless fibration.
    # Riemann--Roch makes B effective because B.F=2 and F is nef.  A vertical
    # effective divisor has numerical class kF: rootlessness means there are
    # no reducible-fibre components.  If B=C+kF has one degree-two horizontal
    # component C, then C^2=-2-4k, which violates C^2>=-2 for an irreducible
    # curve unless k=0.  If B=S1+S2+kF has two degree-one components, they
    # are sections and B^2=-2 gives S1.S2+2k=1.  Hence
    # B.S1=-2+S1.S2+k=-1-k<0, contradicting the certified nonnegative
    # pairing with every section.  Thus B is irreducible; adjunction then
    # makes its arithmetic (and geometric) genus zero.
    vertical_degree_two_square = -2 - 4
    assert vertical_degree_two_square == -6
    two_section_intersection_at_k_zero = 1
    assert -2 + two_section_intersection_at_k_zero == -1

    arguments.orbits_output.parent.mkdir(parents=True, exist_ok=True)
    with arguments.orbits_output.open("w") as stream:
        stream.write(
            "orbit_mask\thex\tmin_norm\tminimal_unoriented_count\t"
            "short_basis_w\tpinned_rank17_w\tlattice_bisection\n"
        )
        for orbit in sorted(representatives):
            short_vector = representatives[orbit]
            pinned_vector = short_vector * change
            assert norm(pinned_vector, pinned) == 10
            stream.write(
                f"{orbit}\t0x{orbit:05x}\t10\t{multiplicities[orbit]}\t"
                f"{entries(short_vector)}\t{entries(pinned_vector)}\t"
                f"2 2 {entries(pinned_vector)}\n"
            )

    multiplicity_histogram = {
        str(count): frequency
        for count, frequency in sorted(Counter(multiplicities.values()).items())
    }
    output = {
        "schema": "elkies-k3.rootless-bisection-orbits.v1",
        "status": "PASS_LATTICE_BISECTION_ORBITS",
        "scope": (
            "Exact Neron--Severi lattice enumeration only; no bisection "
            "equation, branch divisor, quadratic-extension hash, collision, "
            "or base-changed height matrix is asserted."
        ),
        "input": {
            "pinned_rank17_gram": display_path(RANK17),
            "pinned_rank17_gram_sha256": digest(RANK17),
            "short_basis_coordinates": display_path(SHORT_COORDS),
            "short_basis_coordinates_sha256": digest(SHORT_COORDS),
            "short_basis_gram": display_path(SHORT_GRAM),
            "short_basis_gram_sha256": digest(SHORT_GRAM),
            "rank": 17,
            "determinant": 948,
        },
        "translation_action": {
            "class": "B_w=((w.M.w-2)/4,2,w)",
            "section_translation": "w -> w + 2x",
            "quotient": "M/2M",
            "fibrewise_inversion": {
                "action_on_bisection_vector": "w -> -w",
                "action_on_translation_orbits": "identity",
                "proof": "-w-w=-2w belongs to 2M",
                "collision_search_consequence": (
                    "Fibrewise inversion creates no second translation orbit; "
                    "it cannot be counted as an independent equal-extension collision."
                ),
            },
            "total_cosets": 1 << dimension,
            "residue_two_mod_four_cosets": len(residue_two),
            "rootless_fibre_component_orbits": "trivial",
        },
        "section_nonnegative_filter": {
            "excluded_norms": [2, 6],
            "excluded_cosets": len(excluded),
            "surviving_translation_orbits": len(candidates),
            "minimum_norm_of_every_survivor": 10,
            "minimum_section_intersection": str(min_pairing),
        },
        "minimal_representatives": {
            "norm": 10,
            "lattice_bisection": [2, 2, "w"],
            "bisection_square": -2,
            "fibre_degree": 2,
            "intersection_with_zero_section": 0,
            "unoriented_norm_ten_pairs": sum(multiplicities.values()),
            "multiplicity_histogram": multiplicity_histogram,
        },
        "geometric_realization_consequence": {
            "hypothesis": (
                "Any K3 realization with this rootless elliptic fibration; "
                "F is nef and every vertical irreducible component is a whole fibre."
            ),
            "effectivity": (
                "K3 Riemann--Roch and B.F=2 force B, rather than -B, to be effective."
            ),
            "irreducibility": (
                "A vertical summand makes a single degree-two component have square "
                "below -2; two section summands force negative pairing with one section."
            ),
            "conclusion": (
                "Every section-nonnegative degree-two (-2)-class is an irreducible "
                "smooth rational bisection on such a realization."
            ),
        },
        "orbits_tsv": display_path(arguments.orbits_output),
        "orbits_tsv_sha256": digest(arguments.orbits_output),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")

    print(
        "R17BISECT|translation_orbits={}|minimum_norm=10|"
        "minimal_pairs={}|fiber_components=trivial|status=PASS_LATTICE_BISECTION_ORBITS".format(
            len(candidates), sum(multiplicities.values())
        ),
        flush=True,
    )
    print(f"R17BISECT|output={arguments.output}|sha256={digest(arguments.output)}", flush=True)


if __name__ == "__main__":
    main()
