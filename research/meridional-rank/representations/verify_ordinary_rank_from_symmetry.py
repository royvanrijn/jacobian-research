#!/usr/bin/env python3
"""Verify canonical triangulations and the filling bounds for the rank results.

Run with Sage 10.9 and SnapPy 3.3.2. In a split installation, set
MRC_SNAPPY_SAGE_PATH to the directory containing the Sage-compatible snappy.
The topological filling identifications and rank deductions are written proofs.
"""
from argparse import ArgumentParser
from hashlib import sha256
from importlib.metadata import version
from itertools import permutations
import json
import os
from pathlib import Path
import random
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SOURCE_FILES = ["all_parameter_hyperbolic_covers_inputs.json",
                "all_parameter_cover_parent.tri", "twist_exterior_symmetry_parent.tri",
                "rank_three_knot_group_inputs.json", "prime_family_cover_homology_inputs.json"]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def product(a, b):
    return tuple(a[i] for i in b)


def inverse(p):
    return tuple(p.index(i) for i in range(4))


def automorphisms(data):
    """Determine a potential simplicial isomorphism from one framed tetrahedron."""
    n = len(data)
    inv = [[inverse(g) for g in gluing] for _, gluing in data]
    count = 0
    for image in range(n):
        for frame in permutations(range(4)):
            images = {0: (image, frame)}
            used = {image}
            frontier = [0]
            conflict = False
            for i in frontier:
                k, perm = images[i]
                for side, j in enumerate(data[i][0]):
                    target_side = perm[side]
                    l = data[k][0][target_side]
                    travel = data[k][1][target_side]
                    induced = product(product(travel, perm), inv[i][side])
                    if j in images:
                        conflict = images[j] != (l, induced)
                    elif l in used:
                        conflict = True
                    else:
                        images[j] = (l, induced)
                        used.add(l)
                        frontier.append(j)
                    if conflict:
                        break
                if conflict:
                    break
            if not conflict:
                assert len(images) == n and len(used) == n
                count += 1
    return count


def interval_data(x):
    return {"low": str(x.lower().exact_rational()), "high": str(x.upper().exact_rational()),
            "display": str(x)}


def run_geometry(snappy, Link, QQ):
    from snappy.geometric_structure.cusp_neighborhood.real_cusp_cross_section import RealCuspCrossSection
    from snappy.verify.hyperbolicity import check_logarithmic_gluing_equations_and_positively_oriented_tets

    def seed(value=1729):
        random.seed(value)
        snappy.set_rand_seed(value)

    def canonical(M, expected):
        # The interval-only path requires all tilts to be strictly negative.
        seed()
        C = M.canonical_retriangulation(verified=True, interval_bits_precs=[100, 200, 400],
                                        exact_bits_prec_and_degrees=[])
        assert C is not None
        count = automorphisms(C._get_tetrahedra_gluing_data())
        assert count == len(C.isomorphisms_to(C)) == expected
        geometric = snappy.Manifold(C)
        shapes = geometric.tetrahedra_shapes('rect', intervals=True, bits_prec=150)
        check_logarithmic_gluing_equations_and_positively_oriented_tets(geometric, shapes)
        section = RealCuspCrossSection.fromManifoldAndShapes(geometric, shapes)
        if geometric.num_cusps() > 1:
            section.normalize_cusps()
        section.compute_tilts()
        upper_tilt = max(face.Tilt.upper() for face in section.mcomplex.Faces)
        assert upper_tilt < 0
        return {"input_signature": M.triangulation_isosig(), "canonical_signature": C.triangulation_isosig(),
                "canonical_tetrahedra": C.num_tetrahedra(), "tetrahedron_frame_trials": 24*C.num_tetrahedra(),
                "automorphisms": count, "max_face_tilt_upper": str(upper_tilt.exact_rational()),
                "raw_canonical_gluings_sha256":sha256(json.dumps(C._get_tetrahedra_gluing_data()).encode()).hexdigest()}

    def systole(M):
        geodesics = M.length_spectrum_alt(count=1, verified=True, bits_prec=150)
        assert geodesics and all(g.length.real().lower() > QQ(1)/2 for g in geodesics)
        return {"systole_interval": interval_data(geodesics[0].length.real()), "word": geodesics[0].word,
                "certified_lower_bound": "1/2", "verified_shortest_count_result": len(geodesics)}

    def lengths(M, cusp, a, b, c, d, const, linear):
        meridian, longitude = M.cusp_translations(verified=True, bits_prec=150)[cusp]
        area = abs(meridian.imag()) * abs(longitude)
        # Use the real projection of a*mu+b*lambda + n*(c*mu+d*lambda).
        constant = (a*meridian.real()+b*longitude) / area.sqrt()
        coefficient = (c*meridian.real()+d*longitude) / area.sqrt()
        assert constant.lower() > const and coefficient.lower() > linear
        return {"cusp":cusp,"slope":[[a,b],[c,d]], "area":interval_data(area),
                "meridian_real":interval_data(meridian.real()),"meridian_imag":interval_data(meridian.imag()),
                "longitude_real":interval_data(longitude), "projected_normalized_constant":interval_data(constant),
                "projected_normalized_step":interval_data(coefficient),
                "rational_lower_bounds":[str(const),str(linear)]}

    covers = json.loads((HERE / SOURCE_FILES[0]).read_text())
    family = covers['family']
    base = family['base_word']; suffix = family['augmentation_suffix']
    assert base == [1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3,2,2,2] and suffix == [4,3,2,2,3,-4]
    assert covers['parent']['parent'] == SOURCE_FILES[1]
    assert covers['parent']['slope_constant'] == [-1,2] and covers['parent']['slope_n_coefficient'] == [0,6]
    seed()
    augmented = Link(braid_closure=base+suffix)
    assert [list(r) for r in augmented.PD_code()] == family['augmented_PD']
    assert len(augmented.sublink([0]).crossings) == 0 and len(augmented.sublink([1]).crossings) == 17
    assert augmented.linking_matrix() == [[0,2],[2,0]]
    parent = augmented.exterior()
    fixed = snappy.Manifold(str(HERE/SOURCE_FILES[2]))
    assert parent.num_cusps()==2 and all(x['is_complete'] for x in parent.cusp_info())
    assert parent._get_tetrahedra_gluing_data() == fixed._get_tetrahedra_gluing_data()
    assert parent._get_cusp_indices_and_peripheral_curve_data() == fixed._get_cusp_indices_and_peripheral_curve_data()
    parents = {"exterior":canonical(parent,1)}
    parents['exterior']['systole'] = systole(parent)
    parents['exterior']['length'] = lengths(parent,0,-1,0,0,3,QQ(1),QQ(5)/2)
    # The vector above is minus the slope (1,-3n). Its norm is unchanged.
    N = snappy.Manifold(str(HERE / SOURCE_FILES[1]))
    assert N.num_cusps()==1 and N.cusp_info(0)['is_complete']
    parents['cover'] = canonical(N,2)
    parents['cover']['systole'] = systole(N)
    parents['cover']['length'] = lengths(N,0,-1,2,0,6,QQ(3),QQ(5))
    cutoff = QQ(101)/10
    assert 1+(QQ(5)/2)*4 > cutoff and 3+5*2 > cutoff
    assert 8/(cutoff**2-58) < QQ(1)/2

    exteriors = []
    for n in range(4):
        seed()
        L = Link(braid_closure=base+[2]*(6*n))
        assert len(L.link_components) == 1 and len(L.crossings) == 17+6*n
        M = L.exterior()
        exteriors.append({"n":n, **canonical(M,2 if n==0 else 1)})

    homology = json.loads((HERE / SOURCE_FILES[4]).read_text())
    marked_group = N.fundamental_group(False,False,False)
    expected_peripheral = homology['peripheral_words']
    assert marked_group.generators() == homology['generators']
    assert marked_group.relators() == homology['relators']
    assert marked_group.peripheral_curves() == [(expected_peripheral['meridian'], expected_peripheral['longitude'])]
    images = dict(zip(homology['generators'], homology['generator_homology_images']))

    def homology_class(word, modulus):
        return sum(images[x.lower()] * (1 if x.islower() else -1) for x in word) % modulus

    exceptional_covers = []
    for n, cutoff, number in [(0, QQ(29)/50, 2), (1, QQ(1)/10, 1)]:
        filled = N.copy()
        filled.dehn_fill((-1,2+6*n),0)
        assert filled.verify_hyperbolicity(bits_prec=200)[0]
        basis = filled.fundamental_group(simplify_presentation=False)
        assert basis.generators() == marked_group.generators()
        assert basis.relators()[:len(marked_group.relators())] == marked_group.relators()
        assert basis.peripheral_curves() == marked_group.peripheral_curves()
        geodesics = filled.length_spectrum_alt(max_len=cutoff,verified=True,bits_prec=200)
        assert len(geodesics)==number and all(g.length.real().upper()<cutoff for g in geodesics)
        cores = [g for g in geodesics if g.core_curve == 0]
        assert len(cores) == 1
        core = cores[0]
        drilled = filled.drill_word(core.word, verified=True,bits_prec=200)
        assert drilled.num_cusps() == 1 and drilled.cusp_info(0)['is_complete']
        assert drilled.isomorphisms_to(N)
        modulus = 117+366*n
        classes = {g.word: homology_class(g.word, modulus) for g in geodesics}
        if n == 0:
            assert core.word == 'bcckdf' and classes == {'bcckdf':56,'b':43}
            assert classes['b'] not in {classes[core.word], -classes[core.word]%modulus}
        exceptional_covers.append({"n":n,"filled_slope":[-1,2+6*n],"cutoff":str(cutoff),
                                   "core_word":core.word,"homology_modulus":modulus,
                                   "short_geodesics":[{"word":g.word,"length":interval_data(g.length.real()),
                                                       "homology":classes[g.word],"is_core":g.core_curve==0}
                                                      for g in geodesics],
                                   "verified_drill_returns_parent_combinatorially":True})

    named = []
    data = json.loads((HERE / SOURCE_FILES[3]).read_text())
    assert len(data['cases'])==32 and len({d['name'] for d in data['cases']})==32
    for x in data['cases']:
        M = snappy.Manifold(x['snappea_triangulation'])
        assert M.triangulation_isosig()==x['triangulation_isosig']
        named.append({"knot":x['name'],**canonical(M,1)})
    return {"parents":parents,"direct_exteriors":exteriors,"exceptional_covers":exceptional_covers,"named_cases":named}


def main():
    cli = ArgumentParser(description=__doc__)
    cli.add_argument('--snappy-path', default=os.environ.get('MRC_SNAPPY_SAGE_PATH'))
    args = cli.parse_args()
    candidates = [args.snappy_path, REPO/'.scratch/mrc-sage-snappy',
                  Path('/tmp/mrc-knot-audit-venv/lib/python3.12/site-packages')]
    for p in candidates:
        if p and Path(p).is_dir():
            sys.path.append(str(p))
    from sage.all import QQ
    from sage.env import SAGE_VERSION
    import snappy
    from spherogram import Link
    assert SAGE_VERSION=='10.9' and snappy.__version__=='3.3.2' and version('spherogram')=='2.4.1'
    package = Path(snappy.__file__).resolve().parent
    modules = ['verify/canonical.py','verify/cusp_translations.py','len_spec/__init__.py',
               'len_spec/geometric_structure.py','drilling/__init__.py',
               'geometric_structure/geodesic/add_core_curves.py',
               'geometric_structure/cusp_neighborhood/real_cusp_cross_section.py']
    output = {"status":"PASS verified symmetry and filling inputs for ordinary ranks",
              "source_sha256":{name:digest(HERE/name) for name in SOURCE_FILES},
              "checker_sha256":digest(Path(__file__)),"sage":SAGE_VERSION,"snappy":snappy.__version__,
              "spherogram":version('spherogram'),"native_kernel_sha256":digest(Path(snappy.SnapPy.__file__)),
              "implementation_sha256":{name:digest(package/name) for name in modules},
              "geometry":run_geometry(snappy, Link, QQ),
              "scope":"Verified finite geometry, exact automorphism enumeration and interval bounds. Topological filling, effective surgery and ordinary rank deductions are in the written proof.",
              "assurance":{"whole_proof_independent_replay":False,"formal_verification":False,"external_specialist_review":False}}
    print(json.dumps(output,indent=2))


if __name__ == '__main__':
    main()
