#!/usr/bin/env python3
"""Exact complete graded cap23 odd first-jet kernel and normalization controls.

Replay uses Python integers/Fraction only. The existing independent curvature
engine supplies polynomial arithmetic; the canonical producer supplies input
polynomials only. No nullspace solver, SymPy, FLINT or floating arithmetic is
used by this verifier.
"""
from fractions import Fraction as Q
from pathlib import Path
from math import isqrt
from itertools import combinations_with_replacement
import argparse, hashlib, importlib.util, json, platform, tempfile

HERE = Path(__file__).resolve().parent
CAP = 23


def load(path, name):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def basis(weights):
    """All degree<=23 monomials of weight w_i-3: a=b+4c+2d+3-w_i."""
    out = []
    for i, weight in enumerate(weights):
        for b in range(CAP + 1):
            for c in range(CAP + 1 - b):
                for d in range(CAP + 1 - b - c):
                    a = b + 4*c + 2*d + 3 - weight
                    if a >= 0 and a+b+c+d <= CAP:
                        out.append((i, (a,b,c,d)))
    return sorted(out)


def determinant_Q(matrix):
    a = [list(row) for row in matrix]
    n = len(a)
    assert all(len(row) == n for row in a)
    result = Q(1)
    for j in range(n):
        pivot = next((i for i in range(j,n) if a[i][j]), None)
        if pivot is None:
            return Q(0)
        if pivot != j:
            a[pivot],a[j] = a[j],a[pivot]
            result = -result
        q = a[j][j]
        result *= q
        for i in range(j+1,n):
            factor = a[i][j]/q
            for k in range(j+1,n):
                a[i][k] -= factor*a[j][k]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    cert = load(HERE/"check_graded_long_cap23.py", "graded_linear_engine")
    alg = load(HERE/"check_long_graph_curvature.py", "independent_polynomial_engine")
    source = HERE.parent/"dixmier-audit/verify_long_pbw.py"
    producer = load(source, "canonical_polynomial_input")
    images = [{tuple(e):Q(c) for e,c in poly.items()}
              for poly in producer.long_images()]  # INPUT ONLY.
    witness_path = HERE/"odd_graded_long_first_jets_witness.json"
    witness = json.loads(witness_path.read_text())
    assert witness["cap"] == CAP
    assert witness["image_order"] == list(cert.LABELS)
    assert witness["source_weights"] == [-1,1,4,2]
    assert witness["image_weights"] == list(cert.WEIGHTS)
    assert witness["h_weight"] == 3
    p = witness["prime"]
    assert p == 1000003 and all(p%d for d in range(2,isqrt(p)+1))
    unknowns = basis(cert.WEIGHTS)
    counts = [sum(i==j for i,e in unknowns) for j in range(4)]
    assert counts == [80,114,140,102] and len(unknowns) == 436
    for poly,weight in zip(images,cert.WEIGHTS):
        assert all(-e[0]+e[1]+4*e[2]+2*e[3] == weight for e in poly)
    columns = cert.linear_columns(alg,images,unknowns)
    rows = sorted(set().union(*(set(c) for c in columns)))
    assert witness["linear_matrix"] == {
        "rows":len(rows), "columns":len(columns),
        "nonzeros":sum(map(len,columns))}
    assert len(rows) == 2494
    rank_det = cert.minor(columns,{},rows,witness["rank_minor"],p)
    assert len(witness["rank_minor"]["rows"]) == 430

    one = {(0,0,0,0):Q(1)}
    variables = [{tuple(int(i==j) for j in range(4)):Q(1)}
                 for i in range(4)]
    X,Qpos,P,Z = variables
    r,T,D,S = images
    add,mul,scale = alg.add,alg.multiply,alg.scale
    power = lambda poly,n: cert.power(alg,poly,n)
    # Independently check all six Poisson frame equations.
    for i,j in cert.PAIRS:
        assert alg.bidifferential(images[i],images[j],1) == \
               alg.bidifferential(variables[i],variables[j],1)
    t = mul(X,Qpos)
    B = add(scale(mul(power(X,2),P),3),
            scale(mul(add(one,scale(t,-3)),Z),2))
    A0 = add(Qpos,scale(mul(X,power(Qpos,2)),3),
             scale(mul(X,B),Q(-1,3)))
    Y = add(mul(r,A0),scale(t,-3))
    H4 = add(power(Y,4),scale(mul(t,power(Y,2)),4),
             scale(mul(t,Y),Q(-16,3)),scale(power(t,2),10))
    correction = add(scale(power(Y,5),9),
        scale(mul(t,power(Y,3)),45),scale(mul(power(t,2),Y),45),
        scale(mul(t,power(Y,2)),-60),scale(mul(t,Y),110),
        scale(power(t,2),-720))
    H6 = add(mul(power(r,4),D),scale(correction,Q(1,25)))
    generators = [t,Y,power(Y,2),mul(power(r,2),T),H4,H6]
    names = ["XQ", "Y=R*A0-3XQ", "Y^2", "R^2*T",
             "Y^4+4tY^2-(16/3)tY+10t^2",
             "R^4*D+(9Y^5+45tY^3+45t^2Y-60tY^2+110tY-720t^2)/25"]
    jets = []
    vectors = []
    degrees = []
    for G in generators:
        assert all(-e[0]+e[1]+4*e[2]+2*e[3] == 0 for e in G)
        jet = [alg.bidifferential(G,f,1) for f in images]
        vector = cert.vector(unknowns,jet)  # Enforces full cap/weight admission.
        assert not cert.product(columns,vector)
        jets.append(jet)
        vectors.append(vector)
        degrees.append([max(map(sum,poly),default=-1) for poly in jet])
        # Recover the Hamiltonian by the full symplectic Jacobian and
        # ordinary Euler integration; weighted Euler has denominator zero.
        def derivative(poly,i):
            return alg.derivative(poly,tuple(int(i==j) for j in range(4)))
        gradient = [add(scale(mul(jet[2],derivative(r,i)),-1),
                        scale(mul(jet[3],derivative(T,i)),-1),
                        mul(jet[0],derivative(D,i)),
                        mul(jet[1],derivative(S,i))) for i in range(4)]
        recovered = {}
        for i in range(4):
            for exponent,c in gradient[i].items():
                e = list(exponent)
                e[i] += 1
                e = tuple(e)
                recovered[e] = recovered.get(e,Q(0)) + c/Q(sum(e))
        recovered = {e:c for e,c in recovered.items() if c}
        assert recovered == G
        assert all(derivative(G,i) == gradient[i] for i in range(4))
    assert degrees == [[-1,15,23,11],[-1,8,16,4],[-1,15,23,11],
                       [-1,-1,18,6],[-1,17,23,13],[12,17,23,20]]
    independence = witness["kernel_independence_minor"]
    indices = independence["unknown_indices"]
    assert len(indices) == len(set(indices)) == 6
    assert independence["unknowns"] == [
        [unknowns[i][0],list(unknowns[i][1])] for i in indices]
    independent_det = determinant_Q([[v[i] for v in vectors] for i in indices])
    assert independent_det == Q(independence["det_Q"]) == Q(-10240,9)

    # Exact relations identify the two uniformly cap-preserving shears.
    assert mul(r,S) == add(scale(t,-6),scale(power(Y,2),-3),scale(Y,-8))
    assert alg.bidifferential(t,P,1) == scale(Qpos,-1)
    assert alg.bidifferential(t,Z,1) == scale(X,-1)
    assert not alg.bidifferential(t,X,1)
    assert not alg.bidifferential(t,Qpos,1)
    RSjet = [alg.bidifferential(mul(r,S),f,1) for f in images]
    assert RSjet == [{},r,scale(S,-1),{}]
    assert jets[3] == [{},{},scale(mul(r,T),-2),scale(power(r,2),-1)]
    # Natural classical exponentials need not be finite polynomial maps.
    assert alg.bidifferential(Y,X,1) == scale(mul(r,power(X,3)),-1)
    assert not alg.bidifferential(Y,r,1)
    assert jets[5][0] == power(r,4)
    # Applying the R^2T target shear to a jet with H6 coefficient nonzero
    # produces an h^2 D term of degree27. The Weyl product h^1 correction
    # {R,T}/2 is zero, so it cannot remove this leading degree.
    left = mul(jets[5][0],T)
    right = mul(r,jets[5][1])
    cross = add(left,right)
    assert max(map(sum,left)) == 27
    assert max(map(sum,right)) == 20
    assert max(map(sum,cross)) == 27
    assert not alg.bidifferential(r,T,1)

    # FULL coupled commutator-order h^3 compatibility. Only the two
    # uniformly bounded first gauges above have been removed.
    coupled_jets = [jets[i] for i in (1,3,4,5)]
    pair_parameters = list(combinations_with_replacement(range(4),2))
    rhses = [cert.lower_rhs(alg,[images],1)]
    for a,b in pair_parameters:
        rhs = {}
        for pi,(i,j) in enumerate(cert.PAIRS):
            value = alg.bidifferential(coupled_jets[a][i],coupled_jets[b][j],1)
            if a != b:
                value = add(value,alg.bidifferential(
                    coupled_jets[b][i],coupled_jets[a][j],1))
            for exponent,c in value.items():
                rhs[(pi,exponent)] = -c
        rhses.append(rhs)
    second_unknowns = cert.basis(CAP,1)
    second_columns = cert.linear_columns(alg,images,second_unknowns)
    second_rows = sorted(set().union(*(set(c) for c in second_columns),
                                     *(set(b) for b in rhses)))
    h3 = witness["h3_compatibility"]
    assert len(second_unknowns) == h3["h2_unknowns"] == 307
    assert len(second_rows) == h3["rows"] == 2005
    even_witness_path = HERE/"graded_long_cap23_witness.json"
    even_witness = json.loads(even_witness_path.read_text())
    second_rank_det = cert.minor(second_columns,rhses[0],second_rows,
                                 even_witness["h3_minor"],p)
    second_kernels = [
        [{},{},scale(power(r,2),-3),{}],
        [{},power(r,4),scale(mul(power(r,3),S),-4),{}]]
    for kernel in second_kernels:
        assert not cert.product(second_columns,cert.vector(second_unknowns,kernel))
    assert second_kernels[0][2] and not second_kernels[0][1]
    assert second_kernels[1][1]
    # Thus the full second coefficient matrix has rational rank305 and
    # exactly these TWO additional homogeneous directions; they are kept.
    functional = {}
    for index,key,c in h3["functional"]:
        key = (key[0],tuple(key[1]))
        assert second_rows[index] == key and key not in functional
        functional[key] = Q(c)
    def apply(poly):
        return sum((c*poly.get(key,Q(0)) for key,c in functional.items()),Q(0))
    assert all(not apply(column) for column in second_columns)
    h3_values = [apply(b) for b in rhses]
    assert h3_values == [Q(0)]*10+[Q(1)]
    # Necessity is d^2=0. Over every characteristic-zero field this is d=0.
    required_indices = [0,1,2,3,5,6,8]
    required_names = ["1","a^2","a*b","a*c","b^2","b*c","c^2"]
    particulars = []
    for item,rhs_index,monomial in zip(
            h3["h2_particulars_on_d_zero"],required_indices,required_names):
        assert item["rhs_index"] == rhs_index and item["monomial"] == monomial
        candidate = [{} for _ in range(4)]
        for label,exponent,c in item["coefficients"]:
            i = cert.LABELS.index(label)
            exponent = tuple(exponent)
            assert exponent not in candidate[i]
            candidate[i][exponent] = Q(c)
        vector = cert.vector(second_unknowns,candidate)
        assert cert.product(second_columns,vector) == rhses[rhs_index]
        assert not candidate[0]
        assert all(e[2]+e[3] == 0 for e in candidate[3])
        particulars.append(candidate)
    assert len(particulars) == 7
    future_v_degree = max(map(sum,mul(power(r,3),particulars[-1][3])))
    assert future_v_degree == 24
    # Deliberate invalid controls; canonical witness data are never mutated.
    bad = [dict(poly) for poly in particulars[0]]
    bad[3][(5,0,0,0)] += 1
    assert cert.product(second_columns,cert.vector(second_unknowns,bad)) != rhses[0]
    bad_functional = dict(functional)
    bad_key = next(iter(rhses[-1]))
    bad_functional[bad_key] = bad_functional.get(bad_key,Q(0))+1
    assert sum((c*rhses[-1].get(key,Q(0))
                for key,c in bad_functional.items()),Q(0)) != 1
    # Source t fixes positional R0 and R1=dR^4; target RS fixes R exactly.
    # Undoing the two allowed first normalizations cannot change R1 or R2.
    assert not alg.bidifferential(t,power(r,4),1)

    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    receipt = {
        "claim":"Complete cap23 odd first kernel and nonlinear h3 compatibility: H6 vanishes; three parameters survive",
        "status":"PASS",
        "environment":{"python":platform.python_version(),
            "arithmetic":"Python integers and fractions.Fraction; no symbolic/nullspace library"},
        "commutator_coefficients":["h^2 homogeneous first kernel","h^3 coupled quadratic compatibility"],
        "basis_counts":counts, "unknowns":436, "rows":len(rows),
        "nonzeros":sum(map(len,columns)), "prime":p,
        "rank_minor_size":430, "rank_minor_det_modp":rank_det,
        "independent_kernel_minor_size":6,
        "independent_kernel_minor_det_Q":str(independent_det),
        "rank_Q":430, "kernel_dimension_Q":6,
        "compact_hamiltonians":names,
        "Hamiltonian_source_degrees":[max(map(sum,G)) for G in generators],
        "jet_degrees_in_R_T_D_S_order":degrees,
        "full_symplectic_jacobian_and_ordinary_Euler_reconstruction":True,
        "cap_preserving_shears":["source XQ","target RS"],
        "R2T_shear_H6_cross_RT_degree":27,
        "normalization_gap":"No subsequent cap-preserving elimination or evenization is proved",
        "h3_compatibility":{
            "normalized_parameters":["a*Y","b*R^2*T","c*H4","d*H6"],
            "h2_unknowns":307,"rows":len(second_rows),
            "second_rank_Q":305,"second_rank_minor_det_modp":second_rank_det,
            "quadratic_pairs":[list(z) for z in pair_parameters],
            "rhs_nonzeros":[len(b) for b in rhses],
            "functional_support":len(functional),
            "all307columns_annihilated":True,
            "exact_functional_values":[str(c) for c in h3_values],
            "necessary_sufficient_over_fields":"d=0",
            "surviving_parameters":["Y","R^2*T","H4"],
            "exact_h2_particular_monomials":required_names,
            "exact_particular_coefficients":sum(
                sum(map(len,z)) for z in particulars),
            "full_h2_family":"particular(a,b,c)+u*K_R3+v*K_R4S",
            "first_R_and_second_R_forced_zero":True,
            "normalized_h2_S_position_only":True,
            "future_v_gauge_c2_D4_degree":future_v_degree,
            "corrupted_particular_rejected":True,
            "corrupted_functional_rejected":True},
        "natural_classical_flows":{
            "Y":"X(tau)=X/(1+2tau*R*X^2)^(1/2), R fixed",
            "H6":"R(tau)=R/(1-3tau*R^3)^(1/3)",
            "scope":"Natural classical Hamiltonian exponentials only; no exclusion of alternative nonlinear gauges"},
        "hashes":{
            "ODD_GRADED_LONG_FIRST_JETS_AND_NORMALIZATION_GAP.md":sha(HERE/"ODD_GRADED_LONG_FIRST_JETS_AND_NORMALIZATION_GAP.md"),
            "check_odd_graded_long_first_jets.py":sha(Path(__file__)),
            "odd_graded_long_first_jets_witness.json":sha(witness_path),
            "check_graded_long_cap23.py":sha(HERE/"check_graded_long_cap23.py"),
            "graded_long_cap23_witness.json":sha(even_witness_path),
            "check_long_graph_curvature.py":sha(HERE/"check_long_graph_curvature.py"),
            "../dixmier-audit/verify_long_pbw.py":sha(source)},
        "coverage_boundary":"Complete homogeneous cap23 first kernel and exact nonlinear h3 compatibility over characteristic-zero fields; H6 forced zero, all three remaining parameters have retained h2 points with two additional homogeneous freedoms. No h4/h5 compatibility, full corrected lift, evenization, nongraded exclusion or DC1/DC2 theorem.",
        "independent_replay":False, "formal_verification":False, "external_review":False}
    output = args.output_dir or Path(tempfile.mkdtemp(prefix="odd-graded-long-first-jets-"))
    output.mkdir(parents=True,exist_ok=True)
    path = output/"odd_graded_long_first_jets_controls.json"
    path.write_text(json.dumps(receipt,indent=2)+"\n")
    print("PASS: odd first rank430/436; h3 forces H6=0 and retains all Y,R2T,H4 parameters")
    print("Minor determinants:",rank_det,"mod",p,"and",independent_det,"over Q")
    print("Receipt:",path)


if __name__ == "__main__":
    main()
