#!/usr/bin/env python3
"""Independent regular-matrix/affine-point replay of the fixed H3 lifts.

Main replay uses ordinary integer6x6 matrices, not the producer's
quadratic-integer operations or root-matrix constructor. PD parsing and
the retained elementary Smith-step checker are shared.
"""

import hashlib
from itertools import permutations
import json
from pathlib import Path

from parabolic_lifts import oriented_diagram
from verify_integer_lifts import nullspace,replay

HERE=Path(__file__).resolve().parent


def compose(a,b):
    return tuple(a[b[i]] for i in range(5))


def inverse(a):
    return tuple(a.index(i) for i in range(5))


def matrix_product(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def expand(matrix,dimension):
    result=[[0]*(2*dimension) for _ in range(2*dimension)]
    for i in range(dimension):
        for j in range(dimension):
            a,b=matrix[dimension*i+j]
            result[2*i][2*j]=a
            result[2*i][2*j+1]=b
            result[2*i+1][2*j]=b
            result[2*i+1][2*j+1]=a+b
    return result


def root_columns(root):
    return [[x for x in pair] for a,b in root for pair in [(a,b),(b,a+b)]]


def build_group(certificate):
    simple=[
        [(-1,0),(0,1),(0,0),(0,0),(1,0),(0,0),(0,0),(0,0),(1,0)],
        [(1,0),(0,0),(0,0),(0,1),(-1,0),(1,0),(0,0),(0,0),(1,0)],
        [(1,0),(0,0),(0,0),(0,0),(1,0),(0,0),(0,0),(1,0),(-1,0)],
    ]
    assert [[list(v) for v in m] for m in simple]==certificate["simple_matrices"]
    generators=[expand(m,3) for m in simple]
    components=[tuple(p) for p in certificate["simple_A5_components"]]
    identity=tuple(range(5))
    unit=[[int(i==j) for j in range(6)] for i in range(6)]
    group={(identity,0):unit}
    queue=[(identity,0)]
    while queue:
        key=queue.pop()
        for p,generator in zip(components,generators):
            assert compose(p,p)==identity and sum(p[i]!=i for i in range(5))==4
            target=(compose(key[0],p),key[1]^1)
            matrix=matrix_product(group[key],generator)
            if target in group:
                assert group[target]==matrix
            else:
                group[target]=matrix
                queue.append(target)
    assert len(group)==len({tuple(x for row in m for x in row) for m in group.values()})==120
    for key,matrix in group.items():
        assert matrix_product(matrix,group[(inverse(key[0]),key[1])])==unit
    return group,generators,components,unit


def point_matrix(matrices,roots,words):
    n=len(matrices)
    output=[]
    for word in words:
        block=[[0]*(2*n) for _ in range(6)]
        for variable in range(2*n):
            point=[0]*6
            for letter in reversed(word):
                arc=abs(letter)-1
                point=[sum(row[k]*point[k] for k in range(6)) for row in matrices[arc]]
                if arc==variable//2:
                    point=[v+roots[arc][i][variable%2] for i,v in enumerate(point)]
            for i,value in enumerate(point):
                block[i][variable]=value
        output.extend(block)
    return output


def expand_rectangular(rows):
    result=[]
    for row in rows:
        result.append([x for a,b in row for x in [a,b]])
        result.append([x for a,b in row for x in [b,a+b]])
    return result


def universal_translations(group,components):
    identity=(tuple(range(5)),0)
    translations={identity:[[0]*6 for _ in range(6)]}
    queue=[identity]
    while queue:
        key=queue.pop()
        for i,p in enumerate(components):
            target=(compose(key[0],p),key[1]^1)
            value=[row[:] for row in translations[key]]
            for j in range(2):
                for row in range(6):
                    value[row][2*i+j]+=group[key][row][2*i+j]
            if target in translations:
                assert translations[target]==value
            else:
                translations[target]=value
                queue.append(target)
    assert len(translations)==120
    return translations


def binary_control(group,translations):
    from itertools import product
    points=list(product(range(2),repeat=6))
    index={p:i for i,p in enumerate(points)}
    actions=[]
    for key,matrix in group.items():
        action=[]
        for point in points:
            image=tuple((sum(row[j]*point[j] for j in range(6))+translations[key][i][0])%2 for i,row in enumerate(matrix))
            action.append(index[image])
        actions.append(tuple(action))
    assert len(set(actions))==120
    unseen=set(range(64))
    sizes=[]
    while unseen:
        orbit={p[min(unseen)] for p in actions}
        assert orbit<=unseen
        unseen-=orbit
        sizes.append(len(orbit))
    assert sorted(sizes)==[12,12,20,20]
    assert not any(all(p[i]==i for p in actions) for i in range(64))
    return sorted(sizes)


def expansion_controls():
    # These exercise the producer's generic field routines separately from
    # the main integer-only replay, including inert, split and ramified primes.
    from quadratic_integer_lifts import integral_expansion,residue_fields,kernel
    checks=[]
    for label,coefficient,p,expected in [
        ("inert prime",(3,0),3,[1]),
        ("split prime one-embedding jump",(-4,1),11,[1,0]),
        ("ramified discriminant",(-1,2),5,[1]),
        ("quadratic unit",(0,1),5,[0]),
    ]:
        scalar=[[coefficient]]
        a,b=coefficient
        expanded=[[a,b],[b,a+b]]
        assert integral_expansion(scalar)==expanded
        dimensions=[len(kernel(scalar,1,f)) for f in residue_fields(p)]
        assert dimensions==expected
        if p==3:
            assert len(nullspace(expanded,2,p))==2*dimensions[0]
        if p==5 and coefficient==(-1,2):
            # Expanded kernel has dimension1, while the residue-field
            # kernel also has dimension1: blindly halving is invalid.
            assert len(nullspace(expanded,2,p))==1
        checks.append(label)
    for p in [2,3]:
        scalar=[[(1,0),(0,1)],[(0,1),(1,1)]]
        field=residue_fields(p)[0]
        basis=kernel(scalar,2,field)
        assert len(basis)==1 and len(nullspace(integral_expansion(scalar),4,p))==2
        for row in scalar:
            total=(0,0)
            for coefficient,value in zip(row,basis[0]):
                total=field.add(total,field.mul(field.reduce(coefficient),value))
            assert total==(0,0)
    return checks


def main():
    report=json.loads((HERE/"icosahedral_lift_results.json").read_text())
    deck=json.loads((HERE/"core_involution_audit.json").read_text())
    prior=json.loads((HERE/"nonabelian_core_results.json").read_text())
    panel=json.loads((HERE/"parabolic_panel.json").read_text())
    for receipt in [report,deck]:
        for path,digest in receipt["input_sha256"].items():
            assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest
    group,simple,components,unit=build_group(report["H3_group"])
    factors=json.loads((HERE/"icosahedral_factorization.json").read_text())
    for path,digest in factors["input_sha256"].items():
        assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest
    translations=universal_translations(group,components)
    factor_by_name={r["name"]:r for r in factors["results"]}
    a5=[p for p in permutations(range(5)) if sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))%2==0]
    results=[]
    assert len(report["results"])==len(deck["results"])==len(panel["cases"])==len(prior["results"])==16
    for result,involution,case,old in zip(report["results"],deck["results"],panel["cases"],prior["results"]):
        assert result["name"]==involution["name"]==case["name"]==old["name"]
        if old["first_nonconstant_colouring"] is None:
            assert result["status"]=="NO_RETAINED_H3_BASE" and old["normalized_colourings"]==1
            continue
        core=[a5[i] for i in old["first_nonconstant_colouring"]["arc_labels"]]
        implementers=[c for c in permutations(range(5)) if all(compose(compose(c,a),inverse(c))==inverse(a) for a in core)]
        assert len(implementers)==1 and list(implementers[0])==involution["implementers"][0]["deck_implementer"]
        keys=[(compose(a,implementers[0]),1) for a in core]
        matrices=[group[k] for k in keys]
        assert matrices==[expand(m,3) for m in result["meridian_matrices"]]
        roots=[root_columns(v) for v in result["root_vectors"]]
        for key,matrix,root in zip(keys,matrices,roots):
            assert matrix_product(matrix,root)==[[-v for v in row] for row in root]
            # This explicit conjugating group element certifies that the
            # root remains primitive in every residue characteristic.
            assert any(
                compose(compose(k[0],p),inverse(k[0]))==key[0]
                and [[g[row][2*i+j] for j in range(2)] for row in range(6)]==root
                for k,g in group.items() for i,p in enumerate(components))
        crossings=oriented_diagram(case["pd"])
        for over,a,b,_ in crossings:
            assert matrix_product(matrix_product(matrices[over],matrices[a]),matrices[over])==matrices[b]
        words=[[o+1,a+1,o+1,b+1] for o,a,b,_ in crossings]
        expanded=point_matrix(matrices,roots,words)
        # Check every retained O coefficient against the independently
        # evaluated integer affine action, then replay the Smith certificate.
        expected=[]
        for row in result["root_matrix"]:
            expected.append([x for a,b in row for x in [a,b]])
            expected.append([x for a,b in row for x in [b,a+b]])
        assert expanded==expected
        factor=factor_by_name[case["name"]]
        pullback=expand_rectangular(factor["pullback_matrix"])
        left=expand_rectangular(factor["evaluation_left_inverse"])
        assert matrix_product(left,pullback)==unit
        assert not any(any(row) for row in matrix_product(expanded,pullback))
        for arc,key in enumerate(keys):
            assert matrix_product(roots[arc],pullback[2*arc:2*arc+2])==translations[key]
        assert len(factor["simple_generator_words_in_meridians"])==3
        for i,word in enumerate(factor["simple_generator_words_in_meridians"]):
            target=(tuple(range(5)),0)
            for letter in word:
                target=(compose(target[0],keys[letter-1][0]),target[1]^1)
            assert target==(components[i],1)
            evaluated=point_matrix(matrices,roots,[word])
            assert matrix_product(evaluated,pullback)==translations[target]
        diagonal=replay(expanded,result["expanded_smith"])
        assert diagonal==[1]*20 and result["generic_dimension"]==3
        assert result["candidate_primes"]==[5]
        scalar5=[[(a+3*b)%5 for a,b in row] for row in result["root_matrix"]]
        basis5=nullspace(scalar5,13,5)
        assert len(basis5)==3
        assert len(result["residue_fields"])==1
        field=result["residue_fields"][0]
        assert (field["prime"],field["tau_root"],field["field_degree"],field["dimension"])==(5,3,1,3)
        basis=[]
        for vector in field["basis"]:
            assert all(b==0 and 0<=a<5 for a,b in vector)
            values=[a for a,b in vector]
            assert all(sum(a*b for a,b in zip(row,values))%5==0 for row in scalar5)
            basis.append(values)
        assert len(nullspace(basis,13,5))==10
        assert result["maximum_dimension_all_fields"]==3 and not result["sharp_four"]
        results.append({"name":case["name"],"expanded_Smith_factors":[1]*20,
                        "generic_dimension":3,"characteristic5_dimension":3,"maximum_all_fields":3,"sharp_four":False})
    assert len(results)==12
    orbits=binary_control(group,translations)
    assert factors["characteristic2_control"]["orbits"]==orbits
    assert factors["characteristic2_control"]["common_fixed_points"]==0
    assert factors["characteristic2_control"]["group_order"]==120
    inputs=["verify_icosahedral_lifts.py","quadratic_integer_lifts.py","audit_icosahedral_lifts.py",
            "icosahedral_lift_plan.json","icosahedral_lift_results.json","audit_core_involutions.py",
            "core_involution_audit.json","nonabelian_core_results.json","parabolic_panel.json",
            "parabolic_lifts.py","verify_integer_lifts.py","icosahedral_factorization.py","icosahedral_factorization.json"]
    receipt={"status":"PASS","fixed_H3_bases":12,"results":results,
             "quadratic_field_controls":expansion_controls(),
             "integral_split_factorizations":12,
             "every_field_fixed_base_affine_image_order_at_most":120,
             "characteristic2_no_fixed_point_control_orbits":orbits,
             "scope":"Independent integer6x6 group and reverse affine-point evaluation, all20 unit Smith factors replayed, characteristic5 scalar elimination and bases checked. PD parser/Smith replay shared; only fixed first-map H3 bases excluded. Uniform all-field reduction is a written argument, not inferred from controls.",
             "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}}
    (HERE/"icosahedral_lift_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"status":"PASS","fixed_H3_bases":12,"all_field_dimensions":[3],"new_sharp_cases":0,
                      "quadratic_field_controls":receipt["quadratic_field_controls"]},indent=2))


if __name__=="__main__":
    main()
