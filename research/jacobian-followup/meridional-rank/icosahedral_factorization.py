#!/usr/bin/env python3
"""Explicit split pullback of universal H3 root lifts on the fixed bases."""

import hashlib
from itertools import product
import json
from pathlib import Path

from quadratic_integer_lifts import (ZERO,ONE,IDENTITY,add,mul,matmul,matvec,h3_model,
                                     permutation_product,permutation_inverse,ResidueField)

HERE=Path(__file__).resolve().parent


def rectmul(a,b):
    result=[]
    for row in a:
        out=[]
        for j in range(len(b[0])):
            value=ZERO
            for i,v in enumerate(row):
                value=add(value,mul(v,b[i][j]))
            out.append(value)
        result.append(out)
    return result


def universal(model,triple):
    identity=(tuple(range(5)),0)
    cocycles={identity:[[ZERO]*3 for _ in range(3)]}
    queue=[identity]
    while queue:
        key=queue.pop()
        for i,p in enumerate(triple):
            target=(permutation_product(key[0],p),key[1]^1)
            block=[row[:] for row in cocycles[key]]
            for j in range(3):
                block[j][i]=add(block[j][i],model[key][3*j+i])
            if target in cocycles:
                assert cocycles[target]==block
            else:
                cocycles[target]=block
                queue.append(target)
    assert len(cocycles)==120
    return cocycles


def char2_control(model,cocycles,simple):
    field=ResidueField(2,None)
    points=list(product(list(product(range(2),repeat=2)),repeat=3))
    index={point:i for i,point in enumerate(points)}
    actions=[]
    for key,matrix in model.items():
        translation=[row[0] for row in cocycles[key]]
        action=[]
        for point in points:
            image=[]
            for i in range(3):
                value=field.reduce(translation[i])
                for j in range(3):
                    value=field.add(value,field.mul(field.reduce(matrix[3*i+j]),point[j]))
                image.append(value)
            action.append(index[tuple(image)])
        actions.append(tuple(action))
    assert len(set(actions))==120
    unseen=set(range(64))
    orbits=[]
    while unseen:
        point=min(unseen)
        orbit={p[point] for p in actions}
        assert orbit<=unseen
        unseen-=orbit
        orbits.append(len(orbit))
    fixed=[i for i in range(64) if all(p[i]==i for p in actions)]
    return {"field":"F4, tau^2=tau+1","simple_translation_parameters":[[1,0],[0,0],[0,0]],
            "points":64,"group_order":120,"orbits":sorted(orbits),"common_fixed_points":len(fixed)}


def main():
    report=json.loads((HERE/"icosahedral_lift_results.json").read_text())
    deck=json.loads((HERE/"core_involution_audit.json").read_text())
    model,root_vectors,meta=h3_model()
    triple=meta["simple_A5_components"]
    cocycles=universal(model,triple)
    results=[]
    for case,involution in zip(report["results"],deck["results"]):
        if case["status"]!="COMPLETE":
            continue
        keys=[(tuple(p),1) for p in involution["implementers"][0]["meridian_A5_components"]]
        roots=[root_vectors[k] for k in keys]
        pullback=[]
        for key,root in zip(keys,roots):
            witnesses=[(g,i) for g,matrix in model.items() for i in range(3)
                       if tuple(matrix[3*j+i] for j in range(3))==root
                       and permutation_product(permutation_product(g[0],triple[i]),permutation_inverse(g[0]))==key[0]]
            g,i=witnesses[0]
            dual=model[(permutation_inverse(g[0]),g[1])][3*i:3*i+3]
            coefficients=rectmul([dual],cocycles[key])[0]
            assert [[mul(v,a) for a in coefficients] for v in root]==cocycles[key]
            pullback.append(coefficients)
        identity=(tuple(range(5)),0)
        words={identity:[]}
        queue=[identity]
        while queue:
            key=queue.pop()
            for i,p in enumerate(keys):
                target=(permutation_product(key[0],p[0]),key[1]^p[1])
                if target not in words:
                    words[target]=words[key]+[i+1]
                    queue.append(target)
        assert len(words)==120
        left=[]
        retained=[]
        for j,p in enumerate(triple):
            word=words[(p,1)]
            retained.append(word)
            prefix=IDENTITY
            row=[[ZERO]*len(keys) for _ in range(3)]
            for letter in word:
                arc=letter-1
                vector=matvec(prefix,roots[arc])
                for i in range(3):
                    row[i][arc]=add(row[i][arc],vector[i])
                prefix=matmul(prefix,model[keys[arc]])
            assert prefix==model[(p,1)]
            left.append(row[j])
        assert rectmul(left,pullback)==[[ONE if i==j else ZERO for j in range(3)] for i in range(3)]
        assert all(all(v==ZERO for v in row) for row in rectmul(case["root_matrix"],pullback))
        results.append({"name":case["name"],"pullback_matrix":pullback,"evaluation_left_inverse":left,
                        "simple_generator_words_in_meridians":retained,"abstract_H3_image_order":120,
                        "every_field_lift_factors_through_H3":True})
    inputs=["icosahedral_factorization.py","quadratic_integer_lifts.py","icosahedral_lift_results.json","core_involution_audit.json"]
    output={"results":results,"characteristic2_control":char2_control(model,cocycles,triple),
            "scope":"Explicit integral split pullback. Combined with the separately checked all-field nullity3, every restricted affine lift of these fixed bases factors through abstract H3, image order at most120. Characteristic2 need not give a common fixed point.",
            "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}}
    (HERE/"icosahedral_factorization.json").write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps({"fixed_bases":len(results),"characteristic2_control":output["characteristic2_control"]},indent=2))


if __name__=="__main__":
    main()
