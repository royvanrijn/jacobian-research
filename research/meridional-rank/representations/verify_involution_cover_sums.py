#!/usr/bin/env python3
"""Finite controls for the normalized involution-cover sum formulas.

Uses the retained literal seven-sheet seed and adds meridional trefoil
amalgams. All meridian squares are imposed first, which is valid after
branched filling; fixed-sheet meridian slopes are then filled explicitly.
The uniform additivity theorem remains a written topological argument.
--census additionally rereads the pinned current named census on all101 seeds.
"""
from fractions import Fraction
from math import gcd
from pathlib import Path
import hashlib
import json
import sys

from verify_general_branched_covers import covering_matrix, determinant
from verify_seven_sheet_covers import diagram, census_bijection, permutations as binary3_permutations
from itertools import combinations

HERE = Path(__file__).resolve().parent


def rank_mod3(matrix):
    rows = [[x % 3 for x in row] for row in matrix]
    columns = len(rows[0])
    position = 0
    for column in range(columns):
        pivot = next((i for i in range(position,len(rows)) if rows[i][column]),None)
        if pivot is None:
            continue
        rows[position],rows[pivot] = rows[pivot],rows[position]
        factor = 1 if rows[position][column]==1 else 2
        rows[position] = [(factor*x)%3 for x in rows[position]]
        for i in range(position+1,len(rows)):
            factor = rows[i][column]
            if factor:
                rows[i] = [(x-factor*y)%3 for x,y in zip(rows[i],rows[position])]
        position += 1
        if position == len(rows):
            break
    return position


def seed_data():
    inputs = json.loads((HERE/'seven_sheet_cover_inputs.json').read_text())
    case = next(c for c in inputs['cases'] if c['name']=='13a_373')
    components = [{v} for v in range(1,27)]
    for row in case['pd']:
        a = next(c for c in components if row[1] in c)
        b = next(c for c in components if row[3] in c)
        if a is not b:
            components.remove(a)
            components.remove(b)
            components.append(a|b)
    components = sorted(sorted(c) for c in components)
    arc = {v:i for i,c in enumerate(components) for v in c}
    action = next(a for a in case['actions'] if a['module']=='natural')
    permutations = [tuple(sum(((r&v).bit_count()%2)<<j for j,r in
                              enumerate(matrix))-1 for v in range(1,8))
                    for matrix in action['meridian_matrices_F2']]
    assert len(permutations) == 13
    relators = [[arc[row[1]]+1,arc[row[0]]+1,
                 arc[row[1]]+1,arc[row[2]]+1] for row in case['pd']]
    # Meridians become involutions in the orbifold quotient. These
    # additional relators are redundant after physical branched filling.
    relators += [[g,g] for g in range(1,14)]
    return permutations,relators


def bipartite_control(a,b):
    g = gcd(a,b)
    left,right = b//g,a//g
    delta = a*b//g
    left_sets = [set() for _ in range(left)]
    right_sets = [set() for _ in range(right)]
    edge = 0
    for i in range(left):
        for j in range(right):
            for _ in range(g):
                pair = {2*edge,2*edge+1}
                left_sets[i] |= pair
                right_sets[j] |= pair
                edge += 1
    assert edge == delta
    assert all(len(s)==2*a for s in left_sets)
    assert all(len(s)==2*b for s in right_sets)
    n = 2*delta+left+right
    for i in range(left):
        left_sets[i].add(2*delta+i)
    for j in range(right):
        right_sets[j].add(2*delta+left+j)
    # Verify the union of active components is connected on every sheet.
    reached = set(left_sets[0])
    while True:
        enlarged = set(reached)
        for s in left_sets+right_sets:
            if reached & s:
                enlarged |= s
        if enlarged == reached:
            break
        reached = enlarged
    assert reached == set(range(n))
    # Arbitrary abstract component ranks control the sphere arithmetic.
    rank_left,rank_right = a+1,b+2
    components_left = left+right  # active left copies + right fixed labels
    components_right = left+right
    cycles = n-delta
    beta = cycles-components_left-components_right+1
    assert beta >= 0
    rank_total = left*rank_left+right*rank_right+beta
    Q_total = n-1+rank_total
    Q_left,Q_right = 2*a+rank_left,2*b+rank_right
    assert Q_total == left*Q_left+right*Q_right-delta
    assert Fraction(Q_total,delta) == Fraction(Q_left,a)+Fraction(Q_right,b)-1
    return {'a':a,'b':b,'copies_left':left,'copies_right':right,
            'common_delta':delta,'degree_with_one_fixed_per_copy':n,
            'sphere_graph_beta':beta}


def rank_modp(matrix,prime):
    rows=[[v%prime for v in row] for row in matrix]
    columns=len(rows[0]);rank=0
    for j in range(columns):
        i=next((i for i in range(rank,len(rows)) if rows[i][j]),None)
        if i is None:continue
        rows[rank],rows[i]=rows[i],rows[rank]
        inverse=pow(rows[rank][j],-1,prime)
        rows[rank]=[v*inverse%prime for v in rows[rank]]
        for i in range(rank+1,len(rows)):
            factor=rows[i][j]
            if factor:rows[i]=[(v-factor*w)%prime for v,w in zip(rows[i],rows[rank])]
        rank+=1
    return rank

def full_filled_chain(images,relators):
    """Full lifted directed edges; fill every generator cycle, including fixed sheets."""
    n=len(images[0]);m=len(images);parent=list(range(n));tree=set()
    assert all(sorted(p)==list(range(n)) and all(p[p[v]]==v for v in range(n)) for p in images)
    def root(v):
        while parent[v]!=v:v=parent[v]
        return v
    for g,p in enumerate(images):
        for v,w in enumerate(p):
            if root(v)!=root(w):parent[root(v)]=root(w);tree.add(g*n+v)
    assert len(tree)==n-1
    chords=[i for i in range(n*m) if i not in tree];rows=[]
    for word in relators:
        for start in range(n):
            row=[0]*(n*m);v=start
            for g in word:row[g*n+v]+=1;v=images[g][v]
            assert v==start
            rows.append([row[j] for j in chords])
    fill_count=0
    for g,p in enumerate(images):
        unused=set(range(n))
        while unused:
            v=start=min(unused);row=[0]*(n*m)
            while True:
                unused.remove(v);row[g*n+v]+=1;v=p[v]
                if v==start:break
            rows.append([row[j] for j in chords]);fill_count+=1
    return rows,len(chords),fill_count

def bound_seed_sources():
    """Rebind source receipts, literal names, four uppers and all101 scalar chains."""
    files=set();seen=set();seeds=[];controls=[]
    specs=[("seven_sheet_cover",3,7,2,1,24),
           ("a6_affine",3,6,2,2,24),("a7_binary",2,7,2,1,3),
           ("sl3_projective_family",3,13,4,2,50)]
    for stem,prime,n,delta,expected_h,count in specs:
        source=stem+"_inputs.json";receipt=stem+"_verification.json"
        data=json.loads((HERE/source).read_text());verified=json.loads((HERE/receipt).read_text())
        assert verified["status"]=="PASS";files.update([source,receipt])
        for name,pin in verified["input_sha256"].items():
            assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==pin,name
            files.add(name)
        reported={r["name"]:r for r in verified.get("results",verified.get("cases",[]))}
        if "--census" in sys.argv:
            from importlib.metadata import version
            import snappy_15_knots
            from spherogram import Link
            provenance=data["census_provenance"];package=Path(snappy_15_knots.__file__).parent
            assert version("snappy_15_knots")==provenance["version"] and version("spherogram")==provenance["spherogram_version"]
            for pin in provenance["source_files"]:
                raw=(package/pin["path"]).read_bytes()
                assert len(raw)==pin["bytes"] and hashlib.sha256(raw).hexdigest()==pin["sha256"]
        cases=data["cases"]
        if stem=="seven_sheet_cover":
            cases=[r for r in cases if next(a for a in r["actions"] if a["module"]=="natural")["expected_H1_F3_dimension"]==1]
            assert verified["sharp_four_count"]==24
        elif stem=="a6_affine":
            assert verified["optional_cover_homology_replayed"] and verified["optional_named_census_reread"]
        elif stem=="a7_binary":
            assert verified["optional_full_cover_homology_checked"] and verified["optional_named_census_reread_count"]==3
        assert len(cases)==count
        for row in cases:
            name=row["name"];assert name not in seen;seen.add(name);claim=reported[name]
            if "--census" in sys.argv:
                link=Link("K"+name.replace("_",""))
                assert [list(r) for r in link.PD_code()]==row["census_pd"]
                assert json.loads(json.dumps(link.DT_code()))==row["census_DT"]
            crossings,_=diagram(row["pd"]);census_bijection(row["pd"],row["census_pd"],row["census_map_witness"])
            known=set(row["seed_arcs"]);assert len(known)==4
            for o,a,b in row["seed_trace"]:
                assert o in known and a in known and b not in known
                assert (o,a,b) in crossings or (o,b,a) in crossings;known.add(b)
            assert known==set(range(13))
            if stem=="seven_sheet_cover":
                images=binary3_permutations(next(a for a in row["actions"] if a["module"]=="natural")["meridian_matrices_F2"])
                assert claim["sharp_four"] and claim["bridge_upper"]==4 and claim["H1_F3_dimensions"]==[1,1]
            elif stem=="a6_affine":
                images=[tuple(p) for p in row["permutations"]]
                assert claim["bridge_upper"]==4 and claim["independent_cover_homology"]["H1_F3_dimension"]==2
            elif stem=="a7_binary":
                images=[tuple(p) for p in row["permutations"]]
                assert claim["bridge_upper"]==4 and claim["optional_full_cover_F2"]["H1_F2_dimension"]==1
            else:
                images=[tuple(p) for p in row["projective_permutations"]]
                assert claim["four_seed_bridge_upper"]==4 and claim["H1_F3_dimension"]==2
            assert len(images)==13 and len(images[0])==n
            assert sum(images[0][v]!=v for v in range(n))==2*delta
            for o,a,b in crossings:assert tuple(images[o][images[a][images[o][v]]] for v in range(n))==images[b]
            relators=[[o,a,o,b] for o,a,b in crossings]
            matrix,columns,fills=full_filled_chain(images,relators)
            rank=rank_modp(matrix,prime);h=columns-rank
            assert h==expected_h and Fraction(n-1+h,delta)==Fraction(7,2)
            seed={"name":name,"source":source,"field_characteristic":prime,"degree":n,"delta":delta,
                  "H1_dimension":h,"full_chain_columns":columns,"full_chain_rank":rank,
                  "meridian_cycle_fill_rows":fills,"scalar_score":"7/2","bridge_upper":4,
                  "trefoil_family":"mr(K#^s T(2,3))=b(K#^s T(2,3))=4+s for every integer s>=0",
                  "finite_score_controls":[]}
            for t in range(4):
                degree=n+delta*t;score=Fraction(degree-1+h,delta)
                assert score==Fraction(7,2)+t and (score.numerator+score.denominator-1)//score.denominator==4+t
                seed["finite_score_controls"].append({"trefoil_factors":t,"constructed_degree":degree,
                        "constructed_H1_dimension":h,"score":str(score),"sharp_lower_and_upper":4+t})
            seeds.append(seed)
        # An actual all-edge amalgam check on one representative from EACH source,
        # using the trefoil three-sheet action needed also in characteristic two.
        row=cases[0]
        if stem=="seven_sheet_cover":images=binary3_permutations(next(a for a in row["actions"] if a["module"]=="natural")["meridian_matrices_F2"])
        else:images=[tuple(p) for p in row.get("permutations",row.get("projective_permutations"))]
        crossings,_=diagram(row["pd"]);words=[[o,a,o,b] for o,a,b in crossings]
        moving=[(v,images[0][v]) for v in range(n) if v<images[0][v]]
        assert len(moving)==delta
        extended=[tuple(list(p)+list(range(n,n+delta))) for p in images]
        b=list(range(n+delta))
        for j,(v,w) in enumerate(moving):b[w],b[n+j]=n+j,w
        extended.append(tuple(b))
        amalgam_words=words+[[0,13]*3]
        matrix,columns,fills=full_filled_chain(extended,amalgam_words)
        h=columns-rank_modp(matrix,prime);assert h==expected_h
        assert Fraction(n+delta-1+h,delta)==Fraction(9,2)
        controls.append({"source":source,"seed":row["name"],"trefoil_factors":1,"field_characteristic":prime,
                         "degree":n+delta,"delta":delta,"H1_dimension":h,"score":"9/2","sharp_lower_and_upper":5})
    assert len(seeds)==101
    return seeds,controls,files

def trefoil_all_field_control():
    # x=(01), y=(12) on three sheets. This presentation follows from
    # the trefoil relation xyx=yxy after valid meridian-square fillings.
    cover=covering_matrix([(1,0,2),(0,2,1)],[[1,2]*3,[1,1],[2,2]],[1])
    columns=len(cover["retained_edges"]);assert cover["degree"]==3 and cover["delta"]==1 and columns==4
    witness=None
    for indices in combinations(range(len(cover["matrix"])),columns):
        value=determinant([cover["matrix"][i] for i in indices])
        if abs(value)==1:witness={"row_indices":list(indices),"determinant":value};break
    assert witness is not None
    assert all(columns-rank_modp(cover["matrix"],p)==0 for p in [2,3,5,7,11])
    return {"degree":3,"delta":1,"matrix_columns":columns,"H1_over_every_field_dimension":0,
            "unit_maximal_minor":witness,"scalar_score_over_every_field":2}

def main():
    scaled = [bipartite_control(a,b) for a in range(1,9) for b in range(1,9)]
    permutations,relations = seed_data()
    sums = []
    for trefoils in range(4):
        images = list(permutations)
        relators = [list(w) for w in relations]
        for _ in range(trefoils):
            images.append(permutations[0])
            b = len(images)
            relators.append([1,b]*3)
            relators.append([b,b])
        data = covering_matrix(images,relators,[1])
        columns = len(data['retained_edges'])
        rank = rank_mod3(data['matrix'])
        dimension = columns-rank
        assert data['degree']==7 and data['delta']==2
        assert dimension == 2*trefoils+1
        bound = (6+dimension+1)//2
        assert bound == 4+trefoils
        sums.append({'seed':'13a_373','trefoil_factors':trefoils,
                     'orbifold_generators':len(images),
                     'degree':7,'delta':2,
                     'matrix_rows':len(data['matrix']),'matrix_columns':columns,
                     'H1_F3_dimension':dimension,'meridional_lower_bound':bound,
                     'bridge_upper_bound':4+trefoils})
    # A single trefoil's double cover is checked from its own presentation.
    double = covering_matrix([(1,0),(1,0)],[[1,2]*3,[1,1],[2,2]],[1])
    assert len(double['retained_edges'])-rank_mod3(double['matrix']) == 1
    # Corrupt the trefoil factor into the unknot: no false +2 Betti gain.
    bad = covering_matrix(permutations+[permutations[0]],
                          relations+[[1,-14],[14,14]],[1])
    assert len(bad['retained_edges'])-rank_mod3(bad['matrix']) == 1
    seeds,new_amalgams,source_files=bound_seed_sources()
    trefoil=trefoil_all_field_control()
    inputs=sorted(source_files|{'verify_involution_cover_sums.py','INVOLUTION_COVER_SUM_FORMULAS.tex',
                              'verify_general_branched_covers.py','verify_seven_sheet_covers.py'})
    receipt = {'status':'PASS','bipartite_scaling_controls':scaled,
               'literal_trefoil_amalgam_controls':sums,
               'all_field_three_sheet_trefoil':trefoil,
               'fractional_seed_count':101,'fractional_seed_families':seeds,
               'current_named_census_reread_count':101 if '--census' in sys.argv else 0,
               'one_trefoil_actual_chain_controls_by_source':new_amalgams,
               'negative_control':'Replacing trefoil by unknot gives no Betti gain',
               'scope':'101 named fractional-score seed families over their specified fixed fields; all101 full scalar homology chains and four-seed upper traces replayed, plus64 finite pair controls, four old literal0..3trefoil chains and four actual three-sheet-trefoil amalgams. The uniform integer parameter and exact normalized additivity use the written proof; no equality of seed suprema with7/2 is asserted.',
               'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}}

    (HERE/'involution_cover_sum_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','pair_matching_controls':64,
                      'trefoil_amalgam_Betti_numbers':[x['H1_F3_dimension'] for x in sums],
                      'sharp_bounds':[x['meridional_lower_bound'] for x in sums],
                      'unknot_negative_control':True,'fractional_seed_families':101,
                      'three_sheet_trefoil_all_fields':True,'actual_new_amalgam_controls':len(new_amalgams)},indent=2))


if __name__=='__main__':
    main()
