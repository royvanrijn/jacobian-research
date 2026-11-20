#!/usr/bin/env python3
"""Direct exact proof checks for one explicit 48-variable cubic map."""
if not __debug__:
    raise SystemExit("Refusing optimized Python: witness checks are required.")
from argparse import ArgumentParser
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import copy,json,time
HERE=Path(__file__).resolve().parent
TOP=F(8803690708443294969,95329784685296)
DIMS=[1,1,1,1,2,3,4,6,9,10,4,3,1,1,1]
def require(test,message):
    if not test:raise ValueError(message)
def validate(data):
    require(data.get("schema_version")==1 and data.get("dimension")==48,"wrong witness schema/dimension")
    grades=data.get("coordinate_grades")
    require(isinstance(grades,list) and len(grades)==48,"wrong coordinate count")
    require(all(type(g)is int and 0<=g<=14 for g in grades),"invalid coordinate grade")
    require(grades==sorted(grades),"coordinate grades not ordered")
    require([grades.count(g) for g in range(15)]==DIMS,"wrong graded dimensions")
    weights=[2*g+1 for g in grades]
    require(all(0<w<=29 for w in weights) and max(weights)==29,"invalid positive weight bound")
    terms=data.get("H_terms");require(isinstance(terms,list) and len(terms)==1213,"wrong cubic monomial count")
    seen=set()
    for term in terms:
        require(isinstance(term,dict) and set(term)=={"output","inputs","coefficient"},"invalid monomial record")
        out=term["output"];inputs=term["inputs"]
        require(type(out)is int and 0<=out<48,"output index outside witness")
        require(isinstance(inputs,list) and len(inputs)==3,"noncubic monomial")
        require(all(type(i)is int and 0<=i<48 for i in inputs),"input index outside witness")
        require(inputs==sorted(inputs),"noncanonical monomial")
        key=(out,tuple(inputs));require(key not in seen,"duplicate monomial");seen.add(key)
        require(isinstance(term["coefficient"],str),"coefficient must be an exact rational string")
        coefficient=F(term["coefficient"]);require(coefficient!=0,"zero monomial coefficient")
        require(grades[out]==1+sum(grades[i] for i in inputs),"term violates strict grading")
        require(all(grades[i]<grades[out] for i in inputs),"term is not strictly triangular")
        require(weights[out]==sum(weights[i] for i in inputs),"term violates positive weights")
    return grades,weights

def add(poly,mon,c):
    value=poly.get(mon,F(0))+c
    if value:poly[mon]=value
    else:poly.pop(mon,None)
def jacobian(terms):
    matrix=defaultdict(dict)
    for term in terms:
        out=term["output"];inputs=term["inputs"];c=F(term["coefficient"])
        for i,multiplicity in Counter(inputs).items():
            mon=list(inputs);mon.remove(i)
            add(matrix[(out,i)],tuple(mon),c*multiplicity)
    return {key:poly for key,poly in matrix.items() if poly}
def multiply(A,B):
    rows=defaultdict(list)
    for (j,k),poly in B.items():rows[j].append((k,poly))
    output=defaultdict(dict);products=0
    for (i,j),left in A.items():
        for k,right in rows[j]:
            target=output[(i,k)]
            for mon,c in left.items():
                for other,d in right.items():
                    add(target,tuple(sorted(mon+other)),c*d);products+=1
    return {key:poly for key,poly in output.items() if poly},products
def point_inverse(data,grades):
    vectors=[{0:F(1)}]
    for grade in range(1,15):
        row={}
        for term in data["H_terms"]:
            if grades[term["output"]]!=grade:continue
            c=F(term["coefficient"])
            for i in term["inputs"]:c*=vectors[grades[i]].get(i,F(0))
            if c:row[term["output"]]=row.get(term["output"],F(0))+c
        vectors.append({i:c for i,c in row.items() if c})
    return vectors
def verify(data):
    grades,weights=validate(data)
    J=jacobian(data["H_terms"]);J2,p2=multiply(J,J);J3,p3=multiply(J2,J)
    require(not J3,"literal Jacobian cube is nonzero")
    point={(i,j):poly.get((0,0),F(0)) for (i,j),poly in J.items()}
    v1={i:point.get((i,0),F(0)) for i in range(48)}
    v2={i:sum((point.get((i,j),F(0))*v1[j] for j in range(48)),F(0)) for i in range(48)}
    v2={i:c for i,c in v2.items() if c}
    require(v2=={2:F(9)},"literal point Jacobian square differs from9e2")
    point1={(i,j):poly.get((1,1),F(0)) for (i,j),poly in J.items()}
    def apply(matrix,vector):
        return {i:sum((matrix.get((i,j),F(0))*c for j,c in vector.items()),F(0))
                for i in range(48)}
    mixed=apply(point,apply(point,apply(point1,{0:F(1)})))
    mixed={i:c for i,c in mixed.items() if c}
    require(mixed=={6:F(27)},"mixed-point cube control differs from27e6")
    require(J.get((1,0))=={(0,0):F(3)} and not J.get((0,1)),"nonsymmetry control failed")
    h2=[term for term in data["H_terms"] if term["output"]==2]
    require(h2==[{"output":2,"inputs":[0,0,1],"coefficient":"3"}],"mixed cubic component control failed")
    vectors=point_inverse(data,grades)
    require(vectors[14]=={47:TOP},"literal inverse top coefficient is wrong or zero")
    return {"dimension":48,"graded_dimensions":DIMS,"cubic_monomials":len(data["H_terms"]),
            "maximum_positive_weight":max(weights),"jacobian_entries":len(J),
            "jacobian_terms":sum(map(len,J.values())),"jacobian_square_entries":len(J2),
            "jacobian_square_terms":sum(map(len,J2.values())),"jacobian_cube_entries":len(J3),
            "exact_square_term_products":p2,"exact_cube_term_products":p3,
            "J_e0_squared_e0":{str(i):str(c) for i,c in v2.items()},
            "J_e0_J_e0_J_e1_e0":{str(i):str(c) for i,c in mixed.items()},
            "nonsymmetric_jacobian":True,"mixed_non_power_linear_component":"H2=3*y0^2*y1",
            "inverse_top_coefficient":str(TOP),"inverse_top_output":47,"inverse_top_monomial":"x0^29",
            "point_inverse_vectors":[{str(i):str(c) for i,c in v.items()} for v in vectors]}

def controls(data):
    # A valid coefficient perturbation, not merely malformed metadata,
    # must be detected by the actual all-variable polynomial cube.
    changed=copy.deepcopy(data)
    changed["H_terms"][0]["coefficient"]=str(F(changed["H_terms"][0]["coefficient"])+1)
    validate(changed)
    J=jacobian(changed["H_terms"]);square,_=multiply(J,J);cube,_=multiply(square,J)
    require(bool(cube),"coefficient perturbation did not break the literal Jacobian cube")
    rejected=0
    variants=[]
    bad=copy.deepcopy(data);bad["coordinate_grades"][0]=-1;variants.append(bad)
    bad=copy.deepcopy(data);bad["H_terms"][0]["inputs"][0]=48;variants.append(bad)
    bad=copy.deepcopy(data);bad["H_terms"][0]["coefficient"]="0";variants.append(bad)
    bad=copy.deepcopy(data);bad["H_terms"][1]=copy.deepcopy(bad["H_terms"][0]);variants.append(bad)
    for bad in variants:
        try:validate(bad)
        except(ValueError,ZeroDivisionError):rejected+=1
        else:raise ValueError("malformed witness accepted")
    sample=jacobian([{"output":2,"inputs":[0,0,1],"coefficient":"1"}])
    require(sample=={(2,0):{(0,1):F(2)},(2,1):{(0,0):F(1)}},"repeated-variable derivative control")
    product,_=multiply({(2,1):{(0,1):F(3)}},{(1,0):{(1,1):F(2)}})
    require(product=={(2,0):{(0,1,1,1):F(6)}},"commutative matrix-product control")
    return {"coefficient_perturbation_nonzero_cube_entries":len(cube),
            "malformed_witnesses_rejected":rejected,"repeated_variable_derivative":"PASS",
            "commutative_matrix_product":"PASS"}

def main():
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("--input",type=Path,default=HERE/"cubic_cube_zero_sharp48.json")
    ap.add_argument("--ledger",type=Path,default=HERE/"cubic_cube_zero_sharp48_source_ledger.json")
    ap.add_argument("--output","--write-receipt",dest="output",type=Path)
    args=ap.parse_args();start=time.time();ledger=json.loads(args.ledger.read_text())
    for item in ledger["inputs"]:
        path=HERE/item["path"]
        require(path.is_file(),"missing pinned source:"+item["path"])
        require(sha256(path.read_bytes()).hexdigest()==item["sha256"],"changed pinned source:"+item["path"])
    pin=next(item for item in ledger["inputs"] if item["role"]=="input")
    require(sha256(args.input.read_bytes()).hexdigest()==pin["sha256"],"selected input hash mismatch")
    data=json.loads(args.input.read_text());facts=verify(data);negative=controls(data)
    result={"status":"PASS","arithmetic":"exact QQ; standard library only","facts":facts,"controls":negative,
            "input_sha256":sha256(args.input.read_bytes()).hexdigest(),
            "checker_sha256":sha256(Path(__file__).read_bytes()).hexdigest(),
            "source_ledger_sha256":sha256(args.ledger.read_bytes()).hexdigest(),
            "seconds":time.time()-start}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":"PASS","dimension":48,"inverse_degree":29,
                      "nilpotency_index":3,"coefficient_corruption_cube_entries":negative["coefficient_perturbation_nonzero_cube_entries"],
                      "seconds":result["seconds"]}))
if __name__=="__main__":main()

