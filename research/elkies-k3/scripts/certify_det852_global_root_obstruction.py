"""Close the sole remaining determinant852 Niemeier case by exact two-vector enumeration."""
from sage.all import *
from pathlib import Path
import json,itertools,hashlib,argparse
R=Path(__file__).resolve().parents[2];P=R/'artifacts/generated-results/elkies-k3-det852-final-anchor-v1'
source=R/'artifacts/generated-results/elkies-k3-det852-root-gate-v1/certificate.json'
catpath=R/'artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json';ancpath=R/'artifacts/generated-results/elkies-k3-niemeier-d5-anchor-orbits.json'
cert=json.loads(source.read_text());assert cert['status']=='PASS_FIFTEEN_ANCHORS_EXCLUDED_ONE_UNRESOLVED'
K0=matrix(ZZ,cert['auxiliary_gram']);J=matrix(ZZ,cert['D5_basis_change']);B=block_diagonal_matrix(J.transpose(),identity_matrix(ZZ,2));K=B*K0*B.transpose();assert K[:5,:5]==matrix(ZZ,CartanMatrix(['D',5]))
Schur=K[5:,5:]-K[5:,:5]*K[:5,:5].inverse()*K[:5,5:];assert Schur==matrix(QQ,[[QQ(51)/4,QQ(3)/4],[QQ(3)/4,QQ(67)/4]])
cat=json.loads(catpath.read_text());anc=json.loads(ancpath.read_text());G=matrix(ZZ,next(x['gram'] for x in cat['rooted_niemeier_lattices'] if x['label']=='2A7_2D5'));anchor=matrix(ZZ,next(x['D5_basis_in_ambient'] for x in anc['anchors'] if x['niemeier']=='2A7_2D5'))
def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def signed_roots(gram):
    minimum_data = pari(gram).qfminim(2)
    representatives = matrix(ZZ, minimum_data[2].sage()).transpose()
    result = [vector(ZZ, row) for row in representatives.rows()]
    return result + [-root for root in result]


def root_components(gram, roots):
    unseen = set(range(len(roots)))
    result = []
    while unseen:
        component = {min(unseen)}
        frontier = list(component)
        unseen.difference_update(component)
        while frontier:
            current = frontier.pop()
            neighbours = {
                index
                for index in unseen
                if roots[current] * gram * roots[index] != 0
            }
            component.update(neighbours)
            unseen.difference_update(neighbours)
            frontier.extend(sorted(neighbours))
        result.append([roots[index] for index in sorted(component)])
    return result


def simple_roots(gram, roots):
    for trial in range(1, 100):
        chamber_vector = vector(
            ZZ,
            [
                (index + 1) ** 2 + trial * (index + 1) + trial**2
                for index in range(gram.nrows())
            ],
        )
        values = [root * gram * chamber_vector for root in roots]
        if all(value != 0 for value in values):
            break
    else:
        raise AssertionError("failed to choose a regular chamber vector")
    positive = [root for root, value in zip(roots, values) if value > 0]
    positive_set = {tuple(map(int, root)) for root in positive}
    result = matrix(
        ZZ,
        [
            root
            for root in positive
            if not any(
                tuple(map(int, root - other)) in positive_set
                for other in positive
            )
        ],
    )
    assert result.rank() == result.nrows()
    return result


def dominant_labels_up_to(cartan_inverse, bound):
    rank = cartan_inverse.nrows()
    current = [ZZ(0)] * rank
    result = []

    def extend(index, norm):
        if index == rank:
            result.append((tuple(map(int, current)), norm))
            return
        cross = sum(
            current[previous] * cartan_inverse[previous, index]
            for previous in range(index)
        )
        coefficient = 0
        while True:
            new_norm = (
                norm
                + 2 * coefficient * cross
                + coefficient * coefficient * cartan_inverse[index, index]
            )
            if new_norm > bound:
                break
            current[index] = coefficient
            extend(index + 1, new_norm)
            coefficient += 1
        current[index] = 0

    extend(0, QQ(0))
    return result



z=G.__pari__().qfminim(2,100000,2);W=matrix(ZZ,z[2]).transpose();assert 2*W.nrows()==int(z[0]);ambient_roots=list(W.rows())+list((-W).rows());allcases=[]
for orientation in [[0,1,2,3,4],[0,1,2,4,3]]:
 A=anchor.matrix_from_rows(orientation);assert A*G*A.transpose()==K[:5,:5]
 orth=[r for r in ambient_roots if not any(r*G*A.transpose())];comps=root_components(G,orth);simples=[simple_roots(G,c) for c in comps];simples.sort(key=lambda M:-M.nrows());D=A.stack(matrix(ZZ,[r for C in simples for r in C.rows()]));inv=(G*D.transpose()).inverse();assert D.rank()==24
 Cs=[C*G*C.transpose() for C in simples];Is=[C.inverse() for C in Cs];assert [C.nrows() for C in Cs]==[7,7,5]
 assert block_diagonal_matrix(K[:5,:5],*Cs)==D*G*D.transpose()
 bound=QQ(51)/4;vbound=QQ(67)/4;choices=[]
 for C,I in zip(Cs,Is):
  allowed=[]
  for lab,norm in dominant_labels_up_to(I,bound):
   zero=[i for i,z in enumerate(lab) if not z];need=sum(C.matrix_from_rows_and_columns(zero,zero).inverse().list()) if zero else 0
   if need<=vbound:allowed.append((lab,norm,need))
  choices.append(allowed)
 sixth=[];raw=0
 for choice in itertools.product(*choices):
  if sum(x[1] for x in choice)!=bound or sum(x[2] for x in choice)>vbound:continue
  raw+=1;labels=[z for c in choice for z in c[0]];u=vector(QQ,list(K.column(5)[:5])+labels)*inv
  if not all(z in ZZ for z in u):continue
  assert u*G*u==14 and u*G*A.transpose()==vector(QQ,K.column(5)[:5])
  sixth.append((u,choice))
 assert raw==52 and len(sixth)==13
 results=[];cache={}
 for u,choice in sixth:
  needs=[c[2] for c in choice];seventh=[]
  for j,(I,uc) in enumerate(zip(Is,choice)):
   b=vbound-sum(needs)+needs[j];den=I.denominator();key=(j,b)
   if key not in cache:
    z=(den*I).__pari__().qfminim(floor(den*b),1000000,2);V=matrix(ZZ,z[2]).transpose();assert 2*V.nrows()==int(z[0]);cache[key]=list(V.rows())+list((-V).rows())
   lab=uc[0];zero=[i for i,z in enumerate(lab) if z==0];assert zero;allowed=[]
   for v in cache[key]:
    if any(v[i]<=0 for i in zero):continue
    allowed.append((v,v*I*v,vector(ZZ,lab)*I*v))
   seventh.append(allowed)
  keys={(n,dot) for v,n,dot in seventh[2]};matches=sum((vbound-n1-n2,QQ(3)/4-d1-d2) in keys for v1,n1,d1 in seventh[0] for v2,n2,d2 in seventh[1]);assert matches==0
  results.append({'sixth_vector':list(map(int,u)),'sixth_component_labels':[list(c[0]) for c in choice],'component_norms':[str(c[1]) for c in choice],'seventh_regular_lower_bounds':[str(c[2]) for c in choice],'seventh_component_candidates':[[{'labels':list(map(int,v)),'norm':str(n),'pairing':str(dot)} for v,n,dot in vs] for vs in seventh],'norm_and_pairing_matches':matches})
 allcases.append({'D5_orientation':orientation,'D5_basis':rows(A),'residual_simple_bases':[rows(C) for C in simples],'component_Cartans':[rows(C) for C in Cs],'pre_integrality_sixth_count':raw,'integral_sixth_count':len(sixth),'sixth_cases':results})
out={'schema':'elkies-k3.det852-global-root-obstruction.v1','status':'PASS_ALL_FRAMES_ROOTFUL','surface_id':'K3-4ff75fec54d01662','NS_determinant':852,'auxiliary_Cartan_gram':rows(K),'remaining_Niemeier':'2A7_2D5','orientation_cases':allcases,'conclusion':'Every frame of the admitted determinant852 NS contains a root. No MW17 fibration exists on this NS. Full saturated rational marking remains proved.','inputs':[{'path':str(p.relative_to(R)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [source,catpath,ancpath]],'theorem_inputs':['Previously certified primitive opposite-form glue and fifteen-anchor exclusion','Complete Niemeier D5 anchor theorem','D5 root-lattice automorphisms are Weyl group extended by its order-two diagram symmetry','Residual Weyl chamber reduction; roots orthogonal to a dominant vector form its zero-node parabolic subsystem','Complete exact positive-definite quadratic-vector enumeration with stored-count equality checks'],'boundary':'Only this actual NS genus; no determinant-wide exclusion, optimal lower MW rank, independent implementation, formal verification, external review or novelty claim.'}
a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args();output=P/'certificate.json'
if args.check:assert json.loads(output.read_text())==out
else:output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('PASS852 globally rootful: both D5 orientations,13 sixth vectors each, no norm-and-pairing-compatible root-avoiding seventh vector.')
