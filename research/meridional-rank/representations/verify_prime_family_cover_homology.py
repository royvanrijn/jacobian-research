#!/usr/bin/env python3
"""Exact integer homology proof; --topology binds unsimplified marked words."""
from pathlib import Path
from copy import deepcopy
import argparse, hashlib, json, sys
HERE=Path(__file__).resolve().parent
ROOT=HERE
INPUT=HERE/'prime_family_cover_homology_inputs.json'
RECEIPT=HERE/'prime_family_cover_homology_verification.json'
def mul(A,B):
 return [[sum(a*b for a,b in zip(row,col))for col in zip(*B)]for row in A]
def exponents(word,generators):
 assert len(set(generators))==len(generators) and all(len(g)==1 and g.islower()for g in generators)
 out=[0]*len(generators);ix={g:i for i,g in enumerate(generators)}
 for c in word:
  assert c.lower()in ix and c.isalpha() and c.isascii()
  out[ix[c.lower()]]+=1 if c.islower()else-1
 return out

def verify(d):
 k=len(d['generators']);A=[exponents(w,d['generators'])for w in d['relators']]
 assert len(A)==k-1 and A==d['relation_exponent_matrix']
 Q=d['unimodular_completion_Q'];V=d['smith_column_change_V'];I=[[int(i==j)for j in range(k)]for i in range(k)]
 assert len(Q)==len(V)==k and all(len(r)==k and all(type(x)is int for x in r)for r in Q+V)
 assert Q==A+[d['completion_last_row']]
 assert exponents(d['completion_generator_word'],d['generators'])==Q[-1]
 assert mul(Q,V)==mul(V,Q)==I
 assert mul(A,V)==d['smith_matrix_A_V']==I[:-1]
 h=[r[-1]for r in V];assert h==d['generator_homology_images']
 # Q and V are explicit inverse integer basis changes: the relations are the
 # first k-1 basis rows. Their quotient is exactly Z, without hidden torsion.
 per=[sum(x*y for x,y in zip(exponents(d['peripheral_words'][w],d['generators']),h))for w in ['meridian','longitude']]
 assert per==d['peripheral_homology_images']==[5,61]
 f=d['all_parameter_filling'];assert f['constant_coefficients']==[-1,2] and f['n_coefficients']==[0,6]
 base=sum(x*y for x,y in zip(f['constant_coefficients'],per));step=sum(x*y for x,y in zip(f['n_coefficients'],per))
 assert(base,step)==(f['homology_image_constant'],f['homology_image_n_coefficient'])==(117,366)
 c=d['determinant_only_positive_control'];n=c['n'];order=base+n*step
 assert n==3 and order==c['order']==1215 and c['prime_cube_divides']is True and order%c['prime']**3==0
 prod=1
 for p,e in c['factorization']:prod*=p**e
 assert prod==order
 return {'integer_unimodular_inverse_both_products_identity':True,'relation_Smith_form':'[I_10 | 0]','H1_parent':'Z','homology_generator_word':d['completion_generator_word'],'generator_images':h,'meridian_image':per[0],'longitude_image':per[1],'all_n_integral_filling_homology':'Z/(117+366n), n>=0','determinant_prime_cube_positive_control':c,'nonunit_8plat_gcd_contradiction':'For any prime p dividing all plat entries, Fox dimension4 would require H1 cover Fp dimension3; a cyclic group has dimension at most1.'}

def native(d):
 import snappy
 assert snappy.__version__==d['presentation_import']['snappy_version']=='3.3.2'
 kernel_sha=hashlib.sha256(Path(snappy.SnapPy.__file__).read_bytes()).hexdigest()
 assert kernel_sha in d['presentation_import']['allowed_kernel_sha256']
 M=snappy.Triangulation(str(ROOT/d['parent']['path']));assert M.num_tetrahedra()==14 and M.num_cusps()==1
 assert all(tuple(p)==(0,0)for p in M.cusp_info('filling'))
 G=M.fundamental_group(False,False,False)
 assert G.generators()==d['generators'] and G.relators()==d['relators']
 assert G.peripheral_curves()==[(d['peripheral_words']['meridian'],d['peripheral_words']['longitude'])]
 return {'fresh_exact_unsimplified_presentation_and_marked_peripheral_words':True,'kernel_sha256':kernel_sha,'snappy_version':snappy.__version__,'numerical_geometry_used':False}

def controls(d):
 tests=[('relator_word',lambda x:x['relators'].__setitem__(0,'aG')),('relation_matrix',lambda x:x['relation_exponent_matrix'][0].__setitem__(0,2)),('integer_inverse',lambda x:x['smith_column_change_V'][0].__setitem__(0,x['smith_column_change_V'][0][0]+1)),('homology_image',lambda x:x['generator_homology_images'].__setitem__(0,31)),('peripheral_orientation',lambda x:x['peripheral_words'].__setitem__('meridian',x['peripheral_words']['meridian'].swapcase())),('parameter_direction',lambda x:x['all_parameter_filling'].__setitem__('n_coefficients',[0,-6])),('false_determinant_only_exclusion',lambda x:x['determinant_only_positive_control'].__setitem__('prime_cube_divides',False))]
 rejected=[]
 for name,mut in tests:
  bad=deepcopy(d);mut(bad)
  try:verify(bad)
  except (AssertionError,KeyError,IndexError,ValueError):rejected.append(name)
  else:raise AssertionError('Accepted corrupt control '+name)
 return rejected

def run_verification(topology=False):
 d=json.loads(INPUT.read_text())
 parent=ROOT/d['parent']['path'];assert hashlib.sha256(parent.read_bytes()).hexdigest()==d['parent']['sha256']
 for name,expected in d['source_sha256'].items():
  assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected,name
 names=[INPUT.name,Path(__file__).name,d['parent']['path']]+list(d['source_sha256'])
 if 'written_proof' in d:names.append(d['written_proof']['path'])
 edition={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()for name in dict.fromkeys(names)}
 if 'written_proof' in d:assert edition[d['written_proof']['path']]==d['written_proof']['sha256']
 result={'status':'PASS exact integral homology certificate','edition_sha256':edition,'integer_verification':verify(d),'rejected_controls':controls(d),'trust_boundary':['Native exact cell presentation/peripheral interpretation is an imported software contract.','Actual branched-cover identification and all-n marked filling slopes use the separate canonical written binding theorem.','Mayer-Vietoris filling quotient and Fox/core homology deduction are written mathematical arguments.'],'assurance':{'formal_verification':False,'external_specialist_review':False,'whole_proof_independent_replay':False}}
 result['mode_replays']={'standard_library':result['integer_verification']}
 if topology:result['mode_replays']['topology']=native(d)
 result['this_invocation']={'topology':topology}
 return result

def save(out):
 old=json.loads(RECEIPT.read_text())if RECEIPT.exists()else{}
 for key in ['internal_pdf_reviews','review_history','edition_replays','root_review','independent_replays']:
  if key in old:out[key]=old[key]
 if old.get('edition_sha256')==out['edition_sha256']:
  out['mode_replays']={**old.get('mode_replays',{}),**out['mode_replays']}
 RECEIPT.write_text(json.dumps(out,indent=2)+'\n')

def main():
 p=argparse.ArgumentParser();p.add_argument('--topology',action='store_true');a=p.parse_args()
 result=run_verification(topology=a.topology);save(result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
