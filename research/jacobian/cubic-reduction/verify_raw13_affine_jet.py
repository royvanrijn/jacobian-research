#!/usr/bin/env python3
"""Independent Fraction translated-jet/minor audit for the entire raw13 affine orbit."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,time
from verify import decode,variable,add,scale,mul,compose,derivative,determinant3,evaluate
B=Path(__file__).resolve().parent;O=B/'raw13-affine-jet';ROOT=B.parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def con(c,n):return {(0,)*n:Q(c)} if c else {}
def total(*ps):
 z={}
 for p in ps:z=add(z,p)
 return z

def pw(p,k,n):
 z=con(1,n)
 for _ in range(k):z=mul(z,p)
 return z

def scalar(p,j,v):
 z={}
 for e,c in p.items():
  f=e[:j]+e[j+1:];z[f]=z.get(f,Q(0))+c*Q(v)**e[j]
 return {e:c for e,c in z.items() if c}
def minor(M,rows,cols,n):
 if len(rows)==3:return determinant3([[M[i][j] for j in cols] for i in rows],n)
 assert len(rows)==len(cols)==2
 return add(mul(M[rows[0]][cols[0]],M[rows[1]][cols[1]]),mul(M[rows[0]][cols[1]],M[rows[1]][cols[0]]),Q(-1))
def main():
 started=time.monotonic();cp=O/'checkpoint.json';C=json.loads(cp.read_text());assert C['status']=='complete' and C['source_sha256']==sha(B/'raw13_affine_jet.py')
 ip=O/'input.json';D=json.loads(ip.read_text());raw=B/'raw13-specialization/input.json';assert C['input_sha256']==sha(ip) and D['source_sha256']==C['source_sha256'] and D['input_sha256']==C['source_input_sha256']==sha(raw) and C['discovery_sha256']==sha(O/'discovery.json')
 disc=json.loads((O/'discovery.json').read_text());assert len(disc['primary_queries'])==C['preliminary_primary_queries']==3 and all(z['result']=='unit ideal [1]' and not z['used_as_proof'] for z in disc['primary_queries']) and C['new_solver_queries']==0
 rawdata=json.loads(raw.read_text());n=3;F=[variable(j,n) for j in range(n)]
 for atom in rawdata['atoms']:F[atom['row']]=add(F[atom['row']],{tuple(atom['exponent']):Q(atom['coefficient'])})
 assert determinant3([[derivative(z,j) for j in range(3)] for z in F],n)==con(1,n) and [sum(sum(e)>=2 for e in z) for z in F]==[2,5,6]
 assert max(sum(e) for z in F for e in z)==7 and [{e:c for e,c in z.items() if sum(e)<=1} for z in F]==[variable(j,n) for j in range(3)]
 weights=[1,-1,-2];assert all(sum(a*b for a,b in zip(e,weights))==w for z,w in zip(F,weights) for e in z)
 points=[[0,1,-8],[-1,1,4]];assert [[evaluate(z,t) for z in F] for t in points]==[[Q(0),Q(1),Q(0)]]*2
 # Independent production route: full6-variable translation and homogeneous collection,
 # rather than the producer's binomial coefficient formula for each matrix entry.
 n=6;xs=[variable(j,n) for j in range(3)];ts=[variable(j,n) for j in range(3,6)];FF=[compose(z,[add(x,t) for x,t in zip(xs,ts)],n) for z in F];matrices={};exponents={};coefficient_entries=0
 for degree in range(2,8):
  rec=D['degrees'][str(degree)];es=sorted({(i,j,degree-i-j) for i in range(degree+1) for j in range(degree-i+1)},reverse=True);assert list(map(list,es))==rec['monomials'];M=[]
  for i in range(3):
   row=[]
   for e in es:
    z={ee[3:]:c for ee,c in FF[i].items() if ee[:3]==e};row.append(z)
    assert z==decode(rec['matrix'][i][len(row)-1],3)
    assert all(sum(a*b for a,b in zip(ee,weights))+sum(a*b for a,b in zip(e,weights))==weights[i] for ee in z)
    coefficient_entries+=1
   M.append(row)
  matrices[degree]=M;exponents[degree]={e:i for i,e in enumerate(es)}
 def getminor(degree,rows,mons,a_value=None):
  M=matrices[degree];cols=[exponents[degree][tuple(e)] for e in mons];z=minor(M,rows,cols,3)
  return scalar(z,0,a_value) if a_value is not None else z
 a,b,c=[variable(j,3) for j in range(3)]
 # Exact global high-degree ranks with all a=0/a!=0 branches explicit.
 assert all(not z for row in matrices[7][:2] for z in row)
 assert matrices[7][2][exponents[7][(3,3,1)]]==con(8,3) and sum(bool(z) for z in matrices[7][2])==1
 for degree in [5,6]:assert not any(matrices[degree][0])
 assert getminor(6,[1,2],[(3,2,1),(2,4,0)])==con(288,3)
 assert getminor(5,[1,2],[(2,3,0),(1,4,0)])==scale(a,1728)
 assert getminor(5,[1,2],[(2,3,0),(2,2,1)],0)==con(432,2)
 assert getminor(4,[0,1,2],[(3,0,1),(1,3,0),(0,4,0)])==scale(pw(a,3,3),-1728)
 assert getminor(4,[0,1,2],[(3,0,1),(2,1,1),(1,3,0)],0)==con(-336,2)
 assert getminor(3,[0,1,2],[(2,1,0),(1,2,0),(1,1,1)],0)==con(-432,2)
 assert getminor(2,[0,1,2],[(2,0,0),(1,0,1),(0,2,0)],0)==scale(variable(0,2),-72)
 # Two exact selected-minor units close every normalized a!=0 point,including all parameter-zero strata.
 units={};terms=0;maxdegree=0
 for key in ['degree2','degree3']:
  degree=int(key[-1]);rec=D['selected_minor_units'][key];M=[[scalar(z,0,1) for z in row] for row in matrices[degree]];out={};local=0
  for cols,erow,mrow in zip(rec['columns'],rec['minors'],rec['multipliers']):
   g=minor(M,[0,1,2],cols,2);assert g==decode(erow,2);m=decode(mrow,2);out=add(out,mul(m,g));local+=len(m);maxdegree=max(maxdegree,max((sum(e) for e in m),default=0))
  assert out==con(1,2);units[key]={'generators':len(rec['columns']),'multiplier_terms':local,'identity_one':True};terms+=local
 assert terms==C['unit_multiplier_terms']==31 and {k:v['multiplier_terms'] for k,v in units.items()}==C['selected_minor_unit_counts']
 # a=b=0 quadratic rank is exactly2; off this axis it is3.
 axis={degree:[[scalar(scalar(z,0,0),0,0) for z in row] for row in M] for degree,M in matrices.items()}
 assert not any(axis[2][0])
 assert minor(axis[2],[1,2],[exponents[2][(1,0,1)],exponents[2][(0,2,0)]],1)==con(24,1)
 rank_axis=[2,3,3,2,2,1];rank_off_axis=[3,3,3,2,2,1];assert sum(rank_axis)==13 and sum(rank_off_axis)==14
 # Every nonzero c-axis translation has exactly18 nonlinear source columns;
 # each is nonzero for every c!=0 by an explicit constant or pure-linear-c entry.
 cols={};new=[];old={tuple(z['exponent']) for z in rawdata['atoms']}
 for degree,M in axis.items():
  for e,j in exponents[degree].items():
   vector=[M[i][j] for i in range(3)]
   if any(vector):
    cols[e]=vector
    assert any(z==con(z.get((0,),0),1) and z for z in vector) or any(len(z)==1 and (1,) in z for z in vector)
    if e not in old:new.append(e)
 assert len(cols)==18 and set(new)=={(3,0,0),(3,2,0),(1,1,0),(2,2,0),(3,3,0)} and set(cols)==old|set(new)
 at_origin={e:[scalar(z,0,0) for z in v] for e,v in cols.items()};live={e:v for e,v in at_origin.items() if any(v)};assert len(live)==13 and all(sum(bool(z) for z in v)==1 for v in live.values()) and {next(i for i,z in enumerate(v) if z) for v in live.values()}=={0,1,2}
 ledger=json.loads((ROOT/'MATH_STATUS.json').read_text());E={e['id']:e for e in ledger['entries']};entry=E['KCUB40'];assert entry['state']=='proved';pins={'KCUB40':{'scope_fingerprint':hashlib.sha256((entry['state']+'\0'+entry['title']+'\0'+entry['scope']).encode()).hexdigest()[:16],'proof_sha256':sha(ROOT/entry['canonical_source'])}}
 receipt={'status':'pass','claim':'KCUB41','checker_sha256':sha(Path(__file__)),'producer_sha256':C['source_sha256'],'arithmetic_sha256':sha(B/'verify.py'),'input_sha256':sha(ip),'checkpoint_sha256':sha(cp),'discovery_sha256':sha(O/'discovery.json'),'complete_homogeneous_degrees':list(range(2,8)),'complete_jet_matrix_entries':coefficient_entries,'independent_full_translation_not_binomial_entry_production':True,'all_weight_covariance_monomials_checked':True,'selected_minor_units':units,'unit_multiplier_terms':terms,'max_multiplier_degree':maxdegree,'high_degree_and_zero_branch_minors_checked':True,'exact_rank_profile_on_c_axis':rank_axis,'exact_rank_profile_off_c_axis':rank_off_axis,'affine_nonlinear_occurrence_lower_bound':13,'off_axis_translation_lower_bound':14,'nonzero_c_axis_translation_source_columns':18,'extra_axis_columns':[list(e) for e in sorted(new)],'equality_source_centre_must_be_origin':True,'equality_source_and_target_linear_parts_must_be_monomial_matrices':True,'universal_rank_invariance_and_UFD_equality_proof_written_not_formalized':True,'imports':pins,'matching_old_upper':{'source_input_sha256':sha(raw),'counts':[2,5,6],'degree':7,'identity_linear_part':True,'full_determinant_one':True,'points':points,'common_image':[0,1,0],'new_witness':False},'corollary_minimum13_in_every_affine_frame_of_every_degree7_diagonal_weight_collision':True,'all_preliminary_solver_answers_unused_as_proof':True,'preliminary_primary_queries':3,'producer_and_checker_solver_queries':0,'goal_complete':False,'no_global_positive_bound_or_unrestricted_profile_change':True,'trust_boundary':'All complete translated coefficient matrices,selected minors,two rational unit identities and existing upper checked exactly and independently. Universal GL-invariance,centre normalization by weights,summed-rank occurrence bound,top-product UFD rigidity and equality argument are written,not formalized or externally reviewed. KCUB40 imported for the full degree7 diagonal-weight collision classification. This closes affine-coordinate sparsification of that entire classified orbit,not nonlinear coordinate changes,unweighted maps or global simplicity optima.','seconds':time.monotonic()-started}
 proof=B/'RAW13_AFFINE_ORBIT_MINIMUM.md'
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (O/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS full raw13 affine orbit:all348 translated jet entries,two exact minor units/31 multiplier terms;rank sums13/14;18 columns on nonzero c-axis;written equality rigidity;global goal open',flush=True)
if __name__=='__main__':main()
