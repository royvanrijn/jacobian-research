#!/usr/bin/env python3
"""Independent full-translation and extremal-jet audit for the degree8 shear family."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,time
from verify import decode,variable,add,scale,mul,compose,derivative,determinant3,evaluate
B=Path(__file__).resolve().parent;O=B/'quadratic-target-affine8';ROOT=B.parents[1]
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
def hom(p,d):return {e:c for e,c in p.items() if sum(e[:3])==d}
def cf(p,e):return {ee[3:]:c for ee,c in p.items() if ee[:3]==tuple(e)}
def main():
 start=time.monotonic();cp=O/'checkpoint.json';C=json.loads(cp.read_text());ip=O/'input.json';D=json.loads(ip.read_text());raw=B/'raw13-specialization/input.json';assert C['status']=='complete' and C['source_sha256']==D['source_sha256']==sha(B/'quadratic_target_affine8.py') and C['input_sha256']==sha(ip) and C['rigidity_sha256']==sha(O/'rigidity.json') and D['source_input_sha256']==C['source_input_sha256']==sha(raw)
 n=3;F=[variable(j,n) for j in range(3)]
 for atom in json.loads(raw.read_text())['atoms']:F[atom['row']]=add(F[atom['row']],{tuple(atom['exponent']):Q(atom['coefficient'])})
 n=8;x,y,z,alpha,beta,gamma,L,M=[variable(j,n) for j in range(n)];HH=[compose(f,[add(x,alpha),add(y,beta),add(z,gamma)],n) for f in F];p0=hom(HH[0],0);P=add(HH[0],p0,Q(-1));square=pw(P,2,n);GG=[HH[0],add(HH[1],mul(L,square)),add(HH[2],mul(M,square))];original=[HH[0],add(HH[1],mul(L,pw(HH[0],2,n))),add(HH[2],mul(M,pw(HH[0],2,n)))]
 assert p0=={(0,0,0)+e:c for e,c in decode(D['first_constant'],5).items()}
 for i,u in enumerate([{},L,M]):assert original[i]==total(GG[i],scale(mul(mul(p0,u),GG[0]),2),scale(mul(u,pw(p0,2,n)),-1))
 assert all(e[2]<=1 and e[0]<=3 for f in HH for e in f)
 matrices={};base={};indices={};entries=0
 for degree in range(2,9):
  rec=D['degrees'][str(degree)];es=sorted({(i,j,degree-i-j) for i in range(degree+1) for j in range(degree-i+1)},reverse=True);assert list(map(list,es))==rec['monomials'];indices[degree]={e:i for i,e in enumerate(es)};matrices[degree]=[[cf(f,e) for e in es] for f in GG];base[degree]=[[cf(f,e) for e in es] for f in HH]
  assert matrices[degree]==[[decode(v,5) for v in row] for row in rec['matrix']] and base[degree]==[[decode(v,5) for v in row] for row in rec['base_matrix']]
  assert [cf(square,e) for e in es]==[decode(v,5) for v in rec['square_row']];entries+=3*len(es)
 assert entries==483
 # Complete transverse c^2 coefficients separate the centered square from every base output.
 aa=variable(0,5)
 for degree,co in [(2,1),(3,6),(4,15),(5,20),(6,15)]:assert cf(square,(degree-2,0,2))==scale(pw(aa,8-degree,5),co)
 # All alpha=0 boundaries,with beta,gamma,L,M still formal.
 q0=scalar(square,3,0);nn=7;xx,yy,zz,bb,cc,ll,mm=[variable(j,nn) for j in range(nn)]
 expected={2:pw(xx,2,nn),3:scale(mul(bb,pw(xx,3,nn)),-6),4:add(mul(add(scale(pw(bb,2,nn),9),scale(cc,-2)),pw(xx,4,nn)),scale(mul(pw(xx,3,nn),yy),-6)),5:total(scale(mul(mul(bb,cc),pw(xx,5,nn)),6),scale(mul(bb,mul(pw(xx,4,nn),yy)),18),scale(mul(pw(xx,4,nn),zz),-2)),6:total(mul(pw(cc,2,nn),pw(xx,6,nn)),scale(mul(cc,mul(pw(xx,5,nn),yy)),6),scale(mul(bb,mul(pw(xx,5,nn),zz)),6),scale(mul(pw(xx,4,nn),pw(yy,2,nn)),9))}
 for degree,e in expected.items():assert hom(q0,degree)==e
 assert hom(scalar(HH[0],3,0),2)==scale(mul(bb,pw(xx,2,nn)),-3)
 axis2=[[scalar(scalar(v,0,0),0,0) for v in row] for row in base[2]]
 assert not any(axis2[0]) and all(not row[indices[2][(2,0,0)]] for row in axis2)
 assert add(mul(axis2[1][indices[2][(1,0,1)]],axis2[2][indices[2][(0,2,0)]]),mul(axis2[1][indices[2][(0,2,0)]],axis2[2][indices[2][(1,0,1)]]),Q(-1))==con(24,3)
 def minor_at_alpha0(degree,mons):
  rows=[[scalar(base[degree][i][indices[degree][tuple(e)]],0,0) for e in mons] for i in range(3)]
  return determinant3(rows,4)
 assert minor_at_alpha0(3,[(2,1,0),(1,2,0),(1,1,1)])==con(-432,4)
 assert minor_at_alpha0(4,[(3,0,1),(2,1,1),(1,3,0)])==con(-336,4)
 assert all(not cf(q0,e) for e in [(2,1,0),(1,2,0),(1,1,1),(3,0,1),(2,1,1),(1,3,0)])
 q7=total(scale(mul(pw(x,5,n),mul(y,z)),6),scale(mul(gamma,mul(pw(x,6,n),z)),2),scale(mul(alpha,mul(pw(x,5,n),pw(z,2,n))),6));q8=mul(pw(x,6,n),pw(z,2,n));h7=scale(mul(pw(x,3,n),mul(pw(y,3,n),z)),8)
 assert hom(square,7)==q7 and hom(square,8)==q8 and hom(HH[2],7)==h7 and not hom(HH[0],7) and not hom(HH[1],7)
 assert [hom(f,7) for f in GG]==[{},mul(L,q7),add(h7,mul(M,q7))] and [hom(f,8) for f in GG]==[{},mul(L,q8),mul(M,q8)]
 # Conditional equality obstructions use a,c as coordinate axes and a free third source row.
 rd=json.loads((O/'rigidity.json').read_text());nr=8;rx,ry,rz,ra,rg,rm,eta,zeta=[variable(j,nr) for j in range(nr)];BB=total(ry,mul(eta,rx),mul(zeta,rz));top=decode(rd['top_cubic'],nr);sub=[rx,BB,rz,ra,{},rg,{},rm]
 assert compose(hom(GG[2],7),sub,nr)==mul(mul(pw(rx,3,nr),rz),top)
 assert cf(top,(1,2,0))==scale(variable(3,5),24) and cf(top,(0,2,1))==scale(variable(4,5),24)
 reduced=scalar(scalar(top,7,0),6,0);assert cf(reduced,(2,1,0))==scale(variable(2,3),6)
 # Axis-centre degree3 forms are reconstructed from the raw map,then restricted linearly.
 nr=6;rx,ry,rz,rg,eta,zeta=[variable(j,nr) for j in range(nr)];BB=total(ry,mul(eta,rx),mul(zeta,rz));translated=[compose(f,[rx,ry,add(rz,rg)],nr) for f in F];actual=[compose(hom(f,3),[rx,BB,rz,rg,eta,zeta],nr) for f in translated];saved=[decode(v,nr) for v in rd['cubic_outputs']];assert actual==saved
 assert cf(saved[0],(2,1,0))==con(-3,3) and cf(saved[1],(1,2,0))==con(24,3) and cf(saved[2],(1,1,1))==con(6,3)
 cc,ee,zzp=[variable(j,3) for j in range(3)];assert cf(saved[0],(2,0,1))==scale(zzp,-3) and cf(saved[0],(3,0,0))==add(scale(cc,-1),scale(ee,-3))
 assert compose(cf(saved[1],(3,0,0)),[scale(ee,-3),ee,{}],3)==scale(pw(ee,2,3),-12)
 # Nineteen source columns at the original centre for every nonzero (L,M).
 n=3;square0=pw(F[0],2,n);old={e for f in F for e in f if sum(e)>=2};assert len(old)==13 and len(square0)==6 and not old&set(square0)
 control=[F[0],add(F[1],square0),F[2]];assert [sum(sum(e)>=2 for e in f) for f in control]==[2,11,6] and determinant3([[derivative(f,j) for j in range(3)] for f in control],3)==con(1,3)
 pts=[[0,1,-8],[-1,1,4]];assert [[evaluate(f,t) for f in control] for t in pts]==[[Q(0),Q(1),Q(0)]]*2
 ledger=json.loads((ROOT/'MATH_STATUS.json').read_text());E={e['id']:e for e in ledger['entries']};pins={}
 for k in ['KCUB41','KCUB42']:
  e=E[k];assert e['state']=='proved';pins[k]={'scope_fingerprint':hashlib.sha256((e['state']+'\0'+e['title']+'\0'+e['scope']).encode()).hexdigest()[:16],'proof_sha256':sha(ROOT/e['canonical_source'])}
 profiles={'zero_shear_axis':[2,3,3,2,2,1,0],'zero_shear_off_axis':[3,3,3,2,2,1,0],'Lzero_Mnonzero_axis':[2,3,3,2,2,1,1],'Lzero_Mnonzero_off_axis':[3,3,3,2,2,1,1],'Lnonzero_axis':[2,3,3,2,2,2,1],'Lnonzero_off_axis':[3,3,3,2,2,2,1]}
 receipt={'status':'pass','claim':'KCUB43','checker_sha256':sha(Path(__file__)),'producer_sha256':C['source_sha256'],'arithmetic_sha256':sha(B/'verify.py'),'input_sha256':sha(ip),'checkpoint_sha256':sha(cp),'rigidity_sha256':sha(O/'rigidity.json'),'complete_shear_jet_matrix_entries':483,'all_base_matrices_and_square_rows_reconstructed':True,'exact_centre_dependent_output_equivalence_checked':True,'five_formal_c2_separators_checked':True,'all_alpha_zero_square_formulas_through6_checked':True,'zero_branch_base_pivots_checked':True,'full_generic_degree7_8_forms_checked':True,'exact_rank_profiles':profiles,'summed_rank_bounds':{k:sum(v) for k,v in profiles.items()},'conditional_top_and_degree3_equality_coefficients_checked':True,'universal_rank_and_UFD_equality_arguments_written_not_formalized':True,'nonmonomial_perfect_square_has_at_least3_terms_written_char0_proof':True,'minimum_nonlinear_shear_occurrences_all_affine_frames_lower':16,'minimum_nonlinear_shear_occurrences_all_affine_frames_upper':19,'sharp_nonlinear_family_minimum_not_established':True,'all_affine_frames_of_complete_degree8_rational_target_family_minimum':13,'equality13_requires_zero_shear_and_KCUB41_monomial_frames':True,'no_raw13_improvement_in_entire_family':True,'matching_existing_degree8_upper':{'degree':8,'counts':[2,11,6],'determinant_one':True,'points':pts,'image':[0,1,0],'new_global_improvement':False},'imports':pins,'solver_queries':0,'candidate_grids':0,'goal_complete':False,'global_positive_bounds_and_unrestricted_profiles_unchanged':True,'trust_boundary':'Complete483 shear coefficients,base matrices,square rows,centre-output equivalence and conditional equality coefficients are checked independently. Rank stability uses the imported KCUB41 base-rank theorem plus written outside-span arguments. Universal GL rank invariance,Laurent-binomial squarefreeness,UFD and equality exhaustion are written,not formal or externally reviewed. KCUB42 imported only for complete degree8 rational-target coverage. Sixteen is a lower bound;19 is a checked upper,not a claimed sharp nonlinear-family minimum.','seconds':time.monotonic()-start}
 proof=B/'QUADRATIC_TARGET_AFFINE_DEGREE8.md'
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (O/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS complete quadratic-target affine family:483 entries,all rank branches and equality obstructions;nonlinear support[16,19],full family minimum13 onlyzero shear;no solver/grid;global goal open',flush=True)
if __name__=='__main__':main()
