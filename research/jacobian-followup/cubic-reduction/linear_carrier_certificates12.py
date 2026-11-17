#!/usr/bin/env python3
"""Produce constant linear-combination certificates for all linear carrier charts."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib
from verify import decode,variable,mul,add,scale
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
DEGREES={3:6,4:6,5:5,6:4,7:4,8:3,9:3,10:3,11:3,12:3}
def span(polys):
 basis={}
 for i,poly in enumerate(polys):
  v=dict(poly);proof={i:Q(1)}
  while v:
   lead=min(v);c=v[lead]
   if lead not in basis:
    basis[lead]=(scale(v,1/c),{j:a/c for j,a in proof.items()});break
   b,bp=basis[lead];v=add(v,b,-c)
   for j,a in bp.items():proof[j]=proof.get(j,Q(0))-c*a
   proof={j:a for j,a in proof.items() if a}
 def solve(target):
  v=dict(target);proof={}
  while v:
   lead=min(v);c=v[lead]
   if lead not in basis:return None
   b,bp=basis[lead];v=add(v,b,-c)
   for j,a in bp.items():proof[j]=proof.get(j,Q(0))+c*a
  return [{'high_index':j,'coefficient':str(a)} for j,a in sorted(proof.items()) if a]
 return solve

def targets(n,np):
 T=[variable(j,np) for j in range(np)];pair=lambda i,j:mul(T[i],T[j]);out=[]
 if n==4:return [('t2',T[2])]
 if n==5:return [('t0*t1',pair(0,1)),('t1*t8',pair(1,8)),('t1*t6',pair(1,6))]
 if n==7:
  out=[('t4',T[4]),('t5',T[5])]
  for a,b in [(0,1),(6,7),(10,11)]:out.extend([(f't2*t{a}',pair(2,a)),(f't2*t{b}+t3*t{a}',add(pair(2,b),pair(3,a))),(f't3*t{b}',pair(3,b))])
  for j in [8,9,12,13]:out.extend([(f't2*t{j}',pair(2,j)),(f't3*t{j}',pair(3,j))])
 return out

def main():
 rows=[];hashes={};families=[]
 for n,D in DEGREES.items():
  folder=HERE/('linear-subspaces-n%02d-d%d'%(n,D));checkpoint=json.loads((folder/'checkpoint.json').read_text());assert checkpoint['status']=='complete' and checkpoint['source_sha256']==sha(HERE/'linear_carrier_subspaces12.py')
  for row in checkpoint['cases']:
   path=folder/('case-%03d.json'%row['case']);model=json.loads(path.read_text());np=model['parameter_count'];polys=[{tuple(a['exponent']):Q(a['coefficient']) for a in t['parameter_polynomial']} for t in model['high_coefficients']];one={(0,)*np:Q(1)}
   rec={'dimension':n,'degree_cap':D,'case':row['case'],'model':str(path.relative_to(HERE)),'model_sha256':sha(path),'kept':row['kept']}
   solve=span(polys);unit=solve(one)
   if unit is not None:rec.update(status='unit-obstruction',unit=unit)
   else:
    rec['status']='surviving-family';rec['generators']=[]
    for label,target in targets(n,np):
     proof=solve(target);assert proof is not None,(n,label)
     rec['generators'].append({'label':label,'polynomial':[{'exponent':list(e),'coefficient':str(c)} for e,c in sorted(target.items())],'linear_combination':proof})
    families.append((n,row['case']))
   rows.append(rec)
 assert len(rows)==512 and families==[(4,1),(5,10),(7,15),(12,0)]
 result={'status':'complete-primary-certificates','source_sha256':sha(Path(__file__)),'model_producer_sha256':sha(HERE/'linear_carrier_subspaces12.py'),'charts':512,'unit_obstructions':sum(r['status']=='unit-obstruction' for r in rows),'families':families,'cases':rows,'scope':'Constant linear combinations of exact high-degree parameter coefficient polynomials; no Groebner output is used as proof.'}
 folder=HERE/'linear-carrier-subspaces12';folder.mkdir(exist_ok=True);(folder/'certificates.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS PRIMARY',result['charts'],result['unit_obstructions'],'units; families',families,flush=True)
if __name__=='__main__':main()
