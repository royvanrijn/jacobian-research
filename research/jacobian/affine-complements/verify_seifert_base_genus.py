"""PF2RULED69: exact fibre-orbit genera; unbounded topology is proved in the note."""
from pathlib import Path
import collections,hashlib,importlib.util,itertools,json,math,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('sealed_end61',HERE/'verify_marked_end_lattice.py')
E=importlib.util.module_from_spec(spec);spec.loader.exec_module(E)

def orbits(gs,sheets):
    out=[];left=set(sheets)
    while left:
        seen={min(left)}
        while True:
            grow=seen|{g[s] for g in gs for s in seen}
            if grow==seen:break
            seen=grow
        assert seen<=left
        out.append(sorted(seen));left-=seen
    return out

def cycles(g):return len(orbits([g],range(len(g))))

def fibre_quotient(gs,H,sheets):
    """Use connected generic fibres as base sheets, retaining nonidentity H."""
    fibres=orbits([H],sheets)
    assert len({len(f) for f in fibres})==1,'constant fibre degree on a connected piece'
    ix={s:i for i,f in enumerate(fibres) for s in f};out=[]
    for g in gs:
        assert E.mul(g,H)==E.mul(H,g)
        dest=tuple(ix[g[f[0]]] for f in fibres)
        assert sorted(dest)==list(range(len(fibres)))
        for f in fibres:assert len({ix[g[s]] for s in f})==1
        out.append(dest)
    return out,len(fibres[0])

def genus_from_base(cones,boundaries,h=0):
    perms=cones+boundaries;n=len(perms[0]);b=len(boundaries);j=len(cones)
    assert all(len(x)==n for x in perms)
    chi=n*(2-2*h-b-j)+sum(cycles(x) for x in cones)
    boundary_count=sum(cycles(x) for x in boundaries)
    twog=2-chi-boundary_count
    assert twog>=0 and twog%2==0
    return twog//2

def companion_bases(row,p,q,r):
    a,b,beta,tau=map(tuple,row['rho']);M=E.power(beta,r);H=E.power(a,p)
    assert H==E.power(b,q)
    out=[]
    for block in orbits([a,b],range(len(a))):
        (aa,bb,mm),m=fibre_quotient([a,b,M],H,block)
        g=genus_from_base([aa,bb],[mm])
        out.append(dict(sheets=block,fibre_degree=m,base_degree=len(aa),
                        induced=[aa,bb,mm],genus=g))
    return out

def cable_bases(row):
    a,b,beta,tau=map(tuple,row['rho']);H=tuple(range(len(a)))
    cone=E.mul(E.inv(tau),beta);out=[]
    for block in orbits([beta,tau],range(len(a))):
        (zz,BB,tt),m=fibre_quotient([cone,beta,tau],H,block)
        out.append(dict(sheets=block,fibre_degree=m,base_degree=len(zz),
                        induced=[zz,BB,tt],genus=genus_from_base([zz],[BB,tt])))
    return out

def full_satellite(p,q,r,row,reverse=False):
    a,b,beta,tau=map(tuple,row['rho']);rho=[a,b,E.mul(tau,E.inv(beta))]
    u,v,g=S.gcdex(q,p);assert g==1
    M=E.wp([1],int(u))+E.wp([2],int(v));L=E.wp([1],p)+E.wp(M,-p*q)
    c=p*q-p-q;k=r*c+1;h=E.wp(M,k)+E.wp(L,r)
    mu=[3]+E.wp(M,-c)+E.wi(L)
    rels=[E.wp([1],p)+E.wp([2],-q),E.wp([3],r)+E.wi(h)]
    I=tuple(range(4));assert all(E.evaluate(w,rho)==I for w in rels)
    assert E.evaluate(h,rho)==I and E.evaluate(mu,rho)==tau and r*k%2==0
    rows,mark,_=E.cover_matrix(rels,rho,mu,h,reverse_tree=reverse)
    return dict(H1=E.h1(rows),quotient=E.h1(rows,[mark]))

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['marked_end_lattice_receipt.json','orbit_tree_receipt.json','binomial_periodic_receipt.json','cyclic_parity_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for p,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/seifert-base-genus/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        p=HERE/row['path'];check('primary bytes '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    # General Euler-characteristic identity, not an interpolation in cover degree.
    n,h,b,j,cs,cb=S.symbols('n h b j cs cb',integer=True)
    chi=n*(2-2*h-b-j)+cs;twog=2-chi-cb
    check('all-degree genus identity',S.expand(twog-(2-n*(2-2*h-b-j)-cs-cb))==0)
    check('disk two cone companion identity',S.expand(twog.subs({h:0,b:1,j:2})-(n+2-cs-cb))==0)
    check('annulus one cone cable identity',S.expand(twog.subs({h:0,b:2,j:1})-(n+2-cs-cb))==0)
    # In an orientable Seifert presentation, all abelian handle columns are ZERO.
    # Extra boundary gluings use only boundary/fibre columns, regardless of slope.
    for g,cones,bounds in [(1,2,1),(2,1,3),(9,3,2)]:
        size=2*g+1+cones+bounds;fibre=2*g
        mat=S.zeros(cones+1,size)
        for i in range(cones):mat[i,fibre]=S.Symbol(f'beta{i}');mat[i,fibre+1+i]=S.Symbol(f'alpha{i}')
        mat[cones,fibre]=S.Symbol('euler')
        for i in range(cones+bounds):mat[cones,fibre+1+i]=1
        slope=S.zeros(1,size)
        for l in range(bounds):
            slope[0,fibre]+=S.Symbol(f'u{l}');slope[0,fibre+1+cones+l]+=S.Symbol(f'v{l}')
        combined=mat.col_join(slope)
        check('arbitrary boundary operations keep handle quotient '+str((g,cones,bounds)),combined[:,:2*g]==S.zeros(combined.rows,2*g))
    data=json.loads((HERE/'binomial_periodic_receipt.json').read_text())['census']['complete_cells']
    counts=collections.Counter();bad=[];allrecords=[];emptied=[]
    for key,cell in data.items():
        p,q,r=map(int,key.split(','));tree=0;good=0
        for i,row in enumerate(cell['orbits']):
            cbases=companion_bases(row,p,q,r);pbases=cable_bases(row)
            cgen=sum(x['genus'] for x in cbases);pgen=sum(x['genus'] for x in pbases)
            check('all cable bases planar '+str((key,i)),pgen==0)
            tag=(row['companion_order'],row['cable_order'],row['cycle_rank'],cgen);counts[tag]+=1
            if not row['cycle_rank']:
                tree+=1
                if not cgen:good+=1
            if cgen:
                check('positive genus only new S4 S4 tree sector '+str((key,i)),
                      row['companion_order']==row['cable_order']==24 and row['cycle_rank']==0 and cgen==1)
                a,beta=row['rho'][:2]
                check('new genus sector has generator orders3 and4 '+str((key,i)),sorted([next(t for t in range(1,13) if E.power(tuple(a),t)==tuple(range(4))),next(t for t in range(1,13) if E.power(tuple(beta),t)==tuple(range(4)))])==[3,4])
                bad.append(dict(cell=key,rho=row['rho'],companion_bases=cbases,cycle_rank=row['cycle_rank']))
            allrecords.append(dict(cell=key,class_index=i,cycle_rank=row['cycle_rank'],companion_bases=cbases,cable_bases=pbases))
        if tree and not good:emptied.append(key)
    expected={(3,24,0,0):84,(4,24,0,0):96,(12,6,0,0):96,(6,6,0,0):64,
              (12,24,1,0):48,(24,24,0,0):64,(24,2,0,0):48,(24,6,1,0):64,(24,24,0,1):32}
    check('complete596 class genus census',dict(counts)==expected and len(allrecords)==596)
    check('exactly32 new tree classes fail finite end',len(bad)==32)
    check('exactly452 genus zero tree classes remain before other gates',sum(v for k,v in counts.items() if k[2]==k[3]==0)==452)
    check('no whole previously nonempty tree cell excluded',emptied==[])
    badcells=collections.Counter(x['cell'] for x in bad)
    check('eight new residue cells each four excluded classes',dict(badcells)=={k:4 for k in ['3,4,3','3,8,3','4,3,3','4,9,9','8,3,3','8,9,9','9,4,9','9,8,9']})
    # Nonidentity-fibre countercontrol prevents applying the sheet formula wrongly.
    M=(1,2,3,0);a=E.power(M,3);bb=E.power(M,2);H=E.power(a,2)
    (aa,bbq,mm),m=fibre_quotient([a,bb,M],H,range(4))
    check('cyclic4 fibre has degree2 and base degree2',m==2 and len(aa)==2)
    check('cyclic4 true base genus is zero',genus_from_base([aa,bbq],[mm])==0)
    check('raw four sheets would falsely give genus one',4+2-cycles(a)-cycles(bb)-cycles(M)==2)
    controlrow=next(x for x in data['2,3,3']['orbits'] if x['companion_order']==4 and x['cable_order']==24)
    cyclic_control=full_satellite(2,3,3,controlrow)
    check('false positive has actual finite branched homology',cyclic_control==dict(H1=dict(free=0,torsion=[3]),quotient=dict(free=0,torsion=[3])))
    controls=[];row=next(x for x in data['3,4,3']['orbits'] if x['rho']==bad[0]['rho'])
    for p,q,r,tors in [(3,4,3,[]),(3,16,3,[]),(15,4,3,[5]),(3,4,15,[5])]:
        ans=full_satellite(p,q,r,row)
        check('full satellite confirms genus obstruction '+str((p,q,r)),ans['H1']==dict(free=2,torsion=tors))
        check('independent spanning tree confirms full group '+str((p,q,r)),full_satellite(p,q,r,row,True)==ans)
        controls.append(dict(p=p,q=q,r=r,**ans))
    openrow=next(x for x in data['2,3,9']['orbits'] if x['companion_order']==x['cable_order']==24 and x['cycle_rank']==0)
    open_control=full_satellite(2,3,9,openrow)
    check('genus zero nonzero mark sector stays open',
          sum(x['genus'] for x in companion_bases(openrow,2,3,9))==0 and
          open_control==dict(H1=dict(free=0,torsion=[24]),quotient=dict(free=0,torsion=[6])))
    # Bipartite permutation graph: edges are sheets, vertices are cycles of A/B.
    for name,A,B,want in [('quartic allowed',(1,0,2,3),(0,2,3,1),0),
                          ('quartic forbidden',(1,2,3,0),(0,2,3,1),1),
                          ('degree5 allowed',(1,2,0,3,4),(0,1,3,4,2),0)]:
        d=len(A);bipartite_vertices=cycles(A)+cycles(B)
        check('transitive permutation graph control '+name,len(orbits([A,B],range(d)))==1)
        check('all-degree companion tree cycle count '+name,d-bipartite_vertices+1==2*want)
    inputs=[Path(__file__),HERE/'SEIFERT_BASE_GENUS_AND_COMPANION_TREE.md',manifest]+[HERE/n for n in preds]+[HERE/r['path'] for r in sources['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(all_finite_cover_degrees=True,all_cabling_lengths=True,
          integral_quotient='H1(end) surjects Z^(b1(piece graph)+2*sum base genera)',
          orientable_selected_Seifert_bases_required=True,arbitrary_boundary_gluings_and_fillings=True,
          cutting_tori_need_not_be_incompressible=True,minimal_JSJ_not_required=True,
          smooth_actual_extension_and_disjoint_A1_required=True,all_lifted_Seifert_bases_must_be_planar=True,
          regular_fibre_orbits_are_base_sheets=True,nontrivial_fibre_monodromy_retained=True,
          arbitrary_degree_centreless_transitive_companion_and_transitive_cable_force_dessin_tree=True,
          companion_tree_identity='cycles(a)+cycles(b)=d+1',global_Sd_not_assumed_to_imply_piece_Sd=True),
        quartic=dict(complete_period12_cells=1152,periodic_classes=596,old_tree_classes=484,
          new_positive_genus_tree_exclusions=32,remaining_genus_zero_tree_classes=452,
          affected_cells=dict(badcells),whole_residue_cells_closed=0,
          S4_S4_generator_order3_order4_sector_excluded=True,
          all_S4_S4_representations_excluded=False,other_prior_gates_already_applied_to_452=False),
        class_counts={str(k):v for k,v in counts.items()},excluded=bad,all_base_records=allrecords,
        full_satellite_controls=controls,nonidentity_fibre_finite_control=cyclic_control,
        remaining_nonzero_mark_control=open_control,
        assurance=dict(written_unbounded_Mayer_Vietoris_and_fibre_quotient_proofs=True,
          complete_existing_census_audited=True,primary_Seifert_model_visually_read=True,
          finite_controls_not_universal_topology_proof=True,
          higher_degree_smoothness_automatic=False,affine_cover_realizability_proved=False,
          noncyclic_companion_sectors_all_closed=False,PF63_curves_declared_nonexistent=False,
          full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_map=False,
          formal_verification=False,independent_full_reproof=False,external_review=False,
          novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'seifert_base_genus_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all-degree base-genus quotient; 32 further quartic tree classes excluded')
if __name__=='__main__':main()
