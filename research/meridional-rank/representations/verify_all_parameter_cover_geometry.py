#!/usr/bin/env python3
"""Rigorous fixed-triangulation complete geometry and uniform slope lengths.
Run under Sage10.9 with SnapPy3.3.2. No family-topology identification is asserted.
The complete cusp is constructed at area11 and checked directly to be embedded.
"""
from pathlib import Path
from copy import deepcopy
import argparse,hashlib,json,platform,sys,time,csv,base64
HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
for extra in [REPO/'.scratch/mrc-sage-snappy',Path('/tmp/mrc-knot-audit-venv/lib/python3.12/site-packages')]:
    if extra.is_dir():sys.path.append(str(extra))
import snappy
from sage.all import RealIntervalField,QQ
from sage.env import SAGE_VERSION
from snappy.geometric_structure.cusp_neighborhood.complex_cusp_cross_section import ComplexCuspCrossSection
from snappy.verify.hyperbolicity import check_logarithmic_gluing_equations_and_positively_oriented_tets

INPUT=HERE/'all_parameter_cover_geometry_inputs.json'
RECEIPT=HERE/'all_parameter_cover_geometry_verification.json'
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
def require(value,message):
    if not value:raise AssertionError(message)
def native(obj):
    return json.loads(json.dumps(obj))
def interval(x):
    return {'lower':str(x.lower().exact_rational()),'upper':str(x.upper().exact_rational()),'display':str(x)}
def cinterval(z):return {'real':interval(z.real()),'imag':interval(z.imag())}
def snapshot(M):
    return {'gluing':native(M._get_tetrahedra_gluing_data()),
            'cusp_indices_and_peripheral_rows':native(M._get_cusp_indices_and_peripheral_curve_data()),
            'logarithmic_gluing_equations':[[int(x)for x in row]for row in M.gluing_equations()]}
def validate_face_pairings(data):
    for t,(neighbors,gluing) in enumerate(data['gluing']):
        require(len(neighbors)==len(gluing)==4,'four faces per tetrahedron')
        for f,(u,p)in enumerate(zip(neighbors,gluing)):
            require(sorted(p)==list(range(4)),'face permutation')
            require(sum(p[i]>p[j]for i in range(4)for j in range(i+1,4))%2==1,'oriented tetrahedron gluings')
            j=p[f];back=data['gluing'][u]
            require(back[0][j]==t,'reciprocal adjacent tetrahedron')
            q=back[1][j]
            require([q[p[i]]for i in range(4)]==list(range(4)),'reciprocal face map')
def validate_edge_orientations(data):
    g=data['gluing'];corners=[(t,v,w)for t in range(len(g))for v in range(4)for w in range(4)if v!=w]
    parent={x:x for x in corners}
    def find(a):
        while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
        return a
    def union(a,b):parent[find(a)]=find(b)
    for t,(neighbors,ps)in enumerate(g):
        for f,p in enumerate(ps):
            for v in range(4):
                for w in range(4):
                    if v!=w and f not in (v,w):union((t,v,w),(neighbors[f],p[v],p[w]))
    require(all(find((t,v,w))!=find((t,w,v))for t,v,w in corners),'no ideal edge identified to itself with reversed orientation')
    count=len({find(a)for a in corners})//2
    require(count==len(g),'one edge class per tetrahedron')
    return {'unoriented_ideal_edges':count,'no_reversed_self_identifications':True}

def validate_torus_cusp_link(data):
    gluing=data['gluing'];n=len(gluing);faces=[(t,v)for t in range(n)for v in range(4)]
    corners=[(t,v,w)for t,v in faces for w in range(4)if w!=v]
    parent={x:x for x in corners}
    def find(a):
        while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
        return a
    def union(a,b):parent[find(a)]=find(b)
    adjacency={x:[]for x in faces};side_pairs=set()
    for t,v in faces:
        neighbors,ps=gluing[t]
        for f in range(4):
            if f==v:continue
            u=neighbors[f];p=ps[f];adjacency[t,v].append((u,p[v]))
            a=(t,v,f);b=(u,p[v],p[f]);require(a!=b,'no unpaired cusp side')
            side_pairs.add(tuple(sorted([a,b])))
            for w in range(4):
                if w not in (v,f):union((t,v,w),(u,p[v],p[w]))
    seen={faces[0]};todo=[faces[0]]
    for a in todo:
        for b in adjacency[a]:
            if b not in seen:seen.add(b);todo.append(b)
    V=len({find(a)for a in corners});E=len(side_pairs);F=len(faces)
    require(len(seen)==F,'one connected cusp link')
    require(2*E==3*F,'every cusp side paired')
    require(V-E+F==0,'cusp Euler characteristic zero')
    # Face pairings are odd in coherent tetrahedron orientations. Thus the
    # connected closed cusp surface is orientable, and chi0 identifies T².
    return {'vertices':V,'edges':E,'faces':F,'Euler_characteristic':0,
            'connected':True,'orientable_by_odd_tetrahedron_face_maps':True,'surface':'torus'}

def validate_peripheral_cycles(data):
    gluing=data['gluing'];indices,rows=data['cusp_indices_and_peripheral_rows'];n=len(gluing)
    require(indices==[[0]*4 for _ in range(n)],'one cusp index')
    require(len(rows)==4*n and all(len(r)==16 for r in rows),'peripheral row dimensions')
    require(all(not any(rows[4*t+2*c+1])for t in range(n)for c in range(2)),'unused handedness rows zero')
    for t,(neighbors,ps)in enumerate(gluing):
        for curve in range(2):
            row=rows[4*t+2*curve]
            for v in range(4):
                require(row[4*v+v]==0 and sum(row[4*v:4*v+4])==0,'integral peripheral triangle flux')
                for f in range(4):
                    if f!=v:
                        u=neighbors[f];p=ps[f];other=rows[4*u+2*curve]
                        require(row[4*v+f]==-other[4*p[v]+p[f]],'peripheral flux across paired side')
    return {'two_exact_integral_closed_cycles':True,'triangle_flux_and_glued_side_conditions':True}

def verify_covolume_index_one(m,l,area):
    require(m.imag()>0 and l>0,'oriented peripheral period basis')
    ratio=l*m.imag()/area
    require(ratio>QQ(99)/100 and ratio<QQ(101)/100,'period covolume index is one')
    # Closed integral torus cycles have covolume equal to an integer index
    # times the cusp area. This positive integer lies in(.99,1.01), hence1.
    return ratio

def lift_intervals(M,z,equations=None):
    CIF=z[0].parent();RF=CIF.real_field();two_pi_i=CIF(0,2*RF.pi())
    logs=[w for x in z for w in [x.log(),(1/(1-x)).log(),((x-1)/x).log()]]
    rows=equations if equations is not None else M.gluing_equations()
    out=[]
    for i,row in enumerate(rows):
        lhs=sum((int(a)*x for a,x in zip(row,logs)),CIF(0))
        target=two_pi_i if i<M.num_tetrahedra()else CIF(0)
        residual=abs(lhs-target)
        require(residual<RF(QQ(1)/1000),'full logarithmic gluing/completeness lift')
        out.append({'row':i,'kind':'edge'if i<M.num_tetrahedra()else'complete_peripheral',
                    'lhs':cinterval(lhs),'residual_absolute':interval(residual)})
    require(len(out)==M.num_tetrahedra()+2,'all edge and both complete cusp equations')
    return out

def cross_section(M,z,area):
    c=ComplexCuspCrossSection.fromManifoldAndShapes(M,z)
    RF=z[0].real().parent();c.normalize_cusps([RF(QQ(area))])
    return c

def embedded(c,strict_margins=True):
    scales=[c.compute_scale_for_std_form(v)for v in c.mcomplex.Vertices]
    gaps=[c._exp_distance_edge(e)for e in c.mcomplex.Edges]
    # These are tested BEFORE any ensure_std_form/ensure_disjoint auto-shrink.
    # Thus the exact normalized area is retained.
    require(all(x>1 for x in scales),'cusp must miss every opposite face')
    require(all(x>1 for x in gaps),'cusp horoballs separated on every ideal edge')
    if strict_margins:
        require(all(x>QQ(6)/5 for x in scales),'recorded standard-form margin')
        require(all(x>QQ(21)/20 for x in gaps),'recorded edge-separation margin')
    return scales,gaps

def uniform_bound(m,l,p,q0,step,base_bound=10,step_bound=19):
    require(abs(p)==1,'every family slope primitive because meridian coefficient is +/-1')
    require(step>0,'family parameter step positive')
    base=p*m+q0*l;increment=step*l
    require(base.real()>base_bound,'claimed affine real base lower bound')
    require(increment.real()>step_bound,'claimed affine real step lower bound')
    return base,increment

def reject(name,call,allowed):
    try:call()
    except Exception as error:
        require(type(error).__name__ in allowed,'unexpected rejection failure '+name+': '+type(error).__name__)
        return {'name':name,'exception':type(error).__name__,'rejected':True}
    raise AssertionError('corrupted control accepted: '+name)

def verify_case(entry,parameters):
    path=(INPUT.parent/entry['path']).resolve();require(sha(path)==entry['sha256'],'triangulation byte identity')
    M=snappy.Manifold(str(path));require(M.num_tetrahedra()==entry['tetrahedra'],'tetrahedron count')
    require(M.num_cusps()==1 and M.is_orientable(),'orientable one cusp')
    require(M.cusp_info(0)['complete?'] and tuple(M.cusp_info(0)['filling'])==(0,0),'unfilled complete cusp')
    actual=snapshot(M);require(actual==entry['snapshot'],'exact face-pairing/peripheral/gluing binding')
    validate_face_pairings(actual)
    edge_orientations=validate_edge_orientations(actual);link=validate_torus_cusp_link(actual);cycles=validate_peripheral_cycles(actual)
    success,z=M.verify_hyperbolicity(bits_prec=parameters['bits_prec'])
    require(success is True,'certified complete hyperbolicity')
    require(all(x.imag()>0 for x in z),'every tetrahedron positive')
    lifts=lift_intervals(M,z)
    c=cross_section(M,z,parameters['cusp_area']);scales,gaps=embedded(c)
    m,l=c.all_normalized_translations()[0]
    require(l>0,'positive real longitude normalization')
    index=verify_covolume_index_one(m,l,c.cusp_areas()[0])
    p,q0,step=parameters['slope_p'],parameters['slope_q0'],parameters['slope_q_step']
    base,increment=uniform_bound(m,l,p,q0,step)
    # Known short longitude control prevents interpreting every nonzero slope
    # as long. Its certified upper length is below6 at this same cusp.
    require(abs(l)<6,'short longitude positive control')
    bad_shapes=list(z);bad_shapes[0]=bad_shapes[0].conjugate()
    bad_rows=[list(row)for row in actual['logarithmic_gluing_equations']];bad_rows[0][0]+=1
    bad_cycles=deepcopy(actual);bad_cycles['cusp_indices_and_peripheral_rows'][1][0][1]+=1
    controls=[
        reject('corrupted integral peripheral flux',lambda:validate_peripheral_cycles(bad_cycles),{'AssertionError'}),
        reject('doubled meridian is not a primitive basis',lambda:verify_covolume_index_one(2*m,l,c.cusp_areas()[0]),{'AssertionError'}),
        reject('negative tetrahedron orientation',lambda:check_logarithmic_gluing_equations_and_positively_oriented_tets(M,bad_shapes),{'ShapePositiveImaginaryPartNumericalVerifyError'}),
        reject('single changed logarithmic gluing coefficient',lambda:lift_intervals(M,z,bad_rows),{'AssertionError'}),
        reject('area100 cusp is not certified embedded',lambda:embedded(cross_section(M,z,100),False),{'AssertionError'}),
        reject('reversed meridian coefficient gives false long slope',lambda:uniform_bound(m,l,1,q0,step),{'AssertionError'}),
        reject('reversed parameter direction violates affine lower bound',lambda:uniform_bound(m,l,p,q0,-step),{'AssertionError'}),
        reject('excessive base bound11 rejected',lambda:uniform_bound(m,l,p,q0,step,11,19),{'AssertionError'})]
    return {'input':{k:entry[k]for k in ['path','sha256','tetrahedra']},'complete_hyperbolicity':True,
            'shape_intervals':[cinterval(x)for x in z], 'all_positive_shapes':True,'edge_topology':edge_orientations,'exact_torus_cusp_link':link,'peripheral_cycles':cycles,'peripheral_covolume_index_interval':interval(index),'peripheral_basis_certified':True,
            'full_logarithmic_lifts':lifts,'exact_cusp_area':parameters['cusp_area'],
            'area_enclosure':interval(c.cusp_areas()[0]),
            'standard_form_scale_intervals':[interval(x)for x in scales],
            'edge_separation_exp_intervals':[interval(x)for x in gaps],
            'standard_form_lower_bound':'6/5','edge_separation_lower_bound':'21/20',
            'auto_shrinking_used':False,'embedded_cusp_certified':True,
            'translations':{'meridian':cinterval(m),'longitude_real_positive':interval(l)},
            'slope_base_translation':cinterval(base),'slope_step_translation':cinterval(increment),
            'slope_n0_length':interval(abs(base)),
            'uniform_length_inequality':'length(-1,2+6n) >= Re(-m+(2+6n)l) > 10+19n > 6 for every integer n>=0',
            'uniform_justification':'The same m,l occur for every n; Re(-m+2l)>10 and6l>19. No sampled-n extrapolation.',
            'negative_controls':controls,'short_longitude_length':interval(abs(l))}

def run_verification():
    """Return the full rigorous geometry result; called only from Sage mode.

    No file mutation or printing occurs here. The master checker can call this
    function after its independent standard-library binding checks.
    """
    started=time.monotonic();data=json.loads(INPUT.read_text())
    require(snappy.__version__==data['software']['SnapPy'],'pinned SnapPy version')
    require(SAGE_VERSION==data['software']['Sage'],'pinned Sage version')
    pkg=Path(snappy.__file__).resolve().parent
    for name,digest in data['trusted_library_source_sha256'].items():require(sha(pkg/name)==digest,'trusted source pin '+name)
    for item in data['source_snapshots']:
        require(sha(INPUT.parent/item['path'])==item['sha256'],'retained official source snapshot '+item['path'])
    distribution_record=pkg.parent/'snappy-3.3.2.dist-info/RECORD'
    require(sha(distribution_record)==data['software']['distribution_RECORD_sha256'],'pinned distribution RECORD')
    verified_distribution_entries=0
    for name,digest,size in csv.reader(distribution_record.read_text().splitlines()):
        if name.startswith('snappy/') and digest.startswith('sha256='):
            path=pkg.parent/name;actual=base64.urlsafe_b64encode(hashlib.sha256(path.read_bytes()).digest()).rstrip(b'=').decode()
            require(actual==digest.split('=',1)[1],'unmodified distribution file '+name)
            verified_distribution_entries+=1
    require(len(data['cases'])==1 and data['cases'][0]['id']=='direct_bound14','single canonical directly bound parent')
    results={entry['id']:verify_case(entry,data['parameters'])for entry in data['cases']}
    return {'status':'PASS rigorous fixed-parent geometry and uniform slope-length certificates',
            'cases':results,'software':{'SnapPy':snappy.__version__,'Sage':SAGE_VERSION,'Python':platform.python_version(),
                                       'snappy_package':str(pkg),'kernel_sha256':sha(snappy.SnapPy.__file__),
                                       'distribution_RECORD_sha256':sha(distribution_record),'verified_distribution_entries':verified_distribution_entries},
            'trusted_library_source_sha256':data['trusted_library_source_sha256'],
            'source_snapshots':data['source_snapshots'],
            'input_sha256':sha(INPUT),'checker_sha256':sha(__file__),'elapsed_seconds':round(time.monotonic()-started,3),
            'scope':{'fixed_triangulation_only':True,'actual_Kn_cover_identification':'separate master-checker and written topology theorem',
                     'all_parameter_length_proof':True,'hyperbolic_filling_theorem_applied':False,
                     'full_formal_verification':False,'external_specialist_review':False},
            'trusted_algorithms':'SnapPy3.3.2 Krawczyk shape certification/full logarithmic lift check, Sage interval arithmetic, exact kernel face-pairing/peripheral transfer, complex cusp development. The redundant residual check is not a proof of exact gluing by itself. Standard-form and edge separation checked directly before any shrinking; no maximal-cusp computation required.'}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--receipt',type=Path,default=RECEIPT);args=parser.parse_args()
    result=run_verification()
    if args.receipt.exists():
        previous=json.loads(args.receipt.read_text())
        for key in ['internal_pdf_reviews','review_history','edition_replays','root_review']:
            if key in previous:result[key]=previous[key]
    temp=args.receipt.with_suffix('.json.tmp');temp.write_text(json.dumps(result,indent=2)+'\n');temp.replace(args.receipt)
    print(json.dumps({'status':result['status'],'cases':list(result['cases']),'elapsed_seconds':result['elapsed_seconds'],'receipt':str(args.receipt)},indent=2))
if __name__=='__main__':main()
