"""Actual closed-braid encoding and unmodified external algorithm comparison."""
from pathlib import Path
import sys,json,importlib.util,itertools,hashlib
P=Path(__file__).resolve().parent
pin=json.loads((P/'sources/recursive-linking-endpoints/author_code_pin.json').read_text())
assert hashlib.sha256((P/pin['path']).read_bytes()).hexdigest()==pin['sha256']
spec=importlib.util.spec_from_file_location('cahn_unmodified',P/pin['path']);C=importlib.util.module_from_spec(spec);spec.loader.exec_module(C)
def diagram(strands,word):
    positions=list(range(strands));events=[[] for i in range(strands)];signs={}
    for crossing,letter in enumerate(word):
        at=abs(letter)-1;assert 0<=at<strands-1
        over,under=(positions[at],positions[at+1]) if letter>0 else (positions[at+1],positions[at])
        events[over].append((crossing,'o'));events[under].append((crossing,'u'))
        positions[at],positions[at+1]=positions[at+1],positions[at];signs[crossing]=1 if letter>0 else -1
    end={s:i for i,s in enumerate(positions)};strand=0;seen=set();walk=[];toparc={};counter=0
    while strand not in seen:
        seen.add(strand);toparc[strand]=counter%len(word)
        walk+=events[strand];counter+=sum(kind=='u' for cross,kind in events[strand]);strand=end[strand]
    assert strand==0 and len(seen)==strands,'not a knot'
    incoming={};overarc={};counter=0
    for cross,kind in walk:
        if kind=='u':incoming[cross]=counter;counter+=1
        else:overarc[cross]=counter%len(word)
    over=[0]*len(word);eps=[0]*len(word)
    for cross,index in incoming.items():over[index]=overarc[cross];eps[index]=signs[cross]
    # One positive Reidemeister-I curl immediately after the last old underpass.
    # It splits old arc0, placing ALL its old overpasses on new arc1.
    if len(over)%2:
        over=[0]+[j+1 for j in over];eps=[1]+eps
        toparc={s:j+1 for s,j in toparc.items()}
    return over,eps,[toparc[s] for s in range(strands)]
def coloring(over,p=3,accept=lambda vec:True):
    n=len(over);mat=[[0]*n for i in range(n)]
    for i,j in enumerate(over):mat[i][j]+=2;mat[i][i]-=1;mat[i][(i+1)%n]-=1
    mat=[[x%p for x in row] for row in mat];pivots=[];at=0
    for j in range(n):
        pivot=next((i for i in range(at,n) if mat[i][j]),None)
        if pivot is None:continue
        mat[at],mat[pivot]=mat[pivot],mat[at];u=pow(mat[at][j],-1,p);mat[at]=[(u*x)%p for x in mat[at]]
        for i in range(n):
            if i!=at and mat[i][j]:
                q=mat[i][j];mat[i]=[(x-q*y)%p for x,y in zip(mat[i],mat[at])]
        pivots.append(j);at+=1
    freecols=sorted(set(range(n))-set(pivots))
    for vals in itertools.product(range(p),repeat=len(freecols)):
        vec=[0]*n
        for col,val in zip(freecols,vals):vec[col]=val
        for i,j in enumerate(pivots):vec[j]=(-sum(mat[i][f]*vec[f] for f in freecols))%p
        if len(set(vec))>1 and accept(vec):
            assert all((2*vec[j]-vec[i]-vec[(i+1)%n])%p==0 for i,j in enumerate(over))
            return vec
    raise AssertionError('no nontrivial coloring')
def cable_trefoil(r,k):
    # r parallel copies of each sigma1 in sigma1^3, then compensate
    # blackboard writhe3: residual ordinary pattern exponent k-3r.
    block=[j+i for i in range(r) for j in range(r,0,-1)]
    delta=list(range(1,r));power=k-3*r
    twist=delta*power if power>=0 else [-x for x in reversed(delta)]*(-power)
    return 2*r,block*3+twist
def run():
    out=[]
    examples=[('torus',(p,q),p,list(range(1,p))*q,p*q//3) for p,q in [(2,3),(2,9),(2,15),(4,3),(8,3)]]
    for r,k in [(3,7),(3,5),(3,2)]:
        strands,word=cable_trefoil(r,k);examples.append(('cable_trefoil',(r,k),strands,word,2*r*r))
    strands,word=cable_trefoil(2,3);examples.append(('birth_C3_trefoil',(2,3),strands,word,2))
    for kind,params,n,word,want in examples:
        over,eps,top=diagram(n,word)
        if kind=='cable_trefoil':
            r,k=params
            accept=lambda vec:len({vec[j] for j in top[:r]})==1 and len({vec[j] for j in top[r:]})==1 and vec[top[0]]!=vec[top[r]]
        else:accept=lambda vec:True
        colors=coloring(over,accept=accept)
        if kind=='birth_C3_trefoil':
            topcolors=[colors[j] for j in top]
            assert (topcolors[0]-topcolors[1])%3==(topcolors[2]-topcolors[3])%3!=0
        result=C.DLNmatrix(3,over,eps,colors)
        assert result[0][1]==result[1][0]==want,(kind,params,result,want)
        out.append(dict(kind=kind,parameters=list(params),strands=n,braid_word=word,overstrands=over,signs=eps,coloring=colors,
            top_arc_indices=top,top_colors=[colors[j] for j in top],companion_and_pattern_coloring_matched=True,
            author_matrix=[[str(x) for x in row] for row in result],mutual=int(want),preferred_diagonals_identified=False))
        print(kind,params,'crossings',len(over),'author mutual',result[0][1])
    # Same C_(3,2)(trefoil), a DIFFERENT coloring selected by the old
    # companion being CYCLIC: its two bundle meridians are equal reflections.
    n,word=cable_trefoil(3,2);over,eps,top=diagram(n,word)
    def cyclic_old(vec):return (vec[top[0]]-vec[top[1]]+vec[top[2]])%3==(vec[top[3]]-vec[top[4]]+vec[top[5]])%3
    other=coloring(over,accept=cyclic_old);matrix=C.DLNmatrix(3,over,eps,other)
    assert matrix[0][1]==matrix[1][0]==2,matrix
    contrast=dict(parameters=[3,2],coloring=other,top_colors=[other[j] for j in top],
        actual_old_companion_image='cyclic C2',actual_S3_birth=[3,2],author_matrix=[[str(x) for x in row] for row in matrix],
        mutual=2,cyclic_loop_coloring_same_knot_mutual=18,coloring_not_chosen_by_linking_output=True)
    print('SAME cable (3,2), old-C2 BIRTH coloring mutual2 versus old-S3 LOOP coloring mutual18')
    higher=[]
    for prime in [5,7]:
        over,eps,top=diagram(2,[1]*prime);colors=coloring(over,p=prime)
        result=C.DLNmatrix(prime,over,eps,colors)
        off=[result[i][j] for i in range(len(result)) for j in range(i)]
        assert all(x==2 for x in off)
        higher.append(dict(p=prime,knot=[2,prime],overstrands=over,signs=eps,coloring=colors,
            author_matrix=[[str(x) for x in row] for row in result],offdiagonal=[int(x) for x in off],preferred_diagonals_identified=False))
    return dict(outcome='PASS',unmodified_author_code=pin,examples=out,different_coloring_control=contrast,higher_degree_controls=higher)
def main():
    result=run()
    (P/'dihedral_author_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
