import ast,time,random,json,sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'.')
import snappy
from wirt import load, wirt_le
todo=json.load(open(sys.argv[1])); out=sys.argv[2]; trials=int(sys.argv[3])
rows,I=load(); S=set(todo); sel=[x for x in rows if x[I['name']] in S]
random.seed(7); res={}; t=time.time()
for x in sel:
    pd=ast.literal_eval(x[I['pd_notation']])
    for trial in range(trials):
        L=snappy.Link(pd); L.backtrack(random.randint(2,40))
        if random.random()<0.5: L.simplify('level')
        p=L.PD_code(); s=wirt_le(p,3)
        if s: res[x[I['name']]]={'pd':[list(c) for c in p],'seeds':list(s)}; break
    json.dump(res,open(out,'w'))
print('w3',len(sel),len(res),round(time.time()-t,1),flush=True)
