from sage.all import QQ, PolynomialRing, Qp, ZZ
import json
import resource
import time

resource.setrlimit(resource.RLIMIT_CPU,(30,35))
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
start=time.process_time()
Z=PolynomialRing(QQ,'z');z=Z.gen()
T=PolynomialRing(Z,'t');t=T.gen()
c=z*z-3;s=4*z*c;beta=c*c;gamma=-3*c*(z*z+1)
F=((t-gamma)**2-3*s*s)**2-s**4*t
Q=(t-beta)**2+4*s*z*(t-beta)+6*s*s*(z*z-1)
assert F==(t-beta)**2*Q
assert Q.discriminant()==-8*s*s*c
p=131
lift=38+p*((-((38**2-3)//p)*pow(76,-1,p))%p)
rows=[]
# Three declared controls only: odd contact, even contact, adjacent even contact.
for zz in [38,lift,lift+p*p]:
    cc=ZZ(zz*zz-3);ss=4*zz*cc;bb=cc*cc;gg=-3*cc*(zz*zz+1)
    a=cc.valuation(p);disc=-8*ss*ss*cc
    P=Qp(p,30);R=PolynomialRing(P,'t');u=R.gen()
    q=(u-bb)**2+4*ss*zz*(u-bb)+6*ss*ss*(zz*zz-1)
    roots=q.roots(multiplicities=False)
    row={'z':zz,'a':int(a),'s':str(ss),'beta':str(bb),'gamma':str(gg),
         'Q_disc_valuation':int(disc.valuation(p)),'rational_branch_roots':len(roots)}
    if len(roots)==2:
        b1,b2=roots;x=lambda b:-2+(b-gg)**2/ss**2
        row['valuations']={
            'beta':int(bb.valuation(p)),'gamma':int(gg.valuation(p)),
            'branch1':int(b1.valuation()),'branch2':int(b2.valuation()),
            'branch_separation':int((b1-b2).valuation()),
            'x_difference':int((x(b1)-x(b2)).valuation()),
            'node_offsets_sum':int((x(b1)+x(b2)-2).valuation()),
            'node_offset1':int((x(b1)-1).valuation())}
        assert row['valuations']['branch_separation']==3*a//2
        assert row['valuations']['x_difference']==a//2
        assert row['valuations']['node_offsets_sum']>a//2
    rows.append(row)
print(json.dumps({'status':'PASS','symbolic_identity':True,'controls':rows,'cpu_seconds':time.process_time()-start},indent=2))
