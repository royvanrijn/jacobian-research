#!/usr/bin/env python3
"""Exact absolute cacti censuses for the degree-fourteen Hurwitz theorem.

Default replay runs both complete degree-thirteen domains, sequentially.
It needs about 12 GiB of RAM and several minutes. --quick runs only the
smaller algorithm controls; it cannot certify the degree-fourteen theorem.
The two C++ programs below are retained literal source, compiled with
assertions enabled into ignored .scratch/. No detached processes are used.
"""
from __future__ import annotations
import argparse, hashlib, json, os, platform, re, subprocess, sys, time
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import comb, factorial
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
NOTE=HERE/"DEGREE_FOURTEEN_HURWITZ_ORBIT.md"

def sha(data):
    if isinstance(data,Path):data=data.read_bytes()
    if isinstance(data,str):data=data.encode()
    return hashlib.sha256(data).hexdigest()

def partitions(n,upper=None):
    if n==0:
        yield ()
        return
    if upper is None:upper=n
    for a in range(min(n,upper),0,-1):
        for tail in partitions(n-a,a):yield (a,)+tail

def exterior_characters(parts):
    n=sum(parts)
    co=[1]
    for d in parts:
        out=co+[0]*d
        for j,a in enumerate(co):out[j+d]-=(-1)**d*a
        co=out
    # Divide the permutation-representation determinant by 1+t.
    quotient=[]
    for k in range(n):quotient.append(co[k]-(quotient[-1]if k else 0))
    assert co[n]==quotient[-1]
    return quotient

def character_bound(n):
    r=(n-1)//2
    rows=[]
    inner=[[Fraction(0)for _ in range(n)]for _ in range(n)]
    for parts in partitions(n):
        mult=Counter(parts)
        centralizer=1
        for d,m in mult.items():centralizer*=d**m*factorial(m)
        assert factorial(n)%centralizer==0
        char=exterior_characters(parts)
        rows.append([list(parts),centralizer,char])
        for i in range(n):
            for j in range(n):inner[i][j]+=Fraction(char[i]*char[j],centralizer)
    assert inner==[[Fraction(i==j)for j in range(n)]for i in range(n)]
    dims=exterior_characters((1,)*n)
    assert dims==[comb(n-1,k)for k in range(n)]
    cycle=exterior_characters((n,))
    assert cycle==[(-1)**k for k in range(n)]
    assert sum(a*a for a in cycle)==n
    dchar=exterior_characters((2,2)+(1,)*(n-4))
    tchar=exterior_characters((3,)+(1,)*(n-3))
    def bk(a,b):return comb(a,b)if 0<=b<=a else 0
    assert dchar==[bk(n-3,k)-2*bk(n-3,k-1)+bk(n-3,k-2)for k in range(n)]
    assert tchar==[bk(n-4,k)+bk(n-4,k-3)for k in range(n)]
    D,T=3*comb(n,4),2*comb(n,3)
    countD=Fraction(D**r,factorial(n))*sum(
        Fraction((-1)**k*dchar[k]**r,dims[k]**(r-1))for k in range(n))
    countT=Fraction(T*D**(r-1),factorial(n))*sum(
        Fraction((-1)**k*tchar[k]*dchar[k]**(r-1),dims[k]**(r-1))for k in range(n))
    assert countD.denominator==countT.denominator==1
    assert countD==n**(r-1)*(r-1)**r
    assert countT==n**(r-1)*(r-1)**(r-1)
    assert countD%n==countT%n==0
    return dict(polynomial_degree=n,r=r,cycle_types=len(rows),
        exterior_inner_product_identity=True,cycle_column_exhausted=True,
        dimensions=dims,DT_characters=dchar,T3_characters=tchar,
        cycle_type_character_table_sha256=sha(json.dumps(rows,separators=(",",":"))),
        literal_DT_count=int(countD),literal_marked_T3_count=int(countT),
        absolute_DT_count=int(countD/n),absolute_marked_T3_count=int(countT/n))

def compose(p,q):return tuple(p[q[i]]for i in range(len(p)))
def inverse(p):
    out=[0]*len(p)
    for i,a in enumerate(p):out[a]=i
    return tuple(out)
def perm(n,cycles):
    out=list(range(n))
    for c in cycles:
        for i,a in enumerate(c):out[a]=c[(i+1)%len(c)]
    return tuple(out)
def conjugate(p,g):return compose(compose(p,g),inverse(p))
def product(t):
    out=tuple(range(len(t[0])))
    for p in t:out=compose(out,p)
    return out
def power(p,k):
    out=tuple(range(len(p)))
    for _ in range(k):out=compose(out,p)
    return out
def hmove(t,j):
    a,b=t[j:j+2]
    return t[:j]+(b,conjugate(inverse(b),a))+t[j+2:]

def seed_controls(n):
    r=(n-1)//2
    cyc=perm(n,[tuple(range(n))])
    d=perm(n,[(0,2),(3,4)])
    seed=tuple(conjugate(power(cyc,2*j),d)for j in range(r))
    assert product(seed)==cyc
    if r>=4:
        e=compose(seed[0],power(compose(seed[0],seed[1]),2))
        assert e==perm(n,[(0,4),(2,3)])
        assert power(compose(e,seed[2]),2)==perm(n,[(0,6,4)])
    raw=(perm(n,[(0,1,2)]),)+tuple(
        perm(n,[(0,2*j+1),(1,2*j+2)])for j in range(1,r))
    pr=product(raw)
    labels=[-1]*n;a=0
    for i in range(n):
        assert labels[a]<0
        labels[a]=i;a=pr[a]
    assert a==0
    t=tuple(conjugate(tuple(labels),g)for g in raw)
    assert product(t)==cyc
    # Literal convention check for the marked first generator Q1^2.
    a,b=t[:2]
    aa=conjugate(b,a)
    bb=conjugate(inverse(aa),b)
    assert hmove(hmove(t,0),0)==(aa,bb)+t[2:]
    return dict(degree=n,DT_seed_cycles=[[[2*j%n,(2*j+2)%n],
          [(2*j+3)%n,(2*j+4)%n]]for j in range(r)],
        DT_seed_full_cycle=True,DT_three_cycle_certificate=r>=4,
        T3_seed=list(map(list,t)),T3_seed_full_cycle=True,
        marked_double_move_convention_verified=True)

def small_independent_count():
    n=7
    D=[]
    for a,b,c,d in combinations(range(n),4):
        for cs in [[(a,b),(c,d)],[(a,c),(b,d)],[(a,d),(b,c)]]:
            D.append(perm(n,cs))
    Ds=set(D);P=perm(n,[tuple(range(n))])
    tuples={(a,b,compose(inverse(compose(a,b)),P))for a in D for b in D
            if compose(inverse(compose(a,b)),P)in Ds}
    assert len(tuples)==392
    # All 392 literal tuples and both positive elementary moves are checked.
    for t in tuples:
        for j in range(2):
            u=hmove(t,j)
            assert u in tuples and product(u)==P
    cyc=perm(n,[tuple(range(n))]);d=perm(n,[(0,2),(3,4)])
    seed=tuple(conjugate(power(cyc,2*j),d)for j in range(3))
    group={tuple(range(n))};queue=list(group)
    for g in queue:
        for a in seed:
            h=compose(g,a)
            if h not in group:group.add(h);queue.append(h)
    assert len(group)==168
    return dict(degree=7,complete_DT_triples=len(tuples),
        all_positive_neighbors_inside_domain=True,exceptional_DT_seed_group_order=len(group))

def run_case(executable,kind,n,capbits,bound,logdir):
    command=[str(executable),str(n),str(capbits)]
    before=time.monotonic()
    proc=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=None,text=True)
    out,_=proc.communicate()
    elapsed=time.monotonic()-before
    if proc.returncode:raise RuntimeError(f"{kind} n{n} exited {proc.returncode}")
    match=re.search(r"N (\d+) D (\d+) reps (\d+) keybits (\d+) absolute_orbit (\d+) seconds ([0-9.e+-]+) literal (\d+)",out)
    assert match,out
    degree,D,reps,bits,orbit,seconds,literal=match.groups()
    degree,D,reps,bits,orbit,literal=map(int,(degree,D,reps,bits,orbit,literal))
    assert degree==n and D==3*comb(n,4)and literal==orbit*n
    expected=bound["absolute_DT_count"if kind=="DT"else"absolute_marked_T3_count"]
    if n==7 and kind=="DT":
        assert orbit==7 and orbit<expected # Retained small exceptional boundary.
    else:assert orbit==expected,(kind,n,orbit,expected)
    rss=re.search(r"maxrss_kib (\d+)",out)
    report=dict(kind=kind,polynomial_degree=n,command=command,
        absolute_orbit=orbit,literal_orbit=literal,absolute_total_bound=expected,
        complete=orbit==expected,DT_class_size=D,first_factor_representatives=reps,
        stored_key_bits=bits,hash_capacity=1<<capbits,
        BFS_seconds=float(seconds),wall_seconds=round(elapsed,6),
        maximum_RSS_KiB=int(rss.group(1))if rss else None,executable_sha256=sha(executable),
        raw_stdout=out)
    (logdir/f"{kind.lower()}_{n}_stdout.txt").write_text(out)
    return report

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quick",action="store_true",
        help="Run n7/9/11 controls only; does not certify the degree14 result.")
    ap.add_argument("--write-receipt",type=Path)
    ap.add_argument("--collect-first-runs",type=Path,
        help="Bind already completed foreground censuses by source/binary hashes; rerun smaller controls only.")
    args=ap.parse_args()
    scratch=ROOT/".scratch/degree14_absolute_census"
    scratch.mkdir(parents=True,exist_ok=True)
    comp=subprocess.check_output(["g++","--version"],text=True).splitlines()[0]
    flags=["-O3","-std=c++17","-Wall","-Wextra","-UNDEBUG"]
    binaries={}
    for kind,body in (("DT",DT_CPP),("T3",T3_CPP)):
        source=scratch/(kind.lower()+".cpp")
        executable=scratch/(kind.lower()+".bin")
        source.write_text(body)
        subprocess.run(["g++",*flags,str(source),"-o",str(executable)],check=True)
        binaries[kind]=executable
    bounds={n:character_bound(n)for n in(7,9,11,13)}
    assert bounds[13]["absolute_DT_count"]==446265625
    assert bounds[13]["absolute_marked_T3_count"]==89253125
    controls=small_independent_count()
    seeds=[seed_controls(n)for n in(7,9,11,13)]
    cases=[]
    for n in(7,9,11):
        for kind in("DT","T3"):
            caps={("DT",7):7,("DT",9):14,("DT",11):22,
                  ("T3",7):8,("T3",9):13,("T3",11):20}
            cases.append(run_case(binaries[kind],kind,n,caps[kind,n],bounds[n],scratch))
    if not args.quick:
        # Sequential high-memory cases: the first process has exited before
        # the second allocates its table. Both foreground handles are joined.
        if args.collect_first_runs:
            retained=json.loads(args.collect_first_runs.read_text())
            assert retained["schema_version"]==1 and retained["complete_foreground_runs"]
            assert [row["kind"]for row in retained["cases"]]==["DT","T3"]
            for row in retained["cases"]:
                kind=row["kind"]
                assert row["CXX_source_sha256"]==sha(DT_CPP if kind=="DT"else T3_CPP)
                assert sha(ROOT/row["original_source_path"])==row["CXX_source_sha256"]
                assert sha(ROOT/row["original_executable_path"])==row["executable_sha256"]
                expected=bounds[13]["absolute_DT_count"if kind=="DT"else"absolute_marked_T3_count"]
                assert row["absolute_orbit"]==row["absolute_total_bound"]==expected
                assert row["literal_orbit"]==13*expected
                assert row["complete"]and row["polynomial_degree"]==13
                match=re.search(r"absolute_orbit ([0-9]+) seconds ([0-9.e+-]+) literal ([0-9]+)",row["raw_stdout"])
                assert match and int(match.group(1))==expected and int(match.group(3))==13*expected
                assert float(match.group(2))==row["BFS_seconds"]
                cases.append(row)
        else:
            for kind,cap in(("DT",30),("T3",28)):
                cases.append(run_case(binaries[kind],kind,13,cap,bounds[13],scratch))
    paths=[NOTE,Path(__file__).resolve(),
        HERE/"ALTERNATING_HURWITZ_CYCLE_FOLDING.md",
        HERE/"UNMARKED_POLYNOMIAL_HURWITZ_REDUCTION.md",
        HERE/"degree_fourteen_hurwitz_census_inputs.json",
        HERE/"sources/jones-zvonkin-braid-cacti-2002.pdf",
        HERE/"sources/jones-primitive-cycle-1209.5169v1.pdf"]
    payload=dict(schema_version=1,result="PASS",full_degree14_census=not args.quick,
        python=platform.python_version(),platform=platform.platform(),
        compiler=comp,compile_flags=flags,CXX_sources_sha256={
            "DT":sha(DT_CPP),"T3":sha(T3_CPP)},checker_sha256=sha(Path(__file__)),
        source_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        character_bounds=list(bounds.values()),small_independent_count=controls,
        seed_controls=seeds,cases=cases,
        large_execution_mode=("quick controls only"if args.quick else
            "retained completed foreground runs, checked by source/binary hashes"if args.collect_first_runs else
            "large engines executed by this wrapper"),
        collected_record_sha256=sha(args.collect_first_runs)if args.collect_first_runs else None,
        scope="Finite n13 polynomial orbit classification and its written n14 transfer. No all-r connectivity claim.",
        formal_verification=False,external_review=False,independent_external_replay=False)
    dest=args.write_receipt
    if dest is None and not args.quick:dest=HERE/"degree_fourteen_hurwitz_verification.json"
    if dest:
        dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_text(json.dumps(payload,indent=2)+"\n")
    print("PASS: "+("all degree13 absolute cacti; degree14 transfer"if not args.quick else "smaller controls only")+"; "+str(dest))

DT_CPP = r'''#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>
using U=uint64_t;
struct DT{unsigned a,b,c,d; U p;};
U mix(U a){a^=a>>30;a*=0xbf58476d1ce4e5b9ULL;a^=a>>27;a*=0x94d049bb133111ebULL;return a^(a>>31);}
int main(int argc,char**argv){
 unsigned n=argc>1?std::stoul(argv[1]):9,r=(n-1)/2,capbits=argc>2?std::stoul(argv[2]):22;
 assert(n>=7&&n<=13&&n%2==1);
 U I=0,P=0;for(unsigned i=0;i<n;i++){I|=U(i)<<(4*i);P|=U((i+1)%n)<<(4*i);}
 auto dtp=[&](unsigned a,unsigned b,unsigned c,unsigned d){U p=I;p^=U(a^b)<<(4*a);p^=U(a^b)<<(4*b);p^=U(c^d)<<(4*c);p^=U(c^d)<<(4*d);return p;};
 auto normkey=[](unsigned a,unsigned b,unsigned c,unsigned d){if(a>b)std::swap(a,b);if(c>d)std::swap(c,d);if(a>c){std::swap(a,c);std::swap(b,d);}return a|(b<<4)|(c<<8)|(d<<12);};
 std::vector<DT>D;
 for(unsigned a=0;a<n;a++)for(unsigned b=a+1;b<n;b++)for(unsigned c=b+1;c<n;c++)for(unsigned d=c+1;d<n;d++){
  for(auto z:{std::array<unsigned,4>{a,b,c,d},std::array<unsigned,4>{a,c,b,d},std::array<unsigned,4>{a,d,b,c}})D.push_back({z[0],z[1],z[2],z[3],dtp(z[0],z[1],z[2],z[3])});
 }
 std::sort(D.begin(),D.end(),[](DT const&a,DT const&b){return a.p<b.p;});unsigned countD=D.size();
 std::vector<uint16_t>lut(65536,65535);
 for(unsigned i=0;i<countD;i++){auto d=D[i];lut[normkey(d.a,d.b,d.c,d.d)]=i;}
 auto gid=[&](U p){U diff=p^I;assert(diff);unsigned a=__builtin_ctzll(diff)/4,b=(p>>(4*a))&15;diff&=~((U(15)<<(4*a))|(U(15)<<(4*b)));assert(diff);unsigned c=__builtin_ctzll(diff)/4,d=(p>>(4*c))&15;unsigned id=lut[normkey(a,b,c,d)];assert(id!=65535&&D[id].p==p);return id;};
 std::vector<uint16_t>dc(countD*countD),rot(n*countD);
 for(unsigned a=0;a<countD;a++)for(unsigned b=0;b<countD;b++){auto d=D[b];auto im=[&](unsigned t){return unsigned((D[a].p>>(4*t))&15);};dc[a*countD+b]=lut[normkey(im(d.a),im(d.b),im(d.c),im(d.d))];assert(dc[a*countD+b]!=65535);}
 for(unsigned k=0;k<n;k++)for(unsigned i=0;i<countD;i++){auto d=D[i];rot[k*countD+i]=lut[normkey((d.a+k)%n,(d.b+k)%n,(d.c+k)%n,(d.d+k)%n)];}
 std::vector<uint16_t>canon(countD),shift(countD),repnum(countD,65535),reps;
 for(unsigned i=0;i<countD;i++){unsigned j=i,k0=0;for(unsigned k=1;k<n;k++)if(rot[k*countD+i]<j){j=rot[k*countD+i];k0=k;}canon[i]=j;shift[i]=k0;}
 for(unsigned i=0;i<countD;i++)if(canon[i]==i){repnum[i]=reps.size();reps.push_back(i);}
 assert(reps.size()*n==countD);
 unsigned bits=0,fbits=0;while((1U<<bits)<countD)bits++;while((1U<<fbits)<reps.size())fbits++;assert(fbits+(r-2)*bits<63);
 U dmask=(U(1)<<bits)-1,fmask=(U(1)<<fbits)-1;
 auto encode=[&](std::array<unsigned,6>t){unsigned k=shift[t[0]];U key=repnum[canon[t[0]]];for(unsigned j=1;j<r-1;j++)key|=U(rot[k*countD+t[j]])<<(fbits+bits*(j-1));return key;};
 auto decode=[&](U key){std::array<unsigned,6>t{};t[0]=reps[key&fmask];for(unsigned j=1;j<r-1;j++)t[j]=(key>>(fbits+bits*(j-1)))&dmask;U prefix=I;
  for(unsigned j=0;j<r-1;j++){auto d=D[t[j]];for(auto ab:{std::array<unsigned,2>{d.a,d.b},std::array<unsigned,2>{d.c,d.d}}){unsigned a=ab[0],b=ab[1];U v=((prefix>>(4*a))^(prefix>>(4*b)))&15;prefix^=(v<<(4*a))|(v<<(4*b));}}
  U inv=0;for(unsigned i=0;i<n;i++)inv|=U(i)<<(4*((prefix>>(4*i))&15));U last=(inv>>4)|((inv&15)<<(4*(n-1)));t[r-1]=gid(last);return t;};
 std::array<unsigned,6>seed{};for(unsigned j=0;j<r;j++){unsigned k=(2*j)%n;seed[j]=lut[normkey(k,(k+2)%n,(k+3)%n,(k+4)%n)];}
 U key=encode(seed);auto check=decode(key);assert(encode(check)==key);
 size_t cap=size_t(1)<<capbits,mask=cap-1;std::vector<U>table(cap,0),queue;queue.reserve(cap/2);
 auto insert=[&](U k){size_t slot=mix(k)&mask;while(table[slot]&&table[slot]!=k+1)slot=(slot+1)&mask;if(table[slot])return false;table[slot]=k+1;queue.push_back(k);assert(queue.size()*10<cap*8);return true;};
 auto start=std::chrono::steady_clock::now();insert(key);
 for(size_t head=0;head<queue.size();head++){
  U at=queue[head];auto t=decode(at);
  for(unsigned j=0;j<r-1;j++){auto u=t;u[j]=t[j+1];u[j+1]=dc[t[j+1]*countD+t[j]];insert(encode(u));}
  if(head&&head%1000000==0){auto sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();std::cerr<<"processed "<<head<<" seen "<<queue.size()<<" seconds "<<sec<<std::endl;}
 }
 auto sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 std::cout<<"N "<<n<<" D "<<countD<<" reps "<<reps.size()<<" keybits "<<fbits+(r-2)*bits<<" absolute_orbit "<<queue.size()<<" seconds "<<sec<<" literal "<<queue.size()*n<<std::endl;
 U a=1,b=1;for(unsigned j=0;j<r-1;j++)a*=n;for(unsigned j=0;j<r;j++)b*=r-1;assert(queue.size()*n==a*b||n==7);
}
'''

T3_CPP = r'''#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>
#include <unordered_map>
#include <sys/resource.h>
using U=uint64_t;
struct DT{unsigned a,b,c,d; U p;};
U mix(U a){a^=a>>30;a*=0xbf58476d1ce4e5b9ULL;a^=a>>27;a*=0x94d049bb133111ebULL;return a^(a>>31);}
int main(int argc,char**argv){
 unsigned n=argc>1?std::stoul(argv[1]):11,r=(n-1)/2,capbits=argc>2?std::stoul(argv[2]):20;
 assert(n>=7&&n<=13&&n%2==1);
 U I=0,P=0;for(unsigned i=0;i<n;i++){I|=U(i)<<(4*i);P|=U((i+1)%n)<<(4*i);}
 auto dtp=[&](unsigned a,unsigned b,unsigned c,unsigned d){U p=I;p^=U(a^b)<<(4*a);p^=U(a^b)<<(4*b);p^=U(c^d)<<(4*c);p^=U(c^d)<<(4*d);return p;};
 auto normkey=[](unsigned a,unsigned b,unsigned c,unsigned d){if(a>b)std::swap(a,b);if(c>d)std::swap(c,d);if(a>c){std::swap(a,c);std::swap(b,d);}return a|(b<<4)|(c<<8)|(d<<12);};
 std::vector<DT>D;
 for(unsigned a=0;a<n;a++)for(unsigned b=a+1;b<n;b++)for(unsigned c=b+1;c<n;c++)for(unsigned d=c+1;d<n;d++){
  for(auto z:{std::array<unsigned,4>{a,b,c,d},std::array<unsigned,4>{a,c,b,d},std::array<unsigned,4>{a,d,b,c}})D.push_back({z[0],z[1],z[2],z[3],dtp(z[0],z[1],z[2],z[3])});
 }
 std::sort(D.begin(),D.end(),[](DT const&a,DT const&b){return a.p<b.p;});unsigned countD=D.size();
 std::vector<uint16_t>lut(65536,65535);
 for(unsigned i=0;i<countD;i++){auto d=D[i];lut[normkey(d.a,d.b,d.c,d.d)]=i;}
 auto gid=[&](U p){U diff=p^I;assert(diff);unsigned a=__builtin_ctzll(diff)/4,b=(p>>(4*a))&15;diff&=~((U(15)<<(4*a))|(U(15)<<(4*b)));assert(diff);unsigned c=__builtin_ctzll(diff)/4,d=(p>>(4*c))&15;unsigned id=lut[normkey(a,b,c,d)];assert(id!=65535&&D[id].p==p);return id;};
 std::vector<uint16_t>dc(countD*countD),rot(n*countD);
 for(unsigned a=0;a<countD;a++)for(unsigned b=0;b<countD;b++){auto d=D[b];auto im=[&](unsigned t){return unsigned((D[a].p>>(4*t))&15);};dc[a*countD+b]=lut[normkey(im(d.a),im(d.b),im(d.c),im(d.d))];assert(dc[a*countD+b]!=65535);}
 for(unsigned k=0;k<n;k++)for(unsigned i=0;i<countD;i++){auto d=D[i];rot[k*countD+i]=lut[normkey((d.a+k)%n,(d.b+k)%n,(d.c+k)%n,(d.d+k)%n)];}
 std::vector<U>M;
 for(unsigned a=0;a<n;a++)for(unsigned b=a+1;b<n;b++)for(unsigned c=b+1;c<n;c++)for(int rev:{0,1}){
 U p=I;unsigned v[3]={a,rev?c:b,rev?b:c};for(unsigned i=0;i<3;i++)p^=U(v[i]^v[(i+1)%3])<<(4*v[i]);M.push_back(p);}
 std::sort(M.begin(),M.end());unsigned countM=M.size();std::unordered_map<U,unsigned>mi;
 for(unsigned i=0;i<countM;i++)mi[M[i]]=i;
 auto inverse=[&](U p){U out=0;for(unsigned i=0;i<n;i++)out|=U(i)<<(4*((p>>(4*i))&15));return out;};
 auto cp=[&](U a,U b){U out=0;for(unsigned i=0;i<n;i++)out|=((a>>(4*((b>>(4*i))&15)))&15)<<(4*i);return out;};
 std::vector<uint16_t>mr(n*countM),pa(countM*countD),pb(countM*countD);
 for(unsigned k=0;k<n;k++)for(unsigned i=0;i<countM;i++){
 U p=0;for(unsigned a=0;a<n;a++)p|=U((((M[i]>>(4*a))&15)+k)%n)<<(4*((a+k)%n));mr[k*countM+i]=mi.at(p);}
 for(unsigned a=0;a<countM;a++)for(unsigned b=0;b<countD;b++){
 U aa=cp(cp(D[b].p,M[a]),D[b].p),bb=cp(cp(inverse(aa),D[b].p),aa);
 pa[a*countD+b]=mi.at(aa);pb[a*countD+b]=gid(bb);
 assert(cp(aa,bb)==cp(M[a],D[b].p));}
 std::vector<uint16_t>canon(countM),repnum(countM,65535),reps;std::vector<std::vector<unsigned>>shifts(countM);
 for(unsigned i=0;i<countM;i++){unsigned j=i;for(unsigned k=1;k<n;k++)j=std::min(j,unsigned(mr[k*countM+i]));canon[i]=j;for(unsigned k=0;k<n;k++)if(mr[k*countM+i]==j)shifts[i].push_back(k);}
 for(unsigned i=0;i<countM;i++)if(canon[i]==i){repnum[i]=reps.size();reps.push_back(i);}
 unsigned bits=0,fbits=0;while((1U<<bits)<countD)bits++;while((1U<<fbits)<reps.size())fbits++;assert(fbits+(r-2)*bits<63);
 U dmask=(U(1)<<bits)-1,fmask=(U(1)<<fbits)-1;
 auto encode=[&](std::array<unsigned,6>t){U best=~U(0);for(unsigned k:shifts[t[0]]){U key=repnum[canon[t[0]]];for(unsigned j=1;j<r-1;j++)key|=U(rot[k*countD+t[j]])<<(fbits+bits*(j-1));best=std::min(best,key);}return best;};
 auto decode=[&](U key){std::array<unsigned,6>t{};t[0]=reps[key&fmask];for(unsigned j=1;j<r-1;j++)t[j]=(key>>(fbits+bits*(j-1)))&dmask;U prefix=M[t[0]];
  for(unsigned j=1;j<r-1;j++){auto d=D[t[j]];for(auto ab:{std::array<unsigned,2>{d.a,d.b},std::array<unsigned,2>{d.c,d.d}}){unsigned a=ab[0],b=ab[1];U v=((prefix>>(4*a))^(prefix>>(4*b)))&15;prefix^=(v<<(4*a))|(v<<(4*b));}}
  U inv=0;for(unsigned i=0;i<n;i++)inv|=U(i)<<(4*((prefix>>(4*i))&15));U last=(inv>>4)|((inv&15)<<(4*(n-1)));t[r-1]=gid(last);return t;};
 std::array<unsigned,6>seed{};std::array<U,6>raw{};raw[0]=I;
 for(unsigned i=0;i<3;i++)raw[0]^=U(i^((i+1)%3))<<(4*i);
 for(unsigned j=1;j<r;j++)raw[j]=dtp(0,2*j+1,1,2*j+2);
 U prod=I;for(unsigned j=0;j<r;j++)prod=cp(prod,raw[j]);
 U sigma=0;unsigned a=0;for(unsigned i=0;i<n;i++){sigma|=U(i)<<(4*a);a=(prod>>(4*a))&15;}assert(a==0&&cp(cp(sigma,prod),inverse(sigma))==P);
 raw[0]=cp(cp(sigma,raw[0]),inverse(sigma));seed[0]=mi.at(raw[0]);
 for(unsigned j=1;j<r;j++)seed[j]=gid(cp(cp(sigma,raw[j]),inverse(sigma)));
 U key=encode(seed);auto check=decode(key);assert(encode(check)==key);
 size_t cap=size_t(1)<<capbits,mask=cap-1;std::vector<U>table(cap,0),queue;queue.reserve(cap/2);
 auto insert=[&](U k){size_t slot=mix(k)&mask;while(table[slot]&&table[slot]!=k+1)slot=(slot+1)&mask;if(table[slot])return false;table[slot]=k+1;queue.push_back(k);assert(queue.size()*10<cap*8);return true;};
 auto start=std::chrono::steady_clock::now();insert(key);
 for(size_t head=0;head<queue.size();head++){
  U at=queue[head];auto t=decode(at);
  for(unsigned j=0;j<r-1;j++){auto u=t;if(j==0){u[0]=pa[t[0]*countD+t[1]];u[1]=pb[t[0]*countD+t[1]];}else{u[j]=t[j+1];u[j+1]=dc[t[j+1]*countD+t[j]];}insert(encode(u));}
  if(head&&head%10000000==0){auto sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();std::cerr<<"processed "<<head<<" seen "<<queue.size()<<" seconds "<<sec<<std::endl;}
 }
 auto sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 std::cout<<"N "<<n<<" D "<<countD<<" reps "<<reps.size()<<" keybits "<<fbits+(r-2)*bits<<" absolute_orbit "<<queue.size()<<" seconds "<<sec<<" literal "<<queue.size()*n<<std::endl;struct rusage usage;getrusage(RUSAGE_SELF,&usage);std::cout<<"maxrss_kib "<<usage.ru_maxrss<<std::endl;
 U countA=1,countB=1;for(unsigned j=0;j<r-1;j++)countA*=n;for(unsigned j=0;j<r-1;j++)countB*=r-1;assert(queue.size()*n==countA*countB||n==7);
}
'''

if __name__=="__main__":
    main()
