#!/usr/bin/env python3
"""Degree-twelve Hurwitz connectivity through complete polynomial censuses.

Default replay compiles the embedded exact C++17 enumeration, checks both
complete marked polynomial orbits, and independently verifies all short
holonomy loops, their literal B6 cables, and the link between coalescences.
Requires g++ and about 1.6 GB RAM. --certificates-only verifies the retained
short words only and does not establish the finite-domain coverage.
"""
from __future__ import annotations
import argparse,ast,copy,hashlib,json,platform,re,shutil,subprocess,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
CERT=Path(__file__).with_name("degree_twelve_hurwitz_certificate.json")
CPP_SOURCE = r"""#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>
using U=uint64_t; constexpr unsigned N=11,D=990; const U VIS=U(1)<<63;
U cp(U a,U b){U c=0;for(unsigned i=0;i<N;i++)c|=((a>>(4*((b>>(4*i))&15)))&15)<<(4*i);return c;}
U iv(U a){U b=0;for(unsigned i=0;i<N;i++)b|=U(i)<<(4*((a>>(4*i))&15));return b;}
U conj(U a,U b){return cp(cp(a,b),iv(a));}
U id(){U p=0;for(unsigned i=0;i<N;i++)p|=U(i)<<(4*i);return p;}
U dt(unsigned a,unsigned b,unsigned c,unsigned d){U p=id();p^=U(a^b)<<(4*a);p^=U(a^b)<<(4*b);p^=U(c^d)<<(4*c);p^=U(c^d)<<(4*d);return p;}
U cyc(unsigned a,unsigned b,unsigned c){U p=id();p^=U(a^b)<<(4*a);p^=U(b^c)<<(4*b);p^=U(c^a)<<(4*c);return p;}
U mix(U a){a^=a>>30;a*=0xbf58476d1ce4e5b9ULL;a^=a>>27;a*=0x94d049bb133111ebULL;return a^(a>>31);}
void printperm(U p){std::cout<<"[";for(unsigned i=0;i<N;i++){if(i)std::cout<<",";std::cout<<((p>>(4*i))&15);}std::cout<<"]";}
void printword(std::vector<int> const&w){std::cout<<"[";for(size_t i=0;i<w.size();i++){if(i)std::cout<<",";std::cout<<w[i];}std::cout<<"]";}
struct DSU{std::array<unsigned,N>p;unsigned count=0;DSU(){for(unsigned i=0;i<N;i++)p[i]=i;}unsigned root(unsigned a){while(p[a]!=a){p[a]=p[p[a]];a=p[a];}return a;}bool join(unsigned a,unsigned b){a=root(a);b=root(b);if(a==b)return false;p[a]=b;return true;}};
int main(int argc,char**argv){
 bool triple=argc>1&&argv[1][0]=='3';
 std::vector<U>C,M;for(unsigned a=0;a<N;a++)for(unsigned b=a+1;b<N;b++)for(unsigned c=b+1;c<N;c++)for(unsigned d=c+1;d<N;d++){C.push_back(dt(a,b,c,d));C.push_back(dt(a,c,b,d));C.push_back(dt(a,d,b,c));}
 std::sort(C.begin(),C.end());assert(C.size()==D);assert(std::adjacent_find(C.begin(),C.end())==C.end());std::unordered_map<U,unsigned>di,mi;
 for(unsigned i=0;i<D;i++)di[C[i]]=i;
 if(triple){for(unsigned a=0;a<N;a++)for(unsigned b=a+1;b<N;b++)for(unsigned c=b+1;c<N;c++){M.push_back(cyc(a,b,c));M.push_back(cyc(a,c,b));}std::sort(M.begin(),M.end());}else M=C;
 assert(M.size()==(triple?330:990));assert(std::adjacent_find(M.begin(),M.end())==M.end());
 for(U p:C){assert(cp(p,p)==id());unsigned support=0;for(unsigned j=0;j<N;j++)support+=((p>>(4*j))&15)!=j;assert(support==4);}
 if(triple)for(U p:M){assert(cp(cp(p,p),p)==id());unsigned support=0;for(unsigned j=0;j<N;j++)support+=((p>>(4*j))&15)!=j;assert(support==3);}
 for(unsigned i=0;i<M.size();i++)mi[M[i]]=i;
 U P=0;for(unsigned i=0;i<N;i++)P|=U((i+1)%N)<<(4*i);
 std::unordered_map<U,std::vector<unsigned>>right,left;
 for(unsigned a=0;a<D;a++)for(unsigned b=0;b<D;b++)right[cp(C[a],C[b])].push_back(a|(b<<10));
 if(triple){for(unsigned a=0;a<M.size();a++)for(unsigned b=0;b<D;b++)left[cp(M[a],C[b])].push_back(a|(b<<10));}else left=right;
 size_t cap=size_t(1)<<(triple?24:26),mask=cap-1;
 std::vector<U>table(cap,0),transport(cap,0);std::vector<uint32_t>parent(cap,0);std::vector<int8_t>parentmove(cap,0);
 U count=0,initial=~U(0);
 auto insert=[&](U key){size_t s=mix(key)&mask;while(table[s]&&table[s]!=key+1)s=(s+1)&mask;assert(!table[s]);table[s]=key+1;initial=std::min(initial,key);++count;};
 for(auto const& item:left){U target=cp(iv(item.first),P);for(unsigned c=0;c<D;c++){auto it=right.find(cp(C[c],target));if(it==right.end())continue;for(unsigned ab:item.second)for(unsigned de:it->second)insert(U(ab)|(U(c)<<20)|(U(de)<<30));}}
 assert(count==(triple?3748096:14992384));
 std::cout<<"enumerated "<<count<<" marker "<<M.size()<<std::endl;
 std::vector<uint16_t>dc(D*D);for(unsigned a=0;a<D;a++)for(unsigned b=0;b<D;b++)dc[a*D+b]=di.at(conj(C[a],C[b]));
 std::vector<uint32_t>plus(M.size()*D),minus(M.size()*D);
 for(unsigned a=0;a<M.size();a++)for(unsigned b=0;b<D;b++){U aa=conj(C[b],M[a]),bb=conj(iv(aa),C[b]);plus[a*D+b]=mi.at(aa)|(di.at(bb)<<10);bb=conj(M[a],C[b]);aa=conj(bb,M[a]);minus[a*D+b]=mi.at(aa)|(di.at(bb)<<10);}
 for(unsigned a=0;a<M.size();a++)for(unsigned b=0;b<D;b++){
  uint32_t x=plus[a*D+b],y=minus[a*D+b];
  assert(minus[(x&1023)*D+(x>>10)]==(a|(b<<10)));
  assert(plus[(y&1023)*D+(y>>10)]==(a|(b<<10)));
  assert(cp(M[a],C[b])==cp(M[x&1023],C[x>>10]));
  assert(cp(M[a],C[b])==cp(M[y&1023],C[y>>10]));
 }
 auto contains=[&](U key){size_t s=mix(key)&mask;while(table[s]&&((table[s]&~VIS)!=key+1))s=(s+1)&mask;return table[s]!=0;};
 U wrong=(initial&~U(1023))|U(((initial&1023)+1)%M.size());assert(!contains(wrong));
 auto locate=[&](U key){size_t s=mix(key)&mask;while(table[s]&&((table[s]&~VIS)!=key+1))s=(s+1)&mask;assert(table[s]);return s;};
 auto path=[&](size_t slot){std::vector<int>w;while(parentmove[slot]){w.push_back(parentmove[slot]);slot=parent[slot];}std::reverse(w.begin(),w.end());return w;};
 U total=0;unsigned components=0;
 while(total<count){
  if(components){for(U slot:table)if(slot&&!(slot&VIS)){initial=slot-1;break;}}
  size_t base=locate(initial);assert(!(table[base]&VIS));table[base]|=VIS;transport[base]=id();parentmove[base]=0;
  U delta=M[initial&1023];std::vector<unsigned>outside;for(unsigned i=0;i<N;i++)if(((delta>>(4*i))&15)==i)outside.push_back(i);
  DSU dsu;unsigned classes=outside.size();std::vector<U>queue;queue.reserve(count);queue.push_back(initial);
  std::cout<<"component "<<components<<" seed ";for(unsigned j=0;j<5;j++){if(j)std::cout<<" ";printperm(j?C[(initial>>(10*j))&1023]:M[initial&1023]);}std::cout<<std::endl;
  unsigned witnesses=0;
  for(size_t head=0;head<queue.size();head++){
   U key=queue[head];size_t slot=locate(key);
   for(unsigned j=0;j<4;j++)for(int sign:{-1,1}){
    U candidate=key,edge=id();
    if(j==0){unsigned a=key&1023,b=(key>>10)&1023;uint32_t out=sign==1?plus[a*D+b]:minus[a*D+b];candidate=(key&~U(1048575))|out;if(classes>1)edge=sign==1?C[b]:C[out>>10];}
    else{unsigned a=(key>>(10*j))&1023,b=(key>>(10*(j+1)))&1023;unsigned aa=sign==1?b:dc[a*D+b],bb=sign==1?dc[b*D+a]:a;candidate=(key&~(U(1048575)<<(10*j)))|(U(aa)<<(10*j))|(U(bb)<<(10*(j+1)));}
    size_t dest=locate(candidate);int code=sign*int(j+1);
    if(!(table[dest]&VIS)){table[dest]|=VIS;queue.push_back(candidate);if(classes>1){transport[dest]=cp(edge,transport[slot]);assert(conj(transport[dest],delta)==M[candidate&1023]);parent[dest]=slot;parentmove[dest]=code;}}
    else if(classes>1){U h=cp(iv(transport[dest]),cp(edge,transport[slot]));assert(conj(h,delta)==delta);bool changed=false;for(unsigned a:outside){unsigned b=(h>>(4*a))&15;assert(std::find(outside.begin(),outside.end(),b)!=outside.end());if(dsu.join(a,b)){--classes;changed=true;}}
     if(changed){auto w=path(slot),v=path(dest);w.push_back(code);for(auto it=v.rbegin();it!=v.rend();it++)w.push_back(-*it);std::cout<<"holonomy ";printperm(h);std::cout<<" codes ";printword(w);std::cout<<" remaining "<<classes<<std::endl;++witnesses;}
    }
   }
  }
  total+=queue.size();std::cout<<"orbit_size "<<queue.size()<<" fiber_orbits "<<classes<<" witnesses "<<witnesses<<std::endl;++components;
 }
 assert(components==1&&total==count);std::cout<<"DONE components "<<components<<" total "<<total<<std::endl;
}

"""

def compose(p,q):return tuple(p[q[i]]for i in range(len(p)))
def inverse(p):
 out=[0]*len(p)
 for i,j in enumerate(p):out[j]=i
 return tuple(out)
def conjugate(a,b):return compose(compose(a,b),inverse(a))
def edges(p):return[(i,j)for i,j in enumerate(p)if i<j and p[j]==i]
def permutation(n,es):
 out=list(range(n))
 for a,b in es:out[a],out[b]=out[b],out[a]
 return tuple(out)
def product(tup):
 out=tuple(range(len(tup[0])))
 for p in tup:out=compose(out,p)
 return out
def cycle_type(p):
 unseen=set(range(len(p)));out=[]
 while unseen:
  v=next(iter(unseen));w=v;cycle=[]
  while w not in cycle:cycle.append(w);unseen.remove(w);w=p[w]
  if len(cycle)>1:out.append(len(cycle))
 return tuple(sorted(out))
def hurwitz(tup,j,sign):
 a,b=tup[j:j+2]
 pair=(b,conjugate(inverse(b),a))if sign==1 else(conjugate(a,b),a)
 return tup[:j]+pair+tup[j+2:]
def polynomial_move(tup,code):
 j=abs(code)-1;sign=1 if code>0 else-1
 if j==0:
  a,b=tup[:2]
  transport=b if sign==1 else conjugate(a,b)
  return hurwitz(hurwitz(tup,0,sign),0,sign),transport
 return hurwitz(tup,j,sign),tuple(range(11))
def split(delta):
 if cycle_type(delta)==(2,2):return edges(delta)
 assert cycle_type(delta)==(3,)
 a=next(i for i in range(11)if delta[i]!=i);b=delta[a];c=delta[b]
 return[(a,b),(b,c)]
def expand_codes(codes):
 result=[]
 for code in codes:
  if abs(code)==1:result.extend([2,1,1,2]if code>0 else[-2,-1,-1,-2])
  else:result.append((abs(code)+1)*(1 if code>0 else-1))
 return result
def verify_words(certificate):
 out={}
 for name,case in certificate["cases"].items():
  base=tuple(tuple(p)for p in case["seed"])
  assert product(base)==tuple(certificate["polynomial_product"])
  marker=base[0];outside=[a for a in range(11)if marker[a]==a]
  parent={a:a for a in outside}
  def root(a):
   while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
   return a
  witnesses=[]
  for witness in case["holonomy_loops"]:
   current=base;g=tuple(range(11))
   for code in witness["codes"]:
    current,transport=polynomial_move(current,code);g=compose(transport,g)
   assert current==base and g==tuple(witness["permutation"])
   assert conjugate(g,marker)==marker
   p,q=split(marker)
   cable=expand_codes(witness["codes"])
   for a in outside:
    aa=g[a];assert aa in parent
    parent[root(a)]=root(aa)
    original=(permutation(12,[(11,a),p]),permutation(12,[(11,a),q]))+tuple(tuple(v)+(11,)for v in base[1:])
    current=original
    for letter in cable:current=hurwitz(current,abs(letter)-1,1 if letter>0 else-1)
    assert current[2:]==original[2:]
    assert compose(current[0],current[1])==compose(original[0],original[1])
    assert set(edges(current[0]))&set(edges(current[1]))=={tuple(sorted((11,aa)))}
   assert len({root(a)for a in outside})==witness["remaining_fiber_orbits"]
   witnesses.append({"codes":witness["codes"],"literal_B6_word":cable,"attachment_permutation":list(g)})
  assert len({root(a)for a in outside})==1
  # The internal B2 action connects every two-transposition splitting.
  first,second=split(marker)
  local=(permutation(11,[first]),permutation(11,[second]))
  found=set();current=local
  while current not in found:
   found.add(current);current=hurwitz(current,0,1)
  assert len(found)==(2 if cycle_type(marker)==(2,2)else 3)
  out[name]={"outside_attachment_count":len(outside),"holonomy_fiber_orbits":1,
             "internal_splitting_count":len(found),"witnesses":witnesses}
 link=certificate["link"];original=tuple(tuple(p)for p in link["tuple"])
 assert all(cycle_type(p)==(2,2)for p in original)
 assert cycle_type(product(original))==(11,)and product(original)[0]==0
 current=original
 for letter in link["braid_word"]:current=hurwitz(current,abs(letter)-1,1 if letter>0 else-1)
 assert product(current)==product(original)
 for tup,positions,expected in[(original,link["initial_auxiliary_blocks"],(2,2)),
                                (current,link["final_auxiliary_blocks"],(3,))]:
  left,right=positions;assert right==left+1
  assert set(edges(tup[left]))&set(edges(tup[right]))=={tuple(link["auxiliary_edge"])}
  assert cycle_type(compose(tup[left],tup[right]))==expected
  assert all(p[0]==0 for i,p in enumerate(tup)if i not in positions)
 out["link"]={"braid_word":link["braid_word"],"product":list(product(original)),
              "initial_coalescence":"double_transposition","final_coalescence":"three_cycle"}
 return out
def parse_output(text):
 lines=text.splitlines();first=lines[0].split()
 assert first[0]=="enumerated"
 seeds=[line for line in lines if line.startswith("component")]
 assert len(seeds)==1
 seed=[ast.literal_eval(p)for p in re.findall(r"\[[^]]*\]",seeds[0])]
 loops=[]
 for line in lines:
  if line.startswith("holonomy"):
   p,w=[ast.literal_eval(x)for x in re.findall(r"\[[^]]*\]",line)]
   loops.append({"permutation":p,"codes":w,"remaining_fiber_orbits":int(line.split()[-1])})
 last=next(line for line in lines if line.startswith("orbit_size")).split()
 return {"raw_tuple_count":int(first[1]),"marker_class_size":int(first[3]),
         "seed":seed,"holonomy_loops":loops,"base_orbit_size":int(last[1]),
         "fiber_orbits":int(last[3]),"components":int(lines[-1].split()[2])}
def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--write-receipt",type=Path)
 parser.add_argument("--certificates-only",action="store_true")
 args=parser.parse_args()
 certificate=json.loads(CERT.read_text())
 words=verify_words(certificate)
 corrupted=copy.deepcopy(certificate)
 corrupted["cases"]["double_transposition"]["holonomy_loops"][0]["permutation"][0]=1
 try:verify_words(corrupted)
 except AssertionError:pass
 else:raise AssertionError("Corrupted holonomy was accepted")
 corrupted=copy.deepcopy(certificate);corrupted["link"]["braid_word"]=[]
 try:verify_words(corrupted)
 except AssertionError:pass
 else:raise AssertionError("Empty coalescence-change word was accepted")
 compiler=None;flags=["-O3","-std=c++17","-UNDEBUG"];coverage=None
 if not args.certificates_only:
  cxx=shutil.which("g++")
  if not cxx:raise RuntimeError("Full coverage replay requires g++")
  compiler=subprocess.run([cxx,"--version"],text=True,capture_output=True,check=True).stdout.splitlines()[0]
  scratch=ROOT/".scratch";scratch.mkdir(exist_ok=True)
  with tempfile.TemporaryDirectory(prefix="degree12-hurwitz-",dir=scratch)as temporary:
   source=Path(temporary)/"census.cpp";binary=Path(temporary)/"census"
   source.write_text(CPP_SOURCE)
   built=subprocess.run([cxx,*flags,str(source),"-o",str(binary)],text=True,capture_output=True)
   if built.returncode:raise RuntimeError(built.stderr)
   coverage={}
   for name,arguments in[("double_transposition",[]),("three_cycle",["3"])]:
    run=subprocess.run([str(binary),*arguments],text=True,capture_output=True)
    if run.returncode:raise RuntimeError(run.stdout+run.stderr)
    actual=parse_output(run.stdout)
    assert actual==certificate["cases"][name]
    coverage[name]={k:actual[k]for k in("raw_tuple_count","base_orbit_size","fiber_orbits","components")}
 report={
  "schema_version":1,"kind":"complete_degree12_Hurwitz_reduction_certificate"if coverage else"short_word_certificates_only",
  "python":platform.python_version(),"compiler":compiler,"compiler_flags":flags,
  "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  "embedded_cpp_sha256":hashlib.sha256(CPP_SOURCE.encode()).hexdigest(),
  "certificate_sha256":hashlib.sha256(CERT.read_bytes()).hexdigest(),
  "cycle_folding_source_sha256":hashlib.sha256((ROOT/"research/arithmetic/ALTERNATING_HURWITZ_CYCLE_FOLDING.md").read_bytes()).hexdigest(),
  "complete_polynomial_coverage":coverage,"independent_word_verification":words,
  "scope":"Complete degree12 absolute full-A12 orbit proof uses the uniform cycle-folding theorem, both exhaustive marked degree11 polynomial censuses, attachment holonomy transitivity, internal split transitivity and the literal B6 coalescence-change word. No uniform higher-degree connectivity or rational-point theorem follows.",
  "coverage_replayed":coverage is not None,"independent_replay":False,
  "external_review":False,"formal_verification":False,"result":"PASS"
 }
 output=json.dumps(report,indent=2,sort_keys=True)+"\n"
 if args.write_receipt:args.write_receipt.write_text(output)
 print(output,end="")
if __name__=="__main__":main()
