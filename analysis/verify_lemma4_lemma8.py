#!/usr/bin/env python3
"""Exact verification of Lemma 4 (hole-graph clique numbers), Lemma 8 (LP identity),
Theorem 6 (64-point extensions) and the values no other script prints.
Floating point occurs only in the vertex-candidate prefilter; every accepting step is
Fraction arithmetic. Self-contained: rebuilds D5, L5, Q5, R5 from the Cohn-Rajagopal
construction rules. Finally the four rebuilt configurations are compared, as sets of
points, with the committed configs/{d5,l5,q5,r5}_40.json, which ties this independent
rebuild to the shipped data. Exit code 1 on any mismatch with EXPECTED."""
import json, os, sys
from fractions import Fraction as F
from itertools import combinations, product
from collections import Counter
import numpy as np, sympy as sp, networkx as nx

EXPECTED = {
 'antipodal': {'D5':40,'L5':24,'Q5':20,'R5':12},
 'vertices':  {'D5':42,'L5':50,'Q5':92,'R5':100},
 'max_hole_ip':{'D5':F(3,4),'L5':F(3,4),'Q5':F(23,20),'R5':F(23,20)},
 'clique':    {'D5':16,'L5':16,'Q5':12,'R5':12},
}
fails=[]
def check(name,got,exp):
    ok = got==exp
    print(f"  {name:28s} computed={got!s:>8}  expected={exp!s:>8}  {'PASS' if ok else 'FAIL'}")
    if not ok: fails.append(name)

# ---------- configurations (coordinates scaled by 10 -> integers) ----------
D5=set()
for i,j in combinations(range(5),2):
    for si in (1,-1):
        for sj in (1,-1):
            v=[0]*5; v[i]=10*si; v[j]=10*sj; D5.add(tuple(v))
D5=sorted(D5)
L5=[p for p in D5 if p[4]!=10]+[s+(10,) for s in product((5,-5),repeat=4) if sum(x<0 for x in s)%2==1]
def reflect(p): s=sum(p); return tuple(F(c)-F(2*s,5) for c in p)
def qmod(X):
    return [tuple(F(c) for c in p) for p in X if sum(p)!=20]+[reflect(p) for p in X if sum(p)==-20]
cfg={'D5':[tuple(F(c) for c in p) for p in D5],'L5':[tuple(F(c) for c in p) for p in L5],'Q5':qmod(D5),'R5':qmod(L5)}
def ip(a,b): return sum(x*y for x,y in zip(a,b))/100     # true inner product


def exact_rank(rows):
    """rank over Q by Fraction Gaussian elimination"""
    M=[r[:] for r in rows]; rk=0; ncol=len(M[0]) if M else 0
    for c in range(ncol):
        piv=next((i for i in range(rk,len(M)) if M[i][c]!=0),None)
        if piv is None: continue
        M[rk],M[piv]=M[piv],M[rk]
        for i in range(len(M)):
            if i!=rk and M[i][c]!=0:
                f=M[i][c]/M[rk][c]; M[i]=[a-f*b for a,b in zip(M[i],M[rk])]
        rk+=1
    return rk

def clique_bb(adj,n):
    """own branch-and-bound (greedy colouring bound), no networkx"""
    best=[0]
    def color_bound(cand):
        order=sorted(cand); colors={}
        for v in order:
            used={colors[u] for u in adj[v] if u in colors}
            c=0
            while c in used: c+=1
            colors[v]=c
        return (max(colors.values())+1) if colors else 0
    def rec(cur,cand):
        if len(cur)>best[0]: best[0]=len(cur)
        if not cand or len(cur)+color_bound(cand)<=best[0]: return
        for v in sorted(cand):
            rec(cur+[v],cand & adj[v]); cand=cand-{v}
    rec([],set(range(n))); return best[0]

holes={}
for n,X in cfg.items():
    print(f"== {n} ==")
    S=set(X); check('antipodal points',sum(1 for a in X if tuple(-c for c in a) in S),EXPECTED['antipodal'][n])
    # vertex enumeration: float prefilter, exact re-verification
    Xf=np.array([[float(c)/10 for c in p] for p in X]); subs=np.array(list(combinations(range(40),5)))
    M=Xf[subs]; ok=np.abs(np.linalg.det(M))>1e-9
    W=np.linalg.solve(M[ok],np.ones((ok.sum(),5,1)))[:,:,0]
    # float prefilter only; from here on every step is exact:
    # rationalize each candidate, dedup as exact tuples, verify feasibility and
    # vertex-ness (rank 5 of the tight constraints) in Fraction arithmetic.
    cand=W[np.all(Xf@W.T<=1+1e-9,axis=0)]
    seen=set(); verts=[]
    for w in cand:
        wq=tuple(F(v).limit_denominator(1000) for v in w)
        if wq in seen: continue
        seen.add(wq)
        vals=[sum(x*c for x,c in zip(p,wq))/10 for p in X]
        assert all(v<=1 for v in vals), 'candidate not in P_X'
        tight=[[F(c) for c in p] for p,v in zip(X,vals) if v==1]
        assert exact_rank(tight)==5, 'candidate is not a vertex'
        verts.append(wq)
    check('polytope vertices',len(verts),EXPECTED['vertices'][n])
    norms=[sum(c*c for c in w) for w in verts]; mx=max(norms)
    check('max |w|^2',mx,F(5,4))
    H=[w for w,nn in zip(verts,norms) if nn==mx]; holes[n]=H
    check('deep holes',len(H),32)
    check('max hole inner product',max(sum(a*b for a,b in zip(u,v)) for u,v in combinations(H,2)),EXPECTED['max_hole_ip'][n])
    G=nx.Graph(); G.add_nodes_from(range(32)); adj={i:set() for i in range(32)}
    for i,j in combinations(range(32),2):
        if sum(a*b for a,b in zip(H[i],H[j]))<=F(1,4): G.add_edge(i,j); adj[i].add(j); adj[j].add(i)
    c1=max(len(c) for c in nx.find_cliques(G)); c2=clique_bb(adj,32)
    check('clique number (networkx)',c1,EXPECTED['clique'][n]); check('clique number (own B&B)',c2,EXPECTED['clique'][n])
    if n=='D5':
        cl=[set(c) for c in nx.find_cliques(G) if len(c)==16]
        par=lambda w: sum(1 for c in w if c<0)%2
        check('D5 16-cliques = parity classes',sorted(sorted(par(H[i]) for i in c) for c in cl)==[[0]*16,[1]*16],True)

# ---------- Lemma 8 ----------
print("== Lemma 8 ==")
t=sp.symbols('t'); G1=t; G2=(5*t**2-1)/4; G3=(7*t**3-3*t)/4
Fp=1+sp.Rational(30,7)*G1+sp.Rational(25,4)*G2+sp.Rational(125,28)*G3
check('identity F = 125/16 (t+3/5)^2 (t-1/5)',sp.expand(Fp-sp.Rational(125,16)*(t+sp.Rational(3,5))**2*(t-sp.Rational(1,5)))==0,True)
check('F(1)',Fp.subs(t,1),16)
even=[v for v in product((1,-1),repeat=5) if sum(x<0 for x in v)%2==0]
check('demihypercube ips in {1/5,-3/5}',{F(sum(a*b for a,b in zip(u,v)),5) for u,v in combinations(even,2)},{F(1,5),F(-3,5)})

# ---------- Theorem 6 ----------
print("== Theorem 6 ==")
for n in ('Q5','R5'):
    X=[tuple(c/10 for c in p) for p in cfg[n]]; H=holes[n]
    G=nx.Graph(); G.add_nodes_from(range(32))
    for i,j in combinations(range(32),2):
        if sum(a*b for a,b in zip(H[i],H[j]))<=F(1,4): G.add_edge(i,j)
    C=next(c for c in nx.find_cliques(G) if len(c)==12)
    pts=[(x,0) for x in X]+[(H[i],1) for i in C]+[(H[i],-1) for i in C]   # sigma * sqrt(3)/2
    ok=all(sum(a*b for a,b in zip(u,v))+(s*r)*F(3,4)<=1 for (u,s),(v,r) in combinations(pts,2))
    nm=all(sum(a*a for a in u)+abs(s)*F(3,4)==2 for u,s in pts)
    check(f'{n} ext64 valid (64 pts, norms, ips)',(len(pts)==64 and ok and nm),True)

# ---------- committed configuration files ----------
print("== configs/*_40.json vs. rebuilt configurations ==")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for n in ('D5', 'L5', 'Q5', 'R5'):
    with open(os.path.join(ROOT, 'configs', f'{n.lower()}_40.json')) as fh:
        data = json.load(fh)
    committed = {tuple(F(c) for c in v) for v in data['vectors']}
    rebuilt = {tuple(c / 10 for c in p) for p in cfg[n]}          # cfg is scaled by 10
    same = (committed == rebuilt and len(committed) == 40 and len(data['vectors']) == 40
            and F(data['norm_squared']) == 2)
    check(f'{n.lower()}_40.json == rebuilt {n} (point sets)', same, True)

print("\nOVERALL:", "FAIL "+str(fails) if fails else "ALL CHECKS PASS")
sys.exit(1 if fails else 0)
