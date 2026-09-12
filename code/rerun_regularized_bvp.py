#!/usr/bin/env python3
from __future__ import annotations
import csv, json, math, platform, sys, time
from pathlib import Path
import numpy as np, scipy, matplotlib
matplotlib.use("Agg")
from scipy.integrate import solve_bvp
from scipy.optimize import brentq
from scipy.special import lambertw
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
DATA=ROOT/'data'
RESULTS=ROOT/'results'
FIG=ROOT/'figures'
RESULTS.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)
ARCHIVED=DATA/'regularized_bvp_results_archived.csv'
EPSILONS=(1e-4,1e-5,1e-6,1e-7,1e-8); XMAX=20.0; TOL=1e-6; BC_TOL=1e-10; MAX_NODES=60000
A0,A1,A2=0.3933180528,1.1065076820,0.0131925717
X0,X1,X2=1.1090768620,1.4615818195,-1.6180901746
K0,K1,K2=0.6898804165,0.3022408770,0.6573620821

def L_eps(e): return float(lambertw(3/(8*e)).real/3)
def asym(e):
    L=L_eps(e); ao1=A0+A1/L; ao2=ao1+A2/L**2; xo1=X0+X1/L; xo2=xo1+X2/L**2; ko1=K0+K1/L; ko2=ko1+K2/L**2
    return dict(L=L,ao1=ao1,ao2=ao2,xo1=xo1,xo2=xo2,ko1=ko1,ko2=ko2,ov=(e*L)**(1/3)*ao2,curv=e**(-1/3)*L**(2/3)*ko2)
def rates(e):
    d=math.sqrt(1-4*e); return math.sqrt(2/(1+d)), math.sqrt((1+d)/(2*e))
def funcs(e):
    def f(x,y):
        U,Up,Upp,Uppp=y; return np.vstack((Up,Upp,Uppp,((1-U*U)*Upp-U)/e))
    def j(x,y):
        U,Up,Upp,Uppp=y; J=np.zeros((4,4,x.size)); J[0,1]=1; J[1,2]=1; J[2,3]=1; J[3,0]=(-2*U*Upp-1)/e; J[3,2]=(1-U*U)/e; return J
    return f,j
def bcs(e):
    ls,lf=rates(e); s=ls+lf; p=ls*lf
    def bc(a,b): return np.array([a[1],a[3],b[2]+s*b[1]+p*b[0],b[3]+s*b[2]+p*b[1]])
    def bj(a,b):
        A=np.zeros((4,4)); B=np.zeros((4,4)); A[0,1]=1; A[1,3]=1; B[2]=[p,s,1,0]; B[3]=[0,p,s,1]; return A,B
    return bc,bj
def mesh(e,xmax):
    L=L_eps(e); delta=e**(1/3)*L**(-1/6); fast=math.sqrt(e); ce=min(.8,max(12*delta,40*fast,.02)); hc=min(fast/5,delta/30,2e-3); nc=max(220,min(3200,int(math.ceil(ce/hc))+1))
    parts=[np.linspace(0,ce,nc)]
    if ce<2: parts.append(np.linspace(ce,2,650))
    if xmax>2: parts.append(np.linspace(2,xmax,550))
    return np.unique(np.concatenate(parts))
def seed(e,x):
    a=asym(e); A=1+a['ov']; rad=A/a['curv']; r=np.sqrt(x*x+rad*rad); U=A*np.exp(rad-r); q=-x/r; qp=-rad*rad/r**3; qpp=3*rad*rad*x/r**5
    return np.vstack((U,U*q,U*(q*q+qp),U*(q**3+3*q*qp+qpp)))
def solve_one(e,xmax=XMAX,guess=None):
    x=mesh(e,xmax); y=seed(e,x) if guess is None else guess(x); f,j=funcs(e); bc,bj=bcs(e); t=time.perf_counter()
    sol=solve_bvp(f,bc,x,y,fun_jac=j,bc_jac=bj,tol=TOL,bc_tol=BC_TOL,max_nodes=MAX_NODES,verbose=0); wall=time.perf_counter()-t
    if sol.status: raise RuntimeError(f'eps={e:g} xmax={xmax:g}: {sol.message}; nodes={sol.x.size}; rms={np.max(sol.rms_residuals):.3e}')
    return sol,wall,float(np.max(sol.rms_residuals)),float(np.max(np.abs(bc(sol.y[:,0],sol.y[:,-1]))))
def root(sol):
    xx=np.unique(np.r_[np.linspace(0,.1,5001),np.linspace(.1,2,2000)]); vv=sol.sol(xx)[0]-1; ids=np.flatnonzero(vv[:-1]*vv[1:]<=0)
    if not ids.size: raise RuntimeError('no U=1 crossing')
    i=int(ids[0]); return float(brentq(lambda z: float(sol.sol(z)[0]-1),xx[i],xx[i+1],xtol=1e-14,rtol=1e-13))
def row(e,sol,wall,rms,bcres,xmax):
    a=asym(e); xc=root(sol); u=float(sol.y[0,0]-1); k=float(sol.y[2,0]); sov=u/(e*a['L'])**(1/3); sx=xc/(e**(1/3)*a['L']**(-1/6)); sk=(-k)/(e**(-1/3)*a['L']**(2/3)); re=lambda z,p:abs(z-p)/abs(z)
    return dict(epsilon=e,xmax=xmax,L_epsilon=a['L'],solver='scipy.integrate.solve_bvp',tol=TOL,bc_tol=BC_TOL,nodes=sol.x.size,max_rms_residual=rms,max_bc_residual=bcres,wall_seconds=wall,u0_minus1=u,x_cross=xc,u_xx0=k,scaled_overshoot=sov,scaled_cross=sx,scaled_curvature=sk,overshoot_O1=a['ao1'],overshoot_O2=a['ao2'],cross_O1=a['xo1'],cross_O2=a['xo2'],curvature_O1=a['ko1'],curvature_O2=a['ko2'],overshoot_relerr_O1=re(sov,a['ao1']),overshoot_relerr_O2=re(sov,a['ao2']),cross_relerr_O1=re(sx,a['xo1']),cross_relerr_O2=re(sx,a['xo2']),curvature_relerr_O1=re(sk,a['ko1']),curvature_relerr_O2=re(sk,a['ko2']))
def write(path,rows):
    with open(path,'w',newline='',encoding='utf-8') as f: w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
def main():
    print(f'Python {sys.version.split()[0]} | NumPy {np.__version__} | SciPy {scipy.__version__}',flush=True); primary=[]; solutions={}
    for e in EPSILONS:
        print(f'SOLVE eps={e:.0e}',flush=True); s,w,r,b=solve_one(e); rr=row(e,s,w,r,b,XMAX); primary.append(rr); solutions[e]=s; print(f" nodes={rr['nodes']} rms={r:.3e} U0-1={rr['u0_minus1']:.13e} x*={rr['x_cross']:.13e} Uxx0={rr['u_xx0']:.13e}",flush=True)
    write(RESULTS/'bvp_rerun.csv',primary)
    compact=[]
    for r in primary:
        compact.append(dict(
            epsilon=r['epsilon'], L_epsilon=r['L_epsilon'],
            overshoot_error_percent=100*r['overshoot_relerr_O2'],
            interface_error_percent=100*r['cross_relerr_O2'],
            curvature_error_percent=100*r['curvature_relerr_O2']))
    write(RESULTS/'numerical_validation_table_rerun.csv',compact)
    arch={}
    with open(ARCHIVED,newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f): arch[float(r['epsilon'])]=r
    comp=[]
    for r in primary:
        a=arch[r['epsilon']]; rd=lambda x,y:abs(x-float(y))/max(abs(float(y)),1e-300)
        comp.append(dict(epsilon=r['epsilon'],u0_minus1_rerun=r['u0_minus1'],u0_minus1_archived=a['u(0)-1'],u0_relative_difference=rd(r['u0_minus1'],a['u(0)-1']),x_cross_rerun=r['x_cross'],x_cross_archived=a['x_cross'],x_cross_relative_difference=rd(r['x_cross'],a['x_cross']),u_xx0_rerun=r['u_xx0'],u_xx0_archived=a['u_xx(0)'],u_xx0_relative_difference=rd(r['u_xx0'],a['u_xx(0)'])))
    write(RESULTS/'bvp_vs_archived.csv',comp)
    e=1e-8; ref=solutions[e]; domains=[]
    for xm in (12.,16.,20.):
        if xm==20.: rr=next(x for x in primary if x['epsilon']==e)
        elif xm<20.:
            s,w,r,b=solve_one(e,xm,guess=lambda x,S=ref:S.sol(x)); rr=row(e,s,w,r,b,xm)
        else:
            ls,lf=rates(e); y20=ref.sol(20.); cs,cf=np.linalg.solve(np.array([[1.,1.],[-ls,-lf]]),y20[:2])
            def eg(x,S=ref,cs=cs,cf=cf,ls=ls,lf=lf):
                x=np.asarray(x); Y=np.empty((4,x.size)); m=x<=20; Y[:,m]=S.sol(x[m]); z=x[~m]-20
                if z.size:
                    es=cs*np.exp(-ls*z); ef=cf*np.exp(-lf*z); Y[0,~m]=es+ef; Y[1,~m]=-ls*es-lf*ef; Y[2,~m]=ls**2*es+lf**2*ef; Y[3,~m]=-ls**3*es-lf**3*ef
                return Y
            s,w,r,b=solve_one(e,xm,guess=eg); rr=row(e,s,w,r,b,xm)
        domains.append({k:rr[k] for k in ('epsilon','xmax','nodes','max_rms_residual','max_bc_residual','u0_minus1','x_cross','u_xx0')})
    R=next(x for x in domains if x['xmax']==20.)
    for d in domains:
        d['u0_rel_to_xmax20']=abs(d['u0_minus1']-R['u0_minus1'])/abs(R['u0_minus1']); d['x_cross_rel_to_xmax20']=abs(d['x_cross']-R['x_cross'])/abs(R['x_cross']); d['u_xx0_rel_to_xmax20']=abs(d['u_xx0']-R['u_xx0'])/abs(R['u_xx0'])
    write(RESULTS/'bvp_domain_refinement.csv',domains)
    inv=np.array([1/r['L_epsilon'] for r in primary]); fig,ax=plt.subplots(figsize=(7,4.7)); ax.plot(inv,[r['scaled_overshoot'] for r in primary],'o',label='overshoot: rerun BVP'); ax.plot(inv,[r['overshoot_O2'] for r in primary],'--',label=r'overshoot: $O(L^{-2})$'); ax.plot(inv,[r['scaled_cross'] for r in primary],'s',label='interface: rerun BVP'); ax.plot(inv,[r['cross_O2'] for r in primary],'-.',label=r'interface: $O(L^{-2})$'); ax.plot(inv,[r['scaled_curvature'] for r in primary],'^',label='curvature: rerun BVP'); ax.plot(inv,[r['curvature_O2'] for r in primary],':',label=r'curvature: $O(L^{-2})$'); ax.set_xlabel(r'$1/L$'); ax.set_ylabel('scaled observable'); ax.legend(frameon=False,fontsize=8); fig.tight_layout(); fig.savefig(FIG/'fig_observables_validation.pdf'); plt.close(fig)
    env=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,solver='scipy.integrate.solve_bvp',solver_tol=TOL,solver_bc_tol=BC_TOL,solver_max_nodes=MAX_NODES,domain_default=XMAX,boundary_conditions="U'(0)=U'''(0)=0 and the two linear stable-tail conditions associated with lambda_s, lambda_f at Xmax",archived_input=ARCHIVED.name)
    (RESULTS/'bvp_environment.json').write_text(json.dumps(env,indent=2),encoding='utf-8')
    print('DONE',flush=True)
    for k in ('u0_relative_difference','x_cross_relative_difference','u_xx0_relative_difference'): print(k,max(float(x[k]) for x in comp),flush=True)
    for d in domains: print('DOMAIN',d['xmax'],d['u0_rel_to_xmax20'],d['x_cross_rel_to_xmax20'],d['u_xx0_rel_to_xmax20'],flush=True)
if __name__=='__main__': main()
