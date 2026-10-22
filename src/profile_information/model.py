"""Weighted interval ridge on a sparse ligand-by-assay grid.

The model contains molecule and assay main effects, molecule x assay interactions,
and optional pair-specific features. It never creates labels for unobserved cells.
A molecular scalar may act differently at different proteins/assays.
"""
import time
import numpy as np
from scipy.optimize import minimize


def unpack(theta, dx, dz, dq):
    at = 1
    intercept = theta[0]
    mol = theta[at:at+dx]; at += dx
    target = theta[at:at+dz]; at += dz
    interaction = theta[at:at+dx*dz].reshape(dx, dz); at += dx*dz
    pair = theta[at:].reshape(dq, dz+1)
    return intercept, mol, target, interaction, pair


def prediction(theta, x, z, q):
    intercept, mol, target, interaction, pair = unpack(theta, x.shape[1], z.shape[1], q.shape[2])
    pred = intercept + (x @ mol)[:,None] + (z @ target)[None,:] + (x @ interaction) @ z.T
    if q.shape[2]:
        pred += np.einsum('mtq,tq->mt', q, np.c_[np.ones(len(z)),z] @ pair.T, optimize=True)
    return pred


def objective(theta, x, z, q, lo, hi, w, alpha):
    pred = prediction(theta, x, z, q)
    residual = pred - np.clip(pred, lo, hi)
    residual = np.where(w > 0, residual, 0.0)
    total = w.sum()
    value = (np.sum(w * residual**2) + alpha*np.dot(theta[1:],theta[1:])) / total
    g = 2*w*residual / total
    parts = [np.array([g.sum()]), x.T @ g.sum(1), z.T @ g.sum(0), (x.T @ g @ z).ravel()]
    if q.shape[2]:
        parts.append((np.einsum('mt,mtq->tq',g,q).T @ np.c_[np.ones(len(z)),z]).ravel())
    grad = np.concatenate(parts)
    grad[1:] += 2*alpha*theta[1:] / total
    return float(value), grad


def fit(x, z, q, lower, upper, weights, alpha, options=None, initial=None):
    x,z,q = map(lambda a:np.asarray(a,dtype=float),(x,z,q))
    lo,hi,w = map(lambda a:np.asarray(a,dtype=float),(lower,upper,weights))
    shape=(len(x),len(z))
    if lo.shape != shape or hi.shape != shape or w.shape != shape or q.shape[:2] != shape:
        raise ValueError('Grid shape mismatch')
    if any(not np.isfinite(a).all() for a in (x,z,q,w)) or (w<0).any() or not (w>0).any():
        raise ValueError('Invalid features or weights')
    lo=np.where(np.isnan(lo),-np.inf,lo); hi=np.where(np.isnan(hi),np.inf,hi)
    if (lo>hi).any() or np.isposinf(lo).any() or np.isneginf(hi).any():
        raise ValueError('Impossible interval')
    if np.any((w>0)&~(np.isfinite(lo)|np.isfinite(hi))):
        raise ValueError('An unobserved cell has positive weight')
    if not np.isfinite(alpha) or alpha<=0:
        raise ValueError('Positive algorithmic ridge penalty required')
    size=1+x.shape[1]+z.shape[1]+x.shape[1]*z.shape[1]+q.shape[2]*(z.shape[1]+1)
    start=time.perf_counter()
    if initial is None:
        initial=np.zeros(size)
        finite=np.zeros_like(lo)
        both=np.isfinite(lo)&np.isfinite(hi);finite[both]=lo[both]/2+hi[both]/2
        only_lo=np.isfinite(lo)&~np.isfinite(hi);finite[only_lo]=lo[only_lo]
        only_hi=np.isfinite(hi)&~np.isfinite(lo);finite[only_hi]=hi[only_hi]
        initial[0]=np.sum(np.where(w>0,finite,0)*w)/w.sum()
    else:
        initial=np.asarray(initial,float)
        if initial.shape!=(size,):raise ValueError('Initial coefficient shape')
    opts=dict(maxiter=600,ftol=1e-10,gtol=1e-6,maxls=40)
    opts.update(options or {})
    result=minimize(objective,initial,args=(x,z,q,lo,hi,w,float(alpha)),jac=True,method='L-BFGS-B',options=opts)
    if not np.isfinite(result.x).all() or not np.isfinite(result.fun):
        raise RuntimeError('Nonfinite fit')
    diag=dict(success=bool(result.success),status=int(result.status),message=str(result.message),
              iterations=int(result.nit),objective=float(result.fun),max_gradient=float(np.max(np.abs(result.jac))),
              parameters=size,observed_cells=int((w>0).sum()),exact_training_cells=int(((w>0)&(lo==hi)).sum()),interval_training_cells=int(((w>0)&(lo!=hi)).sum()),identity_weight=float(w.sum()),alpha=float(alpha),
              wall_seconds=time.perf_counter()-start,objective_definition='(weighted squared interval distance + alpha*nonintercept ridge)/sum_weights')
    # Nonconvergence is a numerical status, never a scientific negative result.
    if not result.success:
        raise RuntimeError('Interval ridge did not converge: '+str(diag))
    return result.x,diag


def interval_loss(pred,lower,upper,weights):
    w=np.asarray(weights,float); observed=(w>0)&np.isfinite(pred)
    if not observed.any():return np.nan
    lo=np.where(np.isnan(lower),-np.inf,lower);hi=np.where(np.isnan(upper),np.inf,upper)
    err=pred-np.clip(pred,lo,hi)
    return float(np.sum(w[observed]*err[observed]**2)/w[observed].sum())
