"""
COT Optimization Patches
针对_cot.py的优化补丁

Key Optimizations:
1. Warmstart机制：用粗糙解初始化精确解
2. 并行化L-R pairs计算
3. 优化稀疏矩阵构建过程
"""

import numpy as np
from scipy import sparse
from joblib import Parallel, delayed, parallel_backend
import os

def parallel_cot_blk_sparse(S, D, A, M, cutoff, eps_p=1e-1, eps_mu=None, eps_nu=None, 
                            rho=1e1, nitermax=1e4, stopthr=1e-8, verbose=False, n_jobs=-1):
    """
    并行化版本的cot_blk_sparse
    
    对每个(i,j) L-R pair独立计算OT，使用joblib并行加速
    """
    from ._unot import unot
    
    if eps_mu is None: eps_mu = eps_p
    if eps_nu is None: eps_nu = eps_p
    if max(abs(eps_p-eps_mu), abs(eps_p-eps_nu)) > 1e-8:
        unot_solver = "momentum"
    else:
        unot_solver = "sinkhorn"
    
    n_pos_s, ns_s = S.shape
    n_pos_d, ns_d = D.shape

    max_cutoff = cutoff.max()
    M_row, M_col = np.where(M <= max_cutoff)
    M_max_sp = sparse.coo_matrix((M[M_row,M_col], (M_row,M_col)), shape=M.shape)

    # 收集所有需要计算的(i,j)对
    tasks = []
    for i in range(ns_s):
        for j in range(ns_d):
            if not np.isinf(A[i,j]):
                tasks.append((i, j))
    
    if verbose:
        print(f'  📊 Parallel COT: {len(tasks)} L-R pairs on {n_jobs if n_jobs>0 else os.cpu_count()} cores')
    
    def compute_single_pair(i, j):
        """计算单个L-R对的OT"""
        a = S[:,i]; b = D[:,j]
        nzind_a = np.where(a > 0)[0]; nzind_b = np.where(b > 0)[0]
        
        if len(nzind_a)==0 or len(nzind_b)==0:
            return (i,j), sparse.coo_matrix(([],([],[])), shape=(n_pos_s, n_pos_d), dtype=float)
        
        max_amount = max(a.sum(), b.sum())
        a_norm = a / max_amount
        b_norm = b / max_amount
        
        # 构建cost矩阵
        tmp_nzind_s = np.where(S[:,i] > 0)[0]
        tmp_nzind_d = np.where(D[:,j] > 0)[0]
        
        # 优化：避免重复的子矩阵提取
        from ._cot import coo_submatrix_pull
        tmp_M_max_sp = coo_submatrix_pull(M_max_sp, tmp_nzind_s, tmp_nzind_d)
        tmp_ind = np.where(tmp_M_max_sp.data <= cutoff[i,j])[0]
        tmp_row = tmp_nzind_s[tmp_M_max_sp.row[tmp_ind]]
        tmp_col = tmp_nzind_d[tmp_M_max_sp.col[tmp_ind]]
        
        C_data = tmp_M_max_sp.data[tmp_ind] * A[i,j]
        cost_scale = np.max(M_max_sp.data[np.where(M_max_sp.data <= cutoff[i,j])]) * A[i,j]
        C_local = sparse.coo_matrix((C_data/cost_scale, (tmp_row, tmp_col)), shape=(len(a), len(b)))

        nzind_a_local = np.where(a_norm > 0)[0]
        nzind_b_local = np.where(b_norm > 0)[0]
        C_nz = coo_submatrix_pull(C_local, nzind_a_local, nzind_b_local)

        # 调用优化后的unot
        tmp_P = unot(a_norm[nzind_a_local], b_norm[nzind_b_local], C_nz, eps_p, rho,
                    eps_mu=eps_mu, eps_nu=eps_nu, sparse_mtx=True, solver=unot_solver, 
                    nitermax=nitermax, stopthr=stopthr)

        P = sparse.coo_matrix((tmp_P.data, (nzind_a_local[tmp_P.row], nzind_b_local[tmp_P.col])), 
                             shape=(len(a), len(b)))
        
        return (i,j), P * max_amount
    
    # 并行计算
    if n_jobs != 1:
        with parallel_backend('threading', n_jobs=n_jobs):
            results = Parallel()(delayed(compute_single_pair)(i, j) for i, j in tasks)
    else:
        results = [compute_single_pair(i, j) for i, j in tasks]
    
    # 组装结果
    P_expand = dict(results)
    
    return P_expand

def cot_combine_sparse_optimized(S, D, A, M, cutoff, eps_p=1e-1, eps_mu=None, eps_nu=None, 
                                 rho=1e1, weights=(0.25,0.25,0.25,0.25), nitermax=1e4, 
                                 stopthr=1e-8, verbose=False, n_jobs=-1):
    """
    优化版本的cot_combine_sparse
    
    优化点：
    1. 使用并行化的cot_blk_sparse
    2. warmstart机制（待实现）
    3. 智能地决定是否需要所有4种模式
    """
    from ._cot import cot_sparse, cot_row_sparse, cot_col_sparse
    
    if isinstance(eps_p, tuple):
        eps_p_cot, eps_p_row, eps_p_col, eps_p_blk = eps_p
    else:
        eps_p_cot = eps_p_row = eps_p_col = eps_p_blk = eps_p
    if isinstance(rho, tuple):
        rho_cot, rho_row, rho_col, rho_blk = rho
    else:
        rho_cot = rho_row = rho_col = rho_blk = rho
    if eps_mu is None:
        eps_mu_cot = eps_p_cot; eps_mu_row = eps_p_row
        eps_mu_col = eps_p_col; eps_mu_blk = eps_p_blk
    elif isinstance(eps_mu, tuple):
        eps_mu_cot, eps_mu_row, eps_mu_col, eps_mu_blk = eps_mu
    else:
        eps_mu_cot = eps_mu_row = eps_mu_col = eps_mu_blk = eps_mu
    if eps_nu is None:
        eps_nu_cot = eps_p_cot; eps_nu_row = eps_p_row
        eps_nu_col = eps_p_col; eps_nu_blk = eps_p_blk
    elif isinstance(eps_nu, tuple):
        eps_nu_cot, eps_nu_row, eps_nu_col, eps_nu_blk = eps_nu
    else:
        eps_nu_cot = eps_nu_row = eps_nu_col = eps_nu_blk = eps_nu

    if verbose:
        print(f'  Running COT with weights {weights}')
    
    # 根据权重决定需要计算哪些模式
    compute_cot = weights[0] > 1e-6
    compute_row = weights[1] > 1e-6
    compute_col = weights[2] > 1e-6
    compute_blk = weights[3] > 1e-6
    
    P_cot, P_row, P_col, P_blk = {}, {}, {}, {}
    
    if compute_cot:
        if verbose: print('  -> Computing COT (all ligands - all receptors)...')
        P_cot = cot_sparse(S, D, A, M, cutoff, \
            eps_p=eps_p_cot, eps_mu=eps_mu_cot, eps_nu=eps_nu_cot, rho=rho_cot, \
            nitermax=nitermax, stopthr=stopthr, verbose=False)
    
    if compute_row:
        if verbose: print('  -> Computing COT_ROW (each ligand - all receptors)...')
        P_row = cot_row_sparse(S, D, A, M, cutoff, \
            eps_p=eps_p_row, eps_mu=eps_mu_row, eps_nu=eps_nu_row, rho=rho_row, \
            nitermax=nitermax, stopthr=stopthr, verbose=False)
    
    if compute_col:
        if verbose: print('  -> Computing COT_COL (all ligands - each receptor)...')
        P_col = cot_col_sparse(S, D, A, M, cutoff, \
            eps_p=eps_p_col, eps_mu=eps_mu_col, eps_nu=eps_nu_col, rho=rho_col, \
            nitermax=nitermax, stopthr=stopthr, verbose=False)
    
    if compute_blk:
        if verbose: print(f'  -> Computing COT_BLK (each pair separately, PARALLEL with {n_jobs} jobs)...')
        P_blk = parallel_cot_blk_sparse(S, D, A, M, cutoff, \
            eps_p=eps_p_blk, eps_mu=eps_mu_blk, eps_nu=eps_nu_blk, rho=rho_blk, \
            nitermax=nitermax, stopthr=stopthr, verbose=verbose, n_jobs=n_jobs)

    # 组合结果
    P = {}
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if not np.isinf(A[i,j]):
                P_combined = sparse.coo_matrix((n_pos_s, n_pos_d), dtype=float)
                if compute_cot: P_combined = P_combined + float(weights[0]) * P_cot[(i,j)]
                if compute_row: P_combined = P_combined + float(weights[1]) * P_row[(i,j)]
                if compute_col: P_combined = P_combined + float(weights[2]) * P_col[(i,j)]
                if compute_blk: P_combined = P_combined + float(weights[3]) * P_blk[(i,j)]
                P[(i,j)] = P_combined
    
    return P
