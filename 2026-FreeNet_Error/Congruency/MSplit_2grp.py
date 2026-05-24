import numpy as np

def run_msplit_with_stats(A, l, weights, max_iter=100, tolerance=1e-6, t_crit=1.96):
    """
    Performs M-split estimation with complete A-posteriori variance 
    and localized t-test evaluation.
    """
    n, u = A.shape
    
    # Initialize Models
    x1 = np.zeros((u, 1))
    x2 = np.ones((u, 1)) * 0.01 
    
    for iteration in range(1, max_iter + 1):
        v1 = A @ x1 - l
        v2 = A @ x2 - l
        
        p1_mod = weights * (v2.flatten() ** 2)
        p2_mod = weights * (v1.flatten() ** 2)
        
        P1 = np.diag(p1_mod)
        P2 = np.diag(p2_mod)
        
        x1_old = x1.copy()
        x2_old = x2.copy()
        
        N1 = A.T @ P1 @ A
        N2 = A.T @ P2 @ A
        
        x1 = np.linalg.pinv(N1) @ A.T @ P1 @ l
        x2 = np.linalg.pinv(N2) @ A.T @ P2 @ l
        
        if np.max(np.abs(x1 - x1_old)) < tolerance and np.max(np.abs(x2 - x2_old)) < tolerance:
            break

    # ==========================================
    # STATISTICAL ANALYSIS BLOCK
    # ==========================================
    # Final residuals
    v1_final = A @ x1 - l
    v2_final = A @ x2 - l
    
    # Degrees of freedom (Observations minus parameters solved)
    # Since we use a pseudo-inverse for free network, we find effective rank
    rank_A = np.linalg.matrix_rank(A)
    df = n - rank_A 
    
    # A-posteriori variance factor (Pooled variance of the split system)
    # We use initial weights matrix 'P' to scale standard residuals
    P = np.diag(weights)
    pooled_residual_sum = (v1_final.T @ P @ v1_final + v2_final.T @ P @ v2_final) / 2.0
    sigma_0_sq = pooled_residual_sum[0, 0] / df
    
    # Covariance matrices for both models
    # Re-derived standard un-modulated structures factored by the a-posteriori variance
    N_standard = A.T @ P @ A
    Q_xx = np.linalg.pinv(N_standard) # Cofactor matrix
    C_xx = sigma_0_sq * Q_xx          # Covariance matrix
    
    # Standard deviations of parameters (square root of diagonal elements)
    std_errors = np.sqrt(np.diag(C_xx)).reshape(-1, 1)
    
    # Compute T-scores (Model 2 is typically our displacement indicator)
    t_scores_m1 = np.abs(x1) / std_errors
    t_scores_m2 = np.abs(x2) / std_errors
    
    return x1, x2, std_errors, t_scores_m1, t_scores_m2, sigma_0_sq, df

# ==========================================
# TEST RUN SETUP
# ==========================================
if __name__ == "__main__":
    # 10 Nodes (4 RBM, 6 TBM)
    num_nodes = 10
    network_schema = [
        (0, 1), (1, 2), (2, 3), (3, 0),
        (0, 4), (4, 5), (5, 1),         
        (1, 6), (6, 7), (7, 2),         
        (2, 8), (8, 9), (9, 3)          
    ]
    
    A = np.zeros((len(network_schema), num_nodes))
    for idx, (frm, to) in enumerate(network_schema):
        A[idx, frm] = -1
        A[idx, to] = 1
        
    # Real Movements (RBM_3 sank -12mm, TBM 7,8,9 sank -45mm)
    true_displacements = np.array([[0, 0, 0, -12, 0, 0, 0, -45, -45, -45]]).T / 1000.0
    l_observed = A @ true_displacements + np.random.normal(0, 0.002, size=(len(network_schema), 1))
    
    # Define Weights
    weights = np.array([1.0/(0.004**2) if f<4 and t<4 else 1.0/(0.008**2) for f, t in network_schema])
    
    # Execute Model
    t_critical_value = 1.96 # 95% Confidence limit
    x1, x2, errors, t1, t2, var_factor, dof = run_msplit_with_stats(A, l_observed, weights, t_crit=t_critical_value)
    
    # Print Statistical Summary
    print(f"--- GLOBAL STATUS ---")
    print(f"Degrees of Freedom (df): {dof}")
    print(f"A-Posteriori Variance Factor (Sigma_0^2): {var_factor:.4f}")
    print(f"Global Network Sigma (Sigma_0): {np.sqrt(var_factor)*1000:.2f} mm\n")
    
    print(f"--- LOCAL POINT ANALYSIS (MODEL 2 DISPLACEMENTS) ---")
    print(f"{'Node':<8} | {'Move (mm)':<10} | {'Std Error (mm)':<15} | {'t-score':<10} | {'Status (95% Confidence)':<22}")
    print("-" * 72)
    for i in range(num_nodes):
        name = f"RBM_{i}" if i < 4 else f"TBM_{i}"
        move = x2[i, 0] * 1000
        err = errors[i, 0] * 1000
        t_val = t2[i, 0]
        status = "❌ SIGNIFICANT SETTLEMENT" if t_val > t_critical_value else "✔ STABLE / NOISE"
        print(f"{name:<8} | {move:<10.2f} | {err:<15.2f} | {t_val:<10.2f} | {status:<22}")
