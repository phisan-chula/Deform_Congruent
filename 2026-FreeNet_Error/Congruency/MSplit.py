import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def perform_msplit_analysis(file_ep0, file_ep1, eps=1e-8, max_iter=200, tol=1e-6, t_crit=1.96):
    print("======================================================================================")
    print("                  GEODETIC DEFORMATION ANALYSIS: TWO-SPLIT M-ESTIMATION               ")
    print("======================================================================================")
    print(f"[INFO] Loading data pools from {file_ep0} and {file_ep1}...")
    
    # 1. Read and align datasets by Station using Pandas
    df0 = pd.read_csv(file_ep0)
    df1 = pd.read_csv(file_ep1)
    df_merged = pd.merge(df0, df1, on='Station', suffixes=('_ep0', '_ep1'))
    stations = df_merged['Station'].values
    
    # 2. Extract height observations and standard deviations
    y0 = df_merged['Adjusted_Height_m_ep0'].values
    y1 = df_merged['Adjusted_Height_m_ep1'].values
    sigma0 = df_merged['StdDev_mm_ep0'].values
    sigma1 = df_merged['StdDev_mm_ep1'].values
    
    # 3. Formulate the stacked observation vector (y) and weights (w)
    y = np.concatenate([y0, y1])
    w0 = 1.0 / (sigma0 ** 2)
    w1 = 1.0 / (sigma1 ** 2)
    w = np.concatenate([w0, w1])
    
    u = len(stations)  # Matrix parameter dimensions
    n = len(y)         # Total pooled observations
    print(f"[INFO] Initialization parameters mapped. Dimension of design matrix A: {n} x {u}")
    
    # 4. Construct Design Matrix A (Stacked Identity blocks)
    A = np.vstack([np.eye(u), np.eye(u)])
    
    # 5. Initialize split parameter vectors X1 and X2 with an artificial disjoint push
    X1 = y0.reshape(-1, 1) - 0.05
    X2 = y1.reshape(-1, 1) + 0.05
    
    # 6. Execute the M-split Iterative Reweighted Least Squares loop
    converged = False
    for iteration in range(1, max_iter + 1):
        v1 = A @ X1 - y.reshape(-1, 1)
        v2 = A @ X2 - y.reshape(-1, 1)
        
        # Calculate cross-weights based on opposing residuals
        p1 = w * (v2.flatten() ** 2 + eps)
        p2 = w * (v1.flatten() ** 2 + eps)
        
        X1_old = X1.copy()
        X2_old = X2.copy()
        
        # Estimate updated parameters using weighted normal equations
        X1 = np.linalg.inv(A.T @ np.diag(p1) @ A) @ A.T @ np.diag(p1) @ y.reshape(-1, 1)
        X2 = np.linalg.inv(A.T @ np.diag(p2) @ A) @ A.T @ np.diag(p2) @ y.reshape(-1, 1)
        
        # Check coordinate stability convergence
        max_delta = max(np.max(np.abs(X1 - X1_old)), np.max(np.abs(X2 - X2_old)))
        print(f"[INFO] Iteration {iteration}: Max Delta X = {max_delta:.6f} m")
        
        if max_delta < tol:
            print(f"[SUCCESS] M-split functional separation achieved in {iteration} iterations.\n")
            converged = True
            break
            
    if not converged:
        print("[WARNING] Maximum iterations reached without meeting the strict convergence tolerance.\n")
        
    # 7. Post-estimation Congruence and Local Statistical Testing
    displacement_mm = (X2.flatten() - X1.flatten()) * 1000
    sigma_diff_mm = np.sqrt(sigma0 ** 2 + sigma1 ** 2)
    t_scores = np.abs(displacement_mm) / sigma_diff_mm
    
    # Establish stability flags based on the 95% critical t-score threshold
    status = ['DISPLACED' if t > t_crit else 'STABLE' for t in t_scores]
    
    # Compile final results into a summary DataFrame
    results_df = pd.DataFrame({
        'Station': stations,
        'Model1_Height_m': X1.flatten(),
        'Model2_Height_m': X2.flatten(),
        'Displacement_mm': displacement_mm,
        'StdDev_Disp_mm': sigma_diff_mm,
        't_score': t_scores,
        'Congruence_Status': status
    })
    
    # 8. Print Clean Console Report
    print("--------------------------------------------------------------------------------------")
    print("                           FINAL ADJUSTMENT MATRIX SUMMARY                            ")
    print("--------------------------------------------------------------------------------------")
    print(f"{'Station':<7} | {'Model 1 (m)':<11} | {'Model 2 (m)':<11} | {'Displacement (mm)':<17} | {'Std.Err (mm)':<12} | {'t-score':<7} | {'Status'}")
    print("--------------------------------------------------------------------------------------")
    for index, row in results_df.iterrows():
        print(f"{row['Station']:<7} | {row['Model1_Height_m']:<11.4f} | {row['Model2_Height_m']:<11.4f} | {row['Displacement_mm']:<17.2f} | {row['StdDev_Disp_mm']:<12.2f} | {row['t_score']:<7.2f} | {row['Congruence_Status']}")
    print("--------------------------------------------------------------------------------------")
    print(f"Critical value threshold (t_95%): {t_crit}")
    
    stable_nodes = results_df[results_df['Congruence_Status'] == 'STABLE']['Station'].tolist()
    displaced_nodes = results_df[results_df['Congruence_Status'] == 'DISPLACED']['Station'].tolist()
    print(f"Global Network Conclusion: Stations {stable_nodes} form a mutually congruent reference block.")
    if displaced_nodes:
        print(f"Station(s) {displaced_nodes} exhibit statistically significant independent local displacement.")
    print("======================================================================================\n")
    
    # 9. Generate Report Visual Plot
    plt.figure(figsize=(9, 5.5))
    colors = ['#ff6b6b' if s == 'DISPLACED' else '#4dadf7' for s in status]
    
    bars = plt.bar(stations, displacement_mm, yerr=sigma_diff_mm, capsize=6, 
                   color=colors, edgecolor='black', alpha=0.85, error_kw={'ecolor': '#495057', 'lw': 1.5})
    
    plt.axhline(0, color='black', linewidth=1.2, linestyle='-')
    plt.axhline(t_crit * sigma_diff_mm.mean(), color='#e03131', linewidth=0.9, linestyle='--', label=f'95% Confidence Threshold (~{t_crit * sigma_diff_mm.mean():.1f}mm)')
    plt.axhline(-t_crit * sigma_diff_mm.mean(), color='#e03131', linewidth=0.9, linestyle='--')
    
    plt.title('Geodetic Subsidence Monitoring: Local Congruence Analysis (M-split)', fontsize=12, fontweight='bold', pad=15)
    plt.xlabel('Monitoring Station / Benchmark', fontsize=10, labelpad=10)
    plt.ylabel('Isolated Vertical Displacement (mm)', fontsize=10, labelpad=10)
    plt.grid(axis='y', linestyle=':', alpha=0.5)
    
    # Add value labels on top of the bars
    for bar in bars:
        yval = bar.get_height()
        va_dir = 'bottom' if yval >= 0 else 'top'
        offset = 5 if yval >= 0 else -15
        plt.text(bar.get_x() + bar.get_width()/2.0, yval + offset, f"{yval:+.1f} mm", ha='center', va=va_dir, fontsize=9, fontweight='bold')

    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig('displacement_plot.png', dpi=300)
    plt.close()
    print("[INFO] Summary plot compiled and exported successfully as 'displacement_plot.png'.")
    
    return results_df

# --- Main Thread Execution ---
if __name__ == "__main__":
    # Execute the processing function using your local source CSV sheets
    analysis_matrix = perform_msplit_analysis('result_ep0.csv', 'result_ep1.csv')
    
    # Output the matrix to a long-term database standard CSV format
    analysis_matrix.to_csv('congruence_analysis_results.csv', index=False)
    print("[INFO] Full data matrix written to 'congruence_analysis_results.csv'. Processing complete.")
