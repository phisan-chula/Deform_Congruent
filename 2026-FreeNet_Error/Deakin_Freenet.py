import numpy as np
import pandas as pd
import lmfit
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from io import StringIO

# =====================================================================
# --- 1. DATA PREPARATION & SPATIAL GEOMETRY ---
# =====================================================================
csv_data = """from,to,dh_m,dist_km
A,X,6.345,1.6
B,X,4.235,2.5
Z,B,3.060,1.0
Z,A,0.920,4.0
A,Y,3.895,1.6
Y,X,2.410,1.25
Z,Y,4.820,2.0"""

df = pd.read_csv(StringIO(csv_data))

# Fully dynamic structural identification of network stations
all_points = sorted(list(set(df['from']).union(set(df['to']))))
point_idx = {name: i for i, name in enumerate(all_points)}

# Planimetric metric coordinates anchored from reference node Y = (5000, 15000)
pos = {
    'Y': (5000.0, 15000.0),
    'A': (3400.0, 15000.0),
    'X': (4400.0, 16000.0),
    'B': (6900.0, 16000.0),
    'Z': (6000.0, 13400.0)
}

# =====================================================================
# --- 2. CORE OPTIMIZATION ENGINES ---
# =====================================================================

def run_pure_freenet_numpy():
    """Free Network Adjustment using Moore-Penrose Pseudo-Inverse."""
    n_obs = len(df)
    n_unk = len(all_points)
    
    A = np.zeros((n_obs, n_unk))
    L = np.zeros(n_obs)
    P = np.diag(1.0 / (df['dist_km'] * 1000.0)) 
    
    for i, row in df.iterrows():
        A[i, point_idx[row['to']]] = 1
        A[i, point_idx[row['from']]] = -1
        L[i] = row['dh_m']
        
    N = A.T @ P @ A
    U = A.T @ P @ L
    
    X = np.linalg.pinv(N) @ U
    
    v = A @ X - L
    dof = n_obs - (n_unk - 1)
    sigma0_sq = (v.T @ P @ v) / dof
    conf_matrix = sigma0_sq * np.linalg.pinv(N)
    std_errors = np.sqrt(np.diag(conf_matrix))
    
    elevations = {pt: X[point_idx[pt]] for pt in all_points}
    std_devs = {pt: std_errors[point_idx[pt]] for pt in all_points}
    trace_Qxx = np.trace(np.linalg.pinv(N)) 
    return elevations, std_devs, trace_Qxx


def run_pure_freenet_lmfit():
    """Free Network Adjustment via lmfit with Option 2 Inner Constraint Regularization."""
    params = lmfit.Parameters()
    for pt in all_points:
        params.add(pt, value=0.000, vary=True)
        
    def residual_function(p, dataframe):
        res = []
        # 1. Physical Observation Residuals
        for _, row in dataframe.iterrows():
            h_from = p[row['from']].value
            h_to = p[row['to']].value
            
            v_i = (h_to - h_from) - row['dh_m']
            sigma_i = np.sqrt(row['dist_km'] * 1000.0) 
            res.append(v_i / sigma_i)
            
        # 2. Inner Constraint (Option 2: Dynamic Weighting)
        param_values = np.fromiter((param.value for param in p.values()), dtype=float)
        centroid_residual = np.sum(param_values - 0.0)  # Center around 0.0 datum
        
        penalty_weight = np.sqrt(len(p)) * 1.0e4
        res.append(penalty_weight * centroid_residual)
        
        return np.array(res)

    out = lmfit.minimize(residual_function, params, args=(df,), method='leastsq')
    print("\n=== LMFIT NATIVE REPORT (CORRECTED STANDARD ERRORS) ===")
    print(lmfit.fit_report(out))
    
    elevations = {}
    std_devs = {}
    for pt in all_points:
        elevations[pt] = out.params[pt].value
        std_devs[pt] = out.params[pt].stderr if out.params[pt].stderr is not None else 0.0
        
    # --- METHOD 2: DIRECT COVARIANCE MATRIX TRACE EXTRACTOR ---
    # Removes legacy manual layout recreation loops and pinv() math completely!
    if out.covar is not None:
        trace_Qxx = np.trace(out.covar) / out.redchi
    else:
        trace_Qxx = 0.0
        
    return elevations, std_devs, trace_Qxx


def run_constrained_matrix(fixed_pt, fixed_val=0.000):
    """Rigidly constrained adjustments anchoring a given node to a datum target."""
    unknown_points = [p for p in all_points if p != fixed_pt]
    sub_idx = {name: i for i, name in enumerate(unknown_points)}
    
    A = np.zeros((len(df), len(unknown_points)))
    L = np.zeros(len(df))
    P = np.diag(1.0 / (df['dist_km'] * 1000.0))
    
    for i, row in df.iterrows():
        if row['to'] in sub_idx:
            A[i, sub_idx[row['to']]] = 1
        else:
            L[i] -= fixed_val
            
        if row['from'] in sub_idx:
            A[i, sub_idx[row['from']]] = -1
        else:
            L[i] += fixed_val
            
        L[i] += row['dh_m']
        
    N = A.T @ P @ A
    X = np.linalg.solve(N, A.T @ P @ L)
    
    v = A @ X - L
    sigma0_sq = (v.T @ P @ v) / (len(df) - len(unknown_points))
    conf_matrix = sigma0_sq * np.linalg.inv(N)
    std_errors = np.sqrt(np.diag(conf_matrix))
    
    elevations = {fixed_pt: fixed_val}
    std_devs = {fixed_pt: 0.0}
    for idx, pt in enumerate(unknown_points):
        elevations[pt] = X[idx]
        std_devs[pt] = std_errors[idx]
        
    trace_Qxx = np.trace(np.linalg.inv(N))
    return elevations, std_devs, trace_Qxx

# =====================================================================
# --- 3. VISUALIZATION SUB-ENGINE ---
# =====================================================================
VERTICAL_SCALE = 55000.0  

def plot_network_sub(ax, title, elevations, std_devs, fixed_pt=None):
    """Generates structural subplots showing layout lines and lean error shapes."""
    x_vals = np.array([elevations[pt] for pt in all_points])
    sum_x = np.sum(x_vals)
    l2_norm = np.linalg.norm(x_vals)
    
    # Draw Network Links
    for _, row in df.iterrows():
        p1, p2 = pos[row['from']], pos[row['to']]
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 'dimgray', linestyle='--', alpha=0.4, zorder=1)
        
    # Draw Text Labels First
    for pt, (px, py) in pos.items():
        is_fixed = (pt == fixed_pt)
        lbl = f"({pt})\n{elevations[pt]:.3f}m\n" + ("FIXED" if is_fixed else f"±{std_devs[pt]*1000:.1f}mm")
        ax.text(px, py + 180, lbl, ha='center', va='center', fontsize=8, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="gainsboro", alpha=0.9, zorder=2))
        
    # Draw Lean Vertical Error Ellipses (Crimson Red + Translucent Layer)
    for pt, (px, py) in pos.items():
        if std_devs[pt] > 0:
            v_height = std_devs[pt] * VERTICAL_SCALE
            h_width = 60.0  
            
            ellipse = Ellipse(xy=(px, py), width=h_width, height=v_height, angle=0,
                              edgecolor='crimson', facecolor='crimson', alpha=0.20, lw=1.2, zorder=5)
            ellipse_outline = Ellipse(xy=(px, py), width=h_width, height=v_height, angle=0,
                                      edgecolor='crimson', facecolor='none', lw=1.0, zorder=6)
            ax.add_patch(ellipse)
            ax.add_patch(ellipse_outline)

    # Draw Node Markers explicitly on the TOPMOST layer (zorder=10)
    for pt, (px, py) in pos.items():
        if pt == fixed_pt:
            ax.scatter(px, py, color='crimson', marker='^', s=160, zorder=10)
        else:
            ax.scatter(px, py, color='forestgreen', marker='o', s=70, zorder=10)
            
    # Render Scale Bar as a Standalone Lean Vertical Reference Error Ellipse
    scale_bar_val_mm = 10.0  
    ref_v_height = (scale_bar_val_mm / 1000.0) * VERTICAL_SCALE
    ref_h_width = 60.0  
    sb_x, sb_y = 7450.0, 16650.0 
    
    ref_ellipse = Ellipse(xy=(sb_x, sb_y), width=ref_h_width, height=ref_v_height, angle=0,
                          edgecolor='crimson', facecolor='crimson', alpha=0.20, lw=1.2, zorder=5, clip_on=False)
    ref_outline = Ellipse(xy=(sb_x, sb_y), width=ref_h_width, height=ref_v_height, angle=0,
                          edgecolor='crimson', facecolor='none', lw=1.0, zorder=6, clip_on=False)
    ax.add_patch(ref_ellipse)
    ax.add_patch(ref_outline)
    
    ax.text(sb_x - 80, sb_y, f"{scale_bar_val_mm} mm\n(Vert. Scale)",
            ha='right', va='center', fontsize=7, fontweight='bold', color='crimson', clip_on=False)
    
    display_sum = 0.000 if abs(sum_x) < 1e-5 else sum_x
    full_title = f"{title}\n$\sum$X = {display_sum:.3f}m | ||X|| = {l2_norm:.3f}m"
    ax.set_title(full_title, fontsize=10, fontweight='bold', color='navy')
    
    ax.set_xlim(2800, 7800)
    ax.set_ylim(12800, 16800)
    ax.axis('off')

# =====================================================================
# --- 4. PIPELINE EXECUTION WRAPPERS ---
# =====================================================================

def generate_comparison_summary():
    """Computes adjustments across all scenarios and prints an absolute Markdown Table."""
    summary_rows = []
    
    # 1. Pure Freenet Scenario via LMFIT Engine (Using clean direct trace calculation)
    elev_lm, _, trace_lm = run_pure_freenet_lmfit()
    x_lm = np.array(list(elev_lm.values()))
    summary_rows.append({
        'Points': len(all_points),
        'Fixed Point': 'none (LMFIT Refactored)',
        'Sum(X)': np.sum(x_lm),
        'L2-Norm': np.linalg.norm(x_lm),
        'Trace(Qxx)': trace_lm
    })
    
    # 2. Pure Freenet Scenario via Moore-Penrose Pseudo-Inverse (NumPy Engine)
    elev_np, _, trace_np = run_pure_freenet_numpy()
    x_np = np.array(list(elev_np.values()))
    summary_rows.append({
        'Points': len(all_points),
        'Fixed Point': 'none (MoorPseudoPenrose)',
        'Sum(X)': np.sum(x_np),
        'L2-Norm': np.linalg.norm(x_np),
        'Trace(Qxx)': trace_np
    })
    
    # 3-7. Rotational Constraints Across All Network Stations sequentially
    for pt in all_points:
        elev_c, _, trace_c = run_constrained_matrix(pt, fixed_val=0.000)
        x_c = np.array([elev_c[p] for p in all_points])
        summary_rows.append({
            'Points': len(all_points) - 1, 
            'Fixed Point': f"'{pt}'",
            'Sum(X)': np.sum(x_c),
            'L2-Norm': np.linalg.norm(x_c),
            'Trace(Qxx)': trace_c
        })
        
    print("\nGEODETIC NETWORK ADJUSTMENT DIAGNOSTIC SUMMARY (MARKDOWN FORMAT):\n")
    print("| Points | Fixed Point | Sum(X) (m) | L2-Norm (m) | Trace(Qxx) |")
    print("| :---:  | :---        | :---:      | :---:       | :---:      |")
    for row in summary_rows:
        clean_sum = 0.000 if abs(row['Sum(X)']) < 1e-5 else row['Sum(X)']
        print(f"| {row['Points']} | {row['Fixed Point']:25} | {clean_sum:10.3f} | {row['L2-Norm']:11.3f} | {row['Trace(Qxx)']:10.5f} |")
    print("\n")


def plot_all_networks():
    """Generates the multi-graph subplot arrangement and saves outputs to storage."""
    fig, axes = plt.subplots(3, 2, figsize=(14, 16))
    axes = axes.flatten()
    
    # Subplot 1: Pure Freenet Mode (LMFIT Engine)
    lm_elev_free, lm_std_free, _ = run_pure_freenet_lmfit()
    plot_network_sub(axes[0], "1. Pure Freenet Mode (LMFIT Engine)", lm_elev_free, lm_std_free, fixed_pt=None)
    
    # Subplots 2 to 6: Rotational Constraint Shifts for all stations (A, B, X, Y, Z) fixed to 0.0m
    for idx, target_fix in enumerate(all_points):
        ax_target = axes[idx + 1] 
        elevs, stds, _ = run_constrained_matrix(target_fix, fixed_val=0.000)
        plot_network_sub(ax_target, f"{idx+2}. Constrained Net: Fixed [{target_fix}=0.0m]", elevs, stds, fixed_pt=target_fix)
        
    plt.suptitle("Comparative Level Network Adjustments", fontsize=14, fontweight='bold', y=0.965)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    
    PLT_SVG = 'FreeNet_Deakin_All.svg'
    PNG_IMG = 'FreeNet_Deakin_All.png'
    for fmt in (PLT_SVG, PNG_IMG):
        print(f'Plotting {fmt} ...')
        plt.savefig(fmt)
    plt.show()

# =====================================================================
# --- 5. CONTROL LOOP CONTROLLER ---
# =====================================================================
if __name__ == "__main__":
    generate_comparison_summary()
    plot_all_networks()
