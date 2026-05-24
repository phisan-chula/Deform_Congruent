import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2

def create_datum_matrix(coords):
    """Constructs the 2D rigid-body datum matrix G (dx, dy, r) centered on the centroid."""
    n = len(coords)
    G = np.zeros((2 * n, 3))
    centroid = np.mean(coords, axis=0)
    centered_coords = coords - centroid
    
    for i in range(n):
        x, y = centered_coords[i]
        G[2 * i, 0] = 1.0      # dx
        G[2 * i + 1, 1] = 1.0  # dy
        G[2 * i, 2] = -y       # rotation effect on X
        G[2 * i + 1, 2] = x    # rotation effect on Y
    return G

def s_transform(x, Q, G, E):
    """Applies S-Transformation using selection matrix E."""
    GTG_inv = np.linalg.inv(G.T @ E @ G)
    I = np.eye(len(x))
    S = I - G @ GTG_inv @ G.T @ E
    return S @ x, S @ Q @ S.T

# ==========================================
# 1. SETUP SIMULATED 10-POINT NETWORK
# ==========================================
np.random.seed(42)
n_points = 10
nominal_coords = np.random.uniform(100, 300, (n_points, 2))

# Define true displacements (Points 3 and 7 move significantly)
true_displacements = np.zeros((n_points, 2))
true_displacements[3] = [0.045, -0.030]  # ~5.4 cm
true_displacements[7] = [-0.050, 0.040]  # ~6.4 cm

x1_obs = nominal_coords.flatten() + np.random.normal(0, 0.002, 2 * n_points)
x2_obs = (nominal_coords + true_displacements).flatten() + np.random.normal(0, 0.002, 2 * n_points)

# Base inner-constraint cofactor matrices
G_initial = create_datum_matrix(nominal_coords)
S_initial = np.eye(2 * n_points) - G_initial @ np.linalg.inv(G_initial.T @ G_initial) @ G_initial.T
Q1 = S_initial @ (np.eye(2 * n_points) * (0.002**2)) @ S_initial.T
Q2 = S_initial @ (np.eye(2 * n_points) * (0.002**2)) @ S_initial.T

# ==========================================
# 2. ITERATIVE PROCESS & PLOTTING LOOP
# ==========================================
alpha = 0.05
critical_value = chi2.ppf(1 - alpha, df=2)
active_datum_points = list(range(n_points))
G = create_datum_matrix(nominal_coords)

iteration = 1
while True:
    # Build E matrix
    E = np.zeros((2 * n_points, 2 * n_points))
    for pt in active_datum_points:
        E[2 * pt, 2 * pt] = 1.0
        E[2 * pt + 1, 2 * pt + 1] = 1.0
        
    # Transform
    x1_trans, _ = s_transform(x1_obs, Q1, G, E)
    x2_trans, Q_trans = s_transform(x2_obs, Q2, G, E)
    d = x2_trans - x1_trans
    Q_dd = Q1 + Q_trans # Combined epoch variances
    
    # Evaluate individual T-statistics
    t_stats = {}
    max_T = -1.0
    worst_point = None
    
    for pt in range(n_points):
        idx = [2 * pt, 2 * pt + 1]
        d_pt = d[idx]
        Q_dd_pt = Q_dd[np.ix_(idx, idx)]
        T = d_pt.T @ np.linalg.inv(Q_dd_pt) @ d_pt
        t_stats[pt] = T
        if pt in active_datum_points and T > max_T:
            max_T = T
            worst_point = pt

    # --- PLOT GENERATION ---
    plt.figure(figsize=(10, 8))
    
    # Plot nominal layout
    plt.scatter(nominal_coords[:, 0], nominal_coords[:, 1], c='lightgray', zorder=1, label='Base Network Structure')
    
    # Scale displacements for visualization visibility (1 meter on plot = 1 millimeter real motion)
    scale_factor = 1000  
    
    for pt in range(n_points):
        dx = d[2 * pt]
        dy = d[2 * pt + 1]
        T_val = t_stats[pt]
        
        # Color nodes based on current status
        if pt == worst_point and max_T > critical_value:
            color = 'crimson'
            marker = 'X'
            size = 130
        elif pt in active_datum_points:
            color = 'forestgreen'
            marker = 'o'
            size = 80
        else:
            color = 'purple'  # Already excluded from datum
            marker = 's'
            size = 60
            
        plt.scatter(nominal_coords[pt, 0], nominal_coords[pt, 1], color=color, marker=marker, s=size, zorder=3)
        
        # Draw exaggerated displacement vectors
        plt.arrow(nominal_coords[pt, 0], nominal_coords[pt, 1], dx * scale_factor, dy * scale_factor,
                  head_width=3, head_length=4, fc=color, ec=color, alpha=0.7, zorder=2)
        
        # Labels
        plt.text(nominal_coords[pt, 0] + 4, nominal_coords[pt, 1] + 4, 
                 f"P{pt}\nT={T_val:.1f}", fontsize=9, fontweight='bold',
                 bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=1))

    # Formatting
    plt.title(f"Iteration {iteration}: S-Transform Datum Analysis\n"
              f"Active Datum Frame: {active_datum_points} | Critical Threshold = {critical_value:.2f}", 
              fontsize=12, pad=15)
    plt.xlabel("Easting (m)")
    plt.ylabel("Northing (m)")
    plt.grid(True, linestyle='--', alpha=0.5)
    
    # Dynamic Legend handles
    from matplotlib.lines import Line2D
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='Active / Stable Reference Node', markerfacecolor='forestgreen', markersize=10),
        plt.Line2D([0], [0], marker='X', color='w', label='Flagged Outlier (To Exclude)', markerfacecolor='crimson', markersize=12),
        plt.Line2D([0], [0], marker='s', color='w', label='Excluded From Datum', markerfacecolor='purple', markersize=9),
        plt.Line2D([0], [0], color='gray', lw=1.5, linestyle='-', label='Displacement Vector (Magnified 1000x)')
    ]
    plt.legend(handles=legend_elements, loc='upper right')
    
    # Save Image
    filename = f"iteration_{iteration}.png"
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {filename} (Worst point: P{worst_point}, T = {max_T:.2f})")
    
    # Termination Logic
    if max_T > critical_value:
        active_datum_points.remove(worst_point)
        iteration += 1
    else:
        break

print("\nProcessing complete. All plots saved to your directory.")
