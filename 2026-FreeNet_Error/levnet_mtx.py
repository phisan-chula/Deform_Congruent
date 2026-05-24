import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO

# --- 1. DATA PREPARATION (CSV) ---
# Simulating the levelling network from your sketch
# Note: Point 'A' is our Fixed Benchmark (Datum)
csv_data = """from,to,dh,dist_km
A,B,2.456,0.5
B,C,1.123,0.4
C,D,-0.567,0.6
D,A,-3.010,0.3
B,D,0.550,0.8"""

df = pd.read_csv(StringIO(csv_data))
H_FIXED = {'A': 10.000}  # Known Elevation of Point A
UNKNOWN_POINTS = ['B', 'C', 'D']
point_idx = {name: i for i, name in enumerate(UNKNOWN_POINTS)}

# --- 2. LEAST SQUARES ADJUSTMENT (LSA) ---
n_obs = len(df)
n_unk = len(UNKNOWN_POINTS)

A = np.zeros((n_obs, n_unk))
L = np.zeros(n_obs)
# Stochastic Model: Weight P = 1 / Distance
P = np.diag(1.0 / df['dist_km'])

for i, row in df.iterrows():
    # Equation: H_to - H_from = dh
    # Handle 'To' point
    if row['to'] in point_idx:
        A[i, point_idx[row['to']]] = 1
    else: # It's a fixed point
        L[i] -= H_FIXED[row['to']]
        
    # Handle 'From' point
    if row['from'] in point_idx:
        A[i, point_idx[row['from']]] = -1
    else: # It's a fixed point
        L[i] += H_FIXED[row['from']]
        
    L[i] += row['dh']

# Normal Equations: (A.T @ P @ A) @ X = (A.T @ P @ L)
N = A.T @ P @ A
U = A.T @ P @ L
X = np.linalg.solve(N, U)

# Statistics
v = A @ X - L
sigma0_sq = (v.T @ P @ v) / (n_obs - n_unk)
conf_matrix = sigma0_sq * np.linalg.inv(N)
std_dev = np.sqrt(np.diag(conf_matrix))

# --- 3. SUMMARY TABLE ---
results = pd.DataFrame({
    'Point': UNKNOWN_POINTS,
    'Height_m': X,
    'StdDev_mm': std_dev * 1000
})
print("--- Adjusted Elevations ---")
print(results.to_string(index=False))

# --- 4. VISUALIZATION ---
# Assigning arbitrary coordinates for plotting the network nodes
pos = {'A': (0, 1), 'B': (1, 2), 'C': (2, 1), 'D': (1, 0)}
all_heights = {**H_FIXED, **dict(zip(UNKNOWN_POINTS, X))}

plt.figure(figsize=(8, 6))

# Plot observations (lines)
for _, row in df.iterrows():
    p1, p2 = pos[row['from']], pos[row['to']]
    plt.plot([p1[0], p2[0]], [p1[1], p2[1]], 'gray', linestyle='--', alpha=0.6)
    # Midpoint for label
    mid = ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2)
    plt.text(mid[0], mid[1], f"{row['dh']}m", color='blue', fontsize=9)

# Plot Points
for pt, (px, py) in pos.items():
    color = 'red' if pt == 'A' else 'green'
    plt.scatter(px, py, color=color, s=200, zorder=5)
    plt.text(px, py+0.1, f"{pt}\n{all_heights[pt]:.3f}m", 
             ha='center', fontweight='bold')

plt.title("Adjusted Levelling Network (Parametric Least Squares)")
plt.axis('off')
plt.show()
