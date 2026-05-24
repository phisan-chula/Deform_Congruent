import numpy as np
import pandas as pd
import lmfit
import matplotlib.pyplot as plt
from io import StringIO

# --- 1. DATA PREPARATION ---
csv_data = """from,to,dh,dist_km
A,B,2.456,0.5
B,C,1.123,0.4
C,D,-0.567,0.6
D,A,-3.010,0.3
B,D,0.550,0.8"""

df = pd.read_csv(StringIO(csv_data))

# --- 2. LMFIT SETUP ---
# Define the parameters (Elevations)
params = lmfit.Parameters()
# Point A is our Fixed Benchmark (vary=False)
params.add('A', value=10.000, vary=False) 
# Points B, C, and D are unknowns to be solved
params.add('B', value=10.000) 
params.add('C', value=10.000)
params.add('D', value=10.000)

def residual_function(params, df):
    """
    Calculates weighted residuals: sqrt(weight) * (Model - Observation)
    Note: lmfit minimizes the sum of squares of the returned array.
    """
    res = []
    for _, row in df.iterrows():
        h_from = params[row['from']].value
        h_to = params[row['to']].value
        
        # Predicted dh = H_to - H_from
        dh_model = h_to - h_from
        
        # Stochastic Model: Weight P = 1/Distance
        # Because lmfit squares the residuals, we use sqrt(P)
        weight = np.sqrt(1.0 / row['dist_km'])
        
        res.append(weight * (dh_model - row['dh']))
    return np.array(res)

# --- 3. RUN LEAST SQUARES ---
out = lmfit.minimize(residual_function, params, args=(df,))

# --- 4. REPORTING & STATISTICS ---
print("--- Fit Statistics ---")
print(lmfit.report_fit(out))

# Convert results to a clean DataFrame
final_elevations = []
for p in out.params:
    final_elevations.append({
        'Point': p,
        'Elevation': out.params[p].value,
        'StdDev_mm': f"{(out.params[p].stderr*1000):.1f}" if out.params[p].stderr else 'FIXED'
    })

results_df = pd.DataFrame(final_elevations)
print("\n--- Adjusted Results ---")
print(results_df.to_string(index=False))

# --- 5. VISUALIZATION ---
pos = {'A': (0, 1), 'B': (1, 2), 'C': (2, 1), 'D': (1, 0)}
plt.figure(figsize=(8, 6))

# Plot observations
for _, row in df.iterrows():
    p1, p2 = pos[row['from']], pos[row['to']]
    plt.plot([p1[0], p2[0]], [p1[1], p2[1]], 'gray', linestyle='--', alpha=0.5)

# Plot Points
for pt, (px, py) in pos.items():
    res_row = results_df[results_df['Point'] == pt].iloc[0]
    color = 'red' if pt == 'A' else 'green'
    plt.scatter(px, py, color=color, s=250, zorder=5)
    plt.text(px, py+0.12, f"{pt}\n{res_row['Elevation']:.3f}m\n±{res_row['StdDev_mm']}mm", 
             ha='center', fontweight='bold', fontsize=9)

plt.title("Network Adjustment (Refactored via lmfit)")
plt.axis('off')
plt.show()
