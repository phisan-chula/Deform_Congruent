import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

try:
    import tomllib  # Native in Python 3.11+
except ImportError:
    import tomli as tomllib  # Fallback backport for Python 3.7 to 3.10

class LevelingNetwork:
    """
    An integrated geodetic engine to parse compact TOML networks,
    track datum topology, and filter effective vs. inactive benchmarks.
    """
    def __init__(self, file_path: str):
        self.file_path = file_path
        self.dfFixed = pd.DataFrame()
        self.dfPos = pd.DataFrame()
        self.dfDiff = pd.DataFrame()
        
        self.load_network_data()

    def load_network_data(self):
        """Ingests all TOML configurations and maps an 'Effective' validation flag."""
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"Target configuration file '{self.file_path}' not found.")
            
        with open(self.file_path, "rb") as f:
            raw_data = tomllib.load(f)
        
        # 1. Parse Differential Leveling Routes First
        diff_lines = raw_data.get("Differential", {}).get("lines", [])
        diff_headers = ["LineNo", "StaFrom", "StaTo", "Diff_m", "Dist_km"]
        self.dfDiff = pd.DataFrame(diff_lines, columns=diff_headers)

        # 2. Gather all active stations from observations
        observed_stations = set(self.dfDiff['StaFrom'].tolist() + self.dfDiff['StaTo'].tolist())

        # 3. Parse Fixed Reference Benchmarks (Keep ALL records)
        fixed_dict = raw_data.get("Fixed", {})
        self.dfFixed = pd.DataFrame(list(fixed_dict.items()), columns=["Station", "Height_m"])
        
        # 4. Dynamically compute the 'Effective' column mapping
        self.dfFixed['Effective'] = self.dfFixed['Station'].isin(observed_stations)

        # 5. Parse Spatial Station Positions
        pos_dict = raw_data.get("Position", {})
        self.dfPos = pd.DataFrame(
            [[k, v[0], v[1]] for k, v in pos_dict.items()], 
            columns=["Station", "X", "Y"]
        )

    def display(self):
        """Audits network structure logic and prints clean Markdown tables."""
        print("### NETWORK INTEGRITY & DATUM AUDIT\n")
        
        # Extract only the benchmarks that are marked as effective
        effective_fixed = set(self.dfFixed[self.dfFixed['Effective']]['Station'].tolist())

        # 1. Display All Fixed Benchmarks with their respective operational status
        df_fixed_display = self.dfFixed.copy()
        df_fixed_display['Effective'] = df_fixed_display['Effective'].apply(
            lambda eff: "**Yes (Active)**" if eff else "No (No use!)"
        )
        print("#### Fixed Benchmarks (`dfFixed`)")
        print(df_fixed_display.to_markdown(index=False), "\n")
        
        # 2. Display Spatial Inventory
        print("#### Station Positions (`dfPos`)")
        print(self.dfPos.to_markdown(index=False), "\n")
        
        # 3. Audit dfDiff and inject anchor markers ONLY for effective benchmarks
        df_diff_display = self.dfDiff.copy()
        df_diff_display['StaFrom'] = df_diff_display['StaFrom'].apply(
            lambda st: f"**{st}***" if st in effective_fixed else st
        )
        df_diff_display['StaTo'] = df_diff_display['StaTo'].apply(
            lambda st: f"**{st}***" if st in effective_fixed else st
        )
        
        print("#### Differential Leveling Lines (`dfDiff`)")
        print("> `*` indicates an Effective Fixed Anchor Base Node\n")
        print(df_diff_display.to_markdown(index=False))

    def plot_network(self):
        """Draws an isometric geodetic schematic map of your observation scheme."""
        if self.dfPos.empty or self.dfDiff.empty:
            print("[Error] Missing spatial arrays. Aborting plotter.")
            return

        coord_map = {row['Station']: (row['X'], row['Y']) for _, row in self.dfPos.iterrows()}
        
        # Visual markers should look at the 'Effective' column status
        effective_fixed = set(self.dfFixed[self.dfFixed['Effective']]['Station'].tolist())
        all_fixed = set(self.dfFixed['Station'].tolist())

        fig, ax = plt.subplots(figsize=(11, 9))

        # --- Draw Observation Arrows and Vector Midpoint Info Boxes ---
        for _, row in self.dfDiff.iterrows():
            frm, to = row['StaFrom'], row['StaTo']
            if frm not in coord_map or to not in coord_map:
                continue

            x1, y1 = coord_map[frm]
            x2, y2 = coord_map[to]

            ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                        arrowprops=dict(arrowstyle="-|>", color='#555555', 
                                       lw=1.6, mutation_scale=16, ls='-'), zorder=1)

            mx, my = (x1 + x2) / 2, (y1 + y2) / 2
            label_text = f"dh: {row['Diff_m']:+.3f} m\nL: {row['Dist_km']:.2f} km"
            
            ax.text(mx, my, label_text, color='black', fontsize=8, ha='center', va='center',
                    fontweight='medium', bbox=dict(boxstyle='round,pad=0.3', 
                    facecolor='#ffffff', edgecolor='#d3d3d3', alpha=0.9, zorder=2))

        # --- Draw Network Station Points ---
        for _, row in self.dfPos.iterrows():
            st, x, y = row['Station'], row['X'], row['Y']

            if st in effective_fixed:
                color, marker, size = '#dc3545', '^', 140  # Active Fixed Benchmark
            elif st in all_fixed:
                color, marker, size = '#6c757d', 's', 90   # Ineffective/Orphan Benchmark (Gray Square)
            else:
                color, marker, size = '#007bff', 'o', 90   # Normal Target Node

            ax.scatter(x, y, color=color, marker=marker, s=size, zorder=3, edgecolors='black', lw=1)
            ax.text(x + 120, y + 120, st, fontsize=11, fontweight='bold', color='#222222', zorder=4)

        # --- Framework Canvas Adjustments ---
        ax.set_title("Geodetic Leveling Network Scheme Map", fontsize=13, fontweight='bold', pad=15)
        ax.set_xlabel("Planimetric Coordinate X (m)", labelpad=8)
        ax.set_ylabel("Planimetric Coordinate Y (m)", labelpad=8)
        ax.grid(True, linestyle=':', color='#cccccc', alpha=0.7)

        legend_handles = [
            Line2D([0], [0], marker='^', color='w', label='Active Fixed Benchmark (BM)', markerfacecolor='#dc3545', markersize=10, markeredgecolor='black'),
            Line2D([0], [0], marker='s', color='w', label='Unused Fixed Benchmark (Orphan)', markerfacecolor='#6c757d', markersize=8, markeredgecolor='black'),
            Line2D([0], [0], marker='o', color='w', label='Leveling Check Point / Node', markerfacecolor='#007bff', markersize=8, markeredgecolor='black'),
            Line2D([0], [0], color='#555555', lw=1.5, label='Leveling Direction (From -> To)')
        ]
        ax.legend(handles=legend_handles, loc='upper left', frameon=True, facecolor='white', framealpha=0.9)
        ax.set_aspect('equal', adjustable='box')
        
        plt.tight_layout()
        plt.show()

def main():
    config_file = "network.toml"
    try:
        network = LevelingNetwork(config_file)
        network.display()
        print("\n[System] Rendering network vector map...")
        network.plot_network()
    except Exception as e:
        print(f"CRITICAL: Engine initialization failed: {e}")

if __name__ == "__main__":
    main()
