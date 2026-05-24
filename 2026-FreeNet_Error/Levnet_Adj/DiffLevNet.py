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
        self.crs = "Unspecified"
        
        self.load_network_data()

    def load_network_data(self):
        """Ingests all TOML configurations and maps an 'Effective' validation flag."""
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"Target configuration file '{self.file_path}' not found.")
            
        with open(self.file_path, "rb") as f:
            raw_data = tomllib.load(f)
        
        # 1. Parse Project Metadata
        meta_sec = raw_data.get("Metadata", {})
        self.crs = meta_sec.get("crs", "Unspecified")

        # 2. Parse Differential Leveling Routes First
        diff_lines = raw_data.get("Differential", {}).get("lines", [])
        diff_headers = ["LineNo", "StaFrom", "StaTo", "Diff_m", "Dist_km"]
        self.dfDiff = pd.DataFrame(diff_lines, columns=diff_headers)

        # 3. Gather all active stations from observations
        observed_stations = set(self.dfDiff['StaFrom'].tolist() + self.dfDiff['StaTo'].tolist())

        # 4. Parse Fixed Reference Benchmarks (Keep ALL records)
        fixed_dict = raw_data.get("Fixed", {})
        self.dfFixed = pd.DataFrame(list(fixed_dict.items()), columns=["Station", "Height_m"])
        
        # 5. Dynamically compute the 'Effective' column mapping
        self.dfFixed['Effective'] = self.dfFixed['Station'].isin(observed_stations)

        # 6. Parse Spatial Station Positions
        pos_dict = raw_data.get("Position", {})
        self.dfPos = pd.DataFrame(
            [[k, v[0], v[1]] for k, v in pos_dict.items()], 
            columns=["Station", "X", "Y"]
        )

    def display(self):
        """Audits network structure logic and prints clean Markdown tables without scientific notation."""
        print("### NETWORK INTEGRITY & DATUM AUDIT\n")
        
        # Cross-reference observed lines to find orphaned benchmarks
        effective_fixed = set(self.dfFixed[self.dfFixed['Effective']]['Station'].tolist())

        # 1. Display All Fixed Benchmarks
        df_fixed_display = self.dfFixed.copy()
        df_fixed_display['Effective'] = df_fixed_display['Effective'].apply(
            lambda eff: "**Yes (Active)**" if eff else "No (No use!)"
        )
        print("#### Fixed Benchmarks (`dfFixed`)")
        print(df_fixed_display.to_markdown(index=False), "\n")
        
        # 2. Display Spatial Inventory (Saves coordinates cleanly as pure text integers)
        df_pos_display = self.dfPos.copy()
        if 'X' in df_pos_display.columns:
            df_pos_display['X'] = df_pos_display['X'].apply(lambda val: f"{val:.0f}" if pd.notnull(val) else "")
        if 'Y' in df_pos_display.columns:
            df_pos_display['Y'] = df_pos_display['Y'].apply(lambda val: f"{val:.0f}" if pd.notnull(val) else "")
            
        print("#### Station Positions (`dfPos`)")
        print(df_pos_display.to_markdown(index=False), "\n")
        
        # 3. Audit dfDiff and inject anchor markers
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


if __name__ == "__main__":
    # Internal component quick-check script logic
    config_file = "DEAKIN/Deakin_Free.toml"
    if os.path.exists(config_file):
        network = LevelingNetwork(config_file)
        network.display()
