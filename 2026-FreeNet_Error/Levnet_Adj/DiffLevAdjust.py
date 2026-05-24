import os
import sys
import argparse
import numpy as np
import pandas as pd
import lmfit
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import StrMethodFormatter  # For absolute integer formatting

# Ingest your custom network parsing class
from DiffLevNet import LevelingNetwork

class LevelNetworkAdjustment:
    """
    An OOP adjustment engine that accepts a parsed LevelingNetwork,
    configures a weighted least-squares problem via lmfit, 
    and handles both Constrained and Free-Net geodetic configurations.
    """
    def __init__(self, network: LevelingNetwork, fallback_height: float = 100.0):
        self.network = network
        self.params = lmfit.Parameters()
        self.results = None
        self.results_df = pd.DataFrame()
        self.fallback_height = fallback_height  # Reference datum center for Free-Net modes
        self.is_free_net = True

    def _determine_datum_mode(self) -> float:
        """Audits fixed inventory to toggle between Free-Net and Constrained modes."""
        if not self.network.dfFixed.empty and 'Effective' in self.network.dfFixed.columns:
            active_fixed = self.network.dfFixed[self.network.dfFixed['Effective']]
            if not active_fixed.empty:
                self.is_free_net = False
                return float(active_fixed['Height_m'].iloc[0])
        self.is_free_net = True
        return self.fallback_height

    def setup_parameters(self):
        """Maps physical topology to optimization parameters."""
        df_diff = self.network.dfDiff
        active_stations = set(df_diff['StaFrom'].tolist() + df_diff['StaTo'].tolist())
        
        fixed_map = dict(zip(self.network.dfFixed['Station'], self.network.dfFixed['Height_m'])) if not self.network.dfFixed.empty else {}
        effective_fixed = set(self.network.dfFixed[self.network.dfFixed['Effective']]['Station']) if not self.network.dfFixed.empty else set()
        
        initial_height = self._determine_datum_mode()

        for station in active_stations:
            if station in effective_fixed:
                # Constrained Control Anchor Base (Do not vary)
                self.params.add(station, value=fixed_map[station], vary=False)
            else:
                # Unknown target node to estimate
                self.params.add(station, value=initial_height, vary=True)

    def _residual_function(self, params: lmfit.Parameters) -> np.ndarray:
        """Calculates structural observation errors scaled by geometric weights."""
        residuals = []
        for _, row in self.network.dfDiff.iterrows():
            h_from = params[row['StaFrom']].value
            h_to = params[row['StaTo']].value
            
            # Theoretical Model: dh = H_to - H_from
            dh_model = h_to - h_from
            
            # Stochastic Model: Weight P = 1 / Distance (km)
            weight = np.sqrt(1.0 / row['Dist_km'])
            residuals.append(weight * (dh_model - row['Diff_m']))
            
        # --- FREE-NET INNER CONSTRAINT ---
        if self.is_free_net:
            centroid_residual = sum(params[st].value - self.fallback_height for st in params)
            residuals.append(1e4 * centroid_residual)  # Rigid constraint penalty
            
        return np.array(residuals)

    def run(self):
        """Executes the optimization process and builds the results matrix."""
        if not self.params:
            self.setup_parameters()
            
        self.results = lmfit.minimize(self._residual_function, self.params)
        self._compile_results()

    def _compile_results(self):
        """Constructs descriptive results data collection framework."""
        compiled_records = []
        effective_fixed = set(self.network.dfFixed[self.network.dfFixed['Effective']]['Station']) if not self.network.dfFixed.empty else set()
        
        for p_name in self.results.params:
            param = self.results.params[p_name]
            
            if p_name in effective_fixed:
                status = 'FIXED (Anchor)'
            elif self.is_free_net:
                status = f"{(param.stderr * 1000):.2f} (Free-Net)"
            else:
                status = f"{(param.stderr * 1000):.2f}"

            compiled_records.append({
                'Station': p_name,
                'Adjusted_Height_m': param.value,
                'StdDev_mm': status
            })
        self.results_df = pd.DataFrame(compiled_records).sort_values(by="Station")

    def report(self):
        """Outputs statistical fit details and processed measurements."""
        if self.results is None:
            print("[Error] No active solutions found. Call run() first.")
            return
            
        mode_str = "FREE NETWORK ADJUSTMENT" if self.is_free_net else "CONSTRAINED ADJUSTMENT"
        
        print("\n" + "="*50)
        print(f"    LMFIT ENGINE REPORT ({mode_str})    ")
        print("="*50)
        print(lmfit.report_fit(self.results))
        
        print("\n#### ADJUSTED GEODETIC ELEVATIONS")
        print(self.results_df.to_markdown(index=False))
        print("\n" + "="*50)

    def Display_Obs(self, k_factor: float = 12.0):
        """
        Computes absolute observation residuals and flags paths that 
        exceed the maximum allowable error threshold: k_factor * sqrt(Dist_km).
        """
        if self.results is None:
            print("[Error] No active solutions found. Run the adjustment engine first.")
            return

        params = self.results.params
        obs_records = []

        for _, row in self.network.dfDiff.iterrows():
            line_no = int(row['LineNo'])
            frm = row['StaFrom']
            to = row['StaTo']
            diff_m = row['Diff_m']
            dist_km = row['Dist_km']

            h_from = params[frm].value
            h_to = params[to].value
            dh_model = h_to - h_from

            residual_mm = (dh_model - diff_m) * 1000
            max_allowable_mm = k_factor * np.sqrt(dist_km)
            
            remark = ""
            residual_str = f"{residual_mm:+.2f}"
            if abs(residual_mm) > max_allowable_mm:
                remark = "⚠️ Outlier"
                residual_str = f"**{residual_mm:+.2f} !!!**"

            obs_records.append({
                'LineNo': line_no,
                'StaFrom': frm,
                'StaTo': to,
                'Diff_m': f"{diff_m:.3f}",
                'Dist_km': f"{dist_km:.2f}",
                'Residual_mm': residual_str,
                'Allowable_Tol_mm': f"\u00b1{max_allowable_mm:.2f}",
                'Remark': remark
            })

        df_obs = pd.DataFrame(obs_records)
        print("\n#### OBSERVATION ADJUSTMENT & TOLERANCE AUDIT")
        print(f"> Tolerance Rule: abs(Residual_mm) > {k_factor} * \u221a(Dist_km)")
        print(df_obs.to_markdown(index=False))
        print("\n" + "="*50)

    def export_result(self, filename: str = "adjusted_stations.csv"):
        """Compiles and exports adjusted network station attributes into a clean CSV."""
        if self.results is None:
            print("[Error] No active solutions found. Run the adjustment engine first.")
            return

        effective_fixed = set(self.network.dfFixed[self.network.dfFixed['Effective']]['Station']) if not self.network.dfFixed.empty else set()
        export_records = []
        params = self.results.params

        for p_name in params:
            param = params[p_name]
            
            if p_name in effective_fixed and not self.is_free_net:
                std_dev = "0.00"
            else:
                std_dev = f"{(param.stderr * 1000):.2f}" if param.stderr is not None else "0.00"

            if self.is_free_net:
                station_type = "FREENET"
            else:
                if p_name in effective_fixed:
                    station_type = "FIXED"
                else:
                    station_type = "ADJUSTED"

            export_records.append({
                'Station': p_name,
                'Adjusted_Height_m': f"{param.value:.4f}",
                'StdDev_mm': std_dev,
                'type': station_type
            })

        df_export = pd.DataFrame(export_records).sort_values(by="Station")
        df_export.to_csv(filename, index=False)
        
        print("\n#### EXPORTED DATA STATIONS LOG SUMMARY")
        print(f"> Successfully transformed results file written to: {filename}\n")
        print(df_export.to_markdown(index=False))
        print("\n" + "="*50)

    def plot_to_geopackage(self, filename: str = "network.gpkg", k_factor: float = 12.0):
        """Exports adjusted stations and observation lines to a GeoPackage using native TOML coordinates."""
        if self.results is None:
            print("[Error] No active solutions found. Run the adjustment engine first.")
            return

        try:
            import geopandas as gpd
            from shapely.geometry import Point, LineString
        except ImportError:
            print("[Error] geopandas or shapely is missing. Cannot export GeoPackage.")
            return

        df_pos_gpkg = self.network.dfPos.copy()
        coord_map = {row['Station']: (row['X'], row['Y']) for _, row in df_pos_gpkg.iterrows()}
        effective_fixed = set(self.network.dfFixed[self.network.dfFixed['Effective']]['Station']) if not self.network.dfFixed.empty else set()
        params = self.results.params
        project_crs = getattr(self.network, 'crs', 'EPSG:32647')

        # Build Stations Layer
        station_records = []
        for _, row in df_pos_gpkg.iterrows():
            st = row['Station']
            if st in params:
                param = params[st]
                height = param.value
                if st in effective_fixed and not self.is_free_net:
                    std_dev = "0.00"
                    st_type = "FIXED"
                else:
                    std_dev = f"{(param.stderr * 1000):.2f}" if param.stderr is not None else "0.00"
                    st_type = "INNER" if self.is_free_net else "ADJUSTED"
                
                station_records.append({
                    'Station': st,
                    'Adjusted_Height_m': float(f"{height:.4f}"),
                    'StdDev_mm': std_dev,
                    'type': st_type,
                    'X': row['X'],
                    'Y': row['Y']
                })
        
        df_stations_gpkg = pd.DataFrame(station_records)
        gdf_stations = gpd.GeoDataFrame(
            df_stations_gpkg,
            geometry=[Point(r['X'], r['Y']) for _, r in df_stations_gpkg.iterrows()],
            crs=project_crs
        ).drop(columns=['X', 'Y'])

        # Build Lines Layer
        line_records = []
        for _, row in self.network.dfDiff.iterrows():
            frm = row['StaFrom']
            to = row['StaTo']
            if frm in coord_map and to in coord_map:
                h_from = params[frm].value if frm in params else self.fallback_height
                h_to = params[to].value if to in params else self.fallback_height
                dh_model = h_to - h_from
                residual_mm = (dh_model - row['Diff_m']) * 1000
                max_allowable_mm = k_factor * np.sqrt(row['Dist_km'])
                remark = "⚠️ Outlier" if abs(residual_mm) > max_allowable_mm else ""
                
                line_records.append({
                    'LineNo': int(row['LineNo']),
                    'StaFrom': frm,
                    'StaTo': to,
                    'Diff_m': float(f"{row['Diff_m']:.4f}"),
                    'Dist_km': float(f"{row['Dist_km']:.4f}"),
                    'Residual_mm': float(f"{residual_mm:.2f}"),
                    'Allowable_Tol_mm': float(f"{max_allowable_mm:.2f}"),
                    'Remark': remark,
                    'geometry': LineString([coord_map[frm], coord_map[to]])
                })
        
        gdf_lines = gpd.GeoDataFrame(line_records, crs=project_crs)

        gdf_stations.to_file(filename, layer="stations", driver="GPKG")
        gdf_lines.to_file(filename, layer="lines", driver="GPKG")
        print(f"[Spatial Engine] Spatial elements exported cleanly to GeoPackage ({project_crs}): {filename}")

    def plot_adjusted_network(self, k_factor: float = 12.0):
        """Renders planimetric coordinate nodes labeled with integer formatting and 90-degree rotated X ticks."""
        if self.network.dfPos.empty or self.results_df.empty or self.results is None:
            print("[Error] Incomplete dataset coordinates or missing results profiles.")
            return

        coord_map = {row['Station']: (row['X'], row['Y']) for _, row in self.network.dfPos.iterrows()}
        res_map = self.results_df.set_index('Station').to_dict(orient='index')
        effective_fixed = set(self.network.dfFixed[self.network.dfFixed['Effective']]['Station'].tolist()) if not self.network.dfFixed.empty else set()
        all_fixed = set(self.network.dfFixed['Station'].tolist()) if not self.network.dfFixed.empty else set()
        params = self.results.params
        project_crs = getattr(self.network, 'crs', 'EPSG:32647')

        fig, ax = plt.subplots(figsize=(11, 9))

        # --- Draw Vector Arrow Lines and Metadata Text Boxes ---
        for _, row in self.network.dfDiff.iterrows():
            frm, to = row['StaFrom'], row['StaTo']
            if frm not in coord_map or to not in coord_map:
                continue
            x1, y1 = coord_map[frm]
            x2, y2 = coord_map[to]
            
            ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                        arrowprops=dict(arrowstyle="-|>", color='#7f8c8d', lw=1.5, mutation_scale=14), zorder=1)

            mx, my = (x1 + x2) / 2, (y1 + y2) / 2
            
            h_from = params[frm].value
            h_to = params[to].value
            dh_model = h_to - h_from
            residual_mm = (dh_model - row['Diff_m']) * 1000
            max_allowable_mm = k_factor * np.sqrt(row['Dist_km'])
            
            if abs(residual_mm) > max_allowable_mm:
                outlier_flag = " ⚠️"
                box_edge_color = '#e74c3c'  
                box_face_color = '#fdf2f2'  
            else:
                outlier_flag = ""
                box_edge_color = '#bdc3c7'
                box_face_color = '#ffffff'
            
            label_text = (
                f"L{int(row['LineNo'])}: {row['Diff_m']:+.3f} m\n"
                f"{row['Dist_km']:.2f} km: {residual_mm:+.2f} mm{outlier_flag}"
            )
            
            ax.text(mx, my, label_text, color='#2c3e50', fontsize=8, ha='center', va='center',
                    fontweight='bold', bbox=dict(boxstyle='round,pad=0.3', 
                    facecolor=box_face_color, edgecolor=box_edge_color, alpha=0.9, zorder=2))

        # --- Draw Map Nodes ---
        for _, row in self.network.dfPos.iterrows():
            st, x, y = row['Station'], row['X'], row['Y']
            
            if st not in res_map:
                if st in all_fixed:
                    ax.scatter(x, y, color='#95a5a6', marker='s', s=100, zorder=3, edgecolors='black')
                    ax.text(x + 120, y + 120, f"{st}\nUnused BM", fontsize=9, color='#7f8c8d', fontweight='bold')
                continue

            station_data = res_map[st]
            height = station_data['Adjusted_Height_m']
            uncertainty = station_data['StdDev_mm']

            if st in effective_fixed:
                color, marker, size = '#e74c3c', '^', 150
                annotation = f"{st}\n{height:.3f} m\n(Anchor)"
            elif self.is_free_net:
                color, marker, size = '#2ecc71', 'o', 100
                annotation = f"{st}\n{height:.3f} m\n\u00b1{uncertainty.split()[0]} mm"
            else:
                color, marker, size = '#3498db', 'o', 100
                annotation = f"{st}\n{height:.3f} m\n\u00b1{uncertainty} mm"

            ax.scatter(x, y, color=color, marker=marker, s=size, zorder=3, edgecolors='black', lw=1)
            ax.text(x + 120, y + 120, annotation, fontsize=10, fontweight='bold', color='#2c3e50', zorder=4)

        title_mode = "Free-Net Mode" if self.is_free_net else "Constrained Mode"
        ax.set_title(f"Adjusted Leveling Network Map ({title_mode}) - {project_crs}", fontsize=14, fontweight='bold', pad=15)
        ax.set_xlabel("Easting X (m)", labelpad=8)
        ax.set_ylabel("Northing Y (m)", labelpad=8)
        
        # --- APPLIED TICK FORMATTING FORMATS ---
        ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))  # X Axis -> Integer
        ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))  # Y Axis -> Integer
        ax.tick_params(axis='x', rotation=90)                        # Rotate X Axis ticks by 90 degrees
        
        ax.grid(True, linestyle=':', color='#bdc3c7', alpha=0.7)

        if self.is_free_net:
            legend_elements = [
                Line2D([0], [0], marker='o', color='w', label='Free-Net Adjusted Node', 
                       markerfacecolor='#2ecc71', markersize=9, markeredgecolor='black'),
                Line2D([0], [0], color='#7f8c8d', lw=1.5, label='Observed Vector Direction')
            ]
        else:
            legend_elements = [
                Line2D([0], [0], marker='^', color='w', label='Active Fixed Anchor', 
                       markerfacecolor='#e74c3c', markersize=11, markeredgecolor='black'),
                Line2D([0], [0], marker='o', color='w', label='Adjusted Station Node', 
                       markerfacecolor='#3498db', markersize=9, markeredgecolor='black'),
                Line2D([0], [0], color='#7f8c8d', lw=1.5, label='Observed Vector Direction')
            ]
            
        ax.legend(handles=legend_elements, loc='upper left', frameon=True, facecolor='white', framealpha=0.9)
        ax.set_aspect('equal', adjustable='box')
        
        plt.tight_layout()
        plt.savefig( 'test.svg' )
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Adjust Geodetic Leveling Networks from TOML specifications using lmfit Least-Squares."
    )
    parser.add_argument(
        "toml_file", 
        nargs="?",
        default="Deakin_Free.toml",
        help="Path targeting the Network TOML Configuration data file (default: Deakin_Free.toml)"
    )
    parser.add_argument(
        "-k", "--k_factor",
        type=float,
        default=12.0,
        help="Multiplier constant for the allowable tolerance threshold in mm (default: 12.0)"
    )
    parser.add_argument(
        "-f", "--fallback_height",
        type=float,
        default=100.0,
        help="Default elevation height used as datum center for Free-Net adjustments (default: 100.0)"
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="adjusted_stations.csv",
        help="Output CSV target filename path for exported station results (default: adjusted_stations.csv)"
    )
    parser.add_argument(
        "-g", "--gpkg",
        type=str,
        default=None,
        help="Optional output GeoPackage target filename path for spatial layers (e.g., network.gpkg)"
    )
    
    args = parser.parse_args()

    try:
        print(f"[System] Ingesting source architecture layout: {args.toml_file}")
        network_data = LevelingNetwork(args.toml_file)
        network_data.display()
        
        print("\n[System] Initializing Least Squares Engine optimization routines...")
        engine = LevelNetworkAdjustment(network_data, fallback_height=args.fallback_height)
        engine.run()
        
        engine.report()
        engine.Display_Obs(k_factor=args.k_factor)
        engine.export_result(filename=args.output)
        
        if args.gpkg:
            engine.plot_to_geopackage(filename=args.gpkg, k_factor=args.k_factor)
        
        print("[System] Generating planimetric adjustment schema visuals...")
        engine.plot_adjusted_network(k_factor=args.k_factor)

    except Exception as e:
        print(f"\nCRITICAL ENGINE FAILURE: {e}")


if __name__ == "__main__":
    main()
