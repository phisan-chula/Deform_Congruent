import os
import sys
import argparse
import numpy as np
import pandas as pd
import scipy.stats as stats

class DeformationAnalyzer:
    """
    An OOP Geodetic Deformation Analysis Engine that ingests free-network 
    adjustments from two epochs, performs automated global congruency testing, 
    and applies S-transformations to localize unstable or displaced stations.
    """
    def __init__(self, ep0_path: str, ep1_path: str, alpha: float = 0.05, df_denom: int = 1000):
        self.ep0_path = ep0_path
        self.ep1_path = ep1_path
        self.alpha = alpha
        self.df_denom = df_denom  
        
        self.stations = []
        self.d_raw = None        
        self.Sigma_d = None      
        self.m = 0               
        
        self.global_history = []
        self.local_history = []
        self.final_stable_datum = []
        self.final_deformed_points = []
        
        self._load_and_align_epochs()

    def _load_and_align_epochs(self):
        """Ingests epoch dataframes and explicitly inner-merges them to align spatial records."""
        if not os.path.exists(self.ep0_path) or not os.path.exists(self.ep1_path):
            raise FileNotFoundError(f"One or both epoch files were not found. Verify your paths.")

        df0 = pd.read_csv(self.ep0_path)
        df1 = pd.read_csv(self.ep1_path)

        # Robust cross-referencing to protect against row ordering or mismatched inventory errors
        df_merged = pd.merge(df0, df1, on="Station", suffixes=("_ep0", "_ep1")).sort_values(by="Station").reset_index(drop=True)
        
        self.stations = df_merged["Station"].tolist()
        self.m = len(self.stations)

        x0 = df_merged["Adjusted_Height_m_ep0"].values
        x1 = df_merged["Adjusted_Height_m_ep1"].values
        self.d_raw = x1 - x0

        # Build diagonal stochastic covariance matrix from individual standard deviations (meters squared)
        var0 = (df_merged["StdDev_mm_ep0"].values / 1000.0) ** 2
        var1 = (df_merged["StdDev_mm_ep1"].values / 1000.0) ** 2
        self.Sigma_d = np.diag(var0 + var1)

    def _get_s_matrix(self, datum_mask: np.ndarray) -> np.ndarray:
        """Computes the 1D S-transformation matrix based on the active datum mask."""
        q = np.sum(datum_mask)
        I = np.eye(self.m)
        G = np.ones((self.m, 1))
        E = np.diag(datum_mask.astype(float))
        B = E @ G  
        
        # Rigorous 1D S-Transform computation: S = I - G * (B^T * G)^-1 * B^T
        S = I - (1.0 / q) * G @ B.T
        return S

    def analyze(self):
        """Executes the iterative step-by-step global and local congruency test routine."""
        # Start by assuming all shared stations are stable components of the datum reference frame
        datum_mask = np.ones(self.m, dtype=bool)

        for step in range(self.m):
            q = np.sum(datum_mask)
            if q < 2:
                print("[Warning] Datum components exhausted down to a single element. Terminating audit.")
                break

            # 1. Compute current S-transformation projections
            S = self._get_s_matrix(datum_mask)
            d_j = S @ self.d_raw
            Sigma_dj = S @ self.Sigma_d @ S.T

            # 2. Extract active datum sub-arrays for the Global Test
            datum_indices = np.where(datum_mask)[0]
            d_e = d_j[datum_indices]
            Sigma_ee = Sigma_dj[np.ix_(datum_indices, datum_indices)]

            # 3. Calculate Global Congruency Test statistics
            Re = d_e.T @ np.linalg.pinv(Sigma_ee) @ d_e
            he = q - 1
            F_global = Re / he
            F_crit_global = stats.f.ppf(1 - self.alpha, he, self.df_denom)

            passed = F_global <= F_crit_global
            verdict_str = "PASSED (Stable)" if passed else "FAILED (Unstable)"

            self.global_history.append({
                "Step": step,
                "Active_Datum_Points": q,
                "Re_Quadratic": round(Re, 4),
                "DoF": he,
                "F_Global": round(F_global, 4),
                "F_Crit": round(F_crit_global, 4),
                "Outcome": verdict_str
            })

            # 4. Handle Branch Conditions: Global Test Passes vs Fails
            if passed:
                # Log final localized data rows for stable stations
                for i in range(self.m):
                    if datum_mask[i]:
                        # FIXED: Swapped bitwise XOR (^) with proper float exponentiation (**2)
                        f_local_val = (d_j[i]**2) / Sigma_dj[i, i] if Sigma_dj[i, i] > 0 else 0.0
                        self.local_history.append({
                            "Step": step,
                            "Station": self.stations[i],
                            "Displacement_mm": f"{d_j[i]*1000:+.2f}",
                            "StdDev_mm": f"{np.sqrt(Sigma_dj[i,i])*1000:.2f}",
                            "F_Local": f"{f_local_val:.4f}",
                            "Verdict": "Stable Reference Base"
                        })
                self.final_stable_datum = [self.stations[i] for i in range(self.m) if datum_mask[i]]
                self.final_deformed_points = [self.stations[i] for i in range(self.m) if not datum_mask[i]]
                break
            else:
                # Global test fails -> locate the maximum local outlier in the active pool
                max_f_local = -1.0
                worst_idx = -1

                # Compile local tracking statistics for all active candidates
                step_local_stats = []
                for i in datum_indices:
                    f_local = (d_j[i] ** 2) / Sigma_dj[i, i]
                    step_local_stats.append((i, f_local, d_j[i], np.sqrt(Sigma_dj[i, i])))
                    if f_local > max_f_local:
                        max_f_local = f_local
                        worst_idx = i

                for i, f_local, disp, sdev in step_local_stats:
                    verdict = "**Deformed / Outlier !!!**" if i == worst_idx else "Suspect / Smear"
                    self.local_history.append({
                        "Step": step,
                        "Station": self.stations[i],
                        "Displacement_mm": f"{disp*1000:+.2f}",
                        "StdDev_mm": f"{sdev*1000:.2f}",
                        "F_Local": f"{f_local:.4f}",
                        "Verdict": verdict
                    })

                # Drop the single worst offending station out of the reference matrix and loop again
                datum_mask[worst_idx] = False

    def display_report(self):
        """Generates clear, formatted Markdown summaries of the entire adjustment pipeline."""
        print("## GEODETIC DEFORMATION ANALYSIS REPORT")
        print(f"* **Significance Level (\u03b1):** {self.alpha}")
        print(f"* **Reference Framework Model:** 1D S-Transformation Iterative Localization\n")
        print("="*65 + "\n")

        print("### 1. ITERATIVE GLOBAL CONGRUENCY ANALYSIS")
        print("> Objective: Evaluate whether the reference point groups are statistically stable.\n")
        print(pd.DataFrame(self.global_history).to_markdown(index=False), "\n")
        print("="*65 + "\n")

        print("### 2. LOCAL POINT CONGRUENCY & LOCALIZATION AUDIT")
        print("> Objective: Compute relative displacement metrics and filter out moving components.\n")
        print(pd.DataFrame(self.local_history).to_markdown(index=False), "\n")
        print("="*65 + "\n")

        print("### 3. FINAL GEODETIC INTERPRETATION SUMMARY")
        print(f"* **Verified Stable Reference Datum Base:** {self.final_stable_datum}")
        print(f"* **Isolated Displaced / Deformed Stations:** {self.final_deformed_points}")


def main():
    parser = argparse.ArgumentParser(description="Multi-Epoch Geodetic Network Deformation Analysis Engine.")
    parser.add_argument("-e0", "--epoch0", type=str, default="DEAKIN/result_ep0.csv", help="Path to Epoch 0 results CSV")
    parser.add_argument("-e1", "--epoch1", type=str, default="DEAKIN/result_ep1.csv", help="Path to Epoch 1 results CSV")
    parser.add_argument("-a", "--alpha", type=float, default=0.05, help="Significance value alpha (default: 0.05)")
    args = parser.parse_args()

    try:
        analyzer = DeformationAnalyzer(args.epoch0, args.epoch1, alpha=args.alpha)
        analyzer.analyze()
        analyzer.display_report()
    except Exception as e:
        print(f"CRITICAL: Deformation pipeline failed: {e}")

if __name__ == "__main__":
    main()
