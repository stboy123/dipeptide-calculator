# Dipeptide Assembly Analyzer
*A tool for calculating dipeptide phase separation*

A comprehensive tool designed to analyze the phase separation and self-assembly behavior of dipeptides. To run the analyzer, you must provide a trajectory file, a system structure/topology file, and activate your desired calculation flags.

## 🚀 Quick Start Guide

**General Command Structure:**
```bash
python main_nice.py -f <trajectory.xtc> -s <topology.tpr> -gro <structure.gro> -select "<molecule_group>" -cutoff_space 0.7 -cutoff_cz 20 [ANALYSIS FLAGS]
```

---

## 💡 Analysis Examples

To help you get started, here are specific examples grouped by the type of analysis.

### 1. Basic Cluster Analysis
* **Calculate the number of clusters:**
  ```bash
  python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -nb number.xvg --fix-pbc
  ```
* **Track the size of the largest cluster over time:**
  ```bash
  python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -sz size.xvg --fix-pbc
  ```
* **Count the number of molecules within clusters:**
  ```bash
  python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -molnumber molecules_in_clusters.xvg --fix-pbc
  ```

### 2. Physical & Structural Properties
* **Calculate liquidity:**
  ```bash
  python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -liquidity liquidity_factor.xvg --fix-pbc
  ```
* **Calculate density:**
  ```bash
  python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -density density.xvg --fix-pbc
  ```
* **Calculate Solvent Accessible Surface Area (SASA):** *(Automatically selects the "Protein" group)*
  ```bash
  python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -sasa sasa_area.xvg --fix-pbc
  ```

### 3. Advanced Thermodynamics & Anisotropy
* **Calculate the Degree of Clustering (DC) and Aggregation Propensity (AP):**
  ```bash
  python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -dc degree_of_clustering.xvg -ap aggregation_propensity.xvg --fix-pbc
  ```
* **Calculate Orientational Order (P2) and Relative Shape Anisotropy (K2):**
  ```bash
  python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -p2 order_parameter.xvg -k2 shape_anisotropy.xvg --fix-pbc
  ```
* **Calculate the Free Energy Surface (FES) at 300K:**
  ```bash
  python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -fes fes_output.xvg -fes_temp 300 --fix-pbc
  ```
* **Analyze Co-assembly Composition Ratio (e.g., GF/QW):**
  ```bash
  python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -comp composition_ratio.xvg --fix-pbc
  ```

### 4. Extracting Structural Snapshots (PDB)
* **Export the structure file of the largest cluster at a specific time (e.g., 2500000 ps):**
  ```bash
  python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -pdb_system system_python.pdb -pdb_time 2500000
  ```

---

## 🛠 Arguments & Parameters

You can use the command `python main_nice.py -h` at any time to view the detailed help menu directly in your terminal.

### 1. Input & Output Files
| Flag | Description |
| :--- | :--- |
| `-f` | Input trajectory file (`.xtc`) |
| `-s` | Input run input file (`.tpr`) |
| `-gro` | Input structure file (`.gro`) |

### 2. Selection & Calculation Parameters
| Flag | Description |
| :--- | :--- |
| `-select` | Molecule selection string (e.g., `"resname GLY PHE"`) |
| `-cutoff_space` | Cutoff distance for clustering in nm (Default is usually `0.7`) |
| `-cutoff_cz` | Minimum **number of molecules** required to form a valid cluster |
| `-n_pep` | Number of residues per molecule (Default: `2`) |
| `-cutoff_multi` | Multiplier for density search radius (Tang's default = `2.0`) |
| `-sf` | Provide selections from files |
| `-selrpos` | Selection reference positions (e.g., `atom`, `res_com`, `mol_com`) |
| `-seltype` | Default selection output positions |
| `--fix-pbc` | Automatically fix periodic boundary conditions (PBC whole) |
| `--gmx-path` | Path to the GROMACS executable (Default: `'gmx'`) |
| `-b` | START time for analysis in ps |
| `-e` | END time for analysis in ps (`-1` means read to the end of the file) |
| `-dt` | Time step in ps (Highly useful for skipping frames to speed up analysis) |

### 3. Statistical Output (`.xvg`)
| Flag | Description |
| :--- | :--- |
| `-nb` | Output the number of valid clusters over time |
| `-sz` | Output the size of the largest cluster over time |
| `-dc` | Output Degree of Clustering (DC) |
| `-ap` | Output Collapse Degree / Aggregation Propensity (AP) |
| `-fe` | Output Fluctuation Extent of clustering degree (FE) |
| `-p2` | Output Orientational Order Parameter (P2) |
| `-k2` | Output K2 Parameter (Relative Shape Anisotropy) |
| `-pdf` | Output Probability Density Function |
| `-fes` | Output Free Energy Surface `FES(d,theta) = -RT ln[P(d,theta)]` |
| `-fes_temp` | Temperature (K) for FES calculation (Default: `300`) |
| `-comp` | Output co-assembly composition ratio (GF/QW ratio) |
| `-molnumber` | Output the total number of molecules within all valid clusters |
| `-liquidity` | Output liquidity factor and fraction of preservation |
| `-density` | Output the density of aggregated and dispersed phases |
| `-sasa` | Output Solvent Accessible Surface Area (Automatically selects "Protein" group) |

### 4. PDB Snapshot Output
| Flag | Description |
| :--- | :--- |
| `-pdb` | Output the structure of the largest cluster (`.pdb`) |
| `-pdb_system` | Output all valid clusters at their original physical coordinates (`.pdb`) |
| `-pdb_time` | Specific time frame to extract the PDB snapshot (in ps) |
| `-all_clusters`| Flag to process or extract data for all valid clusters |

---

