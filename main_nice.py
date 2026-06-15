import argparse
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
from scipy.spatial import cKDTree
import os
import warnings
import sys
import numpy as np

try:
    from algorithms import (
        GromacsProcessor, ClusteringAnalyzer, 
        OrientationalOrderAnalyzer, ShapeAnisotropyAnalyzer, PDFAnalyzer, FESAnalyzer,
        write_xvg, sync_masses_from_tpr
    )
except ImportError as e:
    print(f"[-] Error: {e}")
    print("[-] Please check if the algorithms.py file has all the required functions and classes.")
    sys.exit(1)

warnings.filterwarnings('ignore')

def parse_arguments():
    description_text = (
        "========================================================================\n"
        "                      Dipeptide Calculation Tool\n"
        "========================================================================"
    )
    parser = argparse.ArgumentParser(
        description=description_text,
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True
    )
    
    group_io = parser.add_argument_group('1. Input/Output Files')
    group_io.add_argument('-f', type=str, required=True, help='Input trajectory file (.xtc)')
    group_io.add_argument('-s', type=str, required=True, help='Input run input file (.tpr)')
    group_io.add_argument('-gro', type=str, required=True, help='Input structure file (.gro)')
    
    group_params = parser.add_argument_group('2. Selection & Calculation Parameters')
    group_params.add_argument('-select', type=str, default="protein", help='Molecule selection string (e.g., "resname GLY PHE")')
    group_params.add_argument('-cutoff_space', type=float, default=0.7, help='Cutoff distance for clustering (nm)')
    group_params.add_argument('-cutoff_cz', type=int, default=20, help='Minimum NUMBER OF MOLECULES to form a valid cluster')
    group_params.add_argument('-n_pep', type=int, default=2, help='Number of residues per molecule (Default: 2)')
    group_params.add_argument('-cutoff_multi', type=float, default=2.0, help="Multiplier for density search radius (Tang's default = 2.0)")
    group_params.add_argument('-sf', type=str, default=None, help='Provide selections from files')
    group_params.add_argument('-selrpos', type=str, default='atom', help='Selection reference positions (atom, res_com, mol_com, etc.)')
    group_params.add_argument('-seltype', type=str, default='atom', help='Default selection output positions')
    group_params.add_argument('--fix-pbc', action='store_true', help='Automatically fix periodic boundary conditions (PBC whole)')
    group_params.add_argument('--gmx-path', type=str, default='gmx', help="Path to the GROMACS executable (Default: 'gmx')")
    
    group_params.add_argument('-b', type=float, default=0.0, help='START time for analysis (ps)')
    group_params.add_argument('-e', type=float, default=-1.0, help='END time for analysis (ps), -1 means to the end of the file')
    group_params.add_argument('-dt', type=float, default=0.0, help='Time step (ps) - Very useful for skipping frames to speed up')

    group_stats = parser.add_argument_group('3. Statistical Output (.xvg)')
    group_stats.add_argument('-nb', type=str, default=None, help='Output number of valid clusters over time')
    group_stats.add_argument('-sz', type=str, default=None, help='Output size of the largest cluster over time')
    group_stats.add_argument('-dc', type=str, default=None, help='Output Degree of Clustering (DC)')
    group_stats.add_argument('-ap', type=str, default=None, help='Output Collapse Degree / Aggregation Propensity (AP)')
    group_stats.add_argument('-fe', type=str, default=None, help='Output Fluctuation Extent of clustering degree (FE)')
    group_stats.add_argument('-p2', type=str, default=None, help='Output Orientational Order Parameter (P2)')
    group_stats.add_argument('-k2', type=str, default=None, help='Output K2 Parameter (Relative Shape Anisotropy)')
    
    group_stats.add_argument('-pdf', type=str, default=None, help='Output Probability Density Function (.xvg)')
    group_stats.add_argument('-fes', type=str, default=None, help='Output Free Energy Surface FES(d,theta) = -RT ln[P(d,theta)] (.xvg)')
    group_stats.add_argument('-fes_temp', type=float, default=300.0, help='Temperature (K) for FES calculation (Default: 300)')
    
    group_stats.add_argument('-comp', type=str, default=None, help='Output co-assembly composition (GF/QW ratio)')
    group_stats.add_argument('-molnumber', type=str, default=None, help='Output total molecules within all valid clusters')
    group_stats.add_argument('-liquidity', type=str, default=None, help='Output liquidity and fraction of preservation')
    group_stats.add_argument('-density', type=str, default=None, help='Output density of aggregated and dispersed phases')
    group_stats.add_argument('-sasa', type=str, default=None, help='Output Solvent Accessible Surface Area (auto-selects "Protein" group)')
    
    group_pdb = parser.add_argument_group('4. PDB Snapshot Output')
    group_pdb.add_argument('-pdb', type=str, default=None, help='Output the largest cluster structure (.pdb)')
    group_pdb.add_argument('-pdb_system', type=str, default=None, help='Output all valid clusters at original physical coordinates (.pdb)')
    group_pdb.add_argument('-pdb_time', type=float, default=-1, help='Time frame to extract the PDB snapshot (ps)')
    group_pdb.add_argument('-all_clusters', action='store_true', help='Export ALL valid clusters (Used in combination with -pdb)')
    
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    if args.sf:
        if os.path.exists(args.sf):
            with open(args.sf, 'r') as file:
                args.select = file.read().strip()
        else:
            print(f"[-] Warning: File '{args.sf}' not found.")

    final_xtc = args.f
    temp_files = []

    if args.fix_pbc:
        gmx = GromacsProcessor(custom_path=args.gmx_path)
        temp_xtc = "tmp_whole_trajectory.xtc"
        if gmx.make_whole(args.s, args.f, temp_xtc):
            final_xtc = temp_xtc
            temp_files.append(temp_xtc)
        else:
            return

    print(f"\n[+] Loading trajectory: {args.gro} + {final_xtc}")
    try:
        u = mda.Universe(args.gro, final_xtc)
        sync_masses_from_tpr(u, args.s, args.gmx_path)
        
        peptide_group = u.select_atoms(args.select)
        residues = peptide_group.residues
        res_names = residues.names
        num_residues = len(residues)
        resindex_to_idx = {res.resindex: i for i, res in enumerate(residues)}
        
        if num_residues == 0:
            print("[-] Error: No molecules found matching the selection.")
            return
            
        num_mols = num_residues // args.n_pep
        print(f"[i] Topology: {args.n_pep} residues/molecule | Total molecules: {num_mols}")

        # =========================================================
        # SYSTEM CONCENTRATION & AVERAGE VOLUME CALCULATION
        # =========================================================
        print("\n[+] Calculating average system volume over the trajectory...")
        v_nm3_list = []
        # Fast iteration to get box dimensions
        for ts in u.trajectory:
            box = ts.dimensions[:3]
            v_nm3_list.append((box[0] / 10.0) * (box[1] / 10.0) * (box[2] / 10.0))
        
        # Rewind trajectory for main analysis
        u.trajectory.rewind()
        
        avg_v_nm3 = np.mean(v_nm3_list) if v_nm3_list else 0.0
        avg_v_ml = avg_v_nm3 * 1e-21
        
        total_mass_g_mol = peptide_group.residues.atoms.total_mass()
        
        if num_mols > 0:
            mass_per_mol = total_mass_g_mol / num_mols
        else:
            mass_per_mol = 0.0
            
        NA = 6.02214076e23
        mol_count = num_mols / NA
        mass_g = mol_count * round(mass_per_mol)
        mass_mg = mass_g * 1000
        
        conc_mg_ml = mass_mg / avg_v_ml if avg_v_ml > 0 else 0.0

        print("========================================================")
        print("                 SYSTEM INFORMATION                     ")
        print("========================================================")
        print(f" Molecule Selection   : {args.select}")
        print(f" Number of Molecules  : {num_mols}")
        print(f" Mass per Molecule    : {mass_per_mol:.0f} g/mol")
        print(f" Average Volume (nm^3): {avg_v_nm3:.3f} nm^3")
        print(f" Average Volume (mL)  : {avg_v_ml:.4e} mL")
        print("--------------------------------------------------------")
        print(f" 👉 Average Concentration : {conc_mg_ml:.3f} mg/mL")
        print("========================================================\n")

    except Exception as e:
        print(f"[-] Error loading trajectory: {e}")
        return

    cluster_calc = ClusteringAnalyzer(args.cutoff_space, args.n_pep)
    pdf_analyzer = PDFAnalyzer(u, args.cutoff_space, args.n_pep)
    fes_analyzer = FESAnalyzer(u, args.cutoff_space, args.n_pep) if args.fes else None
    
    times, cluster_counts, max_sizes, dc_data, ap_data = [], [], [], [], []
    comp_data, molnumber_data, liquidity_data, density_data = [], [], [], []
    
    prev_cluster = set()
    prev_contacts = set()
    pdb_extracted = False 

    print("[+] Starting trajectory scanning...")
    next_time = args.b
    
    for ts in u.trajectory:
        if ts.time < args.b - 0.001: continue
        if args.e >= 0 and ts.time > args.e + 0.001: break
        if args.dt > 0 and ts.time < next_time - 0.001: continue
            
        times.append(ts.time)
        
        all_clusters, _ = cluster_calc.calculate(peptide_group, num_residues, resindex_to_idx, u.dimensions)
        all_clusters_mols = []
        for c in all_clusters:
            mols = set(res_idx // args.n_pep for res_idx in c)
            all_clusters_mols.append(mols)

        qualified_clusters_mols = [mols for mols in all_clusters_mols if len(mols) >= args.cutoff_cz]
        cluster_counts.append(len(qualified_clusters_mols))
        molnumber_data.append(sum(len(mols) for mols in qualified_clusters_mols))

        if args.comp and qualified_clusters_mols:
            largest_c_mols = max(qualified_clusters_mols, key=len)
            n_gf_mol = sum(1 for m in largest_c_mols if res_names[m * args.n_pep] in ['GLY', 'PHE'])
            n_qw_mol = sum(1 for m in largest_c_mols if res_names[m * args.n_pep] in ['GLN', 'TRP'])
            total_mol = n_gf_mol + n_qw_mol
            if total_mol > 0:
                comp_data.append([n_gf_mol, n_qw_mol, (n_gf_mol/total_mol)*100, (n_qw_mol/total_mol)*100])
            else:
                comp_data.append([0, 0, 0.0, 0.0])
        elif args.comp:
            comp_data.append([0, 0, 0.0, 0.0])

        current_max_size = max([len(mols) for mols in all_clusters_mols]) if all_clusters_mols else 0
        max_sizes.append(current_max_size)
        
        if args.dc:
            dc_data.append(current_max_size / num_mols if num_mols > 0 else 0.0)
            
        if args.ap:
            total_mols_in_clusters = sum(len(mols) for mols in qualified_clusters_mols)
            ap_data.append(total_mols_in_clusters / num_mols if num_mols > 0 else 0.0)

        buggy_cluster_this_frame = set()
        for mols in all_clusters_mols:
            mols_list = sorted(list(mols))
            if len(mols_list) > 1:
                for m in mols_list[1:]:
                    buggy_cluster_this_frame.add(m)

        if args.liquidity:
            current_contacts = set()
            if qualified_clusters_mols:
                agg_res_indices = []
                for mols in qualified_clusters_mols:
                    for m in mols:
                        for i in range(args.n_pep):
                            agg_res_indices.append(m * args.n_pep + i)
                agg_atoms = residues[agg_res_indices].atoms
                cutoff_A = args.cutoff_space * 10.0
                pairs = capped_distance(agg_atoms.positions, agg_atoms.positions, max_cutoff=cutoff_A, box=u.dimensions, return_distances=False)
                if len(pairs) > 0:
                    mol_ids = agg_atoms.resindices // args.n_pep
                    mol_pairs = mol_ids[pairs]
                    valid_pairs = np.sort(mol_pairs, axis=1)
                    unique_pairs = np.unique(valid_pairs, axis=0)
                    current_contacts = set(map(tuple, unique_pairs))

            if ts.frame == 0 or len(times) == 1:
                val0, val1, val2, val3, val4 = 0.0, len(buggy_cluster_this_frame) / num_mols, 0.0, len(buggy_cluster_this_frame) / num_mols, 0.0
            else:
                if len(prev_cluster) == 0:
                    val0, val1, val2, val3 = 0.0, len(buggy_cluster_this_frame) / num_mols, 0.0, len(buggy_cluster_this_frame) / num_mols
                else:
                    intersect = buggy_cluster_this_frame.intersection(prev_cluster)
                    val0 = len(intersect) / len(prev_cluster)
                    val1 = (len(buggy_cluster_this_frame) - len(intersect)) / num_mols
                    val2 = (len(prev_cluster) - len(intersect)) / num_mols
                    val3 = len(buggy_cluster_this_frame) / num_mols

                if len(prev_contacts) == 0:
                    val4 = 0.0
                else:
                    intersect_contacts = current_contacts.intersection(prev_contacts)
                    val4 = len(intersect_contacts) / len(prev_contacts)

            prev_cluster = buggy_cluster_this_frame
            prev_contacts = current_contacts
            liquidity_data.append([val0, val1, val2, val3, val4])

        if args.density:
            agg_res_indices_den = []
            for m in buggy_cluster_this_frame:
                for i in range(args.n_pep):
                    agg_res_indices_den.append(m * args.n_pep + i)
            if len(agg_res_indices_den) > 0:
                aggr_atoms = residues[agg_res_indices_den].atoms
                disp_atoms = peptide_group.subtract(aggr_atoms)
            else:
                aggr_atoms = u.select_atoms("none")
                disp_atoms = peptide_group
                
            mass_aggr = aggr_atoms.total_mass() if len(aggr_atoms) > 0 else 0.0
            mass_disp = disp_atoms.total_mass() if len(disp_atoms) > 0 else 0.0
            box_dims = u.dimensions[:3]
            x_grid, y_grid, z_grid = np.arange(0, box_dims[0], 1.0), np.arange(0, box_dims[1], 1.0), np.arange(0, box_dims[2], 1.0)
            X, Y, Z = np.meshgrid(x_grid, y_grid, z_grid, indexing='ij')
            grid_points = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T
            total_grid_points = len(grid_points)
            radius_multi_A = args.cutoff_space * 10.0 * args.cutoff_multi
            
            if len(aggr_atoms) > 0:
                tree_aggr = cKDTree(aggr_atoms.positions % box_dims, boxsize=box_dims)
                distances_aggr, _ = tree_aggr.query(grid_points, k=1, distance_upper_bound=radius_multi_A)
                num_aggr_points = np.sum(distances_aggr <= radius_multi_A)
            else:
                num_aggr_points = 0
                
            num_disp_points = total_grid_points - num_aggr_points
            vol_aggr_nm3, vol_disp_nm3 = num_aggr_points / 1000.0, num_disp_points / 1000.0
            den_aggr = (mass_aggr / vol_aggr_nm3 / 0.602) if vol_aggr_nm3 > 0 else 0.0
            den_disp = (mass_disp / vol_disp_nm3 / 0.602) if vol_disp_nm3 > 0 else 0.0
            density_data.append([den_aggr, den_disp])

        if (args.pdb or args.pdb_system) and not pdb_extracted and abs(ts.time - args.pdb_time) <= 0.1:
            print(f"\n[+] CAPTURED TARGET FRAME {ts.time} ps! Processing PDB output...")
            if all_clusters_mols:
                largest_c_mols = max(all_clusters_mols, key=len)
                agg_res_indices_largest = []
                for m in largest_c_mols:
                    for i in range(args.n_pep):
                        agg_res_indices_largest.append(m * args.n_pep + i)
                largest_cluster_atoms = residues[agg_res_indices_largest].atoms
                box_dims = u.dimensions[:3]
                cluster_com = largest_cluster_atoms.center_of_geometry(pbc=True)
                u.atoms.translate((box_dims / 2.0) - cluster_com)
                u.atoms.positions = u.atoms.positions % box_dims
                
                valid_clusters_atoms = None
                if qualified_clusters_mols:
                    all_agg_res_indices = []
                    for mols in qualified_clusters_mols:
                        for m in mols:
                            for i in range(args.n_pep):
                                all_agg_res_indices.append(m * args.n_pep + i)
                    valid_clusters_atoms = residues[all_agg_res_indices].atoms

                if args.pdb:
                    if args.all_clusters:
                        if valid_clusters_atoms:
                            valid_clusters_atoms.write(args.pdb)
                        else:
                            print(f"  -> WARNING: No clusters met cutoff_cz to export all_clusters!")
                    else:
                        largest_cluster_atoms.write(args.pdb)
                        
                if args.pdb_system:
                    if valid_clusters_atoms:
                        valid_clusters_atoms.write(args.pdb_system)
                    else:
                        print(f"  -> WARNING: No clusters met cutoff_cz to export pdb_system!")
            pdb_extracted = True

        if args.dt > 0: next_time += args.dt
        if len(times) % 100 == 0:
            print(f" Time: {ts.time:10.1f} ps | Valid Clusters: {len(qualified_clusters_mols):3d} | Max Size: {current_max_size:3d}")

    print("\n[+] Exporting analysis results...")
    write_xvg(args.nb, "Number of Clusters", "Count", times, cluster_counts)
    write_xvg(args.sz, "Size of Largest Cluster", "Molecules", times, max_sizes)
    if args.dc: write_xvg(args.dc, "Degree of Clustering (DC)", "DC", times, dc_data)
    if args.ap: write_xvg(args.ap, "Collapse Degree / Aggregation Propensity (AP)", "AP", times, ap_data)
    write_xvg(args.molnumber, "Total Molecules in Clusters", "Molecules", times, molnumber_data)
    if args.comp: write_xvg(args.comp, "Cluster Composition", "Value", times, comp_data, ["GF_Count", "QW_Count", "GF_Percent", "QW_Percent"])
    if args.liquidity:
        legs = ["Fraction of preservation", "Fraction of cluster growth", "Fraction of cluster shrink", "Fraction of aggregation", "Preservation of contacts"]
        write_xvg(args.liquidity, "Fraction of preservation", "Fraction", times, liquidity_data, legs, is_liquidity=True)
    if args.density: write_xvg(args.density, "Density of Phases", "Density (g/cm^3)", times, density_data, ["Aggregated Phase", "Dispersed Phase"])
    if args.sasa: GromacsProcessor(custom_path=args.gmx_path).calculate_sasa(args.s, final_xtc, args.sasa, selection="Protein")

    if args.fe:
        dc_array = [sz / num_mols if num_mols > 0 else 0.0 for sz in max_sizes]
        if len(dc_array) > 0:
            fe_value = np.std(dc_array)
            write_xvg(args.fe, "Fluctuation Extent (FE)", "FE (RMSF of Clustering Degree)", times, [fe_value] * len(times))

    if args.p2:
        pdb_target = args.pdb_system if args.pdb_system else args.pdb
        if not pdb_target or not os.path.exists(pdb_target):
            print(f"[-] Error: Flag -p2 requires exporting a PDB file from the largest Cluster.")
        else:
            try:
                Q_matrix, eigenvalues, P2_val, N_total, RESIDUE_N, RESIDUE_C = OrientationalOrderAnalyzer(pdb_target).calculate_p2()
                if Q_matrix is not None:
                    with open(args.p2, 'w') as f:
                        f.write(f"# Orientational Order Parameter (P2) = {P2_val:.4f}\n")
            except Exception as e: print(f"[-] Error during P2 analysis: {e}")

    if args.k2:
        pdb_target = args.pdb_system if args.pdb_system else args.pdb
        if not pdb_target or not os.path.exists(pdb_target):
            print(f"[-] Error: Flag -k2 requires exporting a PDB file from the largest Cluster.")
        else:
            try:
                S_matrix, eigenvalues, k2_val = ShapeAnisotropyAnalyzer(pdb_target).calculate_k2()
                if S_matrix is not None:
                    with open(args.k2, 'w') as f:
                        f.write(f"# Relative Shape Anisotropy (K2) = {k2_val:.4f}\n")
            except Exception as e: print(f"[-] Error during K2 analysis: {e}")

    # =========================================================
    # 11. CALCULATE PDF OF METRICS
    # =========================================================
    if hasattr(args, 'pdf') and args.pdf and pdf_analyzer:
        print(f"\n[+] Preparing for PDF analysis...")
        try:
            sel_a_str, sel_b_str, sel_a_desc, sel_b_desc = pdf_analyzer.interactive_group_selection()
            
            while True:
                mode_input = input("\nDo you want to calculate Contact Number (type 'cn') or Shortest Distance (type 'd')? ").lower().strip()
                if mode_input in ['cn', 'contact', 'contact number']:
                    calc_mode = 'contact'
                    break
                elif mode_input in ['d', 'distance']:
                    calc_mode = 'distance'
                    break
                else:
                    print("[-] Invalid selection. Please type 'cn' or 'd'.")

            x_range_val = None
            if calc_mode == 'distance':
                while True:
                    try:
                        max_r_input = input("[-] Enter the Max value of the X-axis for the Distance plot (positive number, default Min is 0): ")
                        max_r = float(max_r_input)
                        if max_r <= 0:
                            print("    Please enter a number greater than 0!")
                            continue
                        x_range_val = [0.0, max_r]
                        break
                    except ValueError:
                        print("    Please enter a valid number!")
            elif calc_mode == 'contact':
                 while True:
                    try:
                        max_c_input = input("[-] Enter the Max value of the X-axis for the Contact plot (Press Enter to auto-calculate): ")
                        if not max_c_input.strip():
                            break
                        max_c = float(max_c_input)
                        x_range_val = [0.0, max_c]
                        break
                    except ValueError:
                        print("    Please enter a valid number!")

            try:
                bw_input = input("[-] Enter the Bandwidth for the plot smoothing (e.g., 0.2, 0.5, 1.0) [Press Enter to use default 0.5]: ").strip()
                bw_val = float(bw_input) if bw_input else 0.5
            except ValueError:
                print("    Input error, automatically using default Bandwidth = 0.5")
                bw_val = 0.5

            pdf_times, metric_data = pdf_analyzer.calculate_metric(sel_a_str, sel_b_str, mode=calc_mode, b=args.b, e=args.e, dt=args.dt)
            
            if metric_data is not None:
                title_desc = f"{sel_a_desc} vs {sel_b_desc} ({calc_mode})"
                pdf_analyzer.generate_pdf(metric_data, pdf_times, args.pdf, title_desc, mode=calc_mode, x_range_val=x_range_val, bw=bw_val)
                
        except Exception as e:
            print(f"[-] Error during PDF analysis: {e}")

    # =========================================================
    # 12. CALCULATE FREE ENERGY SURFACE (FES)
    # =========================================================
    if hasattr(args, 'fes') and args.fes and fes_analyzer:
        print(f"\n[+] Preparing for FES (Free Energy Surface) analysis...")
        try:
            # Reuse the interactive group selection mechanism from PDFAnalyzer for synchronization
            sel_a_str, sel_b_str, sel_a_desc, sel_b_desc = pdf_analyzer.interactive_group_selection()
            
            title_desc = f"{sel_a_desc} vs {sel_b_desc}"
            
            # Collect d and theta data from trajectory
            d_data, angle_data = fes_analyzer.calculate_fes_trajectory(
                sel_a_str, sel_b_str, b=args.b, e=args.e, dt=args.dt
            )
            
            if d_data is not None and angle_data is not None:
                fes_analyzer.generate_fes(
                    d_data, angle_data, output_file=args.fes, 
                    title_desc=title_desc, T=args.fes_temp
                )
        except Exception as e:
            print(f"[-] Error during FES analysis: {e}")

    for f in temp_files:
        if os.path.exists(f): os.remove(f)
    print("\n[✔] All tasks completed successfully!")

if __name__ == "__main__":
    main()
