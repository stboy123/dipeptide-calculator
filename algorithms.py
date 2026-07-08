import os
import subprocess
import numpy as np
import MDAnalysis as mda
from scipy.spatial import cKDTree
from MDAnalysis.lib.distances import capped_distance
import multiprocessing
import re

# =========================================================================
# MULTIPROCESSING WORKER FOR FES
# =========================================================================
def _fes_worker_task(task_args):
    gro_file, xtc_file, sel_a_str, sel_b_str, frame_indices, n_pep = task_args
    import warnings
    warnings.filterwarnings('ignore')

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"

    is_pipi = False
    sel_check = (sel_a_str + " " + sel_b_str).lower()
    if "sc" in sel_check or "not name bb" in sel_check:
        is_pipi = True

    u = mda.Universe(gro_file, xtc_file)
    protein_res = u.select_atoms("protein").residues
    num_mols = len(protein_res) // n_pep
    
    mols_a = []
    mols_b = []
    mol_ids_a = []
    mol_ids_b = []
    
    for i in range(num_mols):
        mol_r = protein_res[i*n_pep : (i+1)*n_pep]
        a_atoms = mol_r.atoms.select_atoms(sel_a_str)
        b_atoms = mol_r.atoms.select_atoms(sel_b_str)
        
        if len(a_atoms) > 0:
            mols_a.append(a_atoms)
            mol_ids_a.append(i)
        if len(b_atoms) > 0:
            mols_b.append(b_atoms)
            mol_ids_b.append(i)

    if len(mols_a) == 0 or len(mols_b) == 0:
        return [], []

    box_dims = u.dimensions[:3]
    box_center = box_dims / 2.0
    
    dist_list = []
    angle_list = []

    def get_orientation_vector(pos, center):
        if len(pos) >= 3:
            n = np.cross(pos[1] - pos[0], pos[2] - pos[0])
            norm = np.linalg.norm(n)
            return n / norm if norm > 1e-6 else np.array([0, 0, 1])
        elif len(pos) >= 2:
            v = pos[1] - pos[0]
            norm = np.linalg.norm(v)
            return v / norm if norm > 1e-6 else np.array([0, 0, 1])
        else:
            v = pos[0] - center
            norm = np.linalg.norm(v)
            return v / norm if norm > 1e-6 else np.array([0, 0, 1])

    for frame in frame_indices:
        ts = u.trajectory[frame] 
        
        cog_a = np.array([m.center_of_geometry() for m in mols_a])
        cog_b = np.array([m.center_of_geometry() for m in mols_b])

        tree = cKDTree(cog_b % box_dims, boxsize=box_dims)

        dists, indices = tree.query(cog_a % box_dims, k=5, workers=1)
        dists = np.atleast_2d(dists)
        indices = np.atleast_2d(indices)

        for i in range(len(mols_a)):
            mol_id_i = mol_ids_a[i]
            best_d = 999.0
            best_idx_b = -1
            
            for k_idx in range(dists.shape[1]):
                idx_b = indices[i, k_idx]
                d = dists[i, k_idx] / 10.0
                mol_id_j = mol_ids_b[idx_b]
                
                if mol_id_i != mol_id_j:
                    best_d = d
                    best_idx_b = idx_b
                    break

            if best_d <= 1.2 and best_idx_b != -1:
                atoms_i = mols_a[i].positions
                atoms_j = mols_b[best_idx_b].positions
                
                vec_i = get_orientation_vector(atoms_i, box_center)
                vec_j = get_orientation_vector(atoms_j, box_center)

                cosine_angle = np.dot(vec_i, vec_j)
                angle = np.degrees(np.arccos(np.clip(cosine_angle, -1.0, 1.0)))
                
                if is_pipi and angle > 90:
                    angle = 180 - angle

                dist_list.append(best_d)
                angle_list.append(angle)

    return dist_list, angle_list

class GromacsProcessor:
    def __init__(self, custom_path='gmx'):
        self.gmx_path = custom_path

    def make_whole(self, tpr_file, xtc_file, output_xtc):
        print(f"\n[+] Fixing PBC for {xtc_file}...")
        cmd = f"echo 0 | {self.gmx_path} trjconv -s {tpr_file} -f {xtc_file} -o {output_xtc} -pbc whole"
        try:
            result = subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if result.returncode == 0:
                print(f"[✔] PBC fixed successfully! Output: {output_xtc}")
                return True
            else:
                return False
        except Exception as e:
            return False

    def calculate_sasa(self, tpr_file, xtc_file, output_xvg, selection="Protein"):
        print(f"\n[+] Calculating SASA for {xtc_file}...")
        cmd = f"echo {selection} | {self.gmx_path} sasa -s {tpr_file} -f {xtc_file} -o {output_xvg}"
        try:
            result = subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if result.returncode == 0:
                print(f"[✔] SASA calculated successfully! Output: {output_xvg}")
                return True
            else:
                return False
        except Exception as e:
            return False

    def process_ap_from_sasa(self, tpr_file, xtc_file, selection="Protein"):
        """
        NEW FEATURE: Runs SASA, reads the output, calculates AP (SASA_0 / SASA_t),
        and returns the times and AP values. Does not interfere with other files.
        """
        temp_sasa_xvg = "temp_sasa_for_ap_calculation.xvg"
        print(f"\n[+] Running background SASA to compute Aggregation Propensity (AP)...")
        cmd = f"echo {selection} | {self.gmx_path} sasa -s {tpr_file} -f {xtc_file} -o {temp_sasa_xvg}"
        
        try:
            result = subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if result.returncode != 0 or not os.path.exists(temp_sasa_xvg):
                print("[-] Failed to calculate SASA for AP calculation.")
                return None
        except Exception as e:
            print(f"[-] Exception while running SASA for AP: {e}")
            return None

        # Parse the temporary SASA file
        times = []
        sasa_values = []
        with open(temp_sasa_xvg, 'r') as f:
            for line in f:
                if line.startswith(('@', '#')): 
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    times.append(float(parts[0]))
                    sasa_values.append(float(parts[1]))

        if not sasa_values:
            print("[-] No valid SASA data extracted.")
            return None

        # Calculate AP: SASA_initial / SASA_current
        sasa_0 = sasa_values[0]
        ap_values = []
        for s in sasa_values:
            if s > 0:
                ap_values.append(sasa_0 / s)
            else:
                ap_values.append(1.0) # Prevent division by zero fallback

        # Clean up temporary file
        try:
            os.remove(temp_sasa_xvg)
        except OSError:
            pass

        return times, ap_values


class ClusteringAnalyzer:
    def __init__(self, cutoff_space, n_pep):
        self.cutoff_space = cutoff_space
        self.n_pep = n_pep

    def calculate(self, atom_group, num_residues, resindex_to_idx, box_dimensions):
        positions = atom_group.positions
        cutoff_A = self.cutoff_space * 10.0
        pairs = capped_distance(positions, positions, max_cutoff=cutoff_A, box=box_dimensions, return_distances=False)
        
        resindices = atom_group.resindices
        mol_ids = resindices // self.n_pep
        num_mols = num_residues // self.n_pep
        adj_list = {i: set() for i in range(num_mols)}
        
        if len(pairs) > 0:
            mol_pairs = mol_ids[pairs]
            valid_pairs = mol_pairs[mol_pairs[:, 0] != mol_pairs[:, 1]]
            for p1, p2 in valid_pairs:
                adj_list[p1].add(p2)
                adj_list[p2].add(p1)
                
        visited = set()
        clusters = []
        
        for i in range(num_mols):
            if i not in visited:
                queue = [i]
                current_cluster = set([i])
                visited.add(i)
                
                while queue:
                    node = queue.pop(0)
                    for neighbor in adj_list[node]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            current_cluster.add(neighbor)
                            queue.append(neighbor)
                            
                cluster_res_indices = []
                for mol_idx in current_cluster:
                    for j in range(self.n_pep):
                        cluster_res_indices.append(mol_idx * self.n_pep + j)
                clusters.append(cluster_res_indices)
                
        return clusters, adj_list

class LiquidityAnalyzer: pass
class DensityAnalyzer: pass

def write_xvg(filename, title, y_label, time_data, value_data, legends=None, is_liquidity=False):
    if filename is None: return
    with open(filename, 'w') as f:
        f.write('# File generated by Loc\'s Analysis Tool\n')
        f.write(f'@    title "{title}"\n')
        f.write(f'@    xaxis  label "Time (ps)"\n')
        f.write(f'@    yaxis  label "{y_label}"\n')
        f.write('@TYPE xy\n')
        
        if is_liquidity:
            f.write('@ view 0.15, 0.15, 0.75, 0.85\n')
            f.write('@ legend on\n')
            f.write('@ legend box on\n')
            f.write('@ legend loctype view\n')
            f.write('@ legend 0.78, 0.8\n')
            f.write('@ legend length 2\n')
            
        if legends:
            for i, leg in enumerate(legends):
                f.write(f'@ s{i} legend "{leg}"\n')
                
        for i in range(len(time_data)):
            val = value_data[i]
            line = f"{time_data[i]:10.3f} "
            if isinstance(val, (list, tuple, np.ndarray)):
                line += " ".join([f"{v:.7f}" for v in val])
            else:
                line += f"{val:12.5f}"
            f.write(line + "\n")
    print(f"[*] Successfully written: {filename}")

def sync_masses_from_tpr(u, tpr_file, gmx_path):
    print(f"\n[+] Synchronizing masses from {tpr_file} to ensure 0% error...")
    try:
        u_tpr = mda.Universe(tpr_file)
        if len(u_tpr.atoms) == len(u.atoms):
            u.atoms.masses = u_tpr.atoms.masses
            print("[✔] Method 1: Loaded exact masses from TPR via MDAnalysis!")
            return
    except Exception: pass

    print(f"[i] MDAnalysis failed. Parsing raw GROMACS topology via `gmx dump`...")
    try:
        cmd = f"{gmx_path} dump -s {tpr_file}"
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:
            moltypes = {}
            molblocks = []
            current_moltype = None
            current_mb_moltype = None
            current_mb_nmol = None

            re_moltype_def = re.compile(r'moltype\s*[\(\[]\s*(\d+)\s*[\)\]]', re.IGNORECASE)
            re_mass = re.compile(r'm=\s*([0-9.eE+-]+)')
            re_molblock_def = re.compile(r'molblock\s*[\(\[]\s*(\d+)\s*[\)\]]', re.IGNORECASE)
            re_mb_type = re.compile(r'(?:molb_type|moltype|type)\s*=\s*(\d+)', re.IGNORECASE)
            re_mb_nmol = re.compile(r'(?:nmol|\#molecules)\s*=\s*(\d+)', re.IGNORECASE)

            for line in result.stdout.split('\n'):
                line = line.strip()
                match_mt = re_moltype_def.search(line)
                if match_mt:
                    current_moltype = int(match_mt.group(1))
                    if current_moltype not in moltypes: moltypes[current_moltype] = []
                    continue

                if current_moltype is not None and 'atom[' in line and 'm=' in line:
                    match_mass = re_mass.search(line)
                    if match_mass: moltypes[current_moltype].append(float(match_mass.group(1)))
                    continue

                match_mb = re_molblock_def.search(line)
                if match_mb:
                    current_moltype = None 
                    current_mb_moltype = None
                    current_mb_nmol = None
                    continue

                if current_moltype is None: 
                    match_type = re_mb_type.search(line)
                    if match_type: current_mb_moltype = int(match_type.group(1))
                    match_nmol = re_mb_nmol.search(line)
                    if match_nmol: current_mb_nmol = int(match_nmol.group(1))

                    if current_mb_moltype is not None and current_mb_nmol is not None:
                        molblocks.append((current_mb_moltype, current_mb_nmol))
                        current_mb_moltype = None
                        current_mb_nmol = None

            exact_masses = []
            for mb_type, mb_nmol in molblocks:
                if mb_type in moltypes: exact_masses.extend(moltypes[mb_type] * mb_nmol)

            if len(exact_masses) == len(u.atoms):
                u.atoms.masses = exact_masses
                print("[✔] Method 2: Successfully parsed and unrolled exact masses from `gmx dump`! (100% accurate)")
                return
    except Exception as e: pass

    print("[!] GMX dump parsing failed. Activating Smart Fallback...")
    is_martini_3 = any(atom.name.upper().startswith('T') for atom in u.atoms)
    if is_martini_3:
        for atom in u.atoms:
            name = atom.name.upper()
            if name.startswith('T'): atom.mass = 36.0
            elif name.startswith('S'): atom.mass = 54.0
            else: atom.mass = 72.0
    else:
        for atom in u.atoms:
            name = atom.name.upper()
            if name.startswith('S') or name.startswith('R'): atom.mass = 45.0
            else: atom.mass = 72.0
    print("[✔] Fallback masses assigned successfully!")

class OrientationalOrderAnalyzer:
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.u = mda.Universe(self.pdb_file)

    def calculate_p2(self, atom_name='BB'):
        res_names_unique = list(dict.fromkeys(self.u.residues.resnames))
        if len(res_names_unique) >= 2:
            RESIDUE_N = res_names_unique[0]
            RESIDUE_C = res_names_unique[1]
        else:
            RESIDUE_N = res_names_unique[0]
            RESIDUE_C = res_names_unique[0]
        
        res_N_group = self.u.select_atoms(f"resname {RESIDUE_N}").residues
        res_C_group = self.u.select_atoms(f"resname {RESIDUE_C}").residues
        N_mol_p2 = len(res_N_group)
        
        if N_mol_p2 == 0 or N_mol_p2 != len(res_C_group):
            raise ValueError("The number of N-term and C-term residues does not match.")

        all_mu_vectors = []
        for ts_p2 in self.u.trajectory:
            for i in range(N_mol_p2):
                atom_N = res_N_group[i].atoms.select_atoms(f"name {atom_name}")
                atom_C = res_C_group[i].atoms.select_atoms(f"name {atom_name}")
                if len(atom_N) >= 1 and len(atom_C) >= 1:
                    vec_mu = atom_C.positions[0] - atom_N.positions[0]
                    norm = np.linalg.norm(vec_mu)
                    if norm > 1e-6:
                        all_mu_vectors.append(vec_mu / norm)
                        
        mu_vectors_array = np.array(all_mu_vectors)
        N_total = len(mu_vectors_array)
        
        if N_total == 0: return None, None, 0.0, 0, RESIDUE_N, RESIDUE_C
            
        Q_sum = np.zeros((3, 3))
        I = np.identity(3)
        for mu in mu_vectors_array:
            mu_muT = np.outer(mu, mu)
            Q_sum += 0.5 * (3 * mu_muT - I)
            
        Q_matrix = Q_sum / N_total
        eigenvalues = np.linalg.eigvalsh(Q_matrix)
        P2_val = np.max(eigenvalues)

        return Q_matrix, eigenvalues, P2_val, N_total, RESIDUE_N, RESIDUE_C

class ShapeAnisotropyAnalyzer:
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.u = mda.Universe(self.pdb_file)

    def calculate_k2(self):
        atoms = self.u.atoms
        if len(atoms) == 0: return None, None, 0.0
        com = atoms.center_of_geometry()
        positions = atoms.positions - com
        N = len(atoms)
        S = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                S[i, j] = np.sum(positions[:, i] * positions[:, j]) / N
        eigenvalues = np.linalg.eigvalsh(S)
        l1, l2, l3 = eigenvalues[0], eigenvalues[1], eigenvalues[2]
        denominator = (l1 + l2 + l3)**2
        if denominator < 1e-8: k2 = 0.0
        else: k2 = 1.0 - 3.0 * ((l1*l2 + l2*l3 + l3*l1) / denominator)
        return S, eigenvalues, k2

class PDFAnalyzer:
    def __init__(self, universe, cutoff_space, n_pep):
        self.u = universe
        self.cutoff_space = cutoff_space
        self.n_pep = n_pep
        self.num_frames_analyzed = 0
        self.num_mols_a_for_pdf = 0
        self.is_same_group = False

    def interactive_group_selection(self):
        options = {
            "1": ("Whole Protein (All)", "protein"),
            "2": ("Main Chain (All BB)", "name BB and protein"),
            "3": ("Side Chain (All SC)", "not name BB and protein"), 
            "4": ("Water", "resname W")
        }
        protein_atoms = self.u.select_atoms("protein")
        res_name_pairs = list(dict.fromkeys(zip(protein_atoms.resnames, protein_atoms.names)))
        idx = 5
        for res, name in res_name_pairs:
            options[str(idx)] = (f"Atom {name} of Residue {res}", f"resname {res} and name {name}")
            idx += 1
            
        print("\n--- SELECT INTERACTION GROUPS FOR ANALYSIS ---")
        for key, (desc, sel) in options.items(): print(f"[{key}] {desc}")
            
        def process_input(user_input):
            parts = user_input.strip().split()
            if not parts: return None, None
            if all(p in options for p in parts):
                if len(parts) == 1: return options[parts[0]][0], options[parts[0]][1]
                else: return "Combined(" + " + ".join([options[p][0].replace("Atom ", "") for p in parts]) + ")", " or ".join([f"({options[p][1]})" for p in parts])
            else: return "Custom", user_input.strip()

        while True:
            choice_a = input("Select group A (Enter numbers, CAN TYPE MULTIPLE NUMBERS e.g., '7 8 9') or custom command: ")
            desc_a, sel_a = process_input(choice_a)
            if desc_a:
                sel_a_desc, sel_a_str = desc_a, sel_a
                break
            print("Invalid selection. Please try again.")

        while True:
            choice_b = input("Select group B (Enter numbers, CAN TYPE MULTIPLE NUMBERS e.g., '7 8 9') or custom command: ")
            desc_b, sel_b = process_input(choice_b)
            if desc_b:
                sel_b_desc, sel_b_str = desc_b, sel_b
                break
            print("Invalid selection. Please try again.")
            
        return sel_a_str, sel_b_str, sel_a_desc, sel_b_desc

    def calculate_metric(self, sel_a_str, sel_b_str, mode='contact', b=0.0, e=-1.0, dt=0.0):
        group_a = self.u.select_atoms(sel_a_str)
        group_b = self.u.select_atoms(sel_b_str)
        
        if len(group_a) == 0 or len(group_b) == 0:
            print("[-] Error: Atom selection does not exist in this structure!")
            return None, None
            
        num_mols = len(self.u.select_atoms("protein").residues) // self.n_pep
        if num_mols == 0: num_mols = 1
        
        # CẤP ID PHÂN TỬ ĐỘC LẬP: Tránh lỗi gộp nhầm Nước vào Peptide
        protein_resnames = np.unique(self.u.select_atoms("protein").resnames)
        
        is_prot_a = np.isin(group_a.resnames, protein_resnames)
        mol_ids_a = np.where(is_prot_a, group_a.resindices // self.n_pep, group_a.resindices + 1000000)
        
        is_prot_b = np.isin(group_b.resnames, protein_resnames)
        mol_ids_b = np.where(is_prot_b, group_b.resindices // self.n_pep, group_b.resindices + 1000000)
            
        times = []
        metric_values = []

        print(f"\n[+] Scanning trajectory to calculate strictly INTER-MOLECULAR {'Contact Number' if mode == 'contact' else 'Distance to Centroids'}...")
        self.u.trajectory.rewind() 
        
        mols_a_atoms = []
        mol_ids_a_unique = []
        for mol_id in np.unique(mol_ids_a):
            mol_mask = mol_ids_a == mol_id
            mols_a_atoms.append(group_a[mol_mask])
            mol_ids_a_unique.append(mol_id)
            
        mols_b_atoms = []
        mol_ids_b_unique = []
        for mol_id in np.unique(mol_ids_b):
            mol_mask = mol_ids_b == mol_id
            mols_b_atoms.append(group_b[mol_mask])
            mol_ids_b_unique.append(mol_id)
            
        mol_ids_a_arr = np.array(mol_ids_a_unique)
        mol_ids_b_arr = np.array(mol_ids_b_unique)

        self.num_mols_a_for_pdf = len(mol_ids_a_unique)
        self.is_same_group = (sel_a_str == sel_b_str)
        self.num_frames_analyzed = 0

        next_time = b
        for ts in self.u.trajectory:
            if ts is None: continue 
            if ts.time < b - 0.001: continue
            if e >= 0 and ts.time > e + 0.001: break
            if dt > 0 and ts.time < next_time - 0.001: continue

            times.append(ts.time)
            self.num_frames_analyzed += 1
            box_dims = self.u.dimensions[:3]
            
            if mode == 'contact':
                cutoff_A = self.cutoff_space * 10.0
                # VẬT LÝ BÀI BÁO: Đo Atom-Atom (Đếm toàn bộ hạt chạm nhau)
                pairs = capped_distance(group_a.positions, group_b.positions, max_cutoff=cutoff_A, box=self.u.dimensions, return_distances=False)
                if len(pairs) > 0:
                    mol_a_pairs = mol_ids_a[pairs[:, 0]]
                    mol_b_pairs = mol_ids_b[pairs[:, 1]]
                    
                    valid_mask = mol_a_pairs != mol_b_pairs
                    total_contacts = np.sum(valid_mask)
                    
                    if self.is_same_group:
                        total_contacts = total_contacts // 2
                else: 
                    total_contacts = 0
                    
                # CHUẨN HÓA MẬT ĐỘ TÁC GIẢ TANG: Chia cho TỔNG số lượng hạt của Group A (len(group_a))
                # Điều này trả về chính xác số lượt tiếp xúc trên MỖI HẠT của Group A
                if len(group_a) > 0:
                    metric_values.append(total_contacts / len(group_a))
                else:
                    metric_values.append(0.0)
                
            elif mode == 'distance':
                pos_a_centroids = np.array([m.center_of_geometry() for m in mols_a_atoms]) % box_dims
                pos_b_centroids = np.array([m.center_of_geometry() for m in mols_b_atoms]) % box_dims

                search_radius_A = 60.0
                pairs, dists = capped_distance(pos_a_centroids, pos_b_centroids, max_cutoff=search_radius_A, box=self.u.dimensions, return_distances=True)

                if len(pairs) > 0:
                    mol_a_pairs = mol_ids_a_arr[pairs[:, 0]]
                    mol_b_pairs = mol_ids_b_arr[pairs[:, 1]]

                    valid_mask = mol_a_pairs != mol_b_pairs

                    if self.is_same_group:
                        unique_mask = pairs[:, 0] < pairs[:, 1]
                        valid_mask = valid_mask & unique_mask

                    dist_nm = dists / 10.0  
                    range_mask = (dist_nm >= 0.4) & (dist_nm <= 6.0)
                    valid_mask = valid_mask & range_mask

                    valid_dists = dist_nm[valid_mask]
                    metric_values.extend(valid_dists)

            if dt > 0: next_time += dt

        return times, np.array(metric_values)

    def generate_pdf(self, metric_data, pdf_times, output_file, title_desc, mode='contact', x_range_val=None, bw=0.05):
        from scipy.stats import gaussian_kde
        import matplotlib.pyplot as plt

        if len(metric_data) == 0: return
        metric_data = metric_data[metric_data > 0]
        if len(metric_data) == 0:
            print("[-] No valid intermolecular data points found.")
            return

        print(f"\n[+] Calculating KDE (Probability per Bin) on {len(metric_data)} pairs with Bandwidth = {bw}...")
        
        if len(metric_data) > 100000:
            kde_data = np.random.choice(metric_data, 100000, replace=False)
        else:
            kde_data = metric_data
            
        if np.std(kde_data) == 0: kde_data = kde_data + np.random.normal(0, 1e-5, len(kde_data))
            
        kde = gaussian_kde(kde_data, bw_method=bw)
        x_label = "Distance (nm)" if mode == 'distance' else "Contact Number"
        
        if x_range_val is not None:
            x_min, x_max = x_range_val[0], x_range_val[1]
        else:
            if mode == 'contact':
                x_min, x_max = np.floor(metric_data.min()) - 0.5, np.ceil(metric_data.max()) + 0.5
            else:
                x_min, x_max = 0.4, 2.0
            
        x_axis = np.linspace(x_min, x_max, 1000)
        smooth_pdf_values = kde(x_axis)

        # Chuyển đổi thành Xác suất trên Bin chuẩn
        if mode == 'distance':
            bin_width = 0.05
            smooth_pdf_values = smooth_pdf_values * bin_width

        with open(output_file, 'w') as f:
            f.write(f"# {mode.capitalize()} & PDF Analysis - {title_desc}\n")
            f.write(f"# Column 1: {x_label}\n")
            f.write(f"# Column 2: Probability per Bin (Bin Width = 0.05 nm)\n")
            f.write(f"@    title \"PDF of {mode.capitalize()}\"\n")
            f.write(f"@    xaxis  label \"{x_label}\"\n")
            f.write(f"@    yaxis  label \"PDF\"\n")
            f.write(f"@TYPE xy\n")
            for x, y in zip(x_axis, smooth_pdf_values):
                f.write(f"{x:12.5f} {y:12.6f}\n")
        print(f"[✔] Exported XVG (X, Y) file for plotting to: {output_file}")
        
        while True:
            ans = input(f"\nDo you want to plot the bell curve (PDF) for {mode} and save the image? (yes/no): ").lower()
            if ans in ['yes', 'y']:
                plt.figure(figsize=(8, 6))
                plt.plot(x_axis, smooth_pdf_values, color='red' if mode=='distance' else 'blue', linewidth=2.5, label=f'PDF (BW={bw})')
                plt.title(f'PDF of {title_desc}')
                plt.xlabel(x_label)
                plt.ylabel('PDF')
                
                plt.xlim(x_min, x_max)
                plt.grid(axis='both', linestyle='--', alpha=0.6)
                plt.legend()
                
                img_file = output_file.replace('.xvg', '.png').replace('.txt', '.png')
                if not img_file.endswith('.png'): img_file += '.png'
                plt.savefig(img_file, dpi=300, bbox_inches='tight')
                print(f"[✔] Saved plot image to: {img_file}")
                plt.show()
                break
            elif ans in ['no', 'n']: break
            else: print("Please type 'yes' or 'no'.")

class FESAnalyzer:
    def __init__(self, universe, cutoff_space, n_pep):
        self.u = universe
        self.cutoff_space = cutoff_space
        self.n_pep = n_pep

    def calculate_fes_trajectory(self, sel_a_str, sel_b_str, xtc_file, gro_file, b=0.0, e=-1.0, dt=0.0, n_threads=4):
        group_a = self.u.select_atoms(sel_a_str)
        group_b = self.u.select_atoms(sel_b_str)
        
        if len(group_a) == 0 or len(group_b) == 0:
            print("[-] FES Error: Atom selection does not exist in the structure!")
            return None, None

        print(f"\n[+] Scanning trajectory to identify target frames (Multi-Processing Setup)...")
        self.u.trajectory.rewind()
        
        target_frames = []
        next_time = b
        
        for ts in self.u.trajectory:
            if ts is None: continue
            if ts.time < b - 0.001: continue
            if e >= 0 and ts.time > e + 0.001: break
            if dt > 0 and ts.time < next_time - 0.001: continue

            target_frames.append(ts.frame)
            if dt > 0: next_time += dt

        if not target_frames:
            print("[-] No valid frames found matching your time range and dt.")
            return None, None

        print(f"[i] Found {len(target_frames)} frames to process. Launching {n_threads} parallel processes...")

        chunks = np.array_split(target_frames, n_threads)
        tasks = []
        for chunk in chunks:
            if len(chunk) > 0:
                tasks.append((gro_file, xtc_file, sel_a_str, sel_b_str, list(chunk), self.n_pep))

        dist_list = []
        angle_list = []

        with multiprocessing.Pool(processes=n_threads) as pool:
            results = pool.map(_fes_worker_task, tasks)

        for d_res, a_res in results:
            dist_list.extend(d_res)
            angle_list.extend(a_res)

        return np.array(dist_list), np.array(angle_list)

    def generate_fes(self, d_data, angle_data, output_file, title_desc, T=300):
        from scipy.stats import gaussian_kde
        import matplotlib.pyplot as plt

        if len(d_data) == 0:
            print("[-] No data available to calculate FES.")
            return

        is_pipi = False
        if len(angle_data) > 0 and np.max(angle_data) <= 90.5:
            is_pipi = True
            
        min_d = 0.0 
        max_d = 0.8 if is_pipi else 1.2
        max_angle = 90.0 if is_pipi else 180.0
        cap_val = 5.0 if is_pipi else 6.0
        vmin_val = -5 if is_pipi else -6

        mask = (d_data >= min_d) & (d_data <= max_d)
        d_filtered = d_data[mask]
        angle_filtered = angle_data[mask]

        if len(d_filtered) < 10:
            print(f"[-] Too few data points satisfying distance <= {max_d} nm to construct FES.")
            return

        print(f"[+] Calculating 2D KDE on {len(d_filtered)} samples (Mode: {'Pi-Pi' if is_pipi else 'Mainchain Inter-molecular'})...")
        values = np.vstack([d_filtered, angle_filtered])
        
        if np.var(d_filtered) == 0:
            values[0, :] += np.random.normal(0, 1e-5, len(d_filtered))
        if np.var(angle_filtered) == 0:
            values[1, :] += np.random.normal(0, 1e-5, len(angle_filtered))
            
        kernel = gaussian_kde(values)

        d_grid, theta_grid = np.mgrid[min_d:max_d:100j, 0.0:max_angle:100j]
        positions = np.vstack([d_grid.ravel(), theta_grid.ravel()])
        
        P = np.reshape(kernel(positions).T, d_grid.shape)

        R = 1.987e-3  
        F = -R * T * np.log(P + 1e-15)  

        F = F - np.min(F)
        F = np.clip(F, 0, cap_val) 
        F = F - cap_val 

        with open(output_file, 'w') as f:
            f.write(f"# Free Energy Surface (FES) - {title_desc}\n")
            f.write(f"# Column 1: Distance (nm)\n")
            f.write(f"# Column 2: Angle (degree)\n")
            f.write(f"# Column 3: Free Energy (kcal/mol)\n")
            f.write(f"@    title \"Free Energy Surface\"\n")
            f.write(f"@    xaxis  label \"Distance (nm)\"\n")
            f.write(f"@    yaxis  label \"Angle (degree)\"\n")
            f.write(f"@TYPE xyz\n")
            
            for i in range(d_grid.shape[0]):
                for j in range(d_grid.shape[1]):
                    f.write(f"{d_grid[i, j]:12.5f} {theta_grid[i, j]:12.5f} {F[i, j]:12.6f}\n")
                f.write("\n")
                
        print(f"[✔] Successfully exported FES matrix file (3 columns) to: {output_file}")

        while True:
            ans = input("\nDo you want to plot the Contourf map (FES) and save the image? (yes/no): ").lower()
            if ans in ['yes', 'y']:
                plt.figure(figsize=(8, 6))
                
                contour = plt.contourf(d_grid, theta_grid, F, levels=60, cmap='jet_r', vmin=vmin_val, vmax=0)
                plt.colorbar(contour, label='Free Energy (kcal/mol)')

                plt.xlabel('Distance (nm)')
                plt.ylabel('Angle (degree)')
                plt.title(f'Free Energy Landscape: {title_desc}')
                plt.xlim(min_d, max_d) 

                if is_pipi:
                    plt.annotate('Parallel', xy=(0.45, 10), color='white', weight='bold')
                    plt.annotate('T-shaped', xy=(0.6, 75), color='white', weight='bold')

                img_file = output_file.replace('.xvg', '.png')
                plt.savefig(img_file, dpi=300, bbox_inches='tight')
                print(f"[✔] Saved FES plot image to: {img_file}")
                plt.show()
                break
            elif ans in ['no', 'n']:
                print("[i] Skipped image plot export step.")
                break
            else:
                print("Please type 'yes' or 'no'.")
