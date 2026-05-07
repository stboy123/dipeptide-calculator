import argparse
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
from scipy.spatial import cKDTree
import os
import warnings
import sys
import numpy as np
import subprocess
import re

try:
    from algorithms import GromacsProcessor, ClusteringAnalyzer, LiquidityAnalyzer, DensityAnalyzer
except ImportError:
    print("[-] Lỗi: Không tìm thấy file algorithms.py trong cùng thư mục.")
    sys.exit(1)

warnings.filterwarnings('ignore')

def parse_arguments():
    parser = argparse.ArgumentParser(description="Dipeptide Assembly Analyzer (Clone 100% Yiming Tang's Code)", add_help=True)
    
    parser.add_argument('-f', type=str, required=True, help='File quỹ đạo (.xtc)')
    parser.add_argument('-s', type=str, required=True, help='File cấu trúc cho GROMACS (.tpr)')
    parser.add_argument('-gro', type=str, required=True, help='File cấu trúc cho Python (.gro)')
    parser.add_argument('-select', type=str, required=True, help='Lệnh chọn phân tử (VD: "resname GLY PHE GLN TRP")')
    
    parser.add_argument('-nb', type=str, default=None, help='Xuất số lượng cụm (number.xvg)')
    parser.add_argument('-sz', type=str, default=None, help='Xuất kích thước cụm lớn nhất (size.xvg)')
    parser.add_argument('-comp', type=str, default=None, help='Xuất tỷ lệ Co-assembly (composition.xvg)')
    parser.add_argument('-molnumber', type=str, default=None, help='Xuất tổng số phân tử trong cụm hợp lệ')
    parser.add_argument('-liquidity', type=str, default=None, help='Xuất độ lỏng (liquidity.xvg)')
    parser.add_argument('-density', type=str, default=None, help='Xuất mật độ (density.xvg)')
    parser.add_argument('-pdb', type=str, default=None, help='Xuất cấu trúc cụm lớn nhất (.pdb)')
    
    parser.add_argument('-cutoff_space', type=float, default=0.7, help='Khoảng cách tương tác (nm)')
    parser.add_argument('-cutoff_multi', type=float, default=2.0, help='Hệ số nhân bán kính dò mật độ (Chuẩn của Tang = 2.0)')
    parser.add_argument('-cutoff_cz', type=int, default=20, help='Số PHÂN TỬ tối thiểu để tính là 1 cụm hợp lệ')
    parser.add_argument('-n_pep', type=int, default=2, help='Số GỐC trong 1 phân tử (Mặc định 2)')
    parser.add_argument('-pdb_time', type=float, default=-1, help='Thời điểm trích xuất file PDB (ps)')
    
    parser.add_argument('--fix-pbc', action='store_true', help='Tự động vá lỗi xé hộp (PBC whole)')
    parser.add_argument('--gmx-path', type=str, default='gmx', help='Đường dẫn đến GROMACS')
    
    return parser.parse_args()

def write_xvg(filename, title, y_label, time_data, value_data, legends=None, is_liquidity=False):
    if filename is None: return
    with open(filename, 'w') as f:
        f.write('# File tạo bởi công cụ phân tích của Lộc (100% Tang match)\n')
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
    print(f"[*] Đã ghi xong: {filename}")

def sync_masses_from_tpr(u, tpr_file, gmx_path):
    print(f"\n[+] Đang đồng bộ khối lượng từ {tpr_file} để đảm bảo sai số 0%...")
    masses = []
    try:
        cmd = f"{gmx_path} dump -s {tpr_file}"
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        in_atoms = False
        for line in result.stdout.split('\n'):
            if 'atom (' in line and '):' in line:
                in_atoms = True
            elif in_atoms and line.strip().startswith('type ('):
                in_atoms = False
                
            if in_atoms and 'm=' in line and 'atom[' in line:
                match = re.search(r'm=\s*([0-9.eE+-]+)', line)
                if match:
                    masses.append(float(match.group(1)))
    except Exception as e:
        print(f"[-] Lỗi chạy gmx dump: {e}")
        
    if len(masses) == len(u.atoms):
        for atom, m in zip(u.atoms, masses):
            atom.mass = m
        print("[✔] Đã nạp thành công khối lượng gốc từ GROMACS. Density sẽ khớp 100%!")
    else:
        print(f"[-] Cảnh báo: TPR có {len(masses)} atoms, GRO có {len(u.atoms)} atoms.")
        print("[i] Đang sử dụng khối lượng mô phỏng Martini 2.0 mặc định thay thế...")
        for atom in u.atoms:
            if any(atom.name.upper().startswith(prefix) for prefix in ['S', 'R']):
                atom.mass = 45.0
            else:
                atom.mass = 72.0

def main():
    args = parse_arguments()
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

    print(f"\n[+] Đang đọc dữ liệu: {args.gro} + {final_xtc}")
    try:
        u = mda.Universe(args.gro, final_xtc)
        
        # Gọi hàm ép khối lượng GROMACS vào Python
        sync_masses_from_tpr(u, args.s, args.gmx_path)
        
        peptide_group = u.select_atoms(args.select)
        residues = peptide_group.residues
        res_names = residues.names
        num_residues = len(residues)
        resindex_to_idx = {res.resindex: i for i, res in enumerate(residues)}
        
        if num_residues == 0:
            print("[-] Lỗi: Không tìm thấy phân tử nào.")
            return
            
        num_mols = num_residues // args.n_pep
        print(f"[i] Cấu trúc: {args.n_pep} gốc/phân tử | Tổng số phân tử: {num_mols}")

    except Exception as e:
        print(f"[-] Lỗi nạp dữ liệu: {e}")
        return

    cluster_calc = ClusteringAnalyzer(args.cutoff_space, args.n_pep)
    
    times, cluster_counts, max_sizes = [], [], []
    comp_data, molnumber_data, liquidity_data, density_data = [], [], [], []
    
    prev_cluster = set()
    prev_contacts = set()

    print("\n[+] Bắt đầu quét quỹ đạo (Quá trình tính Density sẽ khá chậm để mô phỏng C++)...")
    for ts in u.trajectory:
        times.append(ts.time)
        
        all_clusters, _ = cluster_calc.calculate(
            peptide_group, num_residues, resindex_to_idx, u.dimensions
        )
        
        all_clusters_mols = []
        for c in all_clusters:
            mols = set(res_idx // args.n_pep for res_idx in c)
            all_clusters_mols.append(mols)

        qualified_clusters_mols = [mols for mols in all_clusters_mols if len(mols) >= args.cutoff_cz]

        # 1. Đếm cụm và Tính tổng phân tử (Bị lọc bởi cutoff)
        cluster_counts.append(len(qualified_clusters_mols))
        molnumber_data.append(sum(len(mols) for mols in qualified_clusters_mols))

        # 2. Co-assembly
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

        # 3. Kích thước lớn nhất
        current_max_size = max([len(mols) for mols in all_clusters_mols]) if all_clusters_mols else 0
        max_sizes.append(current_max_size)

        # BẮT CHƯỚC LỖI TANG: Biến buggy_cluster_this_frame bỏ rơi phân tử đầu tiên của mọi cụm
        buggy_cluster_this_frame = set()
        for mols in all_clusters_mols:
            mols_list = sorted(list(mols))
            if len(mols_list) > 1:
                for m in mols_list[1:]:
                    buggy_cluster_this_frame.add(m)

        # 4. ĐỘ LỎNG & ĐIỂM TIẾP XÚC
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
                pairs = capped_distance(
                    agg_atoms.positions, agg_atoms.positions,
                    max_cutoff=cutoff_A, box=u.dimensions, return_distances=False
                )
                
                if len(pairs) > 0:
                    mol_ids = agg_atoms.resindices // args.n_pep
                    mol_pairs = mol_ids[pairs]
                    valid_pairs = np.sort(mol_pairs, axis=1)
                    unique_pairs = np.unique(valid_pairs, axis=0)
                    current_contacts = set(map(tuple, unique_pairs))

            if ts.frame == 0:
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

        # 5. TÍNH TOÁN MẬT ĐỘ (BẢN CLONE HOÀN HẢO TỪ C++)
        if args.density:
            # Bước 1: Trích xuất các index của tập hợp lỗi
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
            
            # Bước 2: Tạo lưới y chang C++ (Mắt lưới ÉP BUỘC là 1.0 Angstrom = 0.1 nm)
            box_dims = u.dimensions[:3]
            x_grid = np.arange(0, box_dims[0], 1.0)
            y_grid = np.arange(0, box_dims[1], 1.0)
            z_grid = np.arange(0, box_dims[2], 1.0)
            
            X, Y, Z = np.meshgrid(x_grid, y_grid, z_grid, indexing='ij')
            grid_points = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T
            total_grid_points = len(grid_points)
            
            radius_multi_A = args.cutoff_space * 10.0 * args.cutoff_multi
            
            # Bước 3: Chỉ quét nhóm Ngưng tụ 
            if len(aggr_atoms) > 0:
                wrapped_aggr = aggr_atoms.positions % box_dims
                tree_aggr = cKDTree(wrapped_aggr, boxsize=box_dims)
                distances_aggr, _ = tree_aggr.query(grid_points, k=1, distance_upper_bound=radius_multi_A)
                num_aggr_points = np.sum(distances_aggr <= radius_multi_A)
            else:
                num_aggr_points = 0
                
            # Bước 4: Thể tích phân tán = Phần còn lại của không gian (Logic "Else" của Tang)
            num_disp_points = total_grid_points - num_aggr_points
                
            vol_aggr_nm3 = num_aggr_points / 1000.0
            vol_disp_nm3 = num_disp_points / 1000.0
            
            den_aggr = (mass_aggr / vol_aggr_nm3 / 0.602) if vol_aggr_nm3 > 0 else 0.0
            den_disp = (mass_disp / vol_disp_nm3 / 0.602) if vol_disp_nm3 > 0 else 0.0
            
            density_data.append([den_aggr, den_disp])

        if ts.frame % 100 == 0:
            print(f" Time: {ts.time:10.1f} ps | Cụm Hợp Lệ: {len(qualified_clusters_mols):3d} | MaxSize: {current_max_size:3d}")

    print("\n[+] Đang xuất file kết quả...")
    write_xvg(args.nb, "Number of Clusters", "Count", times, cluster_counts)
    write_xvg(args.sz, "Size of Largest Cluster", "Molecules", times, max_sizes)
    write_xvg(args.molnumber, "Total Molecules in Clusters", "Molecules", times, molnumber_data)
    
    if args.comp:
        write_xvg(args.comp, "Cluster Composition", "Value", times, comp_data, ["GF_Count", "QW_Count", "GF_Percent", "QW_Percent"])

    if args.liquidity:
        legs = [
            "Fraction of preservation",
            "Fraction of cluster growth",
            "Fraction of cluster shrink",
            "Fraction of aggregation",
            "Preservation of contacts"
        ]
        write_xvg(args.liquidity, "Fraction of preservation", "Fraction", times, liquidity_data, legs, is_liquidity=True)
        
    if args.density:
        write_xvg(args.density, "Density of Phases", "Density (g/cm^3)", times, density_data, ["Aggregated Phase", "Dispersed Phase"])

    for f in temp_files:
        if os.path.exists(f): os.remove(f)
    print("\n[✔] Mọi công việc đã xong!")

if __name__ == "__main__":
    main()
