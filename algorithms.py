import subprocess
import shutil
import networkx as nx
import numpy as np
from MDAnalysis.analysis import distances
from scipy.spatial import cKDTree

class GromacsProcessor:
    """Chuyên gia điều khiển GROMACS tự động tìm đường dẫn"""
    
    def __init__(self, custom_path=None):
        self.gmx_path = custom_path or shutil.which("gmx")
        
        if not self.gmx_path:
            raise FileNotFoundError("Không tìm thấy GROMACS trên hệ thống. Hãy đảm bảo bạn đã source GMXRC.")
        else:
            print(f"[i] Đã tìm thấy GROMACS tại: {self.gmx_path}")

    def make_whole(self, s_file, f_file, out_file):
        """Tự động gọi lệnh trjconv để vá hộp mô phỏng"""
        print(f"[+] Đang vá lỗi hộp tuần hoàn (PBC) cho {f_file}...")
        cmd = [self.gmx_path, "trjconv", "-s", s_file, "-f", f_file, "-o", out_file, "-pbc", "whole"]
        
        try:
            subprocess.run(cmd, input="0\n", text=True, capture_output=True, check=True)
            print(f"[*] Đã vá thành công: {out_file}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"[-] Lỗi vá hộp (PBC): {e.stderr}")
            return False

class ClusteringAnalyzer:
    """Class chuyên xử lý Phân cụm đồ thị"""
    def __init__(self, cutoff_nm, cutoff_cz):
        self.cutoff_A = cutoff_nm * 10.0
        self.cutoff_cz = cutoff_cz

    def calculate(self, peptide_group, num_mols, resindex_to_idx, box_dimensions):
        dist_matrix = distances.contact_matrix(peptide_group.positions, cutoff=self.cutoff_A, box=box_dimensions)
        G = nx.Graph()
        G.add_nodes_from(range(num_mols))
        
        rows, cols = dist_matrix.nonzero()
        for r, c in zip(rows, cols):
            res_i = peptide_group[r].resindex
            res_j = peptide_group[c].resindex
            if res_i != res_j:
                G.add_edge(resindex_to_idx[res_i], resindex_to_idx[res_j])

        all_clusters = list(nx.connected_components(G))
        valid_clusters = [c for c in all_clusters if len(c) >= self.cutoff_cz]
        
        cluster_this_frame = set()
        for c in valid_clusters:
            cluster_this_frame.update(c)
            
        return valid_clusters, cluster_this_frame

class LiquidityAnalyzer:
    """Class chuyên tính toán Độ lỏng"""
    def __init__(self, total_mols):
        self.total_mols = total_mols
        self.cluster_last_frame = set()

    def calculate(self, cluster_this_frame):
        if not self.cluster_last_frame:
            frac_agg = len(cluster_this_frame) / self.total_mols if self.total_mols > 0 else 0
            result = [0.0, frac_agg, 0.0, frac_agg]
        else:
            intersection = cluster_this_frame.intersection(self.cluster_last_frame)
            preservation = len(intersection) / len(self.cluster_last_frame) if self.cluster_last_frame else 0
            growth = (len(cluster_this_frame) - len(intersection)) / self.total_mols
            shrink = (len(self.cluster_last_frame) - len(intersection)) / self.total_mols
            fraction_agg = len(cluster_this_frame) / self.total_mols
            result = [preservation, growth, shrink, fraction_agg]
            
        self.cluster_last_frame = cluster_this_frame
        return result

class DensityAnalyzer:
    """Class chuyên dò lưới 3D để tính Mật độ"""
    def __init__(self, cutoff_nm, cutoff_multi):
        self.probe_radius = cutoff_nm * 10.0 * cutoff_multi

    def calculate(self, valid_clusters, molecules, peptide_group, box_dimensions):
        if not valid_clusters:
            return [0.0, 0.0]
            
        ag_idx = []
        for c in valid_clusters:
            ag_idx.extend(list(c))
            
        ag_atoms = molecules[ag_idx].atoms
        disp_atoms = peptide_group - ag_atoms
        
        ag_mass = ag_atoms.masses.sum()
        disp_mass = disp_atoms.masses.sum()
        
        grid_step = 1.0
        x = np.arange(0, box_dimensions[0], grid_step)
        y = np.arange(0, box_dimensions[1], grid_step)
        z = np.arange(0, box_dimensions[2], grid_step)
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        grid_points = np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]).T
        
        tree = cKDTree(ag_atoms.positions, boxsize=box_dimensions[:3])
        dists, _ = tree.query(grid_points, k=1, distance_upper_bound=self.probe_radius)
        
        ag_grid_count = np.sum(dists <= self.probe_radius)
        disp_grid_count = len(grid_points) - ag_grid_count
        
        ag_volume = ag_grid_count / 1000.0
        disp_volume = disp_grid_count / 1000.0
        
        den_ag = (ag_mass / ag_volume / 0.602) if ag_volume > 0 else 0
        den_disp = (disp_mass / disp_volume / 0.602) if disp_volume > 0 else 0
        
        return [den_ag, den_disp]
