# Dipeptide Assembly Analyzer
*Clone 100% Yiming Tang's Code*

Công cụ phân tích sự tự lắp ghép của Dipeptide. Để chạy công cụ, bạn cần cung cấp file quỹ đạo, cấu trúc hệ thống và kích hoạt các cờ tính toán mong muốn.

## 🚀 Hướng dẫn sử dụng nhanh

Cấu trúc lệnh tổng quát:
```bash
python main_nice.py -f <trajectory.xtc> -s <topology.tpr> -gro <structure.gro> -select "<nhóm phân tử>" -cutoff_space 0.7 -cutoff_cz 20 [CÁC CỜ PHÂN TÍCH]
```

### Các ví dụ phân tích cụ thể:

1. **Tính số lượng cụm (Number of clusters):**
   ```bash
   python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -nb number.xvg --fix-pbc
   ```

2. **Kích thước cụm lớn nhất theo thời gian (Size of largest cluster):**
   ```bash
   python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -sz size.xvg --fix-pbc
   ```

3. **Số lượng phân tử nằm trong các cụm:**
   ```bash
   python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -molnumber molecules_in_clusters.xvg --fix-pbc
   ```

4. **Tính độ lỏng (Liquidity):**
   ```bash
   python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -liquidity liquidity_factor.xvg --fix-pbc
   ```

5. **Tính mật độ (Density):**
   ```bash
   python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -density density.xvg --fix-pbc
   ```
6. **Xuất file cấu trúc của cụm lớn nhất (pdb):
   ```bash
   python main1.py -f output_reduced.xtc -s md_0_1_515.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -pdb_system system_python.pdb -pdb_time 2500000
   ```
---

## 🛠 Giải thích các tham số (Arguments)

Sử dụng lệnh `python main_nice.py -h` để xem trợ giúp chi tiết.


| Tham số | Mô tả |
| :--- | :--- |
| `-f` | File quỹ đạo (`.xtc`) |
| `-s` | File cấu trúc cho GROMACS (`.tpr`) |
| `-gro` | File cấu trúc định dạng `.gro` |
| `-select` | Lệnh chọn phân tử (Ví dụ: `"resname GLY PHE"`) |
| `-nb` | Xuất số lượng cụm ra file `.xvg` |
| `-sz` | Xuất kích thước cụm lớn nhất ra file `.xvg` |
| `-comp` | Xuất tỷ lệ Co-assembly (`composition.xvg`) |
| `-molnumber` | Xuất tổng số phân tử nằm trong các cụm hợp lệ |
| `-liquidity` | Xuất chỉ số độ lỏng |
| `-density` | Xuất dữ liệu mật độ |
| `-pdb` | Xuất file cấu trúc của cụm lớn nhất (`.pdb`) |
| `-pdb_system ` | Xuất file cấu trúc của cụm phân tử với điều kiện cutoff_cz (`.pdb`) |
| `-cutoff_space` | Khoảng cách tương tác (nm). Mặc định thường là 0.7 |
| `-cutoff_multi` | Hệ số nhân bán kính dò mật độ (Chuẩn của Tang = 2.0) |
| `-cutoff_cz` | Số phân tử tối thiểu để được tính là 1 cụm hợp lệ |
| `-n_pep` | Số gốc (residues) trong 1 phân tử (Mặc định: 2) |
| `-pdb_time` | Thời điểm muốn trích xuất file PDB (ps) |
| `--fix-pbc` | Tự động xử lý lỗi xé hộp (Periodic Boundary Conditions) |

---

