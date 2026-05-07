Để chạy công cụ, bạn cần cung cấp file quỹ đạo, cấu trúc hệ thống và kích hoạt các cờ tính toán mong muốn:
python main1.py -f <trajectory.xtc> -s <topology.tpr> -gro <structure.gro> -select "<nhóm phân tử>" -cutoff_space 0.7 -cutoff_cz 20 [CÁC CỜ PHÂN TÍCH]

number of clusters
python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -nb number.xvg --fix-pbc

Size of the largest cluster as a function of time
python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -sz size.xvg --fix-pbc

Number of molecules that are included in clusters
python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -molnumber molecules_in_clusters.xvg --fix-pbc

liquidity
python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -liquidity liquidity_factor.xvg --fix-pbc

density
python main_nice.py -f output_reduced.xtc -s md_0_1.tpr -gro md_production.gro -select "resname GLY PHE" -cutoff_space 0.7 -cutoff_cz 20 -density density.xvg --fix-pbc

================================================= python main_nice.py -h ====================================================================================================
usage: main_nice.py [-h] -f F -s S -gro GRO -select SELECT [-nb NB] [-sz SZ]
                    [-comp COMP] [-molnumber MOLNUMBER] [-liquidity LIQUIDITY]
                    [-density DENSITY] [-pdb PDB] [-cutoff_space CUTOFF_SPACE]
                    [-cutoff_multi CUTOFF_MULTI] [-cutoff_cz CUTOFF_CZ]
                    [-n_pep N_PEP] [-pdb_time PDB_TIME] [--fix-pbc]
                    [--gmx-path GMX_PATH]

Dipeptide Assembly Analyzer (Clone 100% Yiming Tang's Code)

optional arguments:
  -h, --help            show this help message and exit
  -f F                  File quỹ đạo (.xtc)
  -s S                  File cấu trúc cho GROMACS (.tpr)
  -gro GRO              File cấu trúc cho Python (.gro)
  -select SELECT        Lệnh chọn phân tử (VD: "resname GLY PHE GLN TRP")
  -nb NB                Xuất số lượng cụm (number.xvg)
  -sz SZ                Xuất kích thước cụm lớn nhất (size.xvg)
  -comp COMP            Xuất tỷ lệ Co-assembly (composition.xvg)
  -molnumber MOLNUMBER  Xuất tổng số phân tử trong cụm hợp lệ
  -liquidity LIQUIDITY  Xuất độ lỏng (liquidity.xvg)
  -density DENSITY      Xuất mật độ (density.xvg)
  -pdb PDB              Xuất cấu trúc cụm lớn nhất (.pdb)
  -cutoff_space CUTOFF_SPACE
                        Khoảng cách tương tác (nm)
  -cutoff_multi CUTOFF_MULTI
                        Hệ số nhân bán kính dò mật độ (Chuẩn của Tang = 2.0)
  -cutoff_cz CUTOFF_CZ  Số PHÂN TỬ tối thiểu để tính là 1 cụm hợp lệ
  -n_pep N_PEP          Số GỐC trong 1 phân tử (Mặc định 2)
  -pdb_time PDB_TIME    Thời điểm trích xuất file PDB (ps)
  --fix-pbc             Tự động vá lỗi xé hộp (PBC whole)

