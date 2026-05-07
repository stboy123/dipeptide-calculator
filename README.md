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
