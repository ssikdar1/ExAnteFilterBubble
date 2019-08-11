

sim:
	nohup python3 -u simulations.py > out_file.out

julia-install:
	julia -e 'import Pkg; Pkg.add("IterTools")'
	julia -e 'import Pkg; Pkg.add("Distances")'
	julia -e 'import Pkg; Pkg.add("Distributions")'

latex:
	/Library/TeX/texbin/pdflatex paper.tex  2>&1 > /dev/null
