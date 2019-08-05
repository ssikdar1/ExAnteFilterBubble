
sim:
	nohup python3 -u simulations.py > out_file.out
latex:
	/Library/TeX/texbin/pdflatex paper.tex  2>&1 > /dev/null
