../user-manual.pdf:
	latex user-manual
	bibtex user-manual
	latex user-manual
	latex user-manual
	dvipdf user-manual.dvi
	mv user-manual.pdf ..

dist-hook: ../user-manual.pdf

EXTRA_DIST=../user-manual.pdf