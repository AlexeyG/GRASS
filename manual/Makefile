# Need some constants
DEST = manual.pdf
AUX = manual.aux manual.log manual.out manual.pdf manual.toc

.PHONY : clean

$(DEST) : manual.tex
	pdflatex manual.tex

clean :
	rm -f $(AUX)
