GRAPHICS_ASY=$(wildcard graphics/*.asy)
GRAPHICS_PDF=$(patsubst graphics/%.asy,graphics/%.pdf,$(GRAPHICS_ASY))
VIEW=evince

graphics/%.pdf : graphics/%.asy
	asy -f pdf -o $@ $<

clean :
	/bin/rm -rf main.pdf main.log $$(find . -name '*~' -o -name '*#' -o -name '.#*')

zip : ../luck.zip


.PHONY: ../luck.zip

../luck.zip :
	/bin/rm -rf ../luck.zip
	cd ..; zip luck.zip $$(find luck -path luck/.git -prune -o -type f)

all : graphs application.tex chi2.tex computation.tex conclusion.tex introduction.tex luck.tex main.tex model.tex multinomial.tex normal.tex proofs.tex randomness.tex stirling.tex summary.tex
	if [ -f main.idx ] ; then makeindex main; fi
	pdflatex main


graphs : $(GRAPHICS_PDF)

pdf : all
	$(VIEW) main.pdf &

qqwing : qqwing.cpp
	$(CXX) -g -o $@ $<


slides : graphs slides.tex
	pdflatex slides

