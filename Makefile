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

all : graphs main.tex introduction.tex normal.tex multinomial.tex computation.tex
	pdflatex main
#	latexmk -pdf main

graphs : $(GRAPHICS_PDF)

pdf : all
	$(VIEW) main.pdf &

qqwing : qqwing.cpp
	$(CXX) -g -o $@ $<
