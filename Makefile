LATEXMK ?= latexmk
LUALATEX_FLAGS ?= -cd -lualatex -interaction=nonstopmode -halt-on-error -file-line-error -shell-escape

MAIN_TEX := main.tex
WILD_TEX := wild/main.tex

GRAPHICS_ASY := $(wildcard graphics/*.asy)
GRAPHICS_PDF := $(patsubst graphics/%.asy,graphics/%.pdf,$(GRAPHICS_ASY))

.PHONY: all main wild graphs clean test test-python test-scilab test-r

all: main wild

main: main.pdf

wild: wild/main.pdf

graphs: $(GRAPHICS_PDF)

graphics/%.pdf: graphics/%.asy
	asy -f pdf -o $@ $<

main.pdf: graphs

%.pdf: %.tex
	$(LATEXMK) $(LUALATEX_FLAGS) $<
	@'$(CURDIR)/tools/check-log.sh' '$(<:.tex=.log)' '$(<:.tex=.pdf)'

clean:
	$(LATEXMK) -C $(MAIN_TEX) || true
	$(LATEXMK) -C $(WILD_TEX) || true
	rm -f $(GRAPHICS_PDF)

# --- code & tests ---

test: test-python test-scilab test-r

test-python:
	cd python && python3 test_luck_ties.py
	cd python && python3 test_luck_extremes.py

test-scilab:
	@echo "scilab tests are interactive functions; run them by exec()'ing"
	@echo "the test_*.sci files in a scilab session."

test-r:
	cd R/tests/testthat && Rscript -e 'testthat::test_dir(".", stop_on_failure=TRUE, stop_on_warning=TRUE)'
