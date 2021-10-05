all: build/main.pdf

.PRECIOUS: build/Abgabe.pdf
# hier Python-Skripte:
plots/bins.pdf plots/e-fkt.pdf plots/rate_20Hz.pdf plots/rate_20Hz_plateau.pdf plots/rate_20Hz_plateau2.pdf: programm.py | build
	TEXINPUTS=$$(pwd): python programm.py

# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/main.pdf: plots/bins.pdf plots/e-fkt.pdf plots/rate_20Hz.pdf plots/rate_20Hz_plateau.pdf plots/rate_20Hz_plateau2.pdf

build/main.pdf: FORCE | build
	  TEXINPUTS=build: \
	  BIBINPUTS=build: \
	  max_print_line=1048576 \
	latexmk -pvc\
	  --lualatex \
	  --output-directory=build \
	  --interaction=errorstopmode \
	main.tex


build:
	mkdir -p build
	mkdir -p plots

clean:
	rm -rf build
	rm -rf plots

FORCE:

.PHONY: all clean
