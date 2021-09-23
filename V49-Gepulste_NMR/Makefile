all: build/main.pdf

.PRECIOUS: build/Abgabe.pdf
# hier Python-Skripte:
plots/df1.pdf plots/df54.pdf plots/df52.pdf plots/df2.pdf plots/df2_2.pdf plots/df58_59.pdf plots/fourier.pdf : programm.py | build
	TEXINPUTS=$$(pwd): python programm.py

# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/main.pdf: plots/df1.pdf plots/df54.pdf plots/df52.pdf plots/df2.pdf plots/df2_2.pdf plots/df58_59.pdf plots/fourier.pdf

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
