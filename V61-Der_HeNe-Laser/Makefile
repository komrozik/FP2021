all: build/main.pdf

.PRECIOUS: build/Abgabe.pdf
# hier Python-Skripte:
plots/leistung.pdf plots/pol.pdf plots/stabilitätsbedingung.pdf  plots/TEM00.pdf plots/TEM10.pdf: plot.py | build
	TEXINPUTS=$$(pwd): python plot.py

# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/main.pdf: plots/leistung.pdf plots/pol.pdf plots/stabilitätsbedingung.pdf plots/TEM00.pdf plots/TEM10.pdf

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
