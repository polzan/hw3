HW_NAME:=hw3

all: clean $(HW_NAME).pdf $(HW_NAME)_code.pdf

$(HW_NAME).pdf: $(HW_NAME)_code.pdf

%.pdf: %.tex
	pdflatex $<
	pdflatex $<

clean:
	rm $(HW_NAME)*.pdf || true
	rm *.aux || true	
	rm *.out || true
	rm matlab/*.pdf || true

.PHONY: clean all
