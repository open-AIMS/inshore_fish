RMD_FILES = $(wildcard *.Rmd)
HTML_FILES = $(patsubst %.Rmd, %.html, $(RMD_FILES))

all: $(HTML_FILES)

%.html: %.Rmd
	echo $<
	echo "library(rmarkdown); render(\"$<\", output_format=\"html_document\", output_options=list(self_contained=TRUE))" | R --no-save --no-restore
