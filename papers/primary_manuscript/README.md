# Compiling the paper
To compile the paper reproducibly, install Gnu Guix and run
```bash
SOURCE_DATE_EPOCH=0 guix time-machine -C channels.scm -- shell -m manifest.scm -- latexmk -auxdir=latex.aux -lualatex Report.tex
```
