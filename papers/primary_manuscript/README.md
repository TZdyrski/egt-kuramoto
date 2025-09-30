# Compiling the paper
To compile the paper reproducibly, install Gnu Guix and run
```bash
guix time-machine -C channels.scm -- shell --pure -m manifest.scm -- latexmk Report.tex
```
Alternatively, to compile the paper without Guix, ensure the `texlive` packages listed in `manifest.scm` are installed
and run `latexmk Report.tex`.
