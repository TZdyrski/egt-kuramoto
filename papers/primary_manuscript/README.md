# Compiling the paper
To compile the paper reproducibly, install Gnu Guix and run
```bash
guix time-machine -C channels.scm -- shell --pure -m manifest.scm -- latexmk Report.tex
```
Alternatively, to compile the paper without Guix, ensure the `texlive` packages listed in `manifest.scm` are installed
and run `latexmk Report.tex`.

# Compiling figures
To compile the figures reproducibly, install Gnu Guix and run
```bash
guix time-machine -C channels.scm -- shell --pure -m manifest.scm -- latexmk -cd tikz/<figure name>.tex
```
Alternatively, to compile the paper without Guix, ensure the `texlive` packages listed in `manifest.scm` are installed
and run `latexmk -cd tikz/<figure name>.tex`.

# Arxiv submission
To generate an arxiv submission, compile all the figures
```bash
guix time-machine -C channels.scm -- shell --pure -m manifest.scm -- latexmk -cd tikz/c-elegans.tex tikz/chimera-states.tex tikz/game-payoffs.tex tikz/game-types.tex tikz/model-setup.tex tikz/phase-diagram.tex tikz/well-mixed.tex
```
Change the standalone package mode to `image` to use the generated figure pdfs:
```bash
sed -i 's/usepackage{standalone}/usepackage[mode=image]{standalone}/' preamble.sty
```
Compile the manuscript and supplemental information
```bash
guix time-machine -C channels.scm -- shell --pure -m manifest.scm -- latexmk Report.tex SupplementaryInformation.tex
```
Move the .bbl files
```bash
mv latex.aux/Report.bbl latex.aux/SupplementaryInformation.bbl .
```
Create a zip of files for arxiv
```bash
zip ArxivUpload.zip Report.tex SupplementaryInformation.tex preamble.sty custom-definitions.tex sn-*  Report.bbl SupplementaryInformation.bbl tikz/preamble.tex tikz/*.pdf
```
