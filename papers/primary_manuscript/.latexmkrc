$aux_dir = './latex.aux/';
# Use lualatex
$pdf_mode = 4;
# Disable postscript and dvi output
$postscript_mode = $dvi_mode = 0;
# Disable non-reproducible metadata (creation date, ModDate, and ID)
&alt_tex_cmds;
$pre_tex_code = '\pdfvariable suppressoptionalinfo \numexpr32+64+512\relax \
                 \pdfvariable objcompresslevel=0\relax';
