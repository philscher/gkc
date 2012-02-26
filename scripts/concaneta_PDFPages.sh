#!/bin/bash
cat <<ENDTXT >__tmp__.tex
\documentclass[landscape]{article}
\usepackage{pdfpages}
\begin{document}
ENDTXT
outnm=$1
count=`echo "$# - 1" |bc`
for f in `seq $count`; do
  shift
  echo "\includepdf[pages=-]{$1}" >>__tmp__.tex
done
echo "\end{document}" >>__tmp__.tex
pdflatex __tmp__.tex
rm __tmp__.{tex,aux,log}
mv __tmp__.pdf $outnm
