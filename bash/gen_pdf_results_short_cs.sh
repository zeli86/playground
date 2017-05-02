#!/bin/bash

###############################################################################
###############################################################################

FILES5=$(find . -name "sol.png")

FILES5=`echo ${FILES5//"./"/}`
FILES5=`echo ${FILES5//"/sol.png"/}`
FILES5=`echo ${FILES5//"_"/":"}`
FILES5=`echo ${FILES5//"/"/":"}`
#FILES5=`echo ${FILES5//","/"."}`
echo $FILES5 > delme.txt

tr -s '\040' '\012' <delme.txt >delme1.txt
sort delme1.txt > delme.txt

readarray sorted < delme.txt

counter=0
divisor=30
let divisorminuseins=divisor-1

rm overview_short.tex
touch overview_short.tex

cwd=$(pwd)
cwd=${cwd//"_"/" "}

echo "\documentclass[reprint,amsmath,amssymb,aps,floatfix]{revtex4}" >> overview_short.tex
echo "\usepackage[utf8x]{inputenc}" >> overview_short.tex
echo "\usepackage{fancyhdr}" >> overview_short.tex
echo "\usepackage{graphicx}" >> overview_short.tex
echo "\usepackage{subfigure}" >> overview_short.tex
echo "\renewcommand*{\thesubfigure}{}" >> overview_short.tex
echo "\pagestyle{fancy}" >> overview_short.tex
echo "\fancyhf{}"  >> overview_short.tex
echo "\fancyhead[C]{$cwd}"  >> overview_short.tex
echo "\fancyfoot[C]{\thepage}"  >> overview_short.tex

echo "\begin{document}" >> overview_short.tex

for i in "${sorted[@]}"
do 
  tmp1=`echo ${i/":"/"_"}`
  #tmp1=`echo ${tmp1/":"/"_"}` #for x_x_x folder names
  tmp1=`echo ${tmp1/":"/"/"}`

  if [[ "$tmp1" =~ "0000" ]] || [[ "$tmp1" =~ "0074" ]]
  then
    echo $tmp1
  else
    continue
  fi

  remainder=`echo "${counter}%${divisor}" | bc`
  #echo "Remainder: $remainder"
   
  if [[ "$remainder" == "0" ]]
  then  
    echo "\begin{figure}" >> overview_short.tex
  fi

  echo "\subfigure[$i]{\includegraphics[width=2.5cm]{$tmp1/sol.png}}" >> overview_short.tex
   
  if [[ "$remainder" == "$divisorminuseins" ]]
  then  
    echo "\end{figure}" >> overview_short.tex
    echo "\setcounter{subfigure}{0}" >> overview_short.tex
  fi 
 
  let counter=counter+1
done
 
if [[ "$remainder" != "$divisorminuseins" ]]
then  
  echo "\end{figure}" >> overview_short.tex
fi

echo "\end{document}" >> overview_short.tex

pdflatex -interaction nonstopmode overview_short.tex 
 
