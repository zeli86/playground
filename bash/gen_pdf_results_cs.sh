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

rm overview.tex
touch overview.tex

echo "\documentclass[reprint,amsmath,amssymb,aps,floatfix]{revtex4}" >> overview.tex
echo "\usepackage[utf8x]{inputenc}" 
echo "\usepackage{graphicx}" >> overview.tex
echo "\usepackage{subfigure}" >> overview.tex
echo "\renewcommand*{\thesubfigure}{}" >> overview.tex
echo "\begin{document}" >> overview.tex
#echo "$title1" >> overview.tex
echo "\newpage" >> overview.tex

for i in "${sorted[@]}"
do 
  tmp1=`echo ${i/":"/"_"}`
  #tmp1=`echo ${tmp1/":"/"_"}` #for x_x_x folder names
  tmp1=`echo ${tmp1/":"/"/"}`
  echo $tmp1
 
  remainder=`echo "${counter}%${divisor}" | bc`
  #echo "Remainder: $remainder"
   
  if [[ "$remainder" == "0" ]]
  then  
    #echo "$case1" >> overview.tex
    echo "\begin{figure}" >> overview.tex
  fi

  echo "\subfigure[$i]{\includegraphics[width=2.5cm]{$tmp1/sol.png}}" >> overview.tex
   
  if [[ "$remainder" == "$divisorminuseins" ]]
  then  
    echo "\end{figure}" >> overview.tex
    echo "\setcounter{subfigure}{0}" >> overview.tex
  fi 
 
  let counter=counter+1
done
 
if [[ "$remainder" != "$divisorminuseins" ]]
then  
  echo "\end{figure}" >> overview.tex
fi

echo "\end{document}" >> overview.tex

pdflatex -interaction nonstopmode overview.tex 
 
