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
divisor=15
divisorminuseins=14

rm main_20.tex
touch main_20.tex

echo "\documentclass[reprint,amsmath,amssymb,aps,floatfix]{revtex4}" >> main_20.tex
echo "\usepackage[utf8x]{inputenc}" 
echo "\usepackage{graphicx}" >> main_20.tex
echo "\usepackage{subfigure}" >> main_20.tex
echo "\renewcommand*{\thesubfigure}{}" >> main_20.tex
echo "\begin{document}" >> main_20.tex
#echo "$title1" >> main_20.tex
echo "\newpage" >> main_20.tex

for i in "${sorted[@]}"
do 
  tmp1=`echo ${i/":"/"_"}`
  tmp1=`echo ${tmp1/":"/"/"}`
  echo $tmp1
 
  remainder=`echo "${counter}%${divisor}" | bc`
  #echo "Remainder: $remainder"
   
  if [[ "$remainder" == "0" ]]
  then  
    #echo "$case1" >> main_20.tex
    echo "\begin{figure}" >> main_20.tex
  fi

  echo "\subfigure[$i]{\includegraphics[width=5.3cm]{$tmp1/sol.png}}" >> main_20.tex
   
  if [[ "$remainder" == "$divisorminuseins" ]]
  then  
    echo "\end{figure}" >> main_20.tex
    echo "\setcounter{subfigure}{0}" >> main_20.tex
  fi 
 
  let counter=counter+1
done
 
if [[ "$remainder" != "$divisorminuseins" ]]
then  
  echo "\end{figure}" >> main_20.tex
fi

echo "\end{document}" >> main_20.tex

pdflatex -interaction nonstopmode main_20.tex 
 
