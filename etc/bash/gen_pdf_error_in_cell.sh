#!/bin/bash

###############################################################################
###############################################################################

FILES5=$(find . -name "error.png")

FILES5=`echo ${FILES5//"./"/}`
FILES5=`echo ${FILES5//"/error.png"/}`
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

rm error_20.tex
touch error_20.tex

echo "\documentclass[reprint,amsmath,amssymb,aps,floatfix]{revtex4}" >> error_20.tex
echo "\usepackage[utf8x]{inputenc}" 
echo "\usepackage{graphicx}" >> error_20.tex
echo "\usepackage{subfigure}" >> error_20.tex
echo "\renewcommand*{\thesubfigure}{}" >> error_20.tex
echo "\begin{document}" >> error_20.tex
#echo "$title1" >> error_20.tex
echo "\newpage" >> error_20.tex

for i in "${sorted[@]}"
do 
  tmp1=`echo ${i/":"/"_"}`
  tmp1=`echo ${tmp1/":"/"/"}`
  echo $tmp1
 
  remainder=`echo "${counter}%${divisor}" | bc`
  #echo "Remainder: $remainder"
   
  if [[ "$remainder" == "0" ]]
  then  
    #echo "$case1" >> error_20.tex
    echo "\begin{figure}" >> error_20.tex
  fi

  echo "\subfigure[$i]{\includegraphics[width=5.3cm]{$tmp1/error.png}}" >> error_20.tex
   
  if [[ "$remainder" == "$divisorminuseins" ]]
  then  
    echo "\end{figure}" >> error_20.tex
    echo "\setcounter{subfigure}{0}" >> error_20.tex
  fi 
 
  let counter=counter+1
done
 
if [[ "$remainder" != "$divisorminuseins" ]]
then  
  echo "\end{figure}" >> error_20.tex
fi

echo "\end{document}" >> error_20.tex

pdflatex -interaction nonstopmode error_20.tex 
 
