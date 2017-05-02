#!/bin/bash

FILES=$(find . -type f -name "results.csv" | sort) 
anzahl=`find . -type f -name "results.csv" | wc -l`

#cwd=$(pwd)
#str1=`basename "$cwd"`
#str2=${str1//_/ } 

counter=1
linewidth=2
#xrange="4:20"
#yrange="0:750"
xrange="*:*"
yrange="*:*"
pointsize="4"

#farben=("E41A1C" "377EB8" "4DAF4A" "984EA3" "FF7F00" "1B9E77" "D95F02" "7570B3" "E7298A" "66A61E" "E6AB02" "A6761D" "666666" "8000FF" "FF0000" "E41A1C" "377EB8" "4DAF4A" "984EA3" "FF7F00" "1B9E77" "D95F02" "7570B3" "E7298A" "66A61E" "E6AB02" "A6761D" "666666" "8000FF" "FF0000")
#farben=("A6CEE3" "1F78B4" "B2DF8A" "33A02C" "FB9A99" "E31A1C" "FDBF6F" "FF7F00")
farben=("000000" "990000" "009900" "000099" "FF3300" "FFCC00" "FF00FF" )
anzahl_farben=${#farben[*]}

#cases=$(ls)
#cases=`find . -maxdepth 1 -type d -name [^\.]\* | sed 's:^\./::'`

rm ??????.eps
rm ??????.plt
rm ??????.tex

for i in $FILES
do
  echo $i

  ordner=`dirname "$i"`
  prefix=${ordner//_}
  prefix=${prefix//.}
  prefix=${prefix///}
  random_str=$(cat /dev/urandom | tr -dc 'a-z' | fold -w 4 | head -n 1)
  
  tmp1="$prefix$random_str"
  
  outputfile="$tmp1.plt"
  outputfile2="$tmp1.png"
  outputfile3="$tmp1.tex"
  touch $outputfile
  
  echo $outputfile3
    
  for((j=0;j<$anzahl_farben;j++)); do echo "set style line $((j+1)) lc rgb '#${farben[$j]}' lw $linewidth ps 1" >> $outputfile; done
  echo "set terminal epslatex input newstyle leveldefault color dashed dl 1 lw 2 rounded clip size 8cm,7.12cm noheader" >> $outputfile
  echo "unset key" >> $outputfile
  echo "set grid" >> $outputfile
  echo "set out \""$outputfile3"\"" >> $outputfile
  echo "set datafile separator \";\"" >> $outputfile
  echo "set xlabel \"$\\\\mu$\"" >> $outputfile
  echo "set ylabel \"$\\\\gamma$\" offset 2,0,0" >> $outputfile
      
  printf "plot \"%s\" u 1:(2*pi*\$2*\$3) w l ls 1 notitle\n" "$i" >> $outputfile 
  
  gnuplot $outputfile
done    

