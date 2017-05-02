#!/bin/bash
FILES=$(find . -type f -name results.csv | sort) 
anzahl=`find . -type f -name "results.csv" | wc -l`

counter=1
linewidth=2
#xrange="*:*"
yrange="*:*"
xrange="2:14"
#yrange="0:750"
pointsize="3"

farben=("E41A1C" "377EB8" "4DAF4A" "984EA3" "FF7F00" "1B9E77" "D95F02" "7570B3" "E7298A" "66A61E" "E6AB02" "A6761D" "666666" "8000FF" "FF0000" "E41A1C" "377EB8" "4DAF4A" "984EA3" "FF7F00" "1B9E77" "D95F02" "7570B3" "E7298A" "66A61E" "E6AB02" "A6761D" "666666" "8000FF" "FF0000")
#echo ${#farben[*]}
#set style line 1 lc rgb '#800000' lt 1 lw 2
#set style line 2 lc rgb '#ff0000' lt 1 lw 2
#set style line 3 lc rgb '#ff4500' lt 1 lw 2
#set style line 4 lc rgb '#ffa500' lt 1 lw 2
#set style line 5 lc rgb '#006400' lt 1 lw 2
#set style line 6 lc rgb '#0000ff' lt 1 lw 2
#set style line 7 lc rgb '#9400d3' lt 1 lw 2

cwd=$(pwd)
str1=`basename "$cwd"`
str2=${str1//_/ } 

rm bif.png
touch bif.plt
echo "set terminal pngcairo enhanced color notransparent crop size 1920,1080 font \"Arial,20\"" > bif.plt
echo "set key inside left top" >> bif.plt
echo "set datafile separator \";\"" >> bif.plt
echo "set xr[$xrange]" >> bif.plt
echo "set yr[$yrange]" >> bif.plt
#echo "set log x" >> bif.plt
#echo "set log y" >> bif.plt
echo "set grid" >> bif.plt
echo "set pointsize $pointsize" >> bif.plt 
echo "set out \"bif.png\"" >> bif.plt
echo "set title \"$str2\"" >> bif.plt
for((i=0;i<${#farben[*]};i++)); do echo "set style line $((i+1)) lc rgb '#${farben[$i]}' lw $linewidth" >> bif.plt; done
#for((i=0;i<${#farben[*]};i++)); do echo "set style line $((i+1)) lc rgb '#${farben[$i]}' pt 1 ps 1.5 lt $((i+1)) lw $linewidth" >> bif.plt; done

for i in $FILES
do 
  ordner=`dirname "$i"`
  ordner2=${ordner:2}
  ordner3=${ordner2//_/,} 
  
  if [ "$counter" = 1 ]
  then
    printf "plot \"%s\" u 1:3 w lp ls %d title \"%s\"" "$i" "$counter" "$ordner3" >> bif.plt 
  else
    printf ", \"%s\" u 1:3 w lp ls %d title \"%s\"" "$i" "$counter" "$ordner3" >> bif.plt
  fi    
 
  let counter=counter+1
done

gnuplot bif.plt
display bif.png
