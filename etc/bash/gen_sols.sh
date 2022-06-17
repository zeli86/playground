#!/bin/bash
FILES=$(find . -type f -name results.csv | sort) 
anzahl=`find . -type f -name "results.csv" | wc -l`

counter=1
linewidth=2
yrange="*:*"
xrange="*:*"
pointsize="2"

farben=( "000000" "E41A1C" "377EB8" "4DAF4A" "984EA3" "FF7F00" "1B9E77" "D95F02" "7570B3" "E7298A" "66A61E" "E6AB02" "A6761D" "666666" "8000FF" "FF0000" "E41A1C" "377EB8" "4DAF4A" "984EA3" "FF7F00" "1B9E77" "D95F02" "7570B3" "E7298A" "66A61E" "E6AB02" "A6761D" "666666" "8000FF" "FF0000")
#echo ${#farben[*]}

cwd=$(pwd)
str1=`basename "$cwd"`
str2=${str1//_/ } 
 
for i in $FILES
do 
  echo $i
  ordner=`dirname "$i"`
  ordner2=${ordner:2}
  ordner3=${ordner2//_/,} 
  
  plt_filename="$ordner2.plt"
  png_filename="$ordner2.png"
  sol1_filename="$ordner2/0000/sol.png"
  sol2_filename="$ordner2/0074/sol.png"

  touch $plt_filename
  
  echo "set terminal pngcairo enhanced color notransparent crop size 1600,1000 font \"Arial,20\"" > "$plt_filename"
  echo "set key inside left top" >> "$plt_filename"
  echo "set datafile separator \";\"" >> "$plt_filename"
  echo "set xr[$xrange]" >> "$plt_filename"
  echo "set yr[$yrange]" >> "$plt_filename"
  echo "set pointsize $pointsize" >> "$plt_filename"
  echo "set out \"$png_filename\"" >> "$plt_filename"
  #echo "set title \"$str2\"" >> "$plt_filename"
  for((j=0;j<${#farben[*]};j++)); do echo "set style line $((j+1)) lc rgb '#${farben[$j]}' lw $linewidth" >> "$plt_filename"; done

  echo "set auto fix" >> "$plt_filename"
  echo "set size sq" >> "$plt_filename"
  echo "set multiplot" >> "$plt_filename"
  
  echo "set size 0.5,1" >> "$plt_filename"
  echo "set origin 0,0" >> "$plt_filename"
  printf "plot \"%s\" u 1:3 w lp ls 1 title \"\"\n" "$i" >> "$plt_filename"
  echo "set size 0.3,0.48" >> "$plt_filename"
  echo "set origin 0.4,0" >> "$plt_filename"
  echo "unset border" >> "$plt_filename"
  echo "unset xtics" >> "$plt_filename"
  echo "unset ytics" >> "$plt_filename"
  printf "plot \"%s\" binary filetype=png with rgbimage notitle\n" "$sol1_filename" >> "$plt_filename" 

  echo "set origin 0.75,0" >> "$plt_filename"
  printf "plot \"%s\" binary filetype=png with rgbimage notitle\n" "$sol2_filename" >> "$plt_filename" 
  
  gnuplot "$plt_filename"
  let counter=counter+1
done
