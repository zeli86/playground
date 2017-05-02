try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

import os,sys,glob,struct,math,subprocess,csv
import vtk

#os.chdir("/home/zeli/tmp/delme1/")
cwd_folder=os.getcwd()

print cwd_folder

with open('results.csv', 'rb') as csvfile:
  tmp = csv.reader(csvfile, delimiter=';', quotechar='|')
  results = list(tmp)

#-------------------------------------------------------------------------
#Durchsucht alles unter "pathname" und macht liste aus gefundenen Datein:  
file_list=[]
sp = subprocess.Popen( 'find . -type f -name "final-*.vtu" | sed \'s:^\./::\'', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
file_list = sp.communicate()[0].split()

file_list.sort();

counter=1

f = open('bif.txt', 'w')

for file in file_list:

    print file 

    current_folder = cwd_folder + '/' + os.path.dirname( file ) + '/';
    os.chdir( current_folder );

    filename = os.path.basename( file )
    filename_noext = os.path.splitext(filename)[0]
   
    #-------------------------------------------------------------
    final_pvtu = XMLUnstructuredGridReader(PointArrayStatus=['Psi_sol','error per cell'], FileName=[filename] )
    pd = final_pvtu.PointData

    #for n in range(pd.GetNumberOfArrays()):
    #  print "Range for array ", pd.GetArray(n).GetName(), " ", pd.GetArray(n).GetRange()
    
    my_range = pd.GetArray("Psi_sol").GetRange()
    a = math.fabs(my_range[0])
    b = math.fabs(my_range[1])
    
    if a > b:
      c = a
    else:
      c = b
 
    f.write("%s %e\n" % (results[counter][0], c)) 
      
    Delete(final_pvtu)
    del final_pvtu
    
    counter += 1

f.close()

