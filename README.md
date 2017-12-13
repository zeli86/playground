
# **ATUS-PRO PLAYGROUND**
## **Introduction**
The PLAYGROUND package is a collection of C++ programs designed to solve the non linear Schrödinger equation, especially the Gross--Pitaevskii equation. This is the active development branch of ATUS-PRO.

## **Requirements and Dependencies**  

* make, cmake
* gcc, g++, gfortran
* Modules (highly recommended)
* Steel Bank Common Lisp (recommended for install script)
* gnuplot
* Paraview (optional)

The following packages are required and will be installed via the install script:

* openmpi 3.0.0 (https://www.open-mpi.org/)
* gnu gsl 2.4 (https://www.gnu.org/software/gsl/)
* nlopt 2.4.2 (http://ab-initio.mit.edu/nlopt) 
* muparser 2.2.5 (http://beltoforion.de/article.php?a=muparser)
* lapack 3.8.0 (http://www.netlib.org/lapack/)
* p4est 2.0 (http://p4est.org/)
* petsc 3.8.2 (http://www.mcs.anl.gov/petsc/) 
* deal.ii 8.5.1 (http://www.dealii.org/) 

The following third party packages are included in the source tree:

* pugixml (https://pugixml.org/)
* String Toolkit (https://www.partow.net/programming/strtk/index.html)
* CXXOPTS (https://github.com/jarro2783/cxxopts)

## **Installation**

It is highly recommended to use our install script which downloads all required packages and installs everything in $HOME/local/ by default. In order to run this script Steel Bank Common Lisp (sbcl) is required.

Further we recommended the use of environmental modules (http://modules.sourceforge.net/) for setting up all paths. Modules should be available via your Linux distribution. The install script also generates module files which are located in $HOME/local/modules. The search paths for the module files needs to be extended by adding your path to $MODULEPATH of your shell.

The binaries are installed in $HOME/bin. Make sure that this folder is added to $PATH of your shell.

### **Install local dependencies**
Make sure that the install script *install_8.5.1* is executable. This scripts requires an active internet connection. 

In the playground folder invoke on the command line the install script via: *./install_8.5.1 -j 4*

Then select 'a' for all packages and press enter. This will take a while depending on your computer speed. The -j option specifies the number of utilized cores. At the end of the installation confirm the the extension of the MODULEPATH environmental variable. This adds the local module files folder. 

* Run: *module avail* (check the available modules) 
* Run: *module load atus-pro-git-meta* (load all required modules)
* Create a build folder below playground: *mkdir Build*
* Change to Build folder: *cd Build*
* Run: *cmake ..*
* Run: *ccmake .* or *cmake-gui .* (Check the CACHE. The base paths of the library files should point to /home/yourname/local/opt/. )
* Build: *make -j 4*
* Test: *make test* or *ctest -V* (verbose output)


   