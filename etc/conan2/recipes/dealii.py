# MIT License
# 
# Copyright (c) 2023 Zelimir Marojevic
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMakeDeps, CMake
from conan.tools.files import get, replace_in_file
import os

required_conan_version = ">=2.0"

class dealii_recipe(ConanFile):
    name = "dealii"
    version = "9.5.1"
    default_user = "atus"
    default_channel = "stable"
    settings = "os", "compiler", "build_type", "arch"
    no_copy_source = True
    options = {
        "shared": [True, False],
        "fpic": [True, False],
        "DEAL_II_WITH_64BIT_INDICES": ["ON", "OFF"],
        "DEAL_II_WITH_COMPLEX_VALUES": ["ON", "OFF"],
    }
    default_options = {
        "shared": True,
        "fpic": True,
        "DEAL_II_WITH_64BIT_INDICES": "OFF",
        "DEAL_II_WITH_COMPLEX_VALUES": "ON",
    }

    def requirements(self):
        for req in self.conan_data["tool_requires"]:
            self.tool_requires(req)
        for req in self.conan_data["requires"][self.name]:
            self.requires(req)

    def source(self):
        get(self, **self.conan_data["sources"][self.name][self.version])
        tmp = os.path.join(self.source_folder, "cmake", "modules", "FindDEAL_II_BOOST.cmake")
        replace_in_file(self, tmp, "set(BOOST_VERSION_REQUIRED 1.59)", "set(BOOST_VERSION_REQUIRED 1.82)")
        replace_in_file(self, tmp, "set(Boost_NO_BOOST_CMAKE ON)", "set(Boost_NO_BOOST_CMAKE ON)\nset(Boost_NO_SYSTEM_PATHS ON)")
        replace_in_file(self, tmp, "find_package(Boost ${BOOST_VERSION_REQUIRED}", "find_package(Boost ${BOOST_VERSION_REQUIRED} EXACT REQUIRED ")

    def generate(self):
        tc = CMakeToolchain(self)

        tc.variables["CMAKE_FIND_DEBUG_MODE"] = "OFF"
        tc.variables["DEAL_II_WITH_MPI"] = "ON"
        tc.variables["MAKE_C_COMPILER"] = "mpicc"
        tc.variables["MAKE_CXX_COMPILER"] = "mpicxx"
        tc.variables["DEAL_II_WITH_UMFPACK"] = "ON"
        tc.variables["DEAL_II_WITH_THREADS"] = "ON"
        tc.variables["DEAL_II_WITH_GINKO"] = "OFF"
        tc.variables["DEAL_II_WITH_ARPACK"] = "OFF"
        tc.variables["DEAL_II_WITH_OPENCASCADE"] = "OFF"
        tc.variables["DEAL_II_WITH_GMSH"] = "OFF"
        tc.variables["DEAL_II_WITH_ASSIMP"] = "OFF"
        tc.variables["DEAL_II_WITH_CGAL"] = "OFF"
        tc.variables["DEAL_II_WITH_CUDA"] = "OFF"
        tc.variables["DEAL_II_WITH_METIS"] = "OFF"
        tc.variables["DEAL_II_WITH_SYMENGINE"] = "OFF"
        tc.variables["DEAL_II_WITH_GINKGO"] = "OFF"
        tc.variables["DEAL_II_WITH_ARBORX"] = "OFF"
        tc.variables["DEAL_II_WITH_ARPACK"] = "OFF"

        tc.variables["DEAL_II_WITH_HDF5"] = "OFF"
        self.output.info( "HDF5_DIR = %s" % self.dependencies["hdf5"].package_folder )
        tc.variables["HDF5_DIR"] = self.dependencies["hdf5"].package_folder

        tc.variables["DEAL_II_WITH_SUNDIALS"] = "ON"
        self.output.info( "SUNDIALS_DIR = %s" % self.dependencies["sundials"].package_folder )
        tc.variables["SUNDIALS_DIR"] = self.dependencies["sundials"].package_folder

        tc.variables["DEAL_II_WITH_VTK"] = "ON"
        self.output.info( "VTK_DIR = %s" % self.dependencies["vtk"].package_folder )
        tc.variables["VTK_DIR"] = self.dependencies["vtk"].package_folder

        tc.variables["DEAL_II_WITH_ADOLC"] = "ON"
        self.output.info( "ADOLC_DIR = %s" % self.dependencies["adol-c"].package_folder )
        tc.variables["ADOLC_DIR"] = self.dependencies["adol-c"].package_folder  

        tc.variables["DEAL_II_WITH_GSL"] = "ON"
        self.output.info( "GSL_DIR = %s" % self.dependencies["gsl"].package_folder )
        tc.variables["GSL_DIR"] = self.dependencies["gsl"].package_folder  

        tc.variables["DEAL_II_WITH_BOOST"] = "ON"
        self.output.info( "Boost_DIR = %s" % self.dependencies["boost"].package_folder )
     
        tc.variables["DEAL_II_WITH_MUPARSER"] = "ON"
        self.output.info( "MUPARSER_DIR = %s" % self.dependencies["muparser"].package_folder )
        tc.variables["MUPARSER_DIR"] = self.dependencies["muparser"].package_folder
        
        tc.variables["DEAL_II_WITH_PETSC"] = "ON"
        self.output.info( "PETSC_DIR = %s" % self.dependencies["petsc"].package_folder )
        tc.variables["PETSC_DIR"] = self.dependencies["petsc"].package_folder
        
        tc.variables["DEAL_II_WITH_SLEPC"] = "ON"
        self.output.info( "SLEPC_DIR = %s" % self.dependencies["slepc"].package_folder )
        tc.variables["SLEPC_DIR"] = self.dependencies["slepc"].package_folder

        tc.variables["DEAL_II_WITH_TRILINOS"] = "OFF"
        self.output.info( "TRILINOS_DIR = %s" % self.dependencies["trilinos"].package_folder )
        tc.variables["TRILINOS_DIR"] = self.dependencies["trilinos"].package_folder

        tc.variables["DEAL_II_WITH_LAPACK"] = "ON"
        self.output.info( "LAPACK_DIR = %s" % self.dependencies["lapack"].package_folder )
        tc.variables["LAPACK_LIBRARIES"] = os.path.join(self.dependencies["lapack"].package_folder, "lib", "liblapack.so")
        tc.variables["BLAS_LIBRARIES"] = os.path.join(self.dependencies["lapack"].package_folder, "lib","libblas.so")
        
        tc.variables["DEAL_II_WITH_P4EST"] = "ON"
        self.output.info( "P4EST_DIR = %s" % self.dependencies["p4est"].package_folder )
        tc.variables["P4EST_DIR"] = self.dependencies["p4est"].package_folder
        
        tc.variables["DEAL_II_WITH_TBB"] = "ON"
        self.output.info( "TBB_DIR = %s" % self.dependencies["onetbb"].package_folder )
        tc.variables["TBB_DIR"] = self.dependencies["onetbb"].package_folder
        
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build(build_tool_args=["-j1"])

    def package(self):
        cmake = CMake(self)
        cmake.install()

