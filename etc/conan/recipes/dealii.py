from conans import ConanFile, CMake
from conans.tools import download, unzip
import os
import shutil

class dealii_conan(ConanFile):
    name = "dealii"
    version = "9.3.0"
    settings = "os", "compiler", "build_type", "arch"
    user = "atus"
    channel = "stable"
    generators = "cmake", "cmake_find_package", "virtualenv", "virtualrunenv"
    no_copy_source = True
    requires = "openmpi/4.1.4@atus/stable", "gsl/2.7.1@atus/stable", "lapack/3.10.1@atus/stable", "muparser/2.3.3@atus/stable", "p4est/2.8.0@atus/stable", "boost/1.76.0@atus/stable", "petsc/3.17.2@atus/stable", "slepc/3.17.1@atus/stable", "trilinos/13.4.0@atus/stable"

    def source(self):
        src_archive = "dealii.tar.gz"
        download("https://github.com/dealii/dealii/releases/download/v9.3.0/dealii-9.3.0.tar.gz", src_archive)
        unzip(src_archive, strip_root=True)
        os.unlink(src_archive)

    def configure_cmake(self):
        cmake = CMake(self)
        cmake.verbose = False
        cmake.parallel = True
            
        cmake.definitions["DEAL_II_WITH_MPI"] = "ON"
        cmake.definitions["MAKE_C_COMPILER"] = "mpicc"
        cmake.definitions["MAKE_CXX_COMPILER"] = "mpicxx"
        cmake.definitions["DEAL_II_WITH_UMFPACK"] = "ON"
        cmake.definitions["DEAL_II_WITH_THREADS"] = "OFF"
        cmake.definitions["DEAL_II_WITH_HDF5"] = "OFF"
        cmake.definitions["DEAL_II_WITH_GINKO"] = "OFF"

        cmake.definitions["DEAL_II_WITH_GSL"] = "ON"
        self.output.info( "GSL_DIR = %s" % self.deps_cpp_info["gsl"].rootpath )
        cmake.definitions["GSL_DIR"] = self.deps_cpp_info["gsl"].rootpath;  
        
        cmake.definitions["DEAL_II_WITH_BOOST"] = "ON"
        self.output.info( "BOOST_DIR = %s" % self.deps_cpp_info["boost"].rootpath )
        cmake.definitions["BOOST_DIR"] = self.deps_cpp_info["boost"].rootpath;
        
        cmake.definitions["DEAL_II_WITH_MUPARSER"] = "ON"
        self.output.info( "MUPARSER_DIR = %s" % self.deps_cpp_info["muparser"].rootpath )
        cmake.definitions["MUPARSER_DIR"] = self.deps_cpp_info["muparser"].rootpath;
        
        cmake.definitions["DEAL_II_WITH_PETSC"] = "ON"
        self.output.info( "PETSC_DIR = %s" % self.deps_cpp_info["petsc"].rootpath )
        cmake.definitions["PETSC_DIR"] = self.deps_cpp_info["petsc"].rootpath;
        
        cmake.definitions["DEAL_II_WITH_SLEPC"] = "ON"
        self.output.info( "SLEPC_DIR = %s" % self.deps_cpp_info["slepc"].rootpath )
        cmake.definitions["SLEPC_DIR"] = self.deps_cpp_info["slepc"].rootpath;

        cmake.definitions["DEAL_II_WITH_TRILINOS"] = "ON"
        self.output.info( "TRILINOS_DIR = %s" % self.deps_cpp_info["trilinos"].rootpath )
        cmake.definitions["TRILINOS_DIR"] = self.deps_cpp_info["trilinos"].rootpath;

        cmake.definitions["DEAL_II_WITH_LAPACK"] = "ON"
        self.output.info( "LAPACK_DIR = %s" % self.deps_cpp_info["lapack"].rootpath )
        cmake.definitions["LAPACK_DIR"] = self.deps_cpp_info["lapack"].rootpath;
        
        cmake.definitions["DEAL_II_WITH_P4EST"] = "ON"
        self.output.info( "P4EST_DIR = %s" % self.deps_cpp_info["p4est"].rootpath )
        cmake.definitions["P4EST_DIR"] = self.deps_cpp_info["p4est"].rootpath;

        cmake.configure()
        return cmake

    def build(self):
        cmake = self.configure_cmake()
        self.output.info( "self.source_folder = %s" % self.source_folder )
        self.output.info( "self.command_line = %s" % cmake.command_line )
        self.output.info( "self.build_config = %s" % cmake.build_config )
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()
        #shutil.rmtree( self.source_folder )
        #shutil.rmtree( self.build_folder )
        
    def package_info(self):
        self.env_info.PATH.append(os.path.join(self.package_folder, "bin"))
        self.env_info.LD_LIBRARY_PATH.append(os.path.join(self.package_folder, "lib"))
