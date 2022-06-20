from conans import ConanFile, CMake
from conans.tools import download, unzip
import os
import shutil

class dealii_conan(ConanFile):
    name = "trilinos"
    version = "13.4.0"
    settings = "os", "compiler", "build_type", "arch"
    user = "atus"
    channel = "stable"
    generators = "cmake", "cmake_find_package", "virtualenv", "virtualrunenv"
    no_copy_source = True
    requires = "openmpi/4.1.4@atus/stable", "lapack/3.10.1@atus/stable", "boost/1.76.0@atus/stable"

    def source(self):
        src_archive = "trilinos.zip"
        download("https://github.com/trilinos/Trilinos/archive/refs/tags/trilinos-release-13-4-0.zip", src_archive)
        unzip(src_archive, strip_root=True)
        os.unlink(src_archive)

    def configure_cmake(self):
        cmake = CMake(self)
        cmake.verbose = False
        cmake.parallel = True        
        cmake.definitions["TPL_ENABLE_MPI"]="ON"
        cmake.definitions["BUILD_SHARED_LIBS"]="ON"
        cmake.definitions["Trilinos_ENABLE_Amesos"]="ON"
        cmake.definitions["Trilinos_ENABLE_Epetra"]="ON"
        cmake.definitions["Trilinos_ENABLE_EpetraExt"]="ON"
        cmake.definitions["Trilinos_ENABLE_Ifpack"]="ON"
        cmake.definitions["Trilinos_ENABLE_AztecOO"]="ON"
        cmake.definitions["Trilinos_ENABLE_Sacado"]="ON"
        cmake.definitions["Trilinos_ENABLE_SEACAS"]="OFF"
        cmake.definitions["Trilinos_ENABLE_Teuchos"]="ON"
        cmake.definitions["Trilinos_ENABLE_MueLu"]="OFF"
        cmake.definitions["Trilinos_ENABLE_Kokkos"]="OFF"
        cmake.definitions["Trilinos_ENABLE_ML"]="ON"
        cmake.definitions["Trilinos_ENABLE_ROL"]="ON"
        cmake.definitions["Trilinos_ENABLE_Tpetra"]="OFF"
        cmake.definitions["Trilinos_ENABLE_COMPLEX_DOUBLE"]="ON"
        cmake.definitions["Trilinos_ENABLE_COMPLEX_FLOAT"]="ON"
        cmake.definitions["Trilinos_ENABLE_Zoltan"]="ON"
        cmake.definitions["Trilinos_VERBOSE_CONFIGURE"]="OFF"
        cmake.definitions["TPL_ENABLE_Boost"]="OFF"
        cmake.definitions["TPL_ENABLE_Netcdf"]="OFF"
        pc_dir = os.path.join(self.deps_cpp_info["lapack"].rootpath, "lib")
        cmake.definitions["BLAS_LIBRARY_NAMES"]="blas"
        cmake.definitions["BLAS_LIBRARY_DIRS"]=pc_dir
        cmake.definitions["LAPACK_LIBRARY_NAMES"]="lapack"
        cmake.definitions["LAPACK_LIBRARY_DIRS"]=pc_dir
        cmake.configure()
        return cmake

    def build(self):
        cmake = self.configure_cmake()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()
        tmp = os.path.dirname(self.build_folder)
        shutil.rmtree( os.path.join(os.path.dirname(tmp), "source" ), ignore_errors=True )
        shutil.rmtree( self.build_folder, ignore_errors=True )
        
    def package_info(self):
        self.env_info.PATH.append(os.path.join(self.package_folder, "bin"))
        self.env_info.LD_LIBRARY_PATH.append(os.path.join(self.package_folder, "lib"))