from conans import ConanFile, CMake
from conans.tools import download, unzip
import os
import shutil

class lapack_conan(ConanFile):
    name = "lapack"
    version = "3.10.1"
    settings = "os", "compiler", "build_type", "arch"
    user = "atus"
    channel = "stable"
    generators = "cmake", "cmake_find_package", "virtualenv", "virtualrunenv"
    no_copy_source = True 

    def source(self):
        src_archive = "lapack.tar.gz"
        download("https://github.com/Reference-LAPACK/lapack/archive/v3.10.1.tar.gz", src_archive)
        unzip(src_archive, strip_root=True)
        os.unlink(src_archive)

    def configure_cmake(self):
        cmake = CMake(self)
        cmake.verbose = False
        cmake.parallel = True
        #cmake.build_folder = "BUILD"
        cmake.definitions["BUILD_SHARED_LIBS"] = "ON"
        cmake.configure()
        return cmake

    def build(self):
        cmake = self.configure_cmake()
        self.output.info( "self.source_folder = %s" % self.source_folder )
        self.output.info( "self.command_line = %s" % cmake.command_line )
        self.output.info( "self.build_config = %s" % cmake.build_config )
        cmake.build()

    def package(self):
        cmake = self.configure_cmake()
        cmake.install()
        tmp = os.path.dirname(self.build_folder)
        shutil.rmtree( os.path.join(os.path.dirname(tmp), "source" ), ignore_errors=True )
        shutil.rmtree( self.build_folder, ignore_errors=True )
        
    def package_info(self):
        self.env_info.PATH.append(os.path.join(self.package_folder, "bin"))
        self.env_info.LD_LIBRARY_PATH.append(os.path.join(self.package_folder, "lib"))

 