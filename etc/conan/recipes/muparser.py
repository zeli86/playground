from conans import ConanFile, CMake
from conans.tools import download, unzip
import os
import shutil

class muparser_conan(ConanFile):
    name = "muparser"
    version = "2.3.3"
    settings = "os", "compiler", "build_type", "arch"
    user = "atus"
    channel = "stable"
    generators = "cmake", "cmake_find_package", "virtualenv", "virtualrunenv"
    no_copy_source = True 

    def source(self):
        src_archive = "muparser.zip"
        download("https://github.com/beltoforion/muparser/archive/refs/tags/v2.3.3-1.zip", src_archive)
        unzip(src_archive, strip_root=True)
        os.unlink(src_archive)

    def configure_cmake(self):
        cmake = CMake(self)
        cmake.verbose = False
        cmake.parallel = True
        cmake.configure()
        return cmake

    def build(self):
        cmake = self.configure_cmake()
        cmake.build()

    def package(self):
        cmake = self.configure_cmake()
        cmake.install()
        tmp = os.path.dirname(self.build_folder)
        shutil.rmtree( os.path.join(os.path.dirname(tmp), "source" ), ignore_errors=True )
        shutil.rmtree( self.build_folder, ignore_errors=True )
        
    def package_info(self):
        self.env_info.LD_LIBRARY_PATH.append(os.path.join(self.package_folder, "lib"))