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
    requires = "openmpi/4.1.4@atus/stable"

    def source(self):
        src_archive = "trilinos.tar.gz"
        download("https://github.com/trilinos/Trilinos/archive/refs/tags/trilinos-release-13-4-0.zip", src_archive)
        unzip(src_archive, strip_root=True)
        os.unlink(src_archive)

    def configure_cmake(self):
        cmake = CMake(self)
        cmake.verbose = False
        cmake.parallel = True
        cmake.definitions["TPL_ENABLE_MPI"] = "ON"
        cmake.definitions["Trilinos_ENABLE_ALL_PACKAGES"] = "ON"
        cmake.definitions["MPI_BASE_DIR"] = self.deps_cpp_info["openmpi"].rootpath
        cmake.configure()
        return cmake

    def build(self):
        cmake = self.configure_cmake()
        self.output.info( "self.source_folder = %s" % self.source_folder )
        self.output.info( "self.command_line = %s" % cmake.command_line )
        self.output.info( "self.build_config = %s" % cmake.build_config )
        cmake.build()