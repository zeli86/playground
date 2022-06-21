from conans import ConanFile, CMake
from conans.tools import download, untargz
import os
import shutil

class iceoryx_conan(ConanFile):
    name = "iceoryx"
    version = "2.0.2"
    settings = "os", "compiler", "build_type", "arch"
    user = "atus"
    channel = "stable"
    generators = "cmake", "cmake_find_package", "virtualrunenv"
    no_copy_source = True 

    def source(self):
        src_archive = "iceoryx.tar.gz"
        download("https://github.com/eclipse-iceoryx/iceoryx/archive/refs/tags/v2.0.2.tar.gz", src_archive)
        untargz(src_archive, strip_root=True)
        os.unlink(src_archive)

    def configure_cmake(self):
        cmake = CMake(self)
        cmake.verbose = False
        cmake.parallel = True
        cmake.definitions["BUILD_SHARED_LIBS"] = "ON"
        cmake.definitions["TOML_CONFIG"] = "ON"
        cmake.definitions["DOWNLOAD_TOML_LIB"] = "ON"
        cmake.configure(source_folder="iceoryx_meta")
        return cmake

    def build(self):
        cmake = self.configure_cmake()
        cmake.build()
        cmake.install()

    def package(self):
        tmp = os.path.dirname(self.build_folder)
        shutil.rmtree( os.path.join(os.path.dirname(tmp), "source" ), ignore_errors=True )
        shutil.rmtree( self.build_folder, ignore_errors=True )

    def package_info(self):
        self.env_info.PATH.append(os.path.join(self.package_folder, "bin"))
        self.env_info.LD_LIBRARY_PATH.append(os.path.join(self.package_folder, "lib"))

