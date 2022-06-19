from conans import ConanFile, CMake
from conans.tools import download, untargz
import os
import shutil

class nlopt_conan(ConanFile):
    name = "nlopt"
    version = "2.7.1"
    settings = "os", "compiler", "build_type", "arch"
    user = "atus"
    channel = "stable"
    generators = "cmake", "virtualrunenv"
    no_copy_source = True 

    def source(self):
        src_archive = "nlopt-2.7.1.tar.gz"
        download("https://github.com/stevengj/nlopt/archive/v2.7.1.tar.gz", src_archive)
        untargz(src_archive, strip_root=True)
        os.unlink(src_archive)

    def configure_cmake(self):
        cmake = CMake(self)
        cmake.verbose = False
        cmake.parallel = True
        cmake.definitions["NLOPT_CXX"] = "ON"
        cmake.definitions["NLOPT_GUILE"] = "OFF"
        cmake.definitions["NLOPT_LINK_PYTHON"] = "OFF"
        cmake.definitions["NLOPT_MATLAB"] = "OFF"
        cmake.definitions["NLOPT_OCTAVE"] = "OFF"
        cmake.definitions["NLOPT_PYTHON"] = "OFF"
        cmake.definitions["NLOPT_SWIG"] = "OFF"
        cmake.configure()
        return cmake

    def build(self):
        cmake = self.configure_cmake()
        cmake.build()
        cmake.test()

    def package(self):
        cmake = self.configure_cmake()
        cmake.install()
        tmp = os.path.dirname(self.build_folder)
        shutil.rmtree( os.path.join(os.path.dirname(tmp), "source" ), ignore_errors=True )
        shutil.rmtree( self.build_folder, ignore_errors=True )

    def package_info(self):
        self.env_info.PATH.append(os.path.join(self.package_folder, "bin"))
        self.env_info.LD_LIBRARY_PATH.append(os.path.join(self.package_folder, "lib"))

