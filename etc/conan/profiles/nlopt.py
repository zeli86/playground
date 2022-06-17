
from conans import ConanFile, CMake
from conans.tools import download, untargz, check_md5, check_sha1, check_sha256
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
        zip_name = "nlopt-2.7.1.tar.gz"
        download("https://github.com/stevengj/nlopt/archive/v2.7.1.tar.gz", zip_name)
        # check_md5(zip_name, "51e11f2c02a36689d6ed655b6fff9ec9")
        # check_sha1(zip_name, "8d87812ce591ced8ce3a022beec1df1c8b2fac87")
        # check_sha256(zip_name, "653f983c30974d292de58444626884bee84a2731989ff5a336b93a0fef168d79")
        untargz(zip_name, strip_root=True)
        os.unlink(zip_name)


    def configure_cmake(self):
        cmake = CMake(self)
        cmake.verbose = False
        cmake.parallel = True
        cmake.definitions["NLOPT_CXX"] = "On"
        cmake.definitions["NLOPT_GUILE"] = "Off"
        cmake.definitions["NLOPT_LINK_PYTHON"] = "Off"
        cmake.definitions["NLOPT_MATLAB"] = "Off"
        cmake.definitions["NLOPT_OCTAVE"] = "Off"
        cmake.definitions["NLOPT_PYTHON"] = "Off"
        cmake.definitions["NLOPT_SWIG"] = "Off"

        cmake.configure()
        return cmake

    def build(self):
        cmake = self.configure_cmake()
        self.output.info( "self.source_folder = %s" % self.source_folder )
        self.output.info( "self.build_folder = %s" % self.build_folder )
        #self.output.info( "self.command_line = %s" % cmake.command_line )
        #self.output.info( "self.build_config = %s" % cmake.build_config )
        cmake.build()
        cmake.test()

    def package(self):
        cmake = self.configure_cmake()
        cmake.install()
        shutil.rmtree( self.source_folder )
        shutil.rmtree( self.build_folder )

    def package_info(self):
        self.env_info.PATH.append(os.path.join(self.package_folder, "bin"))
        self.env_info.LD_LIBRARY_PATH.append(os.path.join(self.package_folder, "lib"))

