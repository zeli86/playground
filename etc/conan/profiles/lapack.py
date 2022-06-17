
from conans import ConanFile, CMake
from conans.tools import download, unzip, check_md5, check_sha1, check_sha256
import os
import shutil

class lapack_conan(ConanFile):
    name = "lapack"
    version = "3.15.1"
    settings = "os", "compiler", "build_type", "arch"
    user = "atus"
    channel = "stable"
    generators = "cmake", "virtualrunenv"
    no_copy_source = True 

    def source(self):
        src_archive = "lapack.tar.gz"
        download("https://github.com/Reference-LAPACK/lapack/archive/v3.15.1.tar.gz", src_archive)
        # check_md5(zip_name, "51e11f2c02a36689d6ed655b6fff9ec9")
        # check_sha1(zip_name, "8d87812ce591ced8ce3a022beec1df1c8b2fac87")
        # check_sha256(zip_name, "653f983c30974d292de58444626884bee84a2731989ff5a336b93a0fef168d79")
        unzip(src_archive, strip_root=True)
        os.unlink(src_archive)

    def configure_cmake(self):
        cmake = CMake(self)
        cmake.verbose = False
        cmake.parallel = True
        cmake.build_folder = "BUILD"
        cmake.definitions["BUILD_SHARED_LIBS"] = "On"
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
        shutil.rmtree( self.source_folder )
        shutil.rmtree( self.build_folder )
        
    def package_info(self):
        self.env_info.PATH.append(os.path.join(self.package_folder, "bin"))
        self.env_info.LD_LIBRARY_PATH.append(os.path.join(self.package_folder, "lib"))

 