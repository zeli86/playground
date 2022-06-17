
from conans import ConanFile, AutoToolsBuildEnvironment
from conans.tools import ftp_download, unzip, check_md5, check_sha1, check_sha256
import os
import shutil

class gsl_conan(ConanFile):
    name = "gsl"
    version = "2.7.1"
    settings = "os", "compiler", "build_type", "arch"
    user = "atus"
    channel = "stable"
    generators = "cmake_find_package"

    def source(self):
        src_archive = "gsl-2.7.1.tar.gz"
        ftp_download("ftp.gnu.org", "gnu/gsl/gsl-2.7.1.tar.gz")
        # check_md5(zip_name, "51e11f2c02a36689d6ed655b6fff9ec9")
        # check_sha1(zip_name, "8d87812ce591ced8ce3a022beec1df1c8b2fac87")
        # check_sha256(zip_name, "653f983c30974d292de58444626884bee84a2731989ff5a336b93a0fef168d79")
        untargz(src_archive, strip_root=True)
        os.unlink(src_archive)

    def build(self):
        autotools = AutoToolsBuildEnvironment(self)
        autotools.fpic = True
        env_build_vars = autotools.vars
        env_build_vars['RCFLAGS'] = '-march=native -O3'        
        autotools.configure()
        autotools.make()
        autotools.install()

    def package(self):
        shutil.rmtree( self.source_folder )
        shutil.rmtree( self.build_folder )

    def package_info(self):
        self.env_info.PATH.append(os.path.join(self.package_folder, "bin"))
        self.env_info.LD_LIBRARY_PATH.append(os.path.join(self.package_folder, "lib"))
