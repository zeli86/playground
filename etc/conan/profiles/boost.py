
from conans import ConanFile, CMake
from conans.tools import download, unzip, check_md5, check_sha1, check_sha256, cpu_count
import os
import shutil

class boost_conan(ConanFile):
    name = "boost"
    version = "1.76.0"
    settings = "os", "compiler", "build_type", "arch"
    user = "atus"
    channel = "stable"
    generators = "cmake", "virtualrunenv"

    def source(self):
        src_archive = "boost_src.tar.gz"
        download("https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.bz2", src_archive)
        # check_md5(src_archive, "51e11f2c02a36689d6ed655b6fff9ec9")
        # check_sha1(src_archive, "8d87812ce591ced8ce3a022beec1df1c8b2fac87")
        check_sha256(src_archive, "f0397ba6e982c4450f27bf32a2a83292aba035b827a5623a14636ea583318c41")
        unzip(src_archive, strip_root=True)
        os.unlink(src_archive)

    @property
    def _toolset(self):
        return "gcc-11"

    @property
    def _build_flags(self):
        flags = []

        flags.append("install")
        if self.settings.build_type == "Debug":
            flags.append("variant=debug")
        else:
            flags.append("variant=release")        
        flags.append("architecture=x86")
        flags.append("address-model=64")
        flags.append("link=shared")
        flags.append("runtime-link=shared")
        flags.append("threading=multi")
        flags.append("toolset=%s" % self._toolset)
        flags.append("-q") # Stop at the first error. No need to continue building.
        flags.append("--disable-icu")
        flags.append("--disable-iconv")
        flags.append("--without-mpi")
        flags.append("--without-python")

        flags.extend([
            "--prefix=%s" % self.package_folder,
            "-j%s" % cpu_count(),
            "--abbreviate-paths"
        ])
        return flags

    def build(self):
        b2_flags = " ".join(self._build_flags)
        self.output.info( "build_flags = %s" % b2_flags )
        self.run("./bootstrap.sh")
        self.run("./b2 %s" % b2_flags )

    def package(self):
        #self.output.info( "self.source_folder = %s" % self.source_folder )
        #self.output.info( "self.build_folder = %s" % self.build_folder )
        tmp = os.path.dirname(self.build_folder)
        shutil.rmtree( os.path.join(os.path.dirname(tmp), "source" ), ignore_errors=True )
        shutil.rmtree( self.build_folder, ignore_errors=True )
    
    def package_info(self):
        self.env_info.PATH.append(os.path.join(self.package_folder, "bin"))
        self.env_info.LD_LIBRARY_PATH.append(os.path.join(self.package_folder, "lib"))

