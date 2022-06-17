
from conans import ConanFile, AutoToolsBuildEnvironment
from conans.tools import download, unzip, check_md5, check_sha1, check_sha256
import os
import shutil

class p4est_conan(ConanFile):
    name = "p4est"
    version = "2.8.0"
    settings = "os", "compiler", "build_type", "arch"
    user = "atus"
    channel = "stable"
    generators = "cmake", "cmake_find_package"
    build_requires = "openmpi/4.1.4@atus/stable"

    def source(self):
        src_archive = "p4est.tar.gz"
        download("https://p4est.github.io/release/p4est-2.8.tar.gz", src_archive)
        unzip(src_archive, strip_root=True)
        os.unlink(src_archive)

    def build(self):
        autotools = AutoToolsBuildEnvironment(self)
        autotools.fpic = True
        env_build_vars = autotools.vars
        env_build_vars['RCFLAGS'] = '-march=native -O3'
        args = []
        args.append("--enable-mpi")
        args.append("--enable-shared")
        args.append("--disable-vtk-binary")
        args.append("--without-blas")
        autotools.configure(args=args)
        autotools.make()
        autotools.install()

    def package(self):
        tmp = os.path.dirname(self.build_folder)
        shutil.rmtree( os.path.join(os.path.dirname(tmp), "source" ), ignore_errors=True )
        shutil.rmtree( self.build_folder, ignore_errors=True )

    def package_info(self):
        self.env_info.PATH.append(os.path.join(self.package_folder, "bin"))
        self.env_info.LD_LIBRARY_PATH.append(os.path.join(self.package_folder, "lib"))
