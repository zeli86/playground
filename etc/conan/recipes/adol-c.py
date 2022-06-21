from conans import ConanFile, AutoToolsBuildEnvironment
from conans.tools import download, unzip
import os
import shutil

class gsl_conan(ConanFile):
    name = "adol-c"
    version = "2.7.2"
    settings = "os", "compiler", "build_type", "arch"
    user = "atus"
    channel = "stable"
    generators = "cmake", "cmake_find_package", "virtualenv", "virtualrunenv"
    requires = "boost/1.76.0@atus/stable"

    def source(self):
        src_archive = "adol-c.tar.gz"
        download("https://github.com/coin-or/ADOL-C/archive/refs/tags/releases/2.7.2.tar.gz", src_archive)
        unzip(src_archive, strip_root=True)
        os.unlink(src_archive)

    def build(self):#
        self.run('pwd')
        autotools = AutoToolsBuildEnvironment(self)
        autotools.fpic = True
        env_build_vars = autotools.vars
        env_build_vars['RCFLAGS'] = '-march=native -O3'
        args = []
        args.append('--with-boost={}'.format(self.deps_cpp_info["boost"].rootpath))
        autotools.configure(args=args)
        autotools.make()
        autotools.install()

    def package(self):
        tmp = os.path.dirname(self.build_folder)
        shutil.rmtree( os.path.join(os.path.dirname(tmp), "source" ), ignore_errors=True )
        shutil.rmtree( self.build_folder, ignore_errors=True )

    def package_info(self):
        self.env_info.LD_LIBRARY_PATH.append(os.path.join(self.package_folder, "lib"))
