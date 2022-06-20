from conans import ConanFile, AutoToolsBuildEnvironment
from conans.tools import download, unzip, check_md5
import os
import shutil

class slepc_conan(ConanFile):
    name = "slepc"
    version = "3.17.1"
    settings = "os", "compiler", "build_type", "arch"
    user = "atus"
    channel = "stable"
    generators = "cmake", "cmake_find_package", "virtualenv", "virtualrunenv"
    requires = "openmpi/4.1.4@atus/stable", "petsc/3.17.2@atus/stable"

    def source(self):
        src_archive = "slepc.tar.gz"
        download("https://slepc.upv.es/download/distrib/slepc-3.17.1.tar.gz", src_archive)
        check_md5(src_archive, "a0ee1d2306b5388b4b94339f82f93550")
        unzip(src_archive, strip_root=True)
        os.unlink(src_archive)
    
    def build(self):
        autotools = AutoToolsBuildEnvironment(self)
        autotools.fpic = True
        env_build_vars = autotools.vars
        #env_build_vars['RCFLAGS'] = '-march=native -O3'
        args = []
        args.append("--with-scalapack=1")
        args.append("--download-arpack")
        args.append("--download-blopex")
        args.append("--download-elpa")
        args.append("--download-evsl")
        args.append("--download-hpddm")
        args.append("--download-primme")
        args.append("--download-trlan")
        autotools.configure(args=args)
        args = []
        args.append('PETSC_DIR={}'.format(self.deps_cpp_info["petsc"].rootpath))
        args.append('SLEPC_DIR={}'.format(self.build_folder))
        autotools.make(args=args)
        autotools.install(args=args)
        
    def package(self):
        tmp = os.path.dirname(self.build_folder)
        shutil.rmtree( os.path.join(os.path.dirname(tmp), "source" ), ignore_errors=True )
        shutil.rmtree( self.build_folder, ignore_errors=True )

    def package_info(self):
        self.env_info.PATH.append(os.path.join(self.package_folder, "bin"))
        self.env_info.LD_LIBRARY_PATH.append(os.path.join(self.package_folder, "lib"))