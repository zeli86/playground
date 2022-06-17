
from conans import ConanFile, AutoToolsBuildEnvironment
from conans.tools import download, untargz
import os
import shutil

class petsc_conan(ConanFile):
    name = "petsc"
    version = "3.17.2"
    settings = "os", "compiler", "build_type", "arch"
    user = "atus"
    channel = "stable"
    generators = "cmake", "virtualenv", "virtualrunenv"
    build_requires = "openmpi/4.1.4@atus/stable"

    def source(self):
        src_archive = "petsc.tar.gz"
        download("http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.17.2.tar.gz", src_archive)
        untargz(src_archive, strip_root=True)
        os.unlink(src_archive)

    def build(self):
        autotools = AutoToolsBuildEnvironment(self)
        autotools.fpic = True
        env_build_vars = autotools.vars
        args = []
        args.append("--with-cc=mpicc")
        args.append("--with-cxx=mpicxx") 
        args.append("--with-fc=mpif90") 
        args.append("--with-shared-libraries")
        args.append("--with-x=0") 
        args.append("--with-debugging=0") 
        args.append("--with-mpi=1") 
        args.append("--download-hypre=yes") 
        args.append("--download-fblaslapack=1") 
        args.append("--download-scalapack") 
        args.append("--download-mumps") 
        args.append("--download-ptscotch")
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
