from conans import ConanFile, tools, AutoToolsBuildEnvironment
from conans.tools import download, unzip, check_sha256
from conans.errors import ConanInvalidConfiguration
import os
import shutil

required_conan_version = ">=1.29.1"

class OpenMPIConan(ConanFile):
    name = "openmpi"
    homepage = "https://www.open-mpi.org"
    url = ""
    topics = ("mpi", "openmpi")
    description = "A High Performance Message Passing Library"
    license = "BSD-3-Clause"
    settings = "os", "arch", "compiler", "build_type"
    version = "4.1.4"
    user = "atus"
    channel = "stable"
    generators = "cmake", "virtualenv", "virtualrunenv"

    _autotools = None

    def configure(self):
        del self.options.fPIC
        del self.settings.compiler.libcxx
        del self.settings.compiler.cppstd
        if self.settings.os == "Windows":
            raise ConanInvalidConfiguration("OpenMPI doesn't support Windows")

    def source(self):
        src_archive = "openmpi_src.tar.gz"
        download("https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.4.tar.bz2", src_archive)
        check_sha256(src_archive, "92912e175fd1234368c8730c03f4996fe5942e7479bb1d10059405e7f2b3930d")
        unzip(src_archive, strip_root=True)
        os.unlink(src_archive)

    def _configure_autotools(self):
        if self._autotools:
            return self._autotools
        self._autotools = AutoToolsBuildEnvironment(self)
        args = []
        if self.settings.build_type == "Debug":
            args.append("--enable-debug")
        args.extend(["--enable-shared", "--disable-static"])
        args.append("--with-pic")
        args.append("--enable-mpi-fortran" )
        self._autotools.configure(args=args)  
        return self._autotools

    def build(self):
        with tools.chdir(self.source_folder):
            autotools = self._configure_autotools()
            self._autotools.make()

    def package(self):
        self.copy(pattern="LICENSE", src=self.source_folder, dst="licenses")
        with tools.chdir(self.source_folder):
             self._autotools.install()
        tools.rmdir(os.path.join(self.package_folder, "lib", "pkgconfig"))
        tmp = os.path.dirname(self.build_folder)
        shutil.rmtree( os.path.join(os.path.dirname(tmp), "source" ), ignore_errors=True )
        shutil.rmtree( self.build_folder, ignore_errors=True )

    def package_info(self):
        self.cpp_info.libs = ['mpi', 'open-rte', 'open-pal']
        if self.settings.os == "Linux":
            self.cpp_info.system_libs = ["dl", "pthread", "rt", "util"]

        self.output.info("Creating MPI_HOME environment variable: {}".format(self.package_folder))
        self.env_info.MPI_HOME = self.package_folder
        self.output.info("Creating OPAL_PREFIX environment variable: {}".format(self.package_folder))
        self.env_info.OPAL_PREFIX = self.package_folder
        mpi_bin = os.path.join(self.package_folder, 'bin')
        self.output.info("Creating MPI_BIN environment variable: {}".format(mpi_bin))
        self.env_info.MPI_BIN = mpi_bin
        self.output.info("Appending PATH environment variable: {}".format(mpi_bin))
        self.env_info.PATH.append(mpi_bin)
        mpi_lib = os.path.join(self.package_folder, 'lib')
        self.output.info("Appending LD_LIBRARY_PATH environment variable: {}".format(mpi_lib))
        self.env_info.LD_LIBRARY_PATH.append(mpi_lib)

