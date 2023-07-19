# MIT License
# 
# Copyright (c) 2023 Zelimir Marojevic
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from conan import ConanFile
from conan.tools.gnu import AutotoolsToolchain, Autotools
from conan.tools.files import get
import os
import shutil

required_conan_version = ">=2.0"

class open_mpi_recipe(ConanFile):
    name = "openmpi"
    version = "5.0.1"
    user = "atus"
    channel = "stable"
    homepage = "https://www.open-mpi.org"
    topics = ("mpi", "openmpi")
    description = "A High Performance Message Passing Library"
    license = "BSD-3-Clause"
    settings = "os", "arch", "compiler", "build_type"
    options = {"shared": [True, False], "fpic": [True, False]}
    default_options = {"shared": True, "fpic": True}
    _autotools = None

    def configure(self):
        if self.settings.os == "Windows":
            raise ConanInvalidConfiguration("OpenMPI doesn't support Windows")

    def source(self):
        get(self, **self.conan_data["sources"][self.name][self.version])

    def generate(self):
        tc = AutotoolsToolchain(self)
        tc.update_configure_args({"--prefix": f"{self.package_folder}", 
                                  "--bindir": None, 
                                  "--sbindir": None,
                                  "--libdir": None, 
                                  "--includedir": None, 
                                  "--oldincludedir": None});
        tc.configure_args.append("--enable-mpi-fortran")
        tc.configure_args.append("--disable-sphinx")
        tc.extra_cflags.append("-march=native");
        tc.generate()

    def build(self):
        autotools = Autotools(self)
        autotools.configure()
        autotools.make()

    def package(self):
        self.output.info("package_folder := %s" % self.package_folder)
        self.output.info("build_folder   := %s" % self.build_folder)
        self.output.info("source_folder  := %s" % self.source_folder)
        autotools = Autotools(self)
        autotools.install(args=["DESTDIR="])
    
    def package_id(self):
        del self.info.settings.compiler
        del self.info.settings.build_type

    def package_info(self):
        self.cpp_info.libs = ['mpi', 'open-rte', 'open-pal']
        if self.settings.os == "Linux":
            self.cpp_info.system_libs = ["dl", "pthread", "rt", "util"]
        self.output.info("Creating MPI_HOME environment variable: {}".format(self.package_folder))
        self.runenv_info.define_path("MPI_HOME", self.package_folder)
        self.output.info("Creating OPAL_PREFIX environment variable: {}".format(self.package_folder))
        self.runenv_info.define_path("OPAL_PREFIX", self.package_folder)
        mpi_bin = os.path.join(self.package_folder, 'bin')
        self.output.info("Creating CC environment variable: {}/mpicc".format(mpi_bin))
        self.runenv_info.append("CC", "{}/mpicc".format(mpi_bin))
        self.output.info("Creating CXX environment variable: {}/mpicxx".format(mpi_bin))
        self.runenv_info.append( "CXX", "{}/mpicxx".format(mpi_bin))
        self.output.info("Creating FC environment variable: {}/mpif90".format(mpi_bin))
        self.runenv_info.append("FC", "{}/mpif90".format(mpi_bin))


