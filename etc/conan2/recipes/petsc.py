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
from conan.tools.microsoft import unix_path
from conan.tools.system.package_manager import Apt, Yum, PacMan, Zypper
import os, shutil, glob

required_conan_version = ">=2.0"

class petsc_recipe(ConanFile):
    name = "petsc"
    version = "3.20.3"
    settings = "os", "compiler", "build_type", "arch"
    user = "atus"
    channel = "stable"
    
    def requirements(self):
        for req in self.conan_data["tool_requires"]:
            self.tool_requires(req)
        for req in self.conan_data["requires"][self.name]:
            self.requires(req)

    def system_requirements(self):
        # depending on the platform or the tools.system.package_manager:tool configuration
        # only one of these will be executed
        Apt(self).install(["flex"], update=True, check=True)
        Apt(self).install(["bison"], update=True, check=True)
        #Yum(self).install(["libglvnd-devel"])
        #PacMan(self).install(["libglvnd"])
        #Zypper(self).install(["Mesa-libGL-devel"])

    def source(self):
        get(self, **self.conan_data["sources"][self.name][self.version])

    def generate(self):
        tc = AutotoolsToolchain(self)
        if self.settings.build_type == "Debug":
            tc.configure_args.append("--with-debugging=1")
        else:
            tc.configure_args.append("--with-debugging=0")
        tc.update_configure_args({"--prefix": f"{self.package_folder}"});
        tc.configure_args.append("--with-mpi=1")
        tc.configure_args.append("--with-cc=mpicc")
        tc.configure_args.append("--with-cxx=mpicxx")
        tc.configure_args.append("--with-fc=mpif90")
        tc.configure_args.append("--with-shared-libraries")
        tc.configure_args.append("--with-x=0")
        tc.configure_args.append("--download-make-shared=1")
        tc.configure_args.append("--download-hypre=1")
        tc.configure_args.append("--download-scalapack=1")
        tc.configure_args.append("--download-mumps=1")
        tc.configure_args.append("--download-tetgen=1")
        tc.configure_args.append("--download-ptscotch=1")
        tc.configure_args.append("--download-suitesparse=1")
        tc.configure_args.append("--download-eigen=1")
        tc.configure_args.append("--download-superlu=1")
        tc.configure_args.append("--download-superlu_dist=1")
        tc.configure_args.append("--download-kokkos=0")
        tc.configure_args.append("--download-kokkos-kernels=0")
        pc_dir = os.path.join(self.dependencies["lapack"].package_folder, "lib", "pkgconfig")
        tc.configure_args.append("--with-blaslapack-pkg-config={}".format(pc_dir))
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


