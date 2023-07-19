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
from conan.tools.system.package_manager import Apt, Yum, PacMan, Zypper
import os, shutil, glob

required_conan_version = ">=2.0"

class slepc_recipe(ConanFile):
    name = "slepc"
    version = "3.20.1"
    user = "atus"
    channel = "stable"
    settings = "os", "compiler", "build_type", "arch"

    def requirements(self):
        for req in self.conan_data["tool_requires"]:
            self.tool_requires(req)
        for req in self.conan_data["requires"][self.name]:
            self.requires(req)

    def system_requirements(self):
        Apt(self).install(["dh-autoreconf"], update=True, check=True)
        #Yum(self).install([""])
        #PacMan(self).install([""])
        #Zypper(self).install([""])

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
        tc.configure_args.append("--with-scalapack=1")
        tc.configure_args.append("--download-arpack")
        tc.configure_args.append("--download-blopex")
        tc.configure_args.append("--download-elpa")
        tc.configure_args.append("--download-evsl")
        tc.configure_args.append("--download-hpddm")
        tc.configure_args.append("--download-primme")
        tc.configure_args.append("--download-trlan")
        env = tc.environment()
        env.define("PETSC_DIR", self.dependencies["petsc"].package_folder)
        self.output.info("petsc package_folder := %s" % self.dependencies["petsc"].package_folder)
        tc.make_args.append('SLEPC_DIR={}'.format(self.build_folder))
        tc.make_args.append('PETSC_DIR={}'.format(self.dependencies["petsc"].package_folder))
        tc.generate(env)

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
