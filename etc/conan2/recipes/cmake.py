# MIT License
# 
# Copyright (c) 2023 Zelimir Marojevic (zelimir.marojevic@gmail.com)
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
from conan.tools.files import get
from conan.tools.system.package_manager import Apt, Yum, PacMan, Zypper
import os

required_conan_version = ">=2.0"

class cmake_recipe(ConanFile):
    name = "cmake"
    version = "3.28.1"
    user = "atus"
    channel = "stable"
    homepage = "https://cmake.org/"
    settings = "os", "compiler", "build_type", "arch"

    def system_requirements(self):
        Apt(self).install(["libssl-dev"], update=True, check=True)
        Zypper(self).install(["libopenssl-3-devel"], update=True, check=True)

    def source(self):
        get(self, **self.conan_data["sources"][self.name][self.version])

    def build(self):
        self.output.info("build_folder := %s" % self.build_folder)
        self.output.info("package_folder := %s" % self.package_folder)
        self.run(f"./bootstrap --prefix={self.package_folder}")
        self.run(f"make -j {os.cpu_count()}")
        self.run("make install")

    def package_id(self):
        del self.info.settings.compiler
        del self.info.settings.build_type

    def package_info(self):
        self.buildenv_info.append_path("PATH", os.path.join(self.package_folder, 'bin'))
        self.runenv_info.append_path("PATH", os.path.join(self.package_folder, 'bin'))

