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
from conan.tools.gnu import AutotoolsToolchain, Autotools
from conan.tools.files import get
import os

required_conan_version = ">=2.0"

class gsl_recipe(ConanFile):
    name = "gsl"
    version = "2.7.1"
    user = "atus"
    channel = "stable"
    homepage = "https://www.gnu.org/software/gsl/"
    topics = ("library", "numerical computing")
    description = "A numerical library for C, C++"
    license = "GPLv3"
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fpic": [True, False]}
    default_options = {"shared": True, "fpic": True }

    def source(self):
        get(self, **self.conan_data["sources"][self.name][self.version])

    def generate(self):
        tc = AutotoolsToolchain(self)
        tc.update_configure_args({"--prefix": f"{self.package_folder}",
                                  "--libdir": f"{os.path.join(self.package_folder, 'lib')}",
                                  "--bindir": None, 
                                  "--sbindir": None,
                                  "--includedir": None, 
                                  "--oldincludedir": None});
        tc.extra_cflags.append("-march=native");
        tc.generate()

    def build(self):
        autotools = Autotools(self)
        autotools.configure()
        autotools.make()

    def package(self):
        self.output.info("package_folder := %s" % self.package_folder)
        autotools = Autotools(self)
        autotools.install(args=["DESTDIR="])

    def package_info(self):
        self.runenv_info.prepend_path("PATH", os.path.join(self.package_folder, 'bin'))
        self.runenv_info.prepend_path("LD_LIBRARY_PATH", os.path.join(self.package_folder, 'lib'))