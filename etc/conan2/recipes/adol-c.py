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

required_conan_version = ">=2.0"

class adolc_recipe(ConanFile):
    name = "adol-c"
    version = "2.7.2"
    user = "atus"
    channel = "stable"
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fpic": [True, False]}
    default_options = {"shared": True, "fpic": True }
    requires = "boost/1.82.0@atus/stable"

    def requirements(self):
        for req in self.conan_data["tool_requires"]:
            self.tool_requires(req)
        for req in self.conan_data["requires"][self.name]:
            self.tool_requires(req)

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
        tc.configure_args.append(f"--with-boost={self.dependencies['boost'].package_folder}")
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
        self.runenv_info.prepend_path("LD_LIBRARY_PATH", os.path.join(self.package_folder, 'lib'))
