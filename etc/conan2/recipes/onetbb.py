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
from conan.tools.cmake import CMakeToolchain, CMakeDeps, CMake
from conan.tools.files import get
import os

required_conan_version = ">=2.0"

class onetbb_recipe(ConanFile):
    name = "onetbb"
    version = "2021.11.0"
    settings = "os", "compiler", "build_type", "arch"
    user = "atus"
    channel = "stable"
    no_copy_source = True
    _cmake = None

    def requirements(self):
        for req in self.conan_data["tool_requires"]:
            self.tool_requires(req)
    
    def source(self):
        get(self, **self.conan_data["sources"][self.name][self.version])

    def generate(self):
        tc = CMakeToolchain(self)
        tc.generate()

    def build(self):
        args = {
            "TBB_STRICT": "OFF",
            "TBB_EXAMPLES": "OFF",
            "TBB_SANITIZE": "",
         }
        self._cmake = CMake(self)
        self._cmake.configure(variables=args)
        self._cmake.build()
        self._cmake.test()

    def package(self):
        self._cmake.install()

    def package_info(self):
        self.runenv_info.prepend_path("LD_LIBRARY_PATH", os.path.join(self.package_folder, 'lib'))
        self.cpp_info.set_property("cmake_find_mode", "none")
        self.cpp_info.builddirs.append(os.path.join(self.package_folder, 'lib', 'cmake'))



