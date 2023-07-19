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

class hdf5_recipe(ConanFile):
    name = "hdf5"
    version = "1.14.3"
    user = "atus"
    channel = "stable"
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fpic": [True, False]}
    default_options = {"shared": True, "fpic": True }
    no_copy_source = True

    def requirements(self):
        for req in self.conan_data["tool_requires"]:
            self.tool_requires(req)

    def source(self):
        get(self, **self.conan_data["sources"][self.name][self.version])

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["HDF5_ENABLE_PARALLEL"] = "OFF"
        tc.variables["HDF5_BUILD_PARALLEL_TOOLS"] = "OFF"
        tc.variables["BUILD_TESTING"] = "OFF"
        tc.variables["HDF5_BUILD_TOOLS"] = "ON"
        tc.variables["HDF5_BUILD_CPP_LIB"] = "ON"
        tc.variables["HDF5_BUILD_FORTRAN"] = "ON"
        tc.variables["HDF5_BUILD_EXAMPLES"] = "OFF"
        tc.variables["HDF5_ENABLE_SZIP_ENCODING"] = "OFF"
        tc.variables["HDF5_ENABLE_SZIP_SUPPORT"] = "OFF"
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure(build_script_folder="hdf5-1.14.3")
        cmake.build(build_tool_args=["-j2"])

    def package(self):
        self.output.info("package_folder := %s" % self.package_folder)
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property("cmake_find_mode", "none")
        self.cpp_info.builddirs.append(os.path.join("lib", "pckconfig"))
        self.cpp_info.builddirs.append(os.path.join("cmake"))
