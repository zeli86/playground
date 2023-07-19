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

class trilinos_recipe(ConanFile):
    name = "trilinos"
    version = "15.0.0"
    user = "atus"
    channel = "stable"
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fpic": [True, False]}
    default_options = {"shared": True, "fpic": True}
    no_copy_source = True

    def requirements(self):
        for req in self.conan_data["tool_requires"]:
            self.tool_requires(req)
        for req in self.conan_data["requires"][self.name]:
            self.requires(req)

    def source(self):
        get(self, **self.conan_data["sources"][self.name][self.version])

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["TPL_ENABLE_MPI"] = "ON"
        tc.variables["BUILD_SHARED_LIBS"] = "ON"
        tc.variables["Trilinos_ENABLE_Amesos"] = "ON"
        tc.variables["Trilinos_ENABLE_Epetra"] = "ON"
        tc.variables["Trilinos_ENABLE_EpetraExt"] = "ON"
        tc.variables["Trilinos_ENABLE_Ifpack"] = "ON"
        tc.variables["Trilinos_ENABLE_AztecOO"] = "ON"
        tc.variables["Trilinos_ENABLE_Sacado"] = "ON"
        tc.variables["Trilinos_ENABLE_SEACAS"] = "OFF"
        tc.variables["Trilinos_ENABLE_Teuchos"] = "ON"
        tc.variables["Trilinos_ENABLE_MueLu"] = "ON"
        tc.variables["Trilinos_ENABLE_Kokkos"] = "ON"
        tc.variables["Trilinos_ENABLE_ML"] = "ON"
        tc.variables["Trilinos_ENABLE_ROL"] = "ON"
        tc.variables["Trilinos_ENABLE_Tpetra"] = "ON"
        tc.variables["Trilinos_ENABLE_COMPLEX_DOUBLE"] = "ON"
        tc.variables["Trilinos_ENABLE_COMPLEX_FLOAT"] = "ON"
        tc.variables["Trilinos_ENABLE_Zoltan"] = "ON"
        tc.variables["Trilinos_VERBOSE_CONFIGURE"] = "OFF"
        tc.variables["TPL_ENABLE_Boost"] = "OFF"
        tc.variables["TPL_ENABLE_Netcdf"] = "OFF"
        pc_dir = os.path.join(self.dependencies["lapack"].package_folder, 'lib')
        tc.variables["BLAS_LIBRARY_NAMES"] = "blas"
        tc.variables["BLAS_LIBRARY_DIRS"] = pc_dir
        tc.variables["LAPACK_LIBRARY_NAMES"] = "lapack"
        tc.variables["LAPACK_LIBRARY_DIRS"] = pc_dir
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build(build_tool_args=["-j2"])

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property("cmake_find_mode", "none")
        self.cpp_info.builddirs.append(os.path.join("lib", "cmake"))
        self.runenv_info.prepend_path("LD_LIBRARY_PATH", os.path.join(self.package_folder, 'lib'))
