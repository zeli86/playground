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
from conan.tools.files import get
from conan.tools.build import build_jobs
import os

required_conan_version = ">=2.0"

class boost_recipe(ConanFile):
    name = "boost"
    version = "1.82.0"
    user = "atus"
    channel = "stable"
    settings = "os", "compiler", "build_type", "arch"
    homepage = "https://www.boost.org/"
    topics = ("library", "C++")
    license = "BOOST"

    def source(self):
        get(self, **self.conan_data["sources"][self.name][self.version])

    # @property
    # def _toolset(self):
    #     return "gcc-12"

    @property
    def _build_flags(self):
        flags = []
        flags.append("install")
        if self.settings.build_type == "Debug":
            flags.append("variant=debug")
        else:
            flags.append("variant=release")
        if 'arm' in self.settings.arch:
            flags.append("architecture=arm")
        else: 
            flags.append("architecture=x86")
        flags.append("address-model=64")
        flags.append("link=shared")
        flags.append("runtime-link=shared")
        flags.append("threading=multi")
        #flags.append("toolset=%s" % self._toolset)
        flags.append("-q") # Stop at the first error. No need to continue building.
        flags.append("--disable-icu")
        flags.append("--disable-iconv")
        flags.append("--without-mpi")
        flags.append("--without-python")
        flags.append("--no-cmake-config")

        flags.extend([
            "--prefix=%s" % self.package_folder,
            "-j%s" % build_jobs(self),
            "--abbreviate-paths"
        ])
        return flags

    def build(self):
        b2_flags = " ".join(self._build_flags)
        self.output.info( "build_flags = %s" % b2_flags )
        self.run("./bootstrap.sh")
        self.run("./b2 %s" % b2_flags )

    def package(self):
        pass

    def package_info(self):
        self.buildenv_info.define("Boost_DIR", self.package_folder)
        self.runenv_info.define("Boost_DIR", self.package_folder)
        self.runenv_info.prepend_path("LD_LIBRARY_PATH", os.path.join(self.package_folder, 'lib'))

