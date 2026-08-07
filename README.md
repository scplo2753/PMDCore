# PMDCore

An experimental C++ implementation of the core PMD analysis workflow for ancient DNA.

This project started as a personal attempt to understand and reimplement the
PMD scoring method described in the PMDtools publication.

The underlying method is based on published research, while this particular
implementation should still be considered experimental. It implements the
core alignment reconstruction, filtering and PMD likelihood workflow, but does
not aim to reproduce every option or behavior of the original PMDtools program.

The project is mainly developed for learning, validation and implementation
experiments, including multithreading and modern C++ design.

## Build Instructions

### Distribution dependencies

The third-party C++ libraries are included as Git submodules. Install the
native compiler/build tools for your distribution before configuring:

| Distribution family | Packages |
| --- | --- |
| Debian / Ubuntu | `build-essential cmake git` |
| Fedora / RHEL / Rocky / AlmaLinux | `gcc-c++ cmake git make` |
| Arch Linux / Manjaro | `base-devel cmake git` |
| openSUSE | `gcc-c++ cmake git make` |

No system Matplot++ or argparse development package is required because the
repository builds its pinned submodules. CMake 3.15 or newer and a compiler
with C++20 support are required.

1. After cloning the repository, initialize and update submodules to fetch external dependencies:

   ```bash
   git submodule init
   git submodule update --recursive
   ```

2. For development, configure and build in a separate Debug directory:

   ```bash
   cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug
   cmake --build build-debug --parallel
   ```

   After compilation, run the program directly from the build directory:

   ```bash
   ./build-debug/PMDCore --help
   ```

3. The project supports statically linking its C++ and vendored dependencies.
   This still uses the platform C library dynamically. Enable it with:

   ```bash
   cmake -S . -B build-static \
     -DCMAKE_BUILD_TYPE=Release \
     -DPMDCORE_STATIC_DEPS=ON
   cmake --build build-static --parallel
   ./build-static/PMDCore --help
   ```

4. For a regular Release build, use:

   ```bash
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build --parallel
   ```

   The executable is then available at `build/PMDCore` and can be run without
   installation or administrator privileges:

   ```bash
   ./build/PMDCore --help
   ```

## Cite

P. Skoglund, B.H. Northoff, M.V. Shunkov, A.P. Derevianko, S. Pääbo, J. Krause, & M. Jakobsson, Separating endogenous ancient DNA from modern day contamination in a Siberian Neandertal, Proc. Natl. Acad. Sci. U.S.A. 111 (6) 2229-2234, https://doi.org/10.1073/pnas.1318934111 (2014).
