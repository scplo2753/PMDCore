# PMDCore

An experimental C++ implementation of the core PMD analysis workflow for ancient DNA.

This is an experimental reimplementation and is not yet intended as a drop-in replacement for every PMDtools workflow.

PMDCore is a C++ reimplementation of the core PMDtools workflow, developed to investigate reproducibility,
performance, and compatibility with modern alignment data.

Output has been compared against the original implementation in 4 validation tests using `--platypus` and `--requirebaseq 30`.
Unit tests cover MD/CIGAR parsing and PMD likelihood calculations.
In an internal benchmark, the same workload decreased from approximately 4 days to 3-5 minutes.

The underlying method is based on published research, while this particular
implementation should still be considered experimental. It implements the
core alignment reconstruction, filtering and PMD likelihood workflow, but does
not aim to reproduce every option or behavior of the original PMDtools program.

The project is mainly developed for learning, validation and implementation
experiments, including multithreading and modern C++ design.

Known behavioral differences from the original PMDtools implementation and
questions that still require validation are documented in
[Open Questions and Compatibility Notes](OPEN_QUESTIONS.md).

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

## Acknowledgements
PMDCore is based on the algorithms and methodology introduced by the original PMDtools project.
I would like to acknowledge the authors and contributors of PMDtools for making their implementation and research publicly available,
which made this reimplementation and validation work possible.

The project also relies on several open-source libraries and tools. Their respective authors and maintainers are acknowledged through the project documentation and license notices.

## AI-assisted development
Generative AI tools, including ChatGPT and the Visual Studio Code Codex extension,
were used during the development of PMDCore as development assistants.

Their use included:

* discussing and evaluating implementation approaches;
* reviewing code and identifying potential implementation issues;
* assisting with the initial design of the thread-pool implementation;
* discussing performance optimization and software architecture;
* assisting with documentation and README drafting and revision.

AI-generated suggestions were not treated as authoritative. Design decisions, generated or suggested code, algorithmic behavior,
and performance-related changes were reviewed, integrated, tested, and validated by the maintainer.

The scientific algorithms and expected behavior of PMDCore are derived from the original PMDtools implementation and its associated publications rather than from generative AI output.

This README was also prepared with AI assistance and subsequently reviewed and edited by the maintainer.

## Cite

P. Skoglund, B.H. Northoff, M.V. Shunkov, A.P. Derevianko, S. Pääbo, J. Krause, & M. Jakobsson, Separating endogenous ancient DNA from modern day contamination in a Siberian Neandertal, Proc. Natl. Acad. Sci. U.S.A. 111 (6) 2229-2234, https://doi.org/10.1073/pnas.1318934111 (2014).
