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
s
## Cite

P. Skoglund, B.H. Northoff, M.V. Shunkov, A.P. Derevianko, S. Pääbo, J. Krause, & M. Jakobsson, Separating endogenous ancient DNA from modern day contamination in a Siberian Neandertal, Proc. Natl. Acad. Sci. U.S.A. 111 (6) 2229-2234, https://doi.org/10.1073/pnas.1318934111 (2014).