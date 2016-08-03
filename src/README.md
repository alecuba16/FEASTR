FEAST
=====

A FEAture Selection Toolbox for C/C++ &amp; R, v1.1.4R.

FEAST provides implementations of common mutual information based filter
feature selection algorithms, and an implementation of RELIEF. All functions
expect discrete inputs (except RELIEF, which does not depend on the MIToolbox),
and they return the selected feature indices. These implementations were
developed to help our research into the similarities between these algorithms,
and our results are presented in the following paper:

 Conditional Likelihood Maximisation: A Unifying Framework for Information Theoretic Feature Selection
 G. Brown, A. Pocock, M.-J.Zhao, M. Lujan
 Journal of Machine Learning Research, 13:27-66 (2012)

If you use these implementations for academic research please cite the paper
above.  All FEAST code is licensed under the BSD 3-Clause License.

Contains implementations of:
   mim, mrmr, mifs, cmim, jmi, disr, cife, icap, condred, cmi, relief, fcbf, betagamma

References for these algorithms are provided in the accompanying feast.bib file
(in BibTeX format).

R Example (using "data" as our feature matrix, and "labels" as the class label vector), empty initialFeatures (id's of given preselected features)):

data<-matrix(1:6, 2)
initialFeatures<-NULL
labels<-array(1,c(2,1))
optionalParam1<-0.0
optionalParam2<-0.0

FSToolboxR("CMIM",5,matrix(1:6, 2),NULL,,0.0,0.0)



The library is written in ANSI C for compatibility with the R C interface
compiler,except FCBF and RELIEF, which are written in MATLAB/OCTAVE
script and must be ported to C.

MIToolbox is required to compile these algorithms, and these implementations
supercede the example implementations given in that package (they have more
robust behaviour when used with unexpected inputs).

MIToolbox can be found at: http://www.github.com/Craigacp/MIToolbox/

The C library expects all matrices in column-major format (i.e. Fortran style).
This is for two reasons, a) R generates Fortran-style arrays, and b)
feature selection iterates over columns rather than rows, unlike most other ML
processes. 

Compilation instructions:
 - R - use the included makefile and import the shared library and use the wrapper FSToolboxR.R
 - Linux C shared library - use the included makefile

Update History
 - 13/07/2016 - v1.1.4R - Modified to be a R external library
 - 12/03/2016 - v1.1.4 - Fixed an issue where Matlab would segfault if all features had zero MI with the label.
 - 12/10/2014 - v1.1.2 - Updated documentation to note that FEAST expects column-major matrices.
 - 11/06/2014 - v1.1.1 - Fixed an issue where MIM wasn't compiled into libFSToolbox.
 - 22/02/2014 - v1.1.0 - Bug fixes in memory allocation, added a C implementation of MIM, moved the selected feature increment into the mex code.
 - 12/02/2013 - v1.0.1 - Bug fix for 32-bit Windows MATLAB's lcc.
 - 08/11/2011 - v1.0.0 - Public Release to complement the JMLR publication.

