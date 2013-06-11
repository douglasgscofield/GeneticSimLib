GeneticSimLib
=============

This is a port of John Conery's Genetic Simulation Library (GSL) to "modern" C++.  Genetic Simulation Library was part of several publications and was described in

Conery JS and M Lynch. 1999.  Genetic Simulation Library.  *Bioinformatics* 15(1): 85-86.
[doi:10.1093/bioinformatics/15.1.85](http://dx.doi.org/10.1093/bioinformatics/15.1.85)

From the original README: 

> GSL is a set of C++ classes for incorporating genetics in population models. Source code can be freely copied and distributed subject to the constraints spelled out in the copyright notice at the end of this file and in the file named COPYRIGHT in this directory.

Version 1.0 was released in 1997, and no updates have been made since that 
time.  This version is still available at John's site, <http://www.csi.uoregon.edu/projects/genetics/GSL>, as well as in the **version-1.0**
folder of this repository.

Since 1997, C++ practices have changed, and GSL wouldn't compile with at least a modern (4.0+) g++.  Since I'd like to use GSL I decided to port it and make it available here, with John's permission.  I initiated the repository with John's original source code, then added my port on top, so you can unroll back to 1997 if you'd like or if I've inadvertantly introduced a bug of some kind.

I've used the name *GeneticSimLib* as the original acronym, *GSL*, is now
closely associated with the [Gnu Scientific Library](http://www.gnu.org/software/gsl/).

Changes introduced by the port, some were required for compilation, some to silence warnings:

* Rename `file.C` to `file.cpp`
* Compile with `CPP` instead of `CC`
* Adding `std::` to names in the `std` namespace
* Removing '`.h`' from standard headers
* Parens around assignments within `if` expressions
* Reimplement some `char*` code using `std::string` (more can be done here)
* Make dtor in `RNG` and `RNGSeed` base classes virtual to ensure derived classes clean everything up
* Remove some unnecessary forward declarations of `ostream`
* Bitwise operations with `&`
* Label `#endif` wrappers


Some caveats:

* I haven't tested this port!
* Some original demos which used the XForms library don't link
* Namespaces are not yet incorporated


List of Directories
-------------------

**GeneticSimLib/**   : The ported source code

**version-1.0/**     : The version 1.0 source code


License
--------
<pre>
Copyright (c) 1997 by the University of Oregon.  
ALL RIGHTS RESERVED. 

Permission to use, copy, and distribute this software in 
its entirety for non-commercial purposes and without fee, 
is hereby granted, provided that the above copyright notice 
and this permission notice appear in all copies and their 
documentation. 

Software developers, consultants, or anyone else who wishes 
to use all or part of the software or its documentation for 
commercial purposes should contact the Technology Transfer 
Office at the University of Oregon to arrange a commercial
license agreement.

This software is provided "as is" without expressed or 
implied warranty of any kind.
</pre>
