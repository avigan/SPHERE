2016-02-15: v1.1
2015-08-23: v1.0

This archive contains the SPHERE/IFS reduction pipeline described in
the paper:

  Vigan et al., 2015, MNRAS, 454, 129

The starting point is the script "data_reduction_IFS.sh", which
contains most of the documentation and a step by step procedure to
reduce an IFS data set. Then the IDL script "data_reduction_IFS.pro"
can be used for the final steps. The pro/ subdirectory contains other
IDL routines that are used in the data reduction.

No example data is included, but public data can be obtained from the
SPHERE science verification period:

http://www.eso.org/sci/activities/vltsv/spheresv.html

The program 60.A-9382 (P.I. Lagrange) contains a very good IFS-H
sequence acquired on beta Pictoris. The data and calibrations
corresponding to this program can be downloaded from the ESO archive,
which contains a SPHERE-specific query form:

http://archive.eso.org/wdb/wdb/eso/sphere/form

Do not hesitate to contact me if you have any problem using the
reduction pipeline.

Arthur Vigan
Laboratoire d'Astrophysique de Marseille
arthur.vigan@lam.fr
2015-08-23
