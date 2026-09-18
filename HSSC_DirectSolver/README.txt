===================================================================================================
==                                                                                               ==
==                                        README                                                 ==
== A. Szymczak                                                                                   ==
===================================================================================================
Copyright (C) 2011-2013 E. Corona, P.G. Martinsson, D. Zorin



Program Name        :   HSSC_DirectSolver
Copyright License   :   GNU General Public License
Download            :    
Contact             :   corona at cims.nyu.edu


See <CMZ2012_paper.pdf> for algorithmic details.
See <INSTALL.txt> for installation instructions.
See <DOCUMENTATION_HSS1D.txt> for HSS1D function descriptions.
See <DOCUMENTATION_HSS2D.txt> for HSS2D function descriptions.
See <QUICK_START.pdf> for implementation and usage details.


HSSC_DirectSolver is a MATLAB program that can be used to:
    
    -   Compute the HSS factorization of a matrix. This can be used to compress matrices with 
            certain low-rank properties (e.g. low-rank off-diagonal blocks).
    -   Perform basic HSS matrix operations, such as sum, matrix-vector apply, invert, scalar 
            multiply, and transpose.
    -   Compute an interpolative decomposition of a matrix.
    -   Solve the linear system arising from the discretization of integral equations on planar 
            curves and surfaces in time linear with respect to the number DOFs.


The code has the following directory tree.

    HSSC_DirectSolver/
        Common/
        HSS1D/
            Test/
        HSS2D/
            Test/

The HSS functions are split into three folders. The folder HSS1D/ contains the basic HSS 
matrix operations (e.g. compression, sum, inversion, etc.) as well as the functions needed for 
solving integral equations on planar curves. The folder HSS2D/ contains the HSS compressed-block 
functions needed for solving integral equations on planar surfaces. The Common/ folder has helper 
functions used for both HSS and compressed-block HSS (HSS1D and HSS2D). These include 
interpolative decomposition and kernel evaluation routines, among others. The Test/ folders 
located inside HSS1D/ and HSS2D/ contain routines used to test the error of the HSS algorithms. 
They also provide good examples on how to run the code.








