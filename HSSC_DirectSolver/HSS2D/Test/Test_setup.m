% Test setup adds the necessary file paths to run other tests and calls the
% routine box_constants, defining the global variable BOX used by the 1D
% HSS functions. 
addpath ../
addpath ../../HSS1D/
addpath ../../Common/
warning off; 
box_constants;
