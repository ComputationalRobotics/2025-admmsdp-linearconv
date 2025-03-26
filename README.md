## Codes and Data for the Paper *"Local Linear Convergence of the Alternating Direction Method of Multipliers for Semidefinite Programming under Strict Complementarity"*

### Installation

There is no need to install anything except MATLAB. 

### Usage 

Please directly run `admmdyn.m` in MATLAB to test `BQP/20-1` data with $\sigma = 100$ and random (standard Gaussian) initial guess. In the end, you will get an image like this:

<p align="center">
  <img src="./admmdyn-data/BQP-r1/20-1/demo.png" width="50%">
</p>

To change to different SDP instances, please change the `prefix1` in line 10, `admmdyn.m` to other tags, such as `theta-12`. The program will automatically load `SDP_data.mat` for that tag. 
