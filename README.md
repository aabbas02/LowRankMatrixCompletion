# Low Rank Matrix Completion
AltGDMin compared with three other benchmark methods for the Low Rank Matrix Completion. 
This repository contains MATLAB code to reproduce results in \url{...}

- MATLAB Parallel Computing Toolbox (https://www.mathworks.com/products/parallel-computing.html) is necesarry to reproduce the results in Figures 1(b), 1(d). The Parallel Computing toolbox is needed for the 
'parfor' command which executes the `for' loop iteration in parallel on different workers https://www.mathworks.com/help/parallel-computing/parfor.html. The code might run even without the Parallel computing toolbox installed, in which case the `parfor' loop is replaced by the standard 'for' loop, but the results will not be reproduced.

