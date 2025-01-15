# iono_current_solver
This repository stores the codes and procedures I used to calculate the ionospheric current based on input FAC.

GITM has been developed in Fortran-90. Original code and copyright by University of Michigan. Please refer to the paper by Ridley, Deng & Tóth (2006) at https://doi.org/10.1016/j.jastp.2006.01.008 and the GITM user manual at https://drive.google.com/file/d/14eLZuaxlNwpKO4sl0Ig7FmygDQa83Cvf/view?usp=sharing

Other useful links:

TACC LS6 User Guide: https://docs.tacc.utexas.edu/hpc/lonestar6/

## HPC Environments & Dependencies:

1. GITM runs on TACC machine (you can use different shells, e.g., bash, csh. I use bash myself). 
2. GITM needs MPI to work.

## Quick Start:

1\. Clone the repository on your TACC Home directory

```shell
git clone https://github.com/yuho-yuho/gitm_default.git
```

$J_{Horizontal}$
