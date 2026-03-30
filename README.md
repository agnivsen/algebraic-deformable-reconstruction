# Convex Solutions to SfT and NRSfM under Algebraic Deformation Models
Code and dataset from our IEEE TPAMI - 2026 article 'Convex Solutions to SfT and NRSfM under Algebraic Deformation Models'

## Abstract

We present nonlinear formulations to Shape-fromTemplate (SfT) and Non-Rigid Structure-from-Motion (NRSfM) faithfully exploiting the isometric, conformal and equiareal deformation models. Existing work uses relaxations such as inextensibility or requires knowing the optic flow field around
the correspondences, an impractical assumption. In contrast, the proposed formulations only require point correspondences and resolve all ambiguities using the notions of maximal depth and maximal isometry heuristics. We propose solution methods using Semi-Definite Programming (SDP) for all formulations. We show that straightforward SDP models conflict with the usual maximal depth heuristic and propose an adapted oppositedepth parameterisation demonstrating a lesser relaxation gap. Experimental results on many real-world benchmark datasets demonstrate superior accuracy over existing methods.

## Dependencies and Setup

The code has the following dependencies:

 - [CVX](https://cvxr.com/) [1]
 - Ruled surface generator [2] _(included)_
 - [Mosek](https://www.mosek.com/) [4] 

Please install CVX following the installation instructions of the library and link CVX and your preferred solver to the code.


## How to Run

The scripts to run are grouped into four categories:

- _sftIsometric.m_: Use this to test **alg-SfT-SqMOD** and **alg-SfT-MOD**
- _sftNonIsometric.m_: Use this to test **alg-SfT-mConf** and **alg-SfT-mEqui**
- _nrsfmIsometric.m_: Use this to test **alg-NRSfM-SqMOD** and **alg-NRSfM-MOD**
- _nrsfmNonIsometric.m_: Use this to test **alg-NRSfM-mConf** and **alg-NRSfM-mEqui**

Instructions to run and configure the scripts are provided therein.

If Mosek is unavailable as solver, please change the '_solver_' variable in the four scripts mentioned above. 

All our solvers are used wrapped in CVX, [please see here](https://web.cvxr.com/cvx/beta/doc/solver.html) for a list of solvers supported by CVX.

## Dataset

All dataset used in our article are pre-existing benchmark dataset.

[Bramante39M](https://github.com/agnivsen/Bramante39M) [3] has been already open-sourced.

## How to Cite

If you find our code and dataset useful, please cite using:

    @article{sengupta2025convex,
      title={Convex Solutions to SfT and NRSfM under Algebraic Deformation Models},
      author={Sengupta, Agniva and Bartoli, Adrien},
      journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
      year={2025},
      publisher={IEEE}
    }

## References

[1]. Grant, M., & Boyd, S. (2014, March). CVX: Matlab software for disciplined convex programming, version 2.1.

[2]. Perriollat, M., & Bartoli, A. (2013). A computational model of bounded developable surfaces with application to image‐based three‐dimensional reconstruction. _Computer Animation and Virtual Worlds_, 24(5), 459-476.

[3]. Bartoli, A., & Sengupta, A. (2025). Camera pose in SfT and NRSfM under isometric and weaker deformation models. _Computer Vision and Image Understanding_, 104488.

[4]. ApS, Mosek. "Mosek optimization toolbox for matlab." User’s Guide and Reference Manual, Version 4.1 (2019): 116.
