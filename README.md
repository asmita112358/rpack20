# Optimal Invariant Kernel in RKHS based statistics

1. There is a fast code to compute the HSIC distance covariance. This code is expected to beat the dcov function in R wrt runtime and complexity.
2. There is a testing method based on an RKHS based statistics. It's meant to be used for goodness of fit tests, tests for independence of two multivariate random sample and Homogeneity tests.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

1. Install Anaconda and click install on Rstudio.
2. Once it's installed, close Anaconda and open it again, this time running as administrator.
3. Select the drop-down list marked "Applications on:"
4. Select rstudio in the list.
5. Additionally, install Mosek Optimization Tools, email the distributors for the liscene and store the liscene in the same folder that contains the Mosek software. (e.g for Windows it is usually the Program Files)
6. Install Rmosek package in rstudio that is launched from Anaconda.
7. Install the CVXR package. 
The CVXR package should now be able to recognize the MOSEK solver. To check, write installed_solvers on the portal and check if MOSEK comes up amongst the solvers list. If it still doesn't please contact me at asmita@stat.tamu.edu.

###Author
Asmita Roy
Some of the theoretical development behind the code have been done by Dr. Xianyang Zhang