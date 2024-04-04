# BigClust - Stochastic Nonsmooth Optimization based Incremental Clustering Software (version 0.1) 

BigClust is a nonsmooth optimization based clustering algorithm for solving the minimum sum-of-squares clustering (MSSC) problem in very large-scale and big data sets. BigClust consist of two different algorithms: an incremental algorithm is used to solve clustering problems globally and at each iteration of this algorithm the stochastic limited memory bundle (SLMB) algorithm is used to solve both the clustering and the auxiliary clustering problems with different starting points. In addition to the k-partition problem, BigClust solves all intermediate l-partition problems where l=1,…,k-1 due to the incremental approach used.

## Files included
* bigclust.f95          
  - Mainprogram for clustering software
* initbigclust.f95  
  - Initialization of clustering parameters and SLMB. Includes modules:
    + initclust - Initialization of parameters for clustering.
    + initslmb - Initialization of SLMB.
* clusteringmod.f95     
  - Subroutines for clustering software.
* slmb.f95              
  - SLMB - Stochastic limited memory bundle method.
* objfun.f95            
  - Computation of the cluster function and subgradients values.
* subpro.f95            
  - Subprograms for SLMB.
* parameters.f95        
  - Parameters. Inludes modules:
    + r_precision - Precision for reals,
    + param - Parameters,
    + exe_time - Execution time.

* Makefile              
  - makefile: requires a Fortran compiler (gfortran) to be installed.


## Installation and usage

To use the code:

1) Modify initbigclust.f95 as needed. The least, select the dataset, give the number of data points, features, and the maximum number of clusters "nclust".
2) Run Makefile (by typing "make"). Makefile uses gfortran as default.
3) Finally, just type "./bigclust".


The algorithm returns a txt-file with clustering function values, Dunn and Davies-Bouldin validity indices and elapsed CPU-times up to nclust clusters.
In addition, separate txt-file with the final cluster centers with nclust clusters and the solutions to all intermediate l-clustering problems with l = 1,...,nclust-1 is returned.

## References:

* BigClust and SLMB:
  - N. Karmitsa, V.-P. Eronen, M.M. Mäkelä, A. Airola, T. Pahikkala, "Stochastic limited memory bundle method for clustering in big data", 2024
* Clustering:
  - A. Bagirov, N. Karmitsa, S Taheri, "[Partitional Clustering via Nonsmooth Optimization](https://link.springer.com/book/10.1007/978-3-030-37826-4)", Springer, 2020.

* LMBM:
  - N. Haarala, K. Miettinen, M.M. Mäkelä, "[Globally Convergent Limited Memory Bundle Method for Large-Scale Nonsmooth Optimization](https://link.springer.com/article/10.1007/s10107-006-0728-2)", Mathematical Programming, Vol. 109, No. 1, pp. 181-205, 2007.
  - M. Haarala, K. Miettinen, M.M. Mäkelä, "[New Limited Memory Bundle Method for Large-Scale Nonsmooth Optimization](https://www.tandfonline.com/doi/abs/10.1080/10556780410001689225)", Optimization Methods and Software, Vol. 19, No. 6, pp. 673-692, 2004.
* Nonsmooth optimization:
  - A. Bagirov, N. Karmitsa, M.M. Mäkelä, "[Introduction to nonsmooth optimization: theory, practice and software](https://link.springer.com/book/10.1007/978-3-319-08114-4)", Springer, 2014.

## Acknowledgements
The work was financially supported by the Research Council of Finland projects (Project No. #345804 and #345805) led by Antti Airola and Tapio Pahikkala.


   
