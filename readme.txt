The codes in this folder can be used to generate all experimental results covered in this article. All 
algorithms (our exact scheme) proposed in this paper are written in the form of functions 
based on R or MATLAB. We also provide the codes to compare the numerical experiment results by 
our exact scheme with discretisation scheme and projection scheme.


We first introduce the code files in this folder:

"rlaptrans.R"  --  R function for calculating inverse Laplace transform, download from Ridout (2009)

"R_Simulation_Quadratic_OU_N.R"  --  calculating and summarising the numerical experiment results by 
	our exact scheme

"M_Discretisation.m"  --  calculating and summarising the numerical experiment results by discretisation 
	scheme

"M_Projection.m"  --  calculating and summarising the numerical experiment results by projection scheme 
	(Giesecke, 2017)

"f_PMF_N.m"  --  MATLAB function for calculating the PMF of point process with quadratic OU intensity

"M_PMF_N.m"  --  calculating the PMF of point process with quadratic OU intensity


"R_Simulation_Quadratic_OU_N.R" includes all algorithms proposed in our paper, and most numerical 
experiment results can be plotted by executing its codes. By adjusting the parameter settings, the 
numerical experiment results in the paper can be obtained by executing the following codes, respectively.

[line 30 to 263]  --  Figure 3 + Figure 5

[line 30 to 253, 265 to 276]  --  Figure 4

[line 281 to 377, 403 to 409]  -- Table 1

[line 282 to 406, 411 to 437]  --  Figure 6

[line 442 to 676, 687 to 694]  --  Figure 1

[line 442 to 676, 688, 696 to 700]  --  Figure 7

[line 442 to 685, 702 to 712]  --  Table 2 + Table 3 + Table 4 + Figure 8 + Figure 9 + Figure 10

[line 442 to 676, 714 to 721]  --  Figure 11

[line 442 to 676, 715, 723 to 727]  --  Figure 12

[line 442 to 676, 729 to 731]  --  Figure 13

[line 442 to 676, 733 to 737]  --  Figure 14


Except for "R_Simulation_Quadratic_OU_N.R", the rest code files correspond to tables and figures as 
follows, respectively.

"M_Discretisation.m"  --  Table 3 + Figure 9

"M_Projection.m"  --  Table 4 + Figure 10

"M_PMF_N.m"  --  Figure 14