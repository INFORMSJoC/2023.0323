[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Exact Simulation of Quadratic Intensity Models

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The codes in this repository are a snapshot of the codes
that were used in the research reported on in the paper 
[Exact Simulation of Quadratic Intensity Models](https://doi.org/10.1287/ijoc.2023.0323) by Y. Qu, A. Dassios, A. Liu and H. Zhao. 


## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0323

https://doi.org/10.1287/ijoc.2023.0323.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{ExactSimulation,
  author =        {Qu, Yan and Dassios, Angelos and Liu, Anxin and Zhao, Hongbiao},
  publisher =     {INFORMS Journal on Computing},
  title =         {Exact Simulation of Quadratic Intensity Models},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0323.cd},
  url =           {https://github.com/INFORMSJoC/2023.0323},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0323},
}  
```

## Description

The codes can be used to generate all experimental results in this paper. All 
algorithms (our exact scheme) proposed in this paper are written in the form of functions 
based on R or MATLAB. We also provide the codes to compare the numerical experiment results by 
our exact scheme with discretisation scheme and projection scheme. 
"R_Simulation_Quadratic_OU_N.R" is the main code, and includes all algorithms proposed in our paper, and most numerical 
experiment results can be plotted by executing this code. 


"R_Simulation_Quadratic_OU_N.R" is the R function for calculating and summarising the numerical experiment results by 
	our exact scheme.

"rlaptrans.R" is the R function for calculating inverse Laplace transform as Ridout (2009), downloaded from https://www.kent.ac.uk/smsas/personal/msr/webfiles/software.html.


"M_Discretisation.m" is the MATLAB code for calculating and summarising the numerical experiment results by discretisation 
	scheme.

"M_Projection.m" is the MATLAB code for calculating and summarising the numerical experiment results by projection scheme 
	(Giesecke Kakavand Mousavi, 2011).

"f_PMF_N.m" is the MATLAB function for calculating the PMF of point process with quadratic OU intensity.

"M_PMF_N.m" is the MATLAB code for calculating the PMF of point process with quadratic OU intensity.


## Replicating

By adjusting the parameter settings, the numerical experiment results in the paper can be obtained by executing the following codes, respectively.
 
Run "R_Simulation_Quadratic_OU_N.R", we get 

[line 30 to 263]  Figure 3 and Figure 5

[line 30 to 253, 265 to 276]  Figure 4

[line 281 to 377, 403 to 409] Table 1

[line 282 to 406, 411 to 437]  Figure 6

[line 442 to 676, 687 to 694]  Figure 1

[line 442 to 676, 688, 696 to 700]  Figure 7

[line 442 to 685, 702 to 712]  Table 2, Table 3, Table 4, Figure 8, Figure 9 and Figure 10

[line 442 to 676, 714 to 721]  Figure 11

[line 442 to 676, 715, 723 to 727]  Figure 12

[line 442 to 676, 729 to 731] Figure 13

[line 442 to 676, 733 to 737]  Figure 14.


Run "M_Discretisation.m", we get  

Table 3 and Figure 9.

Run "M_Projection.m", we get   

Table 4 and Figure 10.

Run "M_PMF_N.m", we get 

Figure 14.
 

