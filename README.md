# PPLS-SVD

Replication of the simulation studies in "On some limitations of probabilistic models for dimension-reduction: illustration in the case of probabilistic formulations of Partial Least Squares", Etievant and Viallon [1].


### Required packages 

```
MASS, pracma, matrixcalc, Matrix, permute, methods, robustbase, xtable, parallel and ggplot2.
```

### Script

Script `Simul_PPLS_limitations.R` allows to replicate the simulations studies proposed by Etievant and Viallon [1]. It relies on functions provided in `helper_function.R`.


### Arguments to specify

* **dir_data_sim1** - directory where to save the true parameters values, generated data and estimates obtained via the EM algorithm devised under the PPLS model proposed by el Bouhaddani et al. [2], for the first simulation study. The directory will be used also to save the comparisons of the weight matrices A and B, as well as the figures summarizing the comparisons.

* **dir_data_sim2** - directory where to save the true parameters values, generated data and estimates obtained via the EM algorithm devised under the PPLS model proposed by el Bouhaddani et al. [2], for the second simulation study. The directory can coincide with **dir_data_sim1**.

* **Ncores** - number of cores to use to run the data generation and parameters estimation simultaneously.


### Instructions to run the script `Simul_PPLS_limitations.R`

* Save scripts `Simul_PPLS_limitations.R` and `helper_function.R` in the same directory.

* Open script `Simul_PPLS_limitations.R` and specify arguments **dir_data_sim1**, **dir_data_sim2** and **Ncores** at the top.

* Run the whole script.

### Details of the steps performed in the script `Simul_PPLS_limitations.R`

##### Data generation and parameters estimation

* Data *(X,M)* are successively generated under the PPLS model proposed by el Bouhaddani et al. [2], which is:

        t ~ N(0,Sigma_T)
        u = t D + e_U,      e_U ~ N(0, sigma_U2 I_r)
        x = t A' + e_X,     e_X ~ N(0, sigma_X2 I_pX)
        m = u B' + e_M,     e_M ~ N(0, sigma_M2 I_pM)
with 
A and B two semi-orthogonal matrices,
Sigma_T and D diagonal matrices with strictly positive diagonal elements,

and under the following misspecified PPLS model:

        t ~ N(0,Sigma_T)
        u = t D + e_U,      e_U ~ N(0, sigma_U2 I_r)
        x = t A' + e_X,     e_X ~ N(0, Psi_X)
        m = u B' + e_M,     e_M ~ N(0, Psi_M)
with Psi_X and Psi_M arbitrary positive semi-definite matrices.

* Parameters estimates returned by the EM algorithm devised under the PPLS model proposed by el Bouhaddani et al. [2] are obtained, under both the correctly specified and misspecified settings.

* True parameters values, generated datasets and estimates are saved under both settings. 

##### Comparisons of the weight matrices

* True weight matrices A and B are compared with (*i*) the estimates returned by the EM alogorithm devised under the PPLS model proposed by el Bouhaddani et al. [2], (*ii*) the estimates obtained via PLS-SVD on *(X,M)*, (*iii*) the estimates obtained via PLS-W2A on *(X,M)*, and (*iv*) the estimates obtained via two distinct PCAs on *X* and *M*, under both settings.

* Estimates returned by the EM alogorithm (devised under the PPLS model proposed by el Bouhaddani et al. [2]) are compared with (*i*) the estimates obtained via PLS-SVD on *(X,M)*, (*ii*) the estimates obtained via PLS-W2A on *(X,M)*, and (*iii*) the estimates obtained via two distinct PCAs on *X* and *M*, under both settings.

* The results of the comparisons are saved and diplayed in 3 figures, as proposed in [1].


### Default arguments used in the script `Simul_PPLS_limitations.R` to replicate the simulation studies in [1]

* **Nsim = 1000** - Number of replicates.
* **pX = 20** - dimension of the observed set of variables x. 
* **pM = 20** - dimension of the observed set of variables m. 
* **r = 3** -  dimension of the sets of latent variables t and u. in the study.
* **RATIO = 0.25** - signal to noise ratio.
* **N = c(50, 250, 500, 10^3, 5000)** - sample size(s). It could be a single value instead of a vector of values. 
* **A** - randomly chosen pX times r semi-orthogonal matrix. Used for the two simulation studies.
* **B** - randomly chosen pM times r semi-orthogonal matrix. Used for the two simulation studies.
* **Sigma_T** - diagonal r times r matrix with digonal elements (exp(−(i−1)/5))_{i = 1, ... r}.
* **D** -  diagonal r times r matrix with digonal elements (1.5exp(3(i−1)/10))_{i = 1, ... r}.
* **sigma_U2 = 5.33**
* **sigma_X2 = 0.4**
* **sigma_M2 = 4**
* **Psi_X** - randomly chosen positive semi-definite pX times pX matrix (in accordance with the value of RATIO).
* **Psi_M** - randomly chosen positive semi-definite pM times pM matrix (in accordance with the value of RATIO).



### References

[1] Etievant Lola and Viallon Vivian. On some limitations of probabilistic models for dimension-reduction: illustration in the case of probabilistic formulations of Partial Least Squares [In press]. 


[2] El Bouhaddani et al. Probabilistic partial least squares model: identifiability, estimation and application (2018). Journal of Multivariate Analysis, 167.

