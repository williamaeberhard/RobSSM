RobSSM: R Code for Robust Estimation of State Space Models
----------------------------------------------------------

These separate R scripts provide functions to simulate data according to the North Sea pollock stock assessment state space model and fit the model according to (Laplace-approximated) maximum likelihood (ML) and robust estimation methods. Estimation relies on the R package Template Model Builder (TMB). Details about the methodology and theoretical properties can be found in Aeberhard et al. (2018).

Updates can be found at https://github.com/williamaeberhard/robssm.

Any requests/comments/bug reports should be sent to william.aeberhard@gmail.com.

### Contents

Files contained in this repository:

* RobSSM_Main.r: main R file that goes through all functions, simulates data, contaminates it, and compares ML and robust outputs TODO: finish
* NP_nst_gen.r: generate data according to the non-stationary (nst) state space model (SSM) used for the assessment of pollock in the North Sea;
* NP_nst_nrcorrect.r: Newton-Raphson correction of the robustified gradient by Monte Carlo approximation of hte expected gradient;
* NP_nst_transfo.r: transformation of model parameters, to feed optimization;
* NP_nst_untransfo.r: un-transformatio of model parameters, after optimization;
* NP_nst_weights.r: robustness weights according to smooth semi-Huber or log-logistic rho function;
* NP_nst.cpp: C++ TMP template for nst SSM, for both ML and robust estimation;
* NP_st_gen.r: generate data according to the non-stationary (nst) state space model (SSM) used for the assessment of pollock in the North Sea;
* NP_st_nrcorrect.r: Newton-Raphson correction of the robustified gradient by Monte Carlo approximation of hte expected gradient;
* NP_st_stdist.r: computes parameters (means and variances) of stationary distribution;
* NP_st_transfo.r: transformation of model parameters, to feed optimization;
* NP_st_untransfo.r: un-transformatio of model parameters, after optimization;
* NP_st_weights.r: robustness weights according to smooth semi-Huber or log-logistic rho function;
* NP_st.cpp: C++ TMP template for st SSM, for both ML and robust estimation;
* this README file.

### Version History

This is RobSSM version 0.1. This is the initial release.

### References

Aeberhard, W. H., Cantoni, E., Field, C. Künsch, H. R., Mills Flemming, J., and Xu, X. (2018) Robust Estimation for General State Space Models. Submitted.



