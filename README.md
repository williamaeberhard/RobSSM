RobSSM: R Code for Robust Estimation for State Space Models
-----------------------------------------------------------

These separate R files provide functions to simulate data according to the North Sea pollock stock assessment state space model (SSM) and fit the model according to (Laplace-approximated) maximum likelihood (ML) and robust estimation methods. The RobSSM_Main.r script guides through all functions and generates plots similar to those in the paper. Estimation relies on the R package Template Model Builder (TMB). Details about the robust methodology and theoretical properties can be found in Aeberhard et al. (2020).

Updates can be found at https://github.com/williamaeberhard/robssm.

Any requests/comments/bug reports should be sent to william.aeberhard@gmail.com.

### Contents

Files contained in this repository:

* RobSSM_Main.r: main R file that goes through all functions, simulates data, contaminates them, and compares ML and robust outputs;
* NP_nst_gen.r: generates data according to the non-stationary (nst) SSM used for the assessment of pollock in the North Sea;
* NP_nst_nrcorrect.r: Newton-Raphson correction of the robustified gradient by Monte Carlo approximation of the expected gradient;
* NP_nst_transfo.r: transformation of model parameters, to feed optimization;
* NP_nst_untransfo.r: un-transformation of model parameters, after optimization;
* NP_nst_weights.r: robustness weights according to smooth semi-Huber or log-logistic rho function;
* NP_nst.cpp: C++ TMP template for nst SSM, for both ML and robust estimation;
* NP_st_gen.r: generate data according to the stationary (st) SSM used for the assessment of pollock in the North Sea;
* NP_st_nrcorrect.r: Newton-Raphson correction of the robustified gradient by Monte Carlo approximation of the expected gradient;
* NP_st_stdist.r: computes parameters (means and variances) of stationary distribution, used for initial conditions;
* NP_st_transfo.r: transformation of model parameters, to feed optimization;
* NP_st_untransfo.r: un-transformatio of model parameters, after optimization;
* NP_st_weights.r: robustness weights according to smooth semi-Huber or log-logistic rho function;
* NP_st.cpp: C++ TMP template for st SSM, for both ML and robust estimation;
* this README file.

### Version History

This is RobSSM version 0.1. This is the initial release.

### References

Aeberhard, W. H., Cantoni, E., Field, C. Künsch, H. R., Mills Flemming, J., and Xu, X. (2020) Robust Estimation for Discrete-Time State Space Models. Scandinavian Journal of Statistics, In Press. DOI: 10.1111/sjos.12482. Preprint: https://arxiv.org/abs/2004.05023
