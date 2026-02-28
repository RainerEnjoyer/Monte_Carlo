# Bootstrap Model Selection in Finite Samples - Resample Size, Screening, and True-Model Selection

This repository contains a project developed in the context of the course "Monte Carlo Simulation Methods" at the University of Trier.

In this simulation study, we use different m-out-of-n bootstrap procedures for linear regression model selection in finite samples.  Using OOB-MSE to assess model performance, we test various resample sizes m as well as several model candidate set sizes K, reduced by screening procedures.
We analyze the effects on true-model selection rates and runtime to assess performance. In our finite sample setting, we find that screening substantially shifts the balance toward runtime gains while reducing selection accuracy.