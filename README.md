# EVTools

Extreme value analysis of data above thresholds, either using a classical model (Generalized Pareto (GP) tail) or on models of the tail based on a large deviation principle: Generalized Weibull (GW), log-GW, and Weibull tails. 

The latter models allow extrapolation over a wide range of probabilities without reliance on questionable assumptions. 

FitTailFromAllData.R is the entry to making estimates with all the above models; it offers the most complete functionality. Other modules with names starting with "Fit" perform more basic tasks for specific models. They may be used to tailor the analysis to specific needs. 