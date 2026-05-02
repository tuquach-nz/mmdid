# mmdid

This code includes implementations in both Stata and R. For empirical applications, we use Stata, while simulations can be conducted in either Stata or R. The file mmdid.ado is the required package for running the empirical application, including generating plots. The commands are as follows.

xi: mmdid depvar indepvar, tv(treatment variable) ev(event time variable) att(heterogenous cluster variable)
mmdid_plot
