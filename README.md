# mmdid

These codes include implementations in both Stata and R for estimating heterogeneous treatment effects using multilevel regression in event-study designs. For empirical applications, we use Stata, while simulations can be conducted in either Stata or R. The file mmdid.ado is the required package for running the empirical application, including generating plots. The commands are as follows.

xi: mmdid depvar indepvar, tv(treatment variable) ev(event time variable) att(heterogenous cluster variable)

mmdid_plot

These codes are used in the manuscript name "Estimating Heterogeneous Treatment Effects through Multilevel Modeling" by Quach, D. T., Wichitaksorn, N., Nguyen, T. K., and Hawley, J. D. (2026).
