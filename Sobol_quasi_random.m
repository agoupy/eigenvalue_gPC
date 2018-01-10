function [MC_SOBOL]=Sobol_quasi_random(N_mc,dim)

%Quasi Monte Carlo for normal parameters:
%Generate a sample for Monte Carlo realizations in order to have a
%convergence in 1/N and not 1/sqrt(N)

normal_mean = 0; normal_std = 1;
convergence_check = 0; % or 1 otherwise

Skip_value = 1e3; Leap_value = 100;
%%% Generate a nD Sobol point set, skip the first Skip_value values, and then retain every Leap_value point:
p = sobolset(dim,'Skip',Skip_value,'Leap',Leap_value);
%%% Use scramble to apply a random linear scramble combined with a random digital shift:
p = scramble(p,'MatousekAffineOwen');
%%% Use net to generate the first N_mc points:
MC_SOBOL = (norminv(net(p,N_mc),normal_mean,normal_std)); % this is the population sample to use

