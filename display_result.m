function display_result(ev_mean_gPC,ev_mean_QMC,ev_var_QMC,ev_var_gPC)
% Display mean and variance for each eigenvalue computed with gPC or QMC

n=size(ev_mean_gPC,1);

fprintf('\n')
disp('Comparison of mean and variance for each eigenvalue')
disp('computed with gPC or quasi Monte-Carlo (qMC)')
disp('-------------------------------------------------------------------------------')
disp('ev n°        Mean qMC                 Mean gPC          Var qMC         Var gPC')
for i=1:n
    fprintf('%3g \t %3.3f + %3.3f i \t %3.3f + %3.3f i \t %3.3f \t \t %3.3f\n',...
        i,real(ev_mean_QMC(i)),imag(ev_mean_QMC(i)),...
        real(ev_mean_gPC(i)),imag(ev_mean_gPC(i)),ev_var_QMC(i),ev_var_gPC(i))

end

end