function [c,q,quantile_MSE,quantile_bias] = FindQuantiles_from_weights(s,Estinfo,field,Basis)

% Quantile analysis
%
% Inputs:   s (vector, elements of which are the side of square)
%           Estinfo (estimated state)
%           field (true spatiotemporal field)
%           Basis (basis functions)
%
% Outputs:  c (cumulative score)
%           q (quantiles)
%           quantile_MSE (mean square error in quantile plot)
%           quantile_bias (bias in quantile plot)


J = length(s);
s_sub = linspace(s(1),s(end),25);
[S1_sub,S2_sub] = meshgrid(s_sub);
[s1,s2] = meshgrid(s);
bias = -2;
N = length(Estinfo);
q = [0.01; 0.05; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 0.95; 0.99];
mul_sigma = [-2.33; -1.645; -1.28; -0.84; -0.525; -0.255; 0; 0.255; 0.525; 0.84;  1.28; 1.645; 2.33];    
c = zeros(length(q),N);
quantile_MSE = zeros(length(q),N);
quantile_bias = zeros(length(q),N);

for i = 1:length(q)
    disp(strcat('Checking quantile ',num2str(q(i))))
    for k = 1:N
       Field_original = @(s1,s2) interp2(S1_sub,S2_sub,field(:,:,k),s1,s2);
       mean_field_est = 0; 
       var_field_est = 0;
       for j = 1:Basis.nx 
            mean_field_est = mean_field_est + Estinfo(k).xestRTS(j)*Basis.phi(:,:,j); 
            var_field_est = var_field_est + Estinfo(k).PRTS(j,j)*Basis.phi(:,:,j); 
       end; 
%        field_lower = exp(bias + mean_field_est + mul_sigma(i)*sqrt(var_field_est));
%        int_original = exp(bias + Field_original(s1,s2));
       field_lower =  mean_field_est + mul_sigma(i)*sqrt(var_field_est);
       int_original = Field_original(s1,s2);
       c(i,k) = sum(sum(field_lower > int_original))/J^2;
       quantile_MSE(i,k) = mean(mean((q(i) - (field_lower > int_original)).^2));
       quantile_bias(i,k) = mean(mean(q(i) - (field_lower > int_original)));
    end
end


