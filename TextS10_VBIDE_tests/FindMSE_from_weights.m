function [Error,Toterror] = FindMSE(s,weights,field,Basis,intensity)

J = length(s);
s_sub = linspace(s(1),s(end),25);
[S1_sub,S2_sub] = meshgrid(s_sub);
[s1,s2] = meshgrid(s);
bias = -2;

N = size(weights,2);
Error = zeros(N,1);
Toterror =zeros(J,J,N);
for i = 1:N
    Field_original = @(s1,s2) interp2(S1_sub,S2_sub,field(:,:,i),s1,s2);
    mean_field_est = 0; 
    for j = 1:Basis.nx 
        mean_field_est = mean_field_est + weights(j,i)*Basis.phi(:,:,j); 
    end; 
    if intensity == 'T'
    Error(i) = sum(sum(((exp(bias + Field_original(s1,s2)) - exp(bias + mean_field_est)).^2)))/J^2;
    Toterror(:,:,i) = ((exp(bias + Field_original(s1,s2)) - exp(bias + mean_field_est)));
    else
       Error(i) = sum(sum(((exp(Field_original(s1,s2)) - exp(mean_field_est)).^2)))/J^2;
       Toterror(:,:,i) = ((exp(Field_original(s1,s2)) - exp(mean_field_est)));
    end
end

