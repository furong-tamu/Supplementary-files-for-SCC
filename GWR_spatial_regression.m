function[beta_hat,MSE]=GWR_spatial_regression(x,y,lon,lat,beta,sim_num,info)

% Description: Perform simulation studies based on the GWR method

% Input ARGUMENTS
% sim_num: number of simulations 
% x:       [sim_num,n,p] array storing the p explanatory variables on n locations
% y:       [sim_num,n] matrix storing the response variables on n locations
% lon:     [n,1] matrix storing the longitudes of n locations
% lat:     [n,1] matrix storing the latitudes of n locations
% beta:    [n,p] matrix storing the true values of regression coefficients
% info:    a structure setting the parameters used in the GWR
%          see function gwr.m for details

% OUTPUT ARGUMENTS
% beta_hat:[sim_num,n,p] array storing the estimates of regression
%          coefficient
% MSE:     mean square error of estimation  


[n,p]=size(beta);
beta_hat=nan(sim_num,n,p);
MSE=nan(sim_num,p);

for t=1:sim_num
    t
    result=gwr(squeeze(y(t,:))',squeeze(x(t,:,:)),lon,lat,info);
    beta_hat(t,:,:)=result.beta;
    MSE(t,:)=mean((squeeze(beta_hat(t,:,:))-beta).^2);
end
