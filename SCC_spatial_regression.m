function[beta_hat,MSE]=SCC_spatial_regression(x,y,lon,lat,beta,sim_num,options)

% Description: Perform simulation studies based on the SCC method

% Input ARGUMENTS
% sim_num: number of simulations 
% x:       [sim_num,n,p] array storing the p explanatory variables on n locations
% y:       [sim_num,n] matrix storing the response variables on n locations
% lon:     [n,1] matrix storing the longitudes of n locations
% lat:     [n,1] matrix storing the latitudes of n locations
% beta:    [n,p] matrix storing the true values of regression coefficients
% options: a structure setting the parameters used in the SCC
%          see function SCC.m for details

% OUTPUT ARGUMENTS
% beta_hat:[sim_num,n,p] array storing the estimates of regression
%          coefficient
% MSE:     mean square error of estimation  


[n,p]=size(beta);

beta_hat=nan(sim_num,n,p);

MSE=nan(sim_num,p);
   

for t=1:sim_num
    t
    [aaa,BIC]=SCC(squeeze(x(t,:,:)),squeeze(y(t,:))',lon,lat,options);
    [~,index]=min(BIC);
    beta_hat(t,:,:)=reshape(aaa(:,index),[n,p]);
    MSE(t,:)=mean((squeeze(beta_hat(t,:,:))-beta).^2);
end

    