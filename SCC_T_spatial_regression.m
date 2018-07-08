function[beta_hat,MSE]=SCC_T_spatial_regression(x,y,lon,lat,beta,phi,sim_num,options)

% Description: Perform simulation studies based on the transformed SCC method

% Input ARGUMENTS
% sim_num: number of simulations 
% x:       [sim_num,n,p] array storing the p explanatory variables on n locations
% y:       [sim_num,n] matrix storing the response variables on n locations
% lon:     [n,1] matrix storing the longitudes of n locations
% lat:     [n,1] matrix storing the latitudes of n locations
% beta:    [n,p] matrix storing the true values of regression coefficients
% phi:     the range parameter of the exponential covariance function for
%          beta_(p+1)
% options: a structure setting the parameters used in the SCC
%          see function SCC.m for details

% OUTPUT ARGUMENTS
% beta_hat:[sim_num,n,p] array storing the estimates of regression
%          coefficient
% MSE:     mean square error of estimation  

[Lc]=SCC_Cholesky_matrix(lon,lat,phi); %Compute the Cholesky matrix

[n,p]=size(beta);

beta_hat=nan(sim_num,n,p);

MSE=nan(sim_num,p);
   

for t=1:sim_num
    t
    [aaa,BIC]=SCC_T(squeeze(x(t,:,:)),squeeze(y(t,:))',lon,lat,Lc,options);
    [~,index]=min(BIC);
    beta_hat(t,:,:)=reshape(aaa(:,index),[n,p]);
    MSE(t,:)=mean((squeeze(beta_hat(t,:,:))-beta).^2);
end

end


% Subroutine SCC_Cholesky_matrix
function[Lc]=SCC_Cholesky_matrix(lon,lat,phi)

n=length(lon);
dx=nan(n,n);
dy=nan(n,n);
 for i=1:n
     for j=1:n
         dx(i,j)=abs(lon(i)-lon(j));
         dy(i,j)=abs(lat(i)-lat(j));
     end
 end
 sigma2_mat=exp(-sqrt(dx.^2+dy.^2)/phi)+eye(n)*0.01; %beta3+epsilon
 [Lc] = chol(sigma2_mat,'lower');
end
%--------------------------------------------------------------------------
    