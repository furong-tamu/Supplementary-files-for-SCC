function[x,y,epsilon]=SCC_T_simulation_spatial_data_prescribe(beta,epsilon_mag,sim_num,coef_linear,lon,lat,phi,phi_ep)

% This funcution is used to generate simulation data for the case with
% hybrid coefficients

% Input:
% beta: [n,p] matrix storing the true regression coefficients with the last
%       column corresponds to the spatially varying intercept
% epsilon_mag: standard deviation (sigma) for spatially correlated noise
% sim_num: the number of simulations
% coef_linear: correlation coeffcient between the first and second covariates
% lon: [n,1] vector containing the longitudes of locations
% lat: [n,1] vector containing the latitudes of locations
% phi: range parameter for the covariates
% phi_ep: range parameter for the spatially correlated noise
% x: [sim_num,n,p] arrary storing the covariates
% y: [sim_num,n] matrix storing the response variable




[n,p]=size(beta);
y=nan(sim_num,n);
x=nan(sim_num,n,p);
epsilon=nan(sim_num,n);

[n,p]=size(beta);



 dx=nan(n,n);
 dy=nan(n,n);
 for i=1:n
     for j=1:n
         dx(i,j)=abs(lon(i)-lon(j));
         dy(i,j)=abs(lat(i)-lat(j));
     end
 end


for t=1:sim_num
    sigma2_mat=exp(-sqrt(dx.^2+dy.^2)/phi_ep);
    epsilon(t,:)=epsilon_mag*mvnrnd(zeros(n,1),sigma2_mat)+0.1*randn(1,n);
  
    for k=1:p
        if phi==0
            x(t,:,k)=randn(n,1);
        else
            sigma2_mat=exp(-sqrt(dx.^2/phi^2+dy.^2/phi^2));
            x(t,:,k)=mvnrnd(zeros(n,1),sigma2_mat);
        end
    end
    % colinearity
    x(t,:,2)=squeeze(x(t,:,1)*coef_linear)+x(t,:,2)*sqrt(1-coef_linear^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y(t,:)=sum(beta.*squeeze(x(t,:,:)),2)'+epsilon(t,:)+0.1*randn(1,n);
end
