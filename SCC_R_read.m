function[beta_hat,MSE]=SCC_R_read(filename,sim_num,beta,x,y)
% Description: This function compute the mean square error of estimation
%              from SCC-KNN and SCC-RNN methods

% INPUT ARGUMENTS
% filename: the file name storing the estimates
% sim_num:  number of simulations 
% x:        [sim_num,n,p] array storing the p explanatory variables on n locations
% y:        [sim_num,n] matrix storing the response variables on n locations
% beta:     [n,p] matrix storing the true values of regression coefficients

% OUTPUT ARGUMENTS
% beta_hat:[sim_num,n,p] array storing the estimates of regression
%          coefficient
% MSE:     mean square error of estimation 

[m,n]=size(beta);
beta_hat=nan(sim_num,m,n);
MSE=nan(sim_num,3);
s=sim_num;

for t=1:s
    filename2=strcat(filename,int2str(t),'.nc');
    beta_e=ncread(filename2,'beta');
    [~,p]=size(beta_e);
    beta_e2=nan(p,m,n);
    for k=1:p
        beta_e2(k,:,:)=reshape(beta_e(:,k),[m,n]);
    end
    [BIC]=SCC_R_admm_BIC(beta_e2,squeeze(x(t,:,:)),y(t,:)');
    [~,index]=min(BIC);
    beta_hat(t,:,:)=squeeze(beta_e2(index,:,:));
    MSE(t,:)=mean((squeeze(beta_hat(t,:,:))-beta).^2); 
end

end

% Subroutine SCC_R_admm_BIC
function[BIC]=SCC_R_admm_BIC(beta_RM,data_i,data_o)
% This functions computes the BIC ar each tuning parameter
[s,n,m]=size(beta_RM);
MSE=nan(1,s);
k=nan(1,s);
BIC=nan(1,s);
for t=1:s
    beta=squeeze(beta_RM(t,:,:));
    MSE(t)=sum((data_o-sum(beta.*data_i,2)).^2);
    
    ppp=0; 
    for i=1:m
        aaa=beta(:,i);
        [aaa]=sort(aaa);
        ds=find(abs(diff(aaa))>1e-7);
        ppp=ppp+length(ds)+1;
    end
    k(t)=ppp;
    BIC(t)=n*log(MSE(t))+k(t)*log(n);
end
end