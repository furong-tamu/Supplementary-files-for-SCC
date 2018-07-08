function[RI]=SCC_RI_loop(beta,beta_hat)
% Description: Compute the Rand index (RI) for the estimates

% INPUT ARGUMENTS
% beta:    [n,p] matrix storing the true value of regression coefficients  
%          of p covariates on n locations 
% beta_hat:[s,n,p] array storing the estimates of beta in s simulations 

% OUTPUT ARGUMENTS
% RI:      [p,1] vector storing the Rand index for each covariate averaged 
%          over s simulations        
[s,~,p]=size(beta_hat);
RI=nan(p,s);
for t=1:s
    RI(:,t)=SCC_RI(beta,squeeze(beta_hat(t,:,:)));
end
RI=mean (RI,2);

end

%Subroutine SCC_RI
function[RI]=SCC_RI(beta,beta_hat)
[n,p]=size(beta);
mask = tril(true(n,n),-1);
RI=nan(p,1);

for j=1:p
    aaa=beta(:,j)*ones(1,n);
    aaa=aaa-aaa';
    aaa(abs(aaa)~=0)=1;
    bbb=beta_hat(:,j)*ones(1,n);
    bbb=bbb-bbb';    
    bbb(abs(bbb)>1e-6)=1;%using a small value to account for computational precision
    bbb(abs(bbb)<1e-6)=0;%using a small value to account for computational precision
    aaa=aaa(mask);
    bbb=bbb(mask);
    st=sum(aaa-bbb==0);
    RI(j)=st/length(aaa);
    
end

end
%--------------------------------------------------------------------------