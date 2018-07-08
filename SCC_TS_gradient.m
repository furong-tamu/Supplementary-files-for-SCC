function[F,betag,latg,depthg]=SCC_TS_gradient(latn,depthn,betan)
% Description: Compute spatial gradient of T-S relationship. 
%              No interpolation is invloved in the compitation of spatial
%              gradient.

% INPUT ARGUMENTS:
% latn:    a vector storing nondimensional coordinate in the horizontal 
%          direction at n locations
% depthn:  a vector storing nondimensional coordinate in the vertical 
%          direction at n locations
% betan:   a vector storing T-S relationship at n locations        

% OUTPUT ARGUMENTS:
% F:       magnitude of spatial gradient of T-S relationship at n locations
% betag:   T-S relationship interpolated onto regular grids using nearest 
%          neighbor interpolation
% latg:    nondimensional coordinate in the horizonal direction for regular
%          grids
% depthg:  nondimensional coordinate in the vertical direction for regular
%          grids

n=length(betan);
F=nan(n,1);
for i=1:n
    [F(i)]=SCC_local_gradient(i,betan,latn,depthn);
end


latg=[-0.5:0.005:0.5];
depthg=[0:0.005:1];
[latg,depthg]=meshgrid(latg,depthg);
latg=latg';
depthg=depthg';
betag=griddata(latn,depthn,betan,latg,depthg,'nearest');

end


% Subroutine SCC_local_gradient
function[F]=SCC_local_gradient(i,beta,latn,depthn)
% This subroutine is used to compute the spatial gradient at location i.

dh=latn-latn(i);
dv=depthn-depthn(i);
dh(i)=nan;
dv(i)=nan;
d=sqrt(dh.^2+dv.^2);

[~,index1]=nanmin(d);
dh1=dh(index1);
dv1=dv(index1);
theta1=atan2(dv1,dh1);


d(index1)=nan;
parr=0;
while parr==0
    [~,index2]=nanmin(d);
    dh2=dh(index2);
    dv2=dv(index2);
    theta2=atan2(dv2,dh2);
    if abs(theta1-theta2)>pi/6&&abs(abs(theta1-theta2)-pi)>pi/6
        parr=1;
    else
        d(index2)=nan;
    end
end

Fs1=(beta(index1)-beta(i))/sqrt(dh1^2+dv1^2);
Fs2=(beta(index2)-beta(i))/sqrt(dh2^2+dv2^2);
F=(Fs1^2+Fs2^2-2*Fs1*Fs2*cos(theta2-theta1))/sin(theta2-theta1)^2;
F=sqrt(F);        

end
%---------------------------------------
