function[]=SCC_R_write(x,y,lon,lat,filename,method)
% Description: This function produces the data file used by the R scripts

% INPUT ARGUMENTS
% filename: the file name storing the data 
% x:        [sim_num,n,p] array storing the p explanatory variables on n locations
% y:        [sim_num,n] matrix storing the response variables on n locations
% lon:      [n,1] matrix storing the longitudes of n locations
% lat:      [n,1] matrix storing the latitudes of n locations
% method:   method==1 for SCC-MST, ==-1 for SCC-RNN, and ==-2 for SCC-KNN


[s,n,p]=size(x);

if method==1

    [H]=SCC_spanning_tree(lon,lat,p,0.1);
    for t=1:s
        X=squeeze(x(t,:,:));  
        Xg=zeros(n,n*p);
        for j=1:p
            Xg(:,(j-1)*n+1:j*n)=diag(X(:,j));
        end
        yg=squeeze(y(t,:))';
        G=Xg/H;
        filename2=strcat(filename,int2str(t),'.nc');
        SCC_netcdf_write_lasso(lon,lat,yg,Xg,G,filename2);
    end
else
    [H]=SCC_genlasso_penalty_matrix(lon,lat,p,4,method);
    for t=1:s
        X=squeeze(x(t,:,:));  
        Xg=zeros(n,n*p);
        for j=1:p
            Xg(:,(j-1)*n+1:j*n)=diag(X(:,j));
        end
        yg=squeeze(y(t,:))';
    
        filename2=strcat(filename,int2str(t),'.nc');
        SCC_netcdf_write_genlasso(lon,lat,yg,Xg,H,filename2);
    end
end

            
end




% Subroutine SCC_genlasso_penalty_matrix
% Compute penalty matrix for SCC-RNN and SCC-KNN 
function[H]=SCC_genlasso_penalty_matrix(lon,lat,p,neighbour_num,method)
n=length(lon);
lon=lon*ones(1,n);
lat=lat*ones(1,n);
d=sqrt((lon-lon').^2+(lat-lat').^2);

q=neighbour_num;

if method==-1
    % Compute penalty matrix for RNN
    
    % changanble critical value
    dc=0.05;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    Hj=zeros(1,n);
    t=0;
    for i=1:n
        di=d(i,:);
        [di,index]=sort(di);
        ccc=find(di<dc,1,'last');
        if isempty(ccc)==1
        
        else         
            for k=2:ccc
                t=t+1;
                Hj(t,index(k))=-1;
                Hj(t,i)=1;
            end
        end
    end
    [aaa,~]=size(Hj);
    H=zeros(aaa*p,n*p);
    for j=1:p
        H((j-1)*aaa+1:j*aaa,(j-1)*n+1:j*n)=Hj;
    end  
elseif method==-2
    % Compute penalty matrix for KNN
    Hj=zeros(q*n,n);
    for i=1:n
        di=d(i,:);
        [~,index]=sort(di);
        index=index(2:q+1);
        for k=1:q
            Hj(q*(i-1)+k,index(k))=-1;
            Hj(q*(i-1)+k,i)=1;
        end
    end
    H=zeros(q*n*p,n*p);
    for j=1:p
        H((j-1)*q*n+1:j*q*n,(j-1)*n+1:j*n)=Hj;
    end    
end

end
%----------------------------------------------

% Subroutine SCC_netcdf_write_genlasso
function[output]=SCC_netcdf_write_genlasso(lon,lat,yg,Xg,H,filename)
[n,p]=size(Xg);
q=length(H(:,1));
ncid=netcdf.create(filename,'64BIT_OFFSET');

dimid_obs=netcdf.defDim(ncid,'obs',n);
dimid_var=netcdf.defDim(ncid,'var',p);
dimid_pen=netcdf.defDim(ncid,'pen',q);

varid_lon=netcdf.defVar(ncid,'lon','double',[dimid_obs]);
varid_lat=netcdf.defVar(ncid,'lat','double',[dimid_obs]);
varid_y=netcdf.defVar(ncid,'y','double',[dimid_obs]);
varid_X=netcdf.defVar(ncid,'X','double',[dimid_obs,dimid_var]);
varid_L=netcdf.defVar(ncid,'L','double',[dimid_pen,dimid_var]);
netcdf.endDef(ncid);

netcdf.putVar(ncid,varid_lon,lon);
netcdf.putVar(ncid,varid_lat,lat);
netcdf.putVar(ncid,varid_y,yg);
netcdf.putVar(ncid,varid_X,Xg);
netcdf.putVar(ncid,varid_L,H);
netcdf.close(ncid);
output='OK';
end
%----------------------------------------------

% Subroutine SCC_netcdf_write_lasso
function[output]=SCC_netcdf_write_lasso(lon,lat,yg,Xg,G,filename)
[n,p]=size(Xg);
ncid=netcdf.create(filename,'64BIT_OFFSET');

dimid_obs=netcdf.defDim(ncid,'obs',n);
dimid_var=netcdf.defDim(ncid,'var',p);


varid_lon=netcdf.defVar(ncid,'lon','double',[dimid_obs]);
varid_lat=netcdf.defVar(ncid,'lat','double',[dimid_obs]);
varid_y=netcdf.defVar(ncid,'y','double',[dimid_obs]);
varid_X=netcdf.defVar(ncid,'X','double',[dimid_obs,dimid_var]);
varid_G=netcdf.defVar(ncid,'G','double',[dimid_obs,dimid_var]);
netcdf.endDef(ncid);

netcdf.putVar(ncid,varid_lon,lon);
netcdf.putVar(ncid,varid_lat,lat);
netcdf.putVar(ncid,varid_y,yg);
netcdf.putVar(ncid,varid_X,Xg);
netcdf.putVar(ncid,varid_G,G);
netcdf.close(ncid);
output='OK';
end
%----------------------------------------------