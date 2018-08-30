function[H]=SCC_spanning_tree(lon,lat,p,dc)
% Compute n*(n-1) matrix constructed from edge set of minimum spanning tree
% INPUT ARGUMNETS:
% lon: [n,1] matrix storing the longitudes of n locations
% lat: [n,1] matrix storing the latitudes of n locations
% p  : number of variables
% dc : a critical distance larger than which the two vertices will not be
%      connected. This value is used only for the reduction of computational
%      burden.

%OUTPUT ARGUMENTS
% H  :  n*(n-1) matrix constructed from edge set of minimum spanning tree

n=length(lon);
lon=lon*ones(1,n);
lat=lat*ones(1,n);

tree_exist=0;
while tree_exist==0
    d=sqrt((lon-lon').^2+(lat-lat').^2);
    d(d>dc)=0;
    d=sparse(tril(d));
    [Tree] = graphminspantree(d);
    [index1_path,index2_path]=find(Tree);
    if isempty(index1_path)==1||length(index1_path)<n-1
        dc=dc*2;
    else
        tree_exist=1;
    end
end
    

H1=zeros(n-1,n);

ppp=sub2ind(size(H1),[1:n-1]',index1_path);
H1(ppp)=1;
ppp=sub2ind(size(H1),[1:n-1]',index2_path);
H1(ppp)=-1;

H=zeros((n-1)*p,n*p);
for j=1:p
    H((j-1)*(n-1)+1:j*(n-1),(j-1)*n+1:j*n)=H1;
end 

for j=1:p
    H(p*(n-1)+j,1+(j-1)*n:j*n)=1/n;
end
    