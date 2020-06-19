function places=find_important_places2(BB_part,energ_fact)

% down_sample_factor=round(N_x/2^depth_im/2^depth_proj);
%[a,b]=find(abs(base)==max(abs(base(:))));


%seam=ceil((length(base)-max([a,b]))/2^depth_im);

%u=zeros(sqrt(size(BB_part,2)),sqrt(size(BB_part,2)));
% n=1:sqrt(size(BB_part,2));
% [Nx,Ny]=meshgrid(n,n);
%u=(abs(Nx)<ceil(n(end)/2+seam)).*(abs(Ny)<ceil(n(end)/2+seam));
% down_sampe_mat=zeros(sqrt(size(BB_part,2)),sqrt(size(BB_part,2)));
% down_sampe_mat(  and((mod(Nx,down_sample_factor)==0), (mod(Ny,down_sample_factor)==0)))=1;

%u=u.*down_sampe_mat;

% u=ones(Nx,Ny);
places=1:size(BB_part,2);

%find(u>0);

energy=sum(abs(BB_part).^2);
maxi=max(energy);
places=places((energy>energ_fact*maxi));





