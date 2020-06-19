function [ BB_part,L ] = calc_B_part_arb_lev( base,N_x,N_y,L5,depth )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


n_xA=1:size(base,2);
n_xA=n_xA-(length(n_xA)-2^depth );
n_yA=1:size(base,1);
n_yA=n_yA-(length(n_yA)-2^depth );

u=zeros(1,L5(1)*L5(2)*N_x*N_y);

for ii=1:length(n_xA);%length(n_xA)
    for jj=1:length(n_yA);%length(n_yA)



start_x=n_xA(jj);
start_y=n_yA(ii);

base_pix=base(jj,ii);

x_places=start_x:2^depth:start_x+(L5(1,2)-1)*2^depth;
y_places=start_y:2^depth:N_y+(L5(1,1)-1)*2^depth;
%y_places=start_y:2^depth:start_y+(L5(1,1)-1)*2^depth;

NN=N_x*N_y;

[X_places,Y_places]=meshgrid(x_places,y_places);
 X_places=X_places.';
 Y_places=Y_places.';


X_places_vec=X_places(:);
Y_places_vec=Y_places(:);
is_pos= and(and((X_places_vec>0),(Y_places_vec>0)),and((X_places_vec<=N_x),(Y_places_vec<=N_y)));
%is_pos=(X_places_vec>0).*(Y_places_vec>0).*(X_places_vec<=N_x).*(Y_places_vec<=N_y);

places=X_places_vec+N_y*(Y_places_vec-1);
n=[1:length(places)]';
places2=places+(n-1)*NN;
%is_pos2=is_pos+(n-1)*NN;


places2_pos=places2(is_pos);

u(places2_pos)=base_pix;

    end
end
BB_part=reshape(u,N_x*N_y,L5(1)*L5(2));
end


