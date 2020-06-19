function [base ] = return_base_arb_lev_even_db(lev_n,depth,wl_name)  %calculates the wavelet reconstruction matrix 
%normally, ord=2;


[LO_D,HI_D,LO_R,HI_R] =wfilters(wl_name);
base_cell{1}=(LO_R).'*(LO_R);
base_cell{2}=(HI_R).'*(LO_R);
base_cell{3}=(LO_R).'*(HI_R);
base_cell{4}=(HI_R).'*(HI_R);

lev_n=lev_n-1;

lev=zeros(1,depth);
for kk=1:depth
    lev(kk)=floor(lev_n/4^(depth-kk));
    lev_n=lev_n-4^(depth-kk)*lev(kk);
end
lev=lev+1;
    

% lev1=ceil(lev_n/4);
% lev2=lev_n-4*(lev1-1);

%%%%%%%
% How can we calculate only one element in B every time?
% 
% We need to have the right building blocks. 
% this is done by convolving the bases + down/up sampling taken into account. 
% everything is done with jumps of 2

base=base_cell{lev(end)};
  for ii=depth-1:-1:1
      base_spread=zeros(size(base,1)*2-1,size(base,2)*2-1);
      for kk=1:size(base,1)
         for uu=1:size(base,2)
            base_spread(kk*2-1,uu*2-1)= base(kk,uu);
         end
      end
      base=conv2(base_spread,base_cell{lev(ii)}); 
  end

%   
%   
% base_A_spread=zeros(7,7);
% for kk=1:4
%     for uu=1:4
%        base_spread(kk*2-1,uu*2-1)= base_cell{lev2}(kk,uu);
%     end
% end
% 
% base=conv2(base_spread,base_cell{lev1});  %this is an example of the basic core. needs to be repeated...
% 
