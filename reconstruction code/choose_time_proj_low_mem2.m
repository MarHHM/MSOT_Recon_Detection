function [ places ] = choose_time_proj_low_mem2( M,over_fact,wl_num )  % for second level db2, returns the places of the relevant projections


%MM=sum(abs(M.').^2);

MM=max(abs(M.'));
ele_num=over_fact*wl_num;
M_s=sort(MM,'descend');
thresh=M_s(ele_num);
places=find(MM>thresh);

% 
% 
% MM2=reshape(MM,L(1,1),size(MM,2)/L(1,1));
% %MM2=sum(MM2.^2);
% 
% MM2=max(MM2);
% 
% ele_num=round(wl_num*over_fact/L(1,1));
% M_s=sort(MM2(:),1,'descend');
% thresh=M_s(ele_num);
% proj_vec=MM2>thresh; 
% AAA=0;
% for kk=1:L(1,2)*4^depth
%     if proj_vec(kk)==1
%          AAA(1+(kk-1)*L(1,1):kk*L(1,1))=1;
%     end
% end
% places=find(AAA==1);

%MM=M(:,(im_ele-1)*size(M,2)/16+1:size(M,2)/16*im_ele).';

% 
% A=sum(abs(MM))>0.5*max(sum(abs(MM)));
% AAA=0;
% for kk=1:L(1,2)*16
%     if sum(A(1+(kk-1)*L(1,1):kk*L(1,1)))>0
%         AAA(1+(kk-1)*L(1,1):kk*L(1,1))=1;
%     end
% end
% places=find(AAA==1);

end

