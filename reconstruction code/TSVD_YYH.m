function e_TSVD=TSVD_YYH(U,D,V,lambda,max_s)
D_diag=diag(D);
B=D_diag-lambda*max_s;
[x,index]=min(B(find(B>0)));
cc=index;
D_diag_i=1./D_diag;
D_diag_i_t=D_diag_i(1:cc);
D_i_t=diag(D_diag_i_t);

U_t=U';
U_t=U_t(1:cc,:);
V=V(:,1:cc);
e_TSVD=V*D_i_t*U_t;
end