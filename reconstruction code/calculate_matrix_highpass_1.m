function FILT = calculate_matrix_highpass_1(n);

pos = linspace(1,n*n,n*n).';
values_diag = (8/9)*ones(n*n,1);
FILT = sparse(pos,pos,values_diag,n*n,n*n);

pos_left = pos-n;
pos_right = pos+n;
pos_down = pos-1;
pos_up = pos+1;
pos_a = pos_left-1;
pos_b = pos_right-1;
pos_c = pos_left+1;
pos_d = pos_right+1;
cond_left_in = pos_left>0;
cond_right_in = pos_right<=n*n;
cond_down_in = floor((pos_down-1)/n)==floor((pos-1)/n);
cond_up_in = floor((pos_up-1)/n)==floor((pos-1)/n);

columns = [pos_left, pos_right, pos_down, pos_up, pos_a, pos_b, pos_c, pos_d];
rows = pos*ones(1,8);
values = (-1/9)*ones(n*n,8);
cond = [(cond_left_in), (cond_right_in), (cond_down_in), (cond_up_in),...
    (cond_left_in&cond_down_in), (cond_right_in&cond_down_in), (cond_left_in&cond_up_in), (cond_right_in&cond_up_in)];

FILT = FILT + sparse(rows(cond),columns(cond),values(cond),n*n,n*n);


