function v = find_t_DOF(x)
    d = Mahalanobis(x); 
    v_vec=3.5:0.5:10;
    cnt=1;
    eX = eprob(d); 
    c = eX.eprob(end:-1:1);
    srt_d = sort(d);
    %default -4
    Pi = logspace(0,-4,1000);
    for cnt=1:length(v_vec)
%         cnt
        v = v_vec(cnt);
        d_theor_t_correct = (v/((v-2)*size(x,1)))*d;
        exc_t = 1-fcdf(sort(d_theor_t_correct),size(x,1),v);
           
        exc_metric = 0;
        for i=1:length(Pi)
            [dummy idx1] = min(abs(exc_t - Pi(i)));
            [dummy idx2] = min(abs(c - Pi(i)));
            exc_metric = exc_metric + abs(srt_d(idx1) - srt_d(idx2));
        end
        exc_m(cnt) = exc_metric;
    end
    [dummy idx_min] = min(exc_m);
    v = v_vec(idx_min);
end