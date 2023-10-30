function y = glottalpulse(dt, t0, t1, t2, A)
    t = 0:dt:t0-dt;
    y = zeros(1, size(t,2));
    k = 1/(1-cos(pi*(t2-t1)/t1));
    idx1 = 0 <= t & t <= t1;
    y(idx1) = 0.5*A.*(1-cos(pi.*t(idx1)./t1));
    idx2 = t1 < t & t <= t2;
    y(idx2) = A.*(k.*cos(pi.*(t(idx2)-t1)./t1)-k+1);
    y(~(idx1 | idx2 )) = 0;
    figure;
    plot(t,y);
end