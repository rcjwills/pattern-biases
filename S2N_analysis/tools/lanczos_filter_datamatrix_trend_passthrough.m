function [X_f] = lanczos_filter_datamatrix_trend_passthrough(X,cutoff)

nt = size(X,1);

X_f = X; t = (1:nt)';
for i = 1:size(X,2)
    p = polyfit(t,X(:,i),1);
    tmp = X(:,i)-p(1)*t-p(2);
    tmp = lanczos([flipud(tmp); tmp; flipud(tmp)],1,cutoff);
    X_f(:,i) = tmp((nt+1):2*nt) + (p(1)*t + p(2));
end