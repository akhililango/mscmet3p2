function [stat, pval] = dieboldmariano(msfe,mafe)

d = msfe - mafe;
dmu = mean(d);
dvar = var(d);
n = size(d,2);
gammad = dvar * autocorr(d,n-1);

% ind = zeros(n,1);
% for t = 1:n
%     if t/
% omz = 2*sum(gammad(ind<=1))-gammad(1); 
% end

stat = dmu/(sqrt(dvar/n));

count = 0;
for t = 1:n
    if ((d-dmu)/sqrt(dvar))>0.95
        count = count + 1;
    end
end

pval = count/n;

end