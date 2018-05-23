function [stat, pval] = dieboldmariano(msfe,mafe,perio)

d = msfe - mafe;
dmu = mean(d);
dvar = var(d);
n = size(d,2);

h = 4;
omz = 0;
% ind = zeros(n,1);
for m = -h:h
    omz = omz + (h+1-abs(m))/(h+1)^2 * perio(0); 
end
omz = omz/(2*pi);

stat = dmu/(sqrt(dvar/n));

count = 0;
for t = 1:n
    if ((d-dmu)/sqrt(dvar))>0.95
        count = count + 1;
    end
end

pval = count/n;

end