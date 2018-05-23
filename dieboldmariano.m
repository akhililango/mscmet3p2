function [stat, pval] = dieboldmariano(d,omegadm,mid)

dmu = mean(d);
n = max(size(d));
dT = n;
periodm = zeros(dT,1);

K = size(omegadm,1);
% for k = 1:K
%     for t = 1:n
%     eperio(:,n) = exp(1i*t*omega(k));
%     end
% end

for k= 1:K
    for t = 1:dT
        periodm(k) = periodm(k) + d(t)*(cos(omegadm(k)*t) - i*sin(omegadm(k)*t)); %exp(-1i*t*omega(k))
    end
    periodm(k) = abs(periodm(k))/dT;
end

h = 4;
% h = ceil(kde(d)*dT);
omz = 0;
% ind = zeros(n,1);
for m = -h:h
    omz = omz + (h+1-abs(m))/(h+1)^2 * periodm(mid+m); 
end

stat = dmu/(sqrt(omz/n));

pval = 2*(1-normcdf(abs(stat)));

end