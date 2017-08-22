function P = fix_covariance( P )
%UNTITLED force symmetry and  positive variance on covariance matrix
P = 0.5 * (P + transpose(P));
P = 0.5*(P+transpose(P));
for i=1:3
    P(i,i) = max(P(i,i),1e-12);
end

end

