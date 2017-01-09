function Ap=rescale_matrix01(A)


maxi=max(max(A));
mini=min(min(A));

a=1/(maxi-mini);
b=1-1/(1-mini/maxi);
Ap=a*A+b*ones(size(A));