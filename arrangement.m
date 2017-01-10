function Jpos=arrangement(J,vmed)

%test arrangement

% transfor vect -> nm,1
[l,c]=size(J);
Jvect=zeros(l*c,1);
for i=1:l
    Jvect(c*(i-1)+1:c*i)=J(i,:);
end
Jvect=sort(Jvect);
quartile=Jvect(round(l*c*3/4));
M=Jvect(end);
m=Jvect(1);

a1=vmed/(quartile-m);
b1=vmed-a1*quartile;

a2=.8/(M-quartile);
b2=vmed-a2*quartile;

Jpos=(J<quartile).*(a1*J+b1)+(J>=quartile).*(a2*J+b2);