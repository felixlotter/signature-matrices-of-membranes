f = (i,j) -> if i<=j then if i==j then 1/2 else 1 else 0;
orb_dim = (d,m,n) -> (
R := QQ[s_(1,1)..s_(m*n,d)];
g = map(QQ,R,apply(m*n*d,t->random(-10000,10000)));
S = genericMatrix(R,s_(1,1),m*n,d);
Am = matrix(table(m,m,f));
An = matrix(table(n,n,f));
A = 4*tensor(Am,An);
I = ideal(transpose(S)*A*S);
rank g(jacobian(I)))

orb_tbl = d -> (matrix table(d,d,(i,j)->orb_dim(d,i+1,j+1)));
