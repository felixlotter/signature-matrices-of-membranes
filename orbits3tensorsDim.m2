orb_dim = (d,m,n) -> (
R = QQ[s_(1,1)..s_(d,m*n)];
g = map(QQ,R,apply(d*m*n,t->random(-10000,10000)));
f = (i,j,k) -> if i<=j and j<=k then (
    if(i==j and j==k) then 1_R else
    if(i==j or j==k or i==k) then 3_R else
    6_R
) else 0;
kron = (i,j) -> n * (i-1) + j;
f2 = (i1,i2,i3,i4,i5,i6) -> f(i1,i3,i5)*f(i2,i4,i6);
--sig = (1,1,1,1,1,1)..(m,n,m,n,m,n);
--sig = apply(sig, i -> (f(i#0,i#1,i#2)));
tf = j->sum({1,1,1,1,1,1}..{m,n,m,n,m,n},i->f2(i#0,i#1,i#2,i#3,i#4,i#5)*s_(j#0,kron(i#0,i#1))*s_(j#1,kron(i#2,i#3))*s_(j#2,kron(i#4,i#5)));
I = ideal(apply((1,1,1)..(d,d,d),tf));
rank(g(jacobian(I))))
