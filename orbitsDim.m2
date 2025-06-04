needsPackage "NumericalImplicitization";
f = (i,j) -> if i<=j then if i==j then 1/2 else 1 else 0;
orb_dim = (d,m,n) -> (
R := CC[s_(1,1)..s_(m*n,d)];
S := genericMatrix(R,s_(1,1),m*n,d);
Am := matrix(table(m,m,f));
An := matrix(table(n,n,f));
A := 4*tensor(Am,An);
m = map(R, CC[x_(1,1)..x_(d,d)], flatten entries (transpose(S)*A*S));
numericalImageDim(m,ideal 0_R)
);

orb_dim(8,4,4)
