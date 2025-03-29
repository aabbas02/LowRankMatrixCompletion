function SD=subspace(U1,U2,n)

%o=(eye(n)-U1*transpose(U1))*U2;
o = U2 - U1*(U1'*U2);
SD=norm(o,"fro");
end