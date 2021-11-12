function [orthonorm_e, orthonorm_a]= orthonorm_GS(a, b)
%ORTHONORM_GS Converts column vectors a, b to set of orthonormal vectors (Gram-Schmidt method)
aTb = a'*b;
aTa = a'*a;
x = aTb/aTa;
l = a*x;
e = b-l;

%check_eTa = e'*a
%check_aTe = a'*e
orthonorm_e = e/norm(e);
orthonorm_a = a/norm(a);
end

