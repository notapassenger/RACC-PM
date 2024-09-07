function[MaxEigvalue, Eigvector] = MaxEigenvector(R)
    [V,D] = eig(R);
    [d,ind] = sort(diag(D));
    Ds = D(ind,ind);
    Vs = V(:,ind);
    a = size(R, 1);
    MaxEigvalue = Ds(a, a);
    Eigvector = Vs(:, end);