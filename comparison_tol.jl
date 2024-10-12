all_tol = [1e-5;1e-4;5e-4;1e-3;5e-3;1e-2];
for j = 1:length(all_tol)
    nfig = j;
    tol = all_tol[j];
    gap(50,tol,nfig);
end
