function a = get_lprior(n,p,M,Mmax)
    a = -log(p) - log(Mmax) - lchoose(n,M);
end