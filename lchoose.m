function a = lchoose(n,k)
    a = 0;
    for ii=1:n
        a = a + log(ii);
    end
    for ii=1:k
        a = a - log(ii);
    end
    for ii=1:(n-k)
        a = a - log(ii);
    end
end