 function logpdf = betalogpdf(x,a,b)
        logpdf = (a-1)*log(x) + (b-1)*log(1-x) - betaln(a,b);

end

