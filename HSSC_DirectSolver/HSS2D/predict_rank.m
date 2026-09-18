function [rank,C] = predict_rank(n,acc)
%
%

p = ceil(-log10(acc));
if p<7
    p=7;
else
    if p>12
        p = 12;
    end
end
m(7:12) = [96 138 184 275 526 487];
b(7:12) = [138 215 300 496 1074 954]; 
r0(7:12) = [86 104 116 130 137 140];

if n<150
    rank = min(n,r0(p)); 
    C    = ceil(m(p)*log10(2)/2); 
else
    rank = min(n,ceil(m(p)*log10(n) - b(p))); 
    C    = ceil(m(p)*log10(2)/2); 
end