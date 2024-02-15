function acc = RemoveMean(acc,s)
if s==0
    acc = acc - mean(acc);
else
    acc = acc - mean(acc(1:10));
end

return
