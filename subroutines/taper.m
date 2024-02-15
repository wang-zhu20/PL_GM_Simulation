function acc = taper(acc,a) %  cosine taper
% 
np = length(acc);
% be = acc(1)*ones(a*np,1);
% en = acc(np)*ones(a*np,1);
% acc = [be;acc;en];
taper =tukeywin(np,a); 
acc = acc.*taper;
return