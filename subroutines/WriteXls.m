function [] = WriteXls(XlsDir,ii,ID,nroll,HPNS,LPNS,HPEW,LPEW,HPUD,LPUD)
%WRITEXLS Summary of this function goes here
%   Detailed explanation goes here
a=[ID,nroll,HPNS,LPNS,HPEW,LPEW,HPUD,LPUD];

s2 = strcat('A',num2str(ii+1),':H',num2str(ii+1));


xlswrite(XlsDir,a,1,s2);
end

