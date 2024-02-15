function [ output_args ] = WriteAT2(AT2Dir,ans1,ans2,ans3,acc,dt)
%WRITEAT2 Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(AT2Dir,'wt');
%fprintf(fid,'%s\r',ans1);
fprintf(fid,'Energy-Matched Spectrum-Matched (EMSM) STRONG MOTION DATA, BY Duruo HUANG and Gang WANG %s\r',date);
fprintf(fid,'%s\r',ans2);
fprintf(fid,'%s\r',ans3);
fprintf(fid,'%d   %g   NPTS,   DT\r',length(acc),dt);

for p=1:length(acc);
    fprintf(fid,'%15.6E',acc(p));
    if mod(p,5)==0;
        fprintf(fid,'\r');
    end
end
fclose(fid);
end

