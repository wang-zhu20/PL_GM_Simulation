%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get predicted wavelet parameter from emperical equation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 %% Note: totalEnergy should *0.01
function [outprm,outpl]=fn_PredictPara_PL(M,Rrup,Vs30,SorD,flag)


outprm = fn_initDBs(1);
%%
load COVDATA.mat
epst_total = mvnrnd(zeros(15,1),COVDATA);

%%%% WZ's regression_pulse %%%%%%%%%%%%%%%%%
epst_p = epst_total(1:3);

LnTp = -2.159 + 0.523.*M - 0.113* log(Vs30)+0.012 * SorD;
LnEx = -3.644 + 0.795.*M + 0.086* log(sqrt(Rrup.*Rrup + 36)) -0.04 *log(Vs30); 
LnEacc = 0.656 + 1.298.*M - 1.282* log(sqrt(Rrup.*Rrup + 36)) - 0.685*log(Vs30)+0.024 * SorD;

if flag == 0
    Tp = exp(LnTp);
    Ex = exp(LnEx);
    Eacc = exp(LnEacc)*10000;
else
    Tp = exp(LnTp+epst_p(3));
    Ex = exp(LnEx+epst_p(2));
    Eacc = exp(LnEacc + epst_p(1))*10000;
end
outpl.Tp = Tp;
outpl.Ex = Ex;
outpl.Eacc = Eacc;
%%%% WZ's regression_residual %%%%%%%%%%%%%%%%%
while true

    epst = epst_total(4:end);
    
    
      Coeff = [-1.577191364	0.729127198	0.120925894	-0.207793424
0.280453825	0.593512142	-0.025546664	-0.351687806
0.148841898	-0.024653927	0.026792068	0.2571908
-0.928606713	0.217150467	0.052532866	0.19804042
-0.730566244	0.056764447	-0.042542368	0.054018767
-2.377120634	0.774204268	0.149220441	-0.165168556
-0.709884331	0.591385586	-0.009302822	-0.266564411
0.429226324	-0.138888213	0.077916051	0.249675023
-1.334281647	0.170671638	0.168037198	0.158267337
-0.791377824	0.020587438	-0.055834256	0.109313498
-3.62503769	0.455240251	-0.887971241	0.032392767
-4.861941284	-0.007398093	-0.980545457	0.012758966];



base_x = [1, M, log(sqrt(Rrup^2+36)),log(Vs30)];
Y = Coeff*base_x';

    if flag==0
        minorEx = exp(Y(1));
        minorSx = exp(Y(2));
        minorEy = exp(Y(3));
        minorSy = exp(Y(4));
        minorRxy = Y(5);      minorRxy =2*normcdf(minorRxy,0,1)-1;
        majorEx = exp(Y(6));
        majorSx = exp(Y(7));
        majorEy = exp(Y(8));
        majorSy = exp(Y(9));
        majorRxy =Y(10);        majorRxy=2*normcdf(majorRxy,0,1)-1;
        totalEnergy = exp(Y(12));
        majorEa  = exp(Y(11));
    else
        minorEx = exp(Y(1)+epst(1));
        minorSx = exp(Y(2)+epst(2));
        minorEy = exp(Y(3)+epst(3));
        minorSy = exp(Y(4)+epst(4));
        minorRxy = Y(5)+epst(5);      minorRxy =2*normcdf(minorRxy,0,1)-1;
        majorEx = exp(Y(6)+epst(6));
        majorSx = exp(Y(7)+epst(7));
        majorEy = exp(Y(8)+epst(8));
        majorSy = exp(Y(9)+epst(9));
        majorRxy =Y(10)+epst(10);        majorRxy=2*normcdf(majorRxy,0,1)-1;
        totalEnergy = exp(Y(11)+epst(11));
        majorEa  = exp(Y(12)+epst(12));
    end
    %%%% Version 1: random is longmornal
    %minorRnd = exp(Y(13));
    mn13=1.25;
    sn13=0.11;
    m13=log(mn13^2/sqrt(sn13^2+mn13^2));
    s13=sqrt(log(sn13^2/mn13^2+1));
    minorRnd=lognrnd(m13,s13);
    
    %% to lxly
    majorVx = majorSx^2;
    majorVy = majorSy^2;
    majorExy = majorRxy*majorSx*majorSy + majorEx*majorEy;
    
    majorElx = log(majorEx^2 / sqrt(majorVx+majorEx^2));
    majorVlx = log(majorVx/majorEx^2 + 1);
    majorSlx = sqrt(majorVlx);
    majorEly = log(majorEy^2 / sqrt(majorVy+majorEy^2));
    majorVly = log(majorVy/majorEy^2 + 1);
    majorSly = sqrt(majorVly);
    
    majorCovlxly = log(majorExy/majorEx/majorEy);
    majorRlxly = majorCovlxly/majorSlx/majorSly;
    
    minorVx = minorSx^2;
    minorVy = minorSy^2;
    minorExy = minorRxy*minorSx*minorSy + minorEx*minorEy;
    
    minorElx = log(minorEx^2 / sqrt(minorVx+minorEx^2));
    minorSlx = sqrt(log(minorVx/minorEx^2 + 1));
    minorVlx = minorSlx^2;
    minorEly = log(minorEy^2 / sqrt(minorVy+minorEy^2));
    minorSly = sqrt(log(minorVy/minorEy^2 + 1));
    minorVly = minorSly^2;
    
    minorCovlxly = log(minorExy/minorEx/minorEy);
    minorRlxly = minorCovlxly/minorSlx/minorSly;
    
    majorLLCov=[majorVlx majorCovlxly; majorCovlxly majorVly];
    [T,err] = cholcov(majorLLCov);
    if err == 0
        break
    end
end
outprm(1) = struct( 'M',M,'hdist',Rrup,'vs30',Vs30,...
                'minorElx',minorElx, 'minorSlx',minorSlx, 'minorVlx',minorVlx, 'minorEly',minorEly, 'minorSly',minorSly, 'minorVly',minorVly,'minorRlxly',minorRlxly, ...
                'majorElx',majorElx, 'majorSlx',majorSlx, 'majorVlx',majorVlx, 'majorEly',majorEly, 'majorSly',majorSly, 'majorVly',majorVly,'majorRlxly',majorRlxly, ...
                'totalEnergy',totalEnergy, 'majorEa',majorEa, 'minorRnd',minorRnd); 
