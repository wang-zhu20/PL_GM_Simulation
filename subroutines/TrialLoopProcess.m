function [acc2,DispT] = TrialLoopProcess(acc,dt,HP,LP,nroll,option)

%%%% loop over the process using different combination of parameters

    %%% apply taper
    acc2 = taper(acc,0.05);
    %%% add zero pad
    [acc2,tpad] = zeropad(acc2,dt,nroll,HP);
    %%% acausal filter
    acc2 = acausal(LP,HP,nroll,acc2,dt); %filter
    %%% remove pad
    acc2 = Removepad(acc2,tpad,dt);
    %if option==1
    %%% baseline correction
        [VelF,DispF,VelT,DispT] = integration(acc2,dt);
        acc2 = LoopBaselinecorrectionUPDATED(acc2,DispT,dt,option); % baseline
    %end
    %%% integrate to displ.
    [VelF,DispF,VelT,DispT] = integration(acc2,dt); %integration
    
