function [f,h] = FourierSpc(acc,dt);
%%%% reinterprete dt

np=length(acc);
% Newdt=0.01;
% Newnp=floor(np*dt/Newdt);
% acc=interp1([1:np]*dt, acc,[1:Newnp]*Newdt);
% dt=Newdt;
% np=Newnp;
% 
% time=[1:np]*dt;
 
%   acc=-bb(ij).*acc;
%%% check if delta_f =< 0.05 Hz as required by Rathje for accuracy
%%% To ensure a stable value of Tm is calculated for recorded strong 
%%% ground motions, motions should contain at least the minimum number 
%%% of points indicated above or should be augmented with zeroes to
%%% attain these minimum values.
%%% delta_f=1/(np*dt)=< 0.05 Hz, therefore, 
%%% requires np*dt>20, i.e., np>20/dt; otherwise, append zeros in the back

%%% input: acc -- time sequence; np -- no. of points in acc; dt--delta time
% a zero-filled section is included at the end of adequate length(25% of
% the record or more) to prevent problem with periodicity.

NFFT = 2^nextpow2(np); 
  
%%% question: which is true???
%Y = fft(acc,NFFT)/np;  %% FFT(X,N) is the N-point FFT, padded with zeros if X has less
                        %% than N points and truncated if it has more
Y = fft(acc,NFFT);      % changed by G. Wang, Fourier coefficient in complex form, unit g-sec
%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%h=abs(Y(1:NFFT/2+1));  % this is Fourier coefficient. Norm of complex number, unit g-sec
h=2*abs(Y(1:NFFT/2+1));  % multiple by factor of 2, cf. Kramer, p.540, eq. (A17)
hh=h.^2;               % Ci^2 in Rathje paper
F=1/dt;
f=F/2*linspace(0,1,NFFT/2+1); % the maximum frequency of FFT is related to delta time
f=f';                 % this is Fourier frequency
 
df=1/(NFFT*dt);
 return