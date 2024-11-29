function [rtxx,ttp,mxf,mif,Td_xx]=Txx(time,force,xx)
[mxf,mxf_index] = max(force);
mif = min(force);

df=(diff(force')./diff(time));
[~, dfmax_index]=max(df);   %max dfdt
t_dfmax = time(dfmax_index);

%Time to peak (TP)
%tp=time(mxf_index)-t_dfmax;
tstim=0;
t_peak = time(mxf_index);
ttp=t_peak-tstim; %only valid if stim starts at t=0 ms.

%TDuration (not used)
T_xxrec=mxf-(xx/100)*(mxf-mif);
t_xx = interp1(force(mxf_index:end),time(mxf_index:end),T_xxrec); %puede dar error
%     y1=force(mxf_index:end);
%     x1=time(mxf_index:end);
%     [y2, index2] = unique(y1);
% t_xx = interp1(y2,x1(index2),T_xxrec); 
Td_xx = t_xx-t_dfmax;

%TP to xx% decay
rtxx=t_xx-t_peak;  % time to xx% relaxation