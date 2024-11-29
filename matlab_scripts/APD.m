function [APDxx,Vpp,dVdt_max,Vr,Vp,t_maxderv]=APD(t,v,xx)

%APD at different repolarization % (xx, eg 90%)
%Other parameters: amplitude (Vpp), max upstroke (dVdt_max), Vrest, Vpeak
% inputs: t & v are vectors of 1 AP. xx is the % of repolarization.
% example: apd90 = APD (time,v,90)
    
    [Vp, Vp_index]=max(v); %Caution with elevated domes
    [Vr, Vr_index]=min(v);
    Vpp=Vp-Vr;
    
    %Derivative of AP
    dVdt=diff(v)./diff(t);
    
    
    [dVdt_max, dVdt_max_index]=max(dVdt);
    t_maxderv=t(dVdt_max_index);
    % figure
    % plot(t(1:end-1),dVdt)
    
    V_xxrep=Vp-(xx/100)*Vpp;
    %interpolation (repolarization)
    if V_xxrep < v(end)
        t_vxx=t(end);
        disp('Error in APD calculation. There is no physiological AP')
        APDxx=nan;
    else
        t_vxx = interp1(v(Vp_index:end),t(Vp_index:end),V_xxrep);
        APDxx=t_vxx-t_maxderv;
    end

end