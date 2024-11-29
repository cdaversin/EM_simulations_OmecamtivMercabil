% MATLAB Implementation of the ToR-ORd model with Land for evaluating OM
% effects

clearvars, clc

%%SETTINGS----------------------------------------------------------------
OM = 0.5;   %uM
hf=1;       %0-N, 1-HF
CL=1000;    %pacing cycle length in ms
%-------------------------------------------------------------------------
%global settings
beats=1000; %%number of beats in the simulation
mech.cell='endo';
options=[];%options for ode solver
tspan=[0 CL]; %[0 CL] variable dt

mech.isacs = 0; % 0 / 1: Isacs
mech.emcoupling = 1; % 0: ORd only; 1: ORd-Land interaction and feedback
    mech.mode = 'intact'; % 'skinned'/'intact'
    mech.lambda = 1; %1.1 max in Isacs; default Land: 1
    mech.dLambda = 0;
mech.calib = 1; %1: Margara calibration in ToR

%% OM model and parameters
load('om_paramsTest1.mat')
conc = OM;
% Method (Hill-EC50)
    kuw_sf=1+om_params(1)./(1+(om_params(2)./conc).^om_params(3));
    kwu_sf=1-om_params(4)./(1+(om_params(5)./conc).^om_params(6));
    kws_sf=1-om_params(7)./(1+(om_params(8)./conc).^om_params(9));

% Alternative (only values)
%kuw_sf=1; kwu_sf=1; kws_sf=1; 

mech.factor1=kuw_sf;
mech.factor2=kwu_sf;
mech.factor3=kws_sf;

%% SS from 1000 beats 
if hf==0
    model=@model_Tor_Land_OM; 
    X0 = [-89.2870281355411	12.1338410842619	12.1341640487601	144.784293578711	144.784252229768	7.29387679569686e-05	6.35372454439132e-05	1.54037051063521	1.54172446030222	0.000719856036450669	0.687510908830619	0.838779899426743	0.838813905144327	0.838674000154056	0.000147486870284105	0.543896845093974	0.302705931989102	0.000918014372963919	0.999644044328252	0.597409578711185	0.000467734022901158	0.999644041779297	0.660801853228359	-4.55711479581949e-36	0.999999993556583	0.940313438791041	0.999999993556567	0.999900071673720	0.999978657109148	0.000449142946108774	0.000765706251797605	0.999999993557539	0.999999993557324	0.240063167496228	0.000166938115757851	8.23356359711029e-25	0.0112591899138282	0.998174009932030	0.000820041243970364	0.000675575913830246	0.000319390552512148	1.09803752534757e-05	-1.71057952089023e-22	0.000157349637623475	0.000236855375284938	0.00813593169007039	0.999368827725695	0	0	0];
elseif hf==1
    model=@model_Tor_Land_HF_OM; 
    X0 = [-88.2223664155536	14.0373924161249	14.0378866607386	142.220154324742	142.220109538740	7.64018486992280e-05	6.12753961759702e-05	0.671578239783226	0.677089627971319	0.000905197735336801	0.652629428829952	0.817402176761659	0.816649863035091	0.814935804849531	0.000180541315521976	0.444619634186746	0.220844702265182	0.000986342950809401	0.999570899255959	0.536382820553116	0.000502564716361747	0.999570926996286	0.585144684111519	-2.23894191842084e-24	0.999999991398650	0.918383307239079	0.999999991398899	0.999707312918462	0.999937021705817	0.000390539962112009	0.000916250382408384	0.999999991395468	0.999999991388783	0.304433850552757	0.000188318263004268	-2.72806105472088e-24	0.00541503658501152	0.997137512375064	0.000883942040440329	0.000778878313131793	0.00115744083222512	4.22231764649591e-05	2.37280619795622e-24	0.000642906793140312	0.000907841093524710	0.0181051597601256	0.997544748064031	0	0	0];
end


for n=1:beats
    [time X]=ode15s(@(t,y)model(t,y,mech,1),tspan,X0);
    X0=X(size(X,1),:);
    n; %output beat number to the screen to monitor runtime progress
end

%%
%calculate and name dependent variables for the final beat in the
%simulation (i.e. currents and fluxes)
for i=1:size(X,1)
    %IsJs=model(time(i),X(i,:),mech,0);
    tension=model(time(i),X(i,:),mech,2);

        Ttot(i)=tension(1);
        Ta(i)=tension(2);
        Tp(i)=tension(3);
end

%% Outputs
apd90 = APD(time,X(:,1),90); %ms
caimax = max(X(:,6)); %mM
[rt50,tp,maxTa,minTa]=Txx(time,Ta,50);
devF=maxTa-minTa;

figure(1), subplot(131), plot(time,X(:,1)), ylabel('Voltage (mV)'), hold on
subplot(132), plot(time,X(:,6)), ylabel('Cai (nM)'), hold on
subplot(133), plot(time,Ta), ylabel('Active Tension (kPa)'), hold on
