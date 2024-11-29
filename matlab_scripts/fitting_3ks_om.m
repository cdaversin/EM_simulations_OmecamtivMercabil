%% Fitting OM model parameters to experimental data

clearvars, clc 

%%% Theoretical curve from experimental data 
% Data: F vs [D] (Shen et al. 2021)
load('drug_points_shen21.mat')
x1data = drug_points(:,1); y1data = drug_points(:,2);
x2data=[0.1 3]; y2data=[0.1 0.5];    %TTP delay increase

mult=-3:0.1:2; c=10.^mult; %0.001:0.001:100; %OM dose [uM]
paramH_0 = [1.5 0.5 1];
[paramH_exp,y1] = hill_eq(paramH_0,x1data,y1data,c);

figure(10), subplot(121),semilogx(x1data,y1data,'kx'), hold on
    semilogx(c,y1,'k'), hold on, box off
    xlabel('[OM] (\muM)'), ylabel('Relative Force increase')
    subplot(122),semilogx(x2data,y2data,'kx'), hold on, box off
    xlabel('[OM] (\muM)'), ylabel('Relative TTP delay increase')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OPTIMIZED Parameters of hill equation for 3 k that best fit exp behavoiur

%[D]:
multx=-2:0.2:2; conc=10.^multx;

% Based on MinObjective function ------------------------------------------

%function = @(param,D) (param(1)) ./ (1+(param(2)./D).^param(3)); %hill equation
% x 3 parameters

Param0= [1.5,0.5,1,... %kuw;
        0.2,0.5,1,...     %kwu;
        0.1,0.5,1,];      %kws;


% Initially:
[MinObj0, Y1_p0, Y2_p0] = ObjFuncT(Param0,0);
figure(10)
subplot(121), semilogx(conc,Y1_p0,'yo-')
subplot(122),semilogx(conc,Y2_p0,'yo-'),


%% Constraints
lowerbound = [1, 0.1, 1.1,...
            0.1, 0.1, 1.1,...
            0.01, 0.1, 1.1];   
upperbound = [10, 1, 2,...
             1, 1, 2,...
             0.5, 1, 2];  


n_param = length(Param0);

%%% GLOBAL Optim
opts = optimoptions(@fmincon,Algorithm="interior-point");

problem = createOptimProblem('fmincon', ...
    x0=Param0,...
    objective=@ObjFuncT,...
    lb = lowerbound,...
    ub = upperbound,...
    options=opts);
tic
gs = GlobalSearch;
[paramG,MinG] = run(gs,problem);
tComp=toc;

%save globaloptim   % to save all variables
om_params = paramG;
save('om_paramsTest1','om_params')  % variables needed for a basic_simulation

%--->optim_SOL:
%9.75193572328972	0.635033554678567	1.87999419950932	0.759427395574438	0.781234296422715	1.51841567375560	0.472150484223696	0.224622359002761	1.13136895214319

[MinObj, Y1_g, Y2_g] = ObjFuncT(paramG,0);
figure(10)
subplot(121), semilogx(conc,Y1_g,'b*'), hold on
subplot(122),semilogx(conc,Y2_g,'b*'), hold on

% Hill equation 
[paramH_calc,yfun1] = hill_eq(paramH_exp,conc,Y1_g);
% --> Sol: 1.46806613841073	0.390590561420242	1.91964305697727
figure(10)
subplot(121),semilogx(conc,yfun1,'b.-')