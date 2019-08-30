

clear all,
% close all,
clc
startup; % resetting of default plot parameters

global Hm TmD MwD pD Tm R MwP pP lambda
global Vdrug v1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 8.314472; % gas constant [J/(K*mol)]

% %% Gordon-Taylor
XGT = [1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0]'; %XIMC mass fraction of drug
TgGT = [43.24 49.15 54.41 61.48 68.02 72.29 77.72 84.50 92.95 102.36 114.30]'; %TgIMC



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Measurements of melting point depression
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 


IMC_solu = [ 
    0.90 156.85 155.90 % add your data here --> Indomethacin
    0.85	152.95	153.44	;
    0.80 149.47	148.13	;
    0.75 139.57	139.35	;
    0.70 132.76 133.81  ];

IMC_solu_25 = [ 
    0.90 157.53 156.67% add your data here --> Indomethacin
    0.85 150.42 153.13	;
    0.80 141.41 142.12	;
    0.75 131.26 130.41	;
    0.70 116.36 115.67  ];

IMC_solu_29 = [ 
    %0.90 157.03 154.64% add your data here --> Indomethacin
    0.85	153.76 153.22;
    0.80 148.94 148.17	;
    0.75 136.74 137.75	;
    0.70 125.47 125.28  ];

IMC_solu_34 = [ 
    0.90 155.75 154.56% add your data here --> Indomethacin
    0.85 149.49 148.24	;
    0.80 147.5 142.51	;
    0.75 131.47 131.29	;
    0.70 130.24 127.65  ];



%%
% fig1 = figure('Position',[20,511,1880,606]);%[left,bottom,width,heigth];
% fit1 = figure;
subplot(1,3,3)
box on,hold on
file = 'IMC_solu_29';   % This chunk switches bewtween the different drugs
% define your own values here and switch for different concentraitons
switch file
    
    case 'IMC_solu'     % data for IND are defined above as matrix - just run the name in terminal !!!   
        xDobs = IMC_solu(:,1)';xDobs = repmat(xDobs,1,2); % was 1,3
        xD = xDobs;
        Tmobs = IMC_solu(:,2:3)';Tmobs = [Tmobs(1,:) Tmobs(2,:)  ]; % had 3 Tmobs
        Tm = Tmobs;
        s = size(IMC_solu,1);
        garn = {'MarkerFaceColor','b','MarkerSize',8,'Marker','square','LineStyle','none','Color','b'};% garniture for plot
        scrsz = get(0,'ScreenSize');
        titel = file;
        h =    [1.1239e-03   1.7093e-02   8.0786e-02   2.3433e-01   1.1239e-03   1.7093e-02   8.0786e-02   2.3433e-01   ];
        % Drug parameters
        Hm = 118.64; % enthalpy of drug [J/g]
%         TmD = 160.11; % melting temperature of drug [K]
        TmD = 160.11;
        MwD = 357.79; % molecular weight of drug [g/mol]
        pD = 1.31; % density of drug [g/cm3]
        % Polymer parameters
        MwP = 118000; % molecular weight of polymer [g/mol]
        pP = 1.08; % density of polymer [g/cm3]
        lambda = (MwP/pP)/(MwD/pD); % molar volume ratio of the polymer and the drug
        
        case 'IMC_solu_25'     % data for IND are defined above as matrix - just run the name in terminal !!!   
        xDobs = IMC_solu_25(:,1)';xDobs = repmat(xDobs,1,2); % was 1,3
        xD = xDobs;
        Tmobs = IMC_solu_25(:,2:3)';Tmobs = [Tmobs(1,:) Tmobs(2,:)  ]; % had 3 Tmobs
        Tm = Tmobs;
        s = size(IMC_solu_25,1);
        garn = {'MarkerFaceColor','b','MarkerSize',8,'Marker','square','LineStyle','none','Color','b'};% garniture for plot
        scrsz = get(0,'ScreenSize');
        titel = file;
        h =    [1.1239e-03   1.7093e-02   8.0786e-02   2.3433e-01   1.1239e-03   1.7093e-02   8.0786e-02   2.3433e-01   ];
        % Drug parameters
        Hm = 118.64; % enthalpy of drug [J/g]
%         TmD = 160.11; % melting temperature of drug [K]
        TmD = 160.11;
        MwD = 357.79; % molecular weight of drug [g/mol]
        pD = 1.31; % density of drug [g/cm3]
        % Polymer parameters
        MwP = 118000; % molecular weight of polymer [g/mol]
        pP = 1.12; % density of polymer [g/cm3]
        lambda = (MwP/pP)/(MwD/pD); % molar volume ratio of the polymer and the drug
        
        case 'IMC_solu_29'     % data for IND are defined above as matrix - just run the name in terminal !!!   
        xDobs = IMC_solu_29(:,1)';xDobs = repmat(xDobs,1,2); % was 1,3
        xD = xDobs;
        Tmobs = IMC_solu_29(:,2:3)';Tmobs = [Tmobs(1,:) Tmobs(2,:)  ]; % had 3 Tmobs
        Tm = Tmobs;
        s = size(IMC_solu_29,1);
        garn = {'MarkerFaceColor','b','MarkerSize',8,'Marker','square','LineStyle','none','Color','b'};% garniture for plot
        scrsz = get(0,'ScreenSize');
        titel = file;
        h =    [1.1239e-03   1.7093e-02   8.0786e-02   2.3433e-01   1.1239e-03   1.7093e-02   8.0786e-02   2.3433e-01   ];
        % Drug parameters
        Hm = 118.64; % enthalpy of drug [J/g]
%         TmD = 160.11; % melting temperature of drug [K]
        TmD = 160.11;
        MwD = 357.79; % molecular weight of drug [g/mol]
        pD = 1.31; % density of drug [g/cm3]
        % Polymer parameters
        MwP = 118000; % molecular weight of polymer [g/mol]
        pP = 1.13; % density of polymer [g/cm3]
        lambda = (MwP/pP)/(MwD/pD); % molar volume ratio of the polymer and the drug
        
        case 'IMC_solu_34'     % data for IND are defined above as matrix - just run the name in terminal !!!   
        xDobs = IMC_solu_34(:,1)';xDobs = repmat(xDobs,1,2); % was 1,3
        xD = xDobs;
        Tmobs = IMC_solu_34(:,2:3)';Tmobs = [Tmobs(1,:) Tmobs(2,:)  ]; % had 3 Tmobs
        Tm = Tmobs;
        s = size(IMC_solu_34,1);
        garn = {'MarkerFaceColor','b','MarkerSize',8,'Marker','square','LineStyle','none','Color','b'};% garniture for plot
        scrsz = get(0,'ScreenSize');
        titel = file;
        h =    [1.1239e-03   1.7093e-02   8.0786e-02   2.3433e-01   1.1239e-03   1.7093e-02   8.0786e-02   2.3433e-01   ];
        % Drug parameters
        Hm = 118.64; % enthalpy of drug [J/g]
%         TmD = 160.11; % melting temperature of drug [K]
        TmD = 160.11;
        MwD = 357.79; % molecular weight of drug [g/mol]
        pD = 1.31; % density of drug [g/cm3]
        % Polymer parameters
        MwP = 118000; % molecular weight of polymer [g/mol]
        pP = 1.14; % density of polymer [g/cm3]
        lambda = (MwP/pP)/(MwD/pD); % molar volume ratio of the polymer and the drug
end

%%
% Volumic fraction of drug Eq. 4
vD = xD/pD./(xD/pD + (1-xD)/pP);  % calculates volume fractions form given weight fractions  

guess = -3; % why do we have to guess ? --> initial starting value, thus random (reproducible)
modelfun = @(Chi,v)(-1)./(R/(Hm*MwD)*(log(v)+ 1-v + Chi*(1-v).^2) - 1./(TmD+273.15));
mdl = fitnlm(vD,Tm+273.15,modelfun,guess);
%mdl = fitnlm(vD,Tm+273.15,modelfun,guess,'weights',1./h);

Chi = mdl.Coefficients.Estimate
residual = mdl.Residuals.Raw;
SE = mdl.Coefficients.SE;
%%
Text = 25:TmD-2; % extrapolated temperature range
Text = cat(2,Text,linspace(TmD-2,TmD,50));% refine resolution around melting temperature of drug

RunFH_MPD(Chi,Text,1:length(Text)); % calculate Vdrug for extrapolated temperature range
Xdrug = Vdrug*pD./((pD-pP)*Vdrug + pP);% Convert back to mass fraction x

df = length(Tm)-1;% Degrees of freedom

% ypred = predict(mdl,Xnew)

% sigma = sqrt(1/df*sum(residual.^2));% sample variance 1/(n-1)*sum((yi-mean(y))^2)
% SE = sigma*sqrt((full(jacobian)'*full(jacobian)))^(-1)   %Seber & Wild Eq. 1.11 sample variance/std.dev(Chi) similar to nlparci(Chi,residual,'jacobian',full(jacobian))
t975 = tinv(0.975,df);% = sqrt(finv(0.95,1,df))

%     CI = nlparci(Chi,residual,'jacobian',full(jacobian));
%     CIR = CI(2)-Chi;SE = CIR/t975
%     [Q,RR] = qr(full(jacobian)) % RR = sqrt(full(jacobian)'*full(jacobian))
%     sigma*abs(RR(1))^(-1)

options = optimset('Display','off','MaxFunEvals',10000,'MaxIter',30,'TolFun',1e-6);
%     i = 1;gl(1)=0.05;gu(1)=0.5;
%     for Tt = Text
%         vDl(i) = fsolve(@(vu)FH(vu,Chi+t975*SE,Tt),gl(i),options); %f(Vdrug,Chi,T)
%         vDu(i) = fsolve(@(vu)FH(vu,Chi-t975*SE,Tt),gu(i),options);
%         i = i+1;
%         gl(i)=vDl(i-1);gu(i)=vDu(i-1);
%     end
%     nlparci(Chi,residual,'jacobian',full(jacobian)) % Chi plus/minus
%     t975*SE*sqrt(df) = Chi plus/minus t975*sigma
i = 1;gl(1)=0.05;gu(1)=0.05;
for Tt = Text
    vDpl(i) = fsolve(@(vu)FH(vu,Chi+t975*SE*sqrt(df)*sqrt(1 + 1/df),Tt),gl(i),options); %f(Vdrug,Chi,T)
    vDpu(i) = fsolve(@(vu)FH(vu,Chi-t975*SE*sqrt(df)*sqrt(1 + 1/df),Tt),gu(i),options); % sigma^2(chi) = SE*sqrt(df) - this looks like a mistake!  t975*SE*sqrt(1+n-1)
    i = i+1;
    gl(i)=vDpl(i-1);gu(i)=vDpu(i-1);
end

% Convert back to mass fraction x
%     xDl = vDl*pD./((pD-pP)*vDl + pP);xDu = vDu*pD./((pD-pP)*vDu + pP);
xDpl = vDpl*pD./((pD-pP)*vDpl + pP);xDpu = vDpu*pD./((pD-pP)*vDpu + pP);
%% plot
%     subplot(1,3,MP),box on,hold on
%K=0.33/0.39;
%yfit = (XGT*43.24 + K*(1-XGT)*114.3)./(XGT + K*(1-XGT));%Gordon-Taylor
%plot(TgGT,XGT,'.','MarkerSize',24,'color',[0.8 0.8 0.8])
%plot(yfit,XGT,'-','color',[0.8 0.8 0.8])
% xlabel('$X_{drug}$','interpreter','latex');
% ylabel('$T_g$ [$^{\circ}$C]','interpreter','latex')

%% prediction interval
plot(Text,xDpu,'linestyle',':','linewidth',1,'color','k');plot(Text,xDpl,'linestyle',':','linewidth',1,'color','k');
%     plot(Text,xDu,'linestyle',':','linewidth',2,'color','k');plot(Text,xDl,'linestyle',':','linewidth',2,'color','k');

%% plot experimental data and fitted values
plot(Text,Xdrug,'k','linewidth',2)
% plot(T(1:s),mean([xD(1:s)',xD(s+1:2*s)',xD(2*s+1:end)'],2),garn{:})
plot(mean(IMC_solu_34(:,2:3),2),IMC_solu_34(:,1),garn{:}) % change the sample for each plot !!! 


% plot(Tm(1:s),xDobs(1:s)',garn{:},'color','g','MarkerFaceColor','g')
% plot(Tm(s+1:2*s),xDobs(s+1:2*s)',garn{:},'color','r','MarkerFaceColor','r')
% plot(Tm(2*s+1:end),xDobs(2*s+1:end)',garn{:},'color','b','MarkerFaceColor','b')

plot(TmD,1,garn{:})
% Vdrug drug determined from Gordon-Taylor Tg = (xD*TgD + K(1-xD)TgP)/(xD + K(1-xD))

axis([25 TmD 0 1]);

% axis tight
xlabel(['{\itT} [',char(176),'C]'],'interpreter','tex');
% ylabel('{\it X_d}','interpreter','tex')
set(gca,'YDir','reverse','XTick',[25 50 75 100 125 150],...
    'TickDir','out');
% title(titel,'FontSize',20,'FontWeight','bold');
% title(['$\chi = $',num2str(Chi,3),' $\pm$ ',num2str(t975*SE*sqrt(df)*sqrt(1+1/df),2),', $X_{d}(25^{\circ}$C) = ',num2str(Xdrug(1),2)...
%     ,' [',num2str(xDpl(1),2),',',num2str(xDpu(1),2),']'],'interpreter','latex')


% Extract axes handles of all subplots from the figure
% axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes')
% Set the axis property to square
% axis(axesHandles,'square')

addpath('D:\OneDrive - Københavns Universitet\Pharmaceutical_Sciences\Master Thesis\Matlab_scripts\Irish Collaboration\')
% export_fig FHnew2.tif -r600 -painters

%     plot(T(1:s),residual(1:s),'ko');plot(T(s+1:2*s),residual(s+1:2*s),'ro');plot(T(2*s+1:end),residual(2*s+1:end),'bo')
%     [h,p,stats]= chi2gof(residual)
%     plot(T(1:s),xD(1:s),garn{:},'color','k');plot(T(s+1:2*s),xD(s+1:2*s),garn{:},'color','r');plot(T(2*s+1:end),xD(2*s+1:end),garn{:},'color','b')
%
% % plot residuals vs. independent variable
% figure,box on,hold on
%
% Tsave = T(1:s);rsave = mean([residual(1:s)',residual(s+1:2*s)',residual(2*s+1:end)'],2)';save(['Tsave',num2str(MP)],'Tsave');save(['rsave',num2str(MP)],'rsave')
%
% plot(T(1:s),mean([residual(1:s)',residual(s+1:2*s)',residual(2*s+1:end)'],2),garn{:})
% plot([min(T) max(T)],[0 0],'k')
% xlabel('$T_{a}$ [$^{\circ}$C]','interpreter','latex');
% ylabel('$X_{IMC}-\hat{X}_{IMC}$','interpreter','latex')
% xlabel(['{\itT}_a [',char(176),'C]'],'interpreter','tex');
% ylabel('Residuals')
% title('Residuals vs. independent variable')
% set(gca,'YDir','reverse');


% set(gcf,'position',[-14.5410 7.7200 50.0740 14.2502])



% export the result as csv for further plotting and analysis

result_ex = [Text; Xdrug; xDpu; xDpl].'


writematrix(result_ex,'solu_0.34_export.csv')% change for each sample !!!!
