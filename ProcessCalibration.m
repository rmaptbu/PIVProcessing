clear;clc;
basepaths={'G:\New\PIV\Calibration\A\';
            'G:\New\PIV\Calibration\N\'};
H=[0.05 0.1 0.15 0.2 0.25 0.3 0.4]; % HF - feed hematocrit
Hall=zeros(1,21);Hall(1:3:end)=H;Hall(2:3:end)=H;Hall(3:3:end)=H; % sort into an array for all cases
Hp=0.99*Hall;                       % corrected for cell packing
dt =1e-3*reshape([ones(1,7);5*ones(1,7);3*ones(1,7)],1,[]);
%% Correct for cell screening
for Agg=1:2;
    clear prof;
    basepath=char(basepaths(Agg));
    cd(basepath);files=struct2cell(dir);files=files(1,3:end);

    X=[0 (8:4:76)/84 1]-0.5;
    XX=(0:4:84)/84-0.5;
    mmPpx=50/84*1e-3;
    for i=1:length(files)
        disp(['Reading Velocities ',char(files(i))])
        filepath=[basepath,char(files(i))]; 
        data=dlmread(filepath);
        [xpx,ypx,upx,vpx]=jvc2mat(data,1,1);            % import pixel data
        [x,y,u,v]=jvc2mat(data,mmPpx,dt(i));          % import mm/s data
        u=u(2:end,:);v=v(2:end,:);
        [residuals,invalid]=normres(u,v,1,0.01);      % normalised residual test
        ur=repmedian(u,invalid)*1.5;                    % replace invalid vectors 
        vr=repmedian(v,invalid)*1.5;                    % and correct for underestimation Poelma2012)
        
        trim=3; % number of data points for extrapolation

        uProfC=mean(ur,2);                % inner bit of profile 
        sProfC=std(ur,[],2)*1;          % get std along the profile                  
        USym=(uProfC+flipud(uProfC))/2; %  build the symmetric profile
        dUSym=sqrt(sProfC.^2+flipud(sProfC).^2);  % add the errors in quadrature
        [USYM,dUSYM,~,~,~,XH]=FitCoshU2(X(2:end-1)',USym,dUSym.^-2);
        USYMN=USYM/max(USYM);
        dUSYMN=dUSYM/max(USYM);
        U3=USYM*USYMN';
        dU3=dUSYM*dUSYMN';
        Umean(i,Agg)=trapz(XH,trapz(XH,U3));
        dUmean(i,Agg)=trapz(XH,trapz(XH,dU3));
    end 
end
UD=Umean/0.05;
dUDu=(Umean+dUmean)/0.05;
dUDl=(Umean-dUmean)/0.05;

HD(1:21,1)=Hall'.*(0.878+0.042*log10(UD(1:21,1)));     % correct for screening - gaehtgens
HD(1:21,2)=Hall'.*(0.901+0.029*log10(UD(1:21,2)));     % correct for screening - gaehtgens

HDl(1:21,1)=Hall'.*(0.878+0.042*log10(dUDl(1:21,1)));     % correct for screening - gaehtgens
HDl(1:21,2)=Hall'.*(0.901+0.029*log10(dUDl(1:21,2)));     % correct for screening - gaehtgens

HDu(1:21,1)=Hall'.*(0.878+0.042*log10(dUDu(1:21,1)));     % correct for screening - gaehtgens
HDu(1:21,2)=Hall'.*(0.901+0.029*log10(dUDu(1:21,2)));     % correct for screening - gaehtgens

(HD-HDl)./HD*100  % Error in velocity makes negligible difference to estimate in x axis.

%% Correct for Fahraeus effect - Pries
d=50;
HT=HD.*(HD+(1-HD).*(1+1.7*exp(-0.415*d)-0.6*exp(-0.011*d)));
clearvars -except HT Hall

%%  Generate IStar Profiles
basepaths={'G:\New\Images\Calibration\A\';
            'G:\New\Images\Calibration\N\'};
Trim=4;
X=(0.5:83.5)/84-0.5;XT=X(Trim+1:end-Trim);
I=200;
o=6;
D=84;
% theory is that instensity AVERAGED OVER TIME at a location is
% proportional to haematocrit. so it is ok to average over time.
for Agg=1:2
    basepath=char(basepaths(Agg));
    Filespath=[basepath,'JustChannel\'];
    cd(Filespath);
    Files=struct2cell(dir);Files=Files(1,3:23);
    for i=1:length(Files)
        [Agg i]
        bgIm=imread([basepath,'PreProcess\Mean\',char(Files(i)),'_mean.tif']);  % import image with bits above and below
        bgOut=bgIm([1:ceil(0.4*D),ceil(1.6*D+1):2*D],:);                % extract bits away from the channel
%         bgImS=imread([basepath,'PreProcess\Mad\',char(Files(i)),'_mad.tif']);  % import image with bits above and below
%         bgOutS=bgImS([1:ceil(0.4*D),ceil(1.6*D+1):2*D],:);                % extract bits away from the channel
        MaxPx(i)=mean(bgOut(:));                                        % use the average for maxpx
        MaxPxS(i)=std(double(bgOut(:)));%/sqrt(numel(bgOut(:)));                % use the average for maxpx
        e(i)=MaxPxS(i)/MaxPx(i);
        cd([Filespath,char(Files(i))]) 
        Images=struct2cell(dir);Images=Images(1,3:end);   
        Int=zeros(84,4*D,I);
        for j=1:I                                   %length(Images)              
            Im=double(imread(char(Images(j))))*4095;     % import image and normalise
            Ic=Im(1:end,3*D+1:7*D);       % crop to smaller region and trim edges
            Int(:,:,j)=Ic;                          % store in one big array        
        end
        IProfMean=mean(Int(Trim+1:end-Trim,:),2);      % Mean Profile
        IProfStd=std(Int(Trim+1:end-Trim,:),[],2);     % Std Profile
        IProfSEM=IProfStd./sqrt(numel(Int(1,:,:)));      % SEM Profile std/root(n)
               
        IProfC=IProfMean/MaxPx(i);                       % profile based only on means
        dIProfC=sqrt( (IProfSEM/MaxPx(i)).^2 + (IProfMean.*MaxPxS(i)/(MaxPx(i).^2)).^2);
        %% for all - was previously i<12 , i.e. for agg only
        if i<25;        [IProfC,dIProfC]=WallFix2(IProfC,dIProfC,X);end  % fix for low hct cases where tiny blurring inverts profile - AGG
        %%
        IProfL=polyval(polyfit(XT(1:Trim)',IProfC(1:Trim),1),X(1:Trim)); % extrapolate LHS
        IProfR=polyval(polyfit(XT(end-Trim:end)',IProfC(end-Trim:end),1),X(end-Trim+1:end));
        IProf=[IProfL'; IProfC; IProfR']/4095;
        dIProfL=polyval(polyfit(XT(1:Trim)',dIProfC(1:Trim),1),X(1:Trim)); % extrapolate LHS
        dIProfR=polyval(polyfit(XT(end-Trim:end)',dIProfC(end-Trim:end),1),X(end-Trim+1:end));
        dIProf=[dIProfL';dIProfC; dIProfR']/4095;
        
        IProf1=coerce1(IProf,0,1);
        dIProf1=dIProf;
        dIProf1(IProf1==1)=0;
        
        IsProf(:,i,Agg)=1-IProf1;
        dIsProf(:,i,Agg)=dIProf1;
    end
end

%% We now have Istar profiles and average haematocrit for each case 

clearvars -except HT IsProf dIsProf X

mI(:,:)=trapz(X,IsProf,1);
sI(:,:)=trapz(X,dIsProf,1);
errorbar(HT,mI,sI,'.')
save('G:\New\MatlabProcessing\PreProcessImages\calib')

clear
load('G:\New\MatlabProcessing\PreProcessImages\calib')

%% 
r=10;
for Agg=1:3 % look at agg only, pbs and combined data sets
    if Agg==3 % Combined agg and pbs data
        IProf=[IsProf(:,:,1) IsProf(:,:,2)];
        dIProf=[dIsProf(:,:,1) dIsProf(:,:,2)];
        Ht=HT(:);
    else
    IProf=IsProf(:,:,Agg);
    dIProf=dIsProf(:,:,Agg);
    Ht=HT(:,Agg);
    end
    clear Iguess1
    for i=1:numel(Ht)
        Iguess1(i)=trapz(X,IProf(:,i));  % guess linear relationship i.e. mean I is mean H
    end
    
    [p1,R]=FitExpCalib(Ht,Iguess1');     % guess f1 to convert 
    a0=0.001*round(p1(1)*1000);  % guess inial a
    b0=0.001*round(p1(2)*1000);  % guess initial b   
    A0(Agg)=a0;
    B0(Agg)=b0;
    
    clear dH;close all
    n=20;
    for res=1:3
        A=1;
        ares=10^-res;amin=a0-n/2*ares;amax=a0+n/2*ares;A=1;
        bres=10^-res;bmin=b0-n/2*bres;bmax=b0+n/2*bres;B=1;
        for a=amin:ares:amax  
            ac=coerce1(a,0,amax);
            for b=bmin:bres:bmax   
                Hg=zeros(numel(Ht),1);
                clear dHg
                for i=1:numel(Ht)                           % calculate predicted haematocrits
                    Hprof=ApplyCalib([ac b],IProf(:,i));
                    dHprof=ApplyCalib([ac b],(IProf(:,i)+dIProf(:,i)))-Hprof; % differences small enough to consider symmetric
                    Hg(i)=(trapz(X,Hprof));
                    dHg(i)=trapz(X,dHprof);
                end
                if RootMS(Hg-Ht)<r
                    figure(res); 
                    errorbar(Ht,(Hg-Ht)./Ht,(dHg')./Ht,(dHg')./Ht,'.') 
                    r=RootMS(Hg-Ht);
                end
                 dH(A,B)=RootMS((Hg-Ht)./dHg');   % weight errors with std
%                  dH(A,B)=RootMS((Hg-Ht));   % 
                
                B=B+1;
            end
            A=A+1;
            B=1;
        end
        [dhf,DH]=min(dH(:));
        [row,col]=ind2sub(size(dH),DH);
        bb=bmin+bres*(col-1);b0=bb;
        aa=amin+ares*(row-1);a0=aa;
        
        DHF(Agg,res)=dhf;
        AAA(Agg,res)=aa;
        BBB(Agg,res)=bb;
        R(Agg,res)=r;
        r=10;
    end
    DHfinal=dhf;
    AA(Agg)=aa;
    BB(Agg)=bb;
end

%% with errors
%AA=[0.7280    0.6590    0.6850]
%BB=[7.5790   10.5000    9.2440]
%% without errors
%AA=[0.7410    0.6660    0.6950]
%BB=[7.5890   11.3400    9.5940]
 clear Hg dh dHg
    %% Therefore select values from combined fit
    col='brm';
XH=(0.5:83.5)/84-0.5;
close all
f=figure(1);
ha = tight_subplot(1,1,[.04 .04],[.15 .05],[.15 .04]);
hold all
errorbar(-1,-1,0.01,'o','markersize',4,'Color',col(1),'MarkerFaceColor',col(1));
errorbar(-1,-1,0.01,'o','markersize',4,'Color',col(2),'MarkerFaceColor',col(2));

            
for t=1:3;  % each set 
    clear Hg dHg
    col='brm';
    a=AA(t);  %
    b=BB(t);  % 
    %% All
    for Agg=1:2
        for i=1:size(HT,1)     % calculate predicted haematocrits
            Hprof=ApplyCalib([a b],IsProf(:,i,Agg));
            dHprof=ApplyCalib([a b],(IsProf(:,i,Agg)-dIsProf(:,i,Agg)))-Hprof; % differences small enough to consider symmetric
            Hg(i,Agg)=trapz(X,Hprof);       % calculate hct
            dHg(i,Agg)=trapz(X,dHprof);   % dhct
        end
    end
    e(:,t)=Ht(:)-Hg(:);
    [E(t),~,S(t)]=Wmean(Ht(:)-Hg(:),dHg(:).^-2,1);
    %% Zero
    a0=A0(t);  %
    b0=B0(t);  % 
    for Agg=1:2
        for i=1:size(HT,1)     % calculate predicted haematocrits
            Hprof=ApplyCalib([a0 b0],IsProf(:,i,Agg));
            dHprof=ApplyCalib([a0 b0],(IsProf(:,i,Agg)-dIsProf(:,i,Agg)))-Hprof; % differences small enough to consider symmetric
            Hg0(i,Agg)=trapz(X,Hprof);       % calculate hct
            dHg0(i,Agg)=trapz(X,dHprof);   % dhct
        end
    end
    e0(:,t)=Ht(:)-Hg(:);
    [E0(t),~,S0(t)]=Wmean(Ht(:)-Hg0(:),dHg0(:).^-2,1);
    %%
   
    I=0:0.01:a;
    Hfit=ApplyCalib([a b],I);
    Hfit0=ApplyCalib([a0 b0],I);
    plot(Hfit0,I,'k','LineWidth',2,'color',0.7*[1 1 1])
    plot(Hfit,I,'k','LineWidth',2);hold all
    for Agg=1:2         % plot
        for i=1:21
            II=trapz(XH,IsProf(:,i,Agg));
            dII=trapz(XH,dIsProf(:,i,Agg));
            h=errorbar(HT(i,Agg),II,1.96*dII,'o','markersize',4,'Color',col(Agg),'MarkerFaceColor',col(Agg));
        errorbar_tick(h,0.005,'units')
        end
    end
    
    
    axis([ 0 0.35 0 0.7])
end 
    
fs=15;
ylabel('I^*','FontSize',fs);
xlabel('H_C','FontSize',fs);

set(gca,'xticklabel',0:0.1:0.3,'xtick',0:0.1:0.3)
set(gca,'yticklabel',0:0.2:0.6,'ytick',0:0.2:0.6)
set(gca,'fontsize',fs)

l=legend('Dex','PBS')
set(l,'location','northwest')
box on

set(f,'Position',[50 50 525 500],'color','w')

savecp=['G:\New\Images\Calibration\FigS1_Calibration.png'];
export_fig f1 savefig.png '-r600' -nocrop
movefile('savefig.png',savecp)

[~,p]=ttest2(e(:,1),e(:,3)) ;  
[~,p]=ttest2(e(:,2),e(:,3));     

a=AA(3);
b=BB(3);

S(3)
S0(3)







    
    
    
    
    
    
    
    
    
    
    