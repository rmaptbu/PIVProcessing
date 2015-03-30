clear;clc;close all;
%
load('G:\New\Mats\Analytical');

X=[-0.5 (8:4:76)/84-0.5 0.5];                                                   % initial vector positions
XR=X(2:end-1);
XT=[-0.5 XR 0.5];
XX=(-25:25)/50;                                                                 % micron positions
XXX=(-5000:5000)/10000;
XH=(0.5:83.5)/84-0.5;
VtH=spline(Xt,VtN,XH);
Vt=normalise(VtH,1,XH);
dt_agg=1e-3*[4*ones(9,1); 2*ones(10,1); ones(10,1); 8*ones(6,1); 4*ones(9,1)];  % dt for agg cases
dt_pbs=1e-3*[4*ones(10,1); 2*ones(8,1); ones(10,1); 8*ones(10,1)];              % dt for pbs cases
mmPpx=50/84*1e-3;                                                               % calibration constant
ro=[1125 1025];
BraNam=['D1';'D2';'O1';'O2';'Pa'];
AggNam='AN';
LigNam='RS';
col='br';
CI=1; % to give one std
%%

basepath='G:\New\PIV\';
for Agg=1:2
    if Agg==1;A='A';dt=dt_agg;                                                  % define dt
    else A='N';dt=dt_pbs; end                                                   % define dt
    aggpath=['G:\New\PIV\',A];                                                  % set path 
    hctpath=['G:\New\HctImages\',A]; 
    hctprofpath=['G:\New\HctProfiles\',A];
    cd(hctprofpath);hctprof=struct2cell(dir);hctprof=hctprof(1,3:end);          % go to the path and get filnames
    cd(aggpath);branches=struct2cell(dir);branches=branches(1,3:end);           % go to the path and get filnames
    for Bra=[5 1 2 3 4]                                                         % for each branch
        branchpath=[aggpath,'\',char(branches(Bra))];                           % select a branch
        cd(branchpath)                                                          % go to the branch
        if Bra==5;Range=15:24;  % flow in, so 1:1.5 diameters
        elseif Bra==3;Range=41:50;           % central so 2.25-2.75    
        else Range=56:65; end     %from 3 to 3.5 diameters                      % this defines the averaging region. For parent, take start, for all others, take end
        for Lig=1:2                                                             % for each illumination method
            k=1;
            if Lig==1,L='Str';Mark='x';Col='k';
            else L='Las';Mark='+';Col='r';end
            lightpath=[branchpath,'\',L];                                       % go to the illumination method
            cd(lightpath)
            cases=struct2cell(dir);cases=cases(1,3:end);                        % list cases
            for Cas=1:length(cases)       % for each flow case
                %% import velocity data and clean up
                data=dlmread(char(cases(Cas)));                                   % import data
%                       mmPpx=1;dt=ones(size(dt));                                  % uncomment this line to analyse px/fr instead of mm/s
                [x,y,u,v]=jvc2mat(data,mmPpx,(dt(Cas)));                          % convert column format to array and calibrate
                u=u(2:end,Range);v=v(2:end,Range);                                % trim top row
                [~,invalid]=normres(u,v,1,0.1);                                   % Normalised median test - strict (normal is 2,0.1
                ur=repmedian(u,invalid)*1.5;                                      % replace invalid u components with median
                vr=repmedian(v,invalid)*1.5;                                      % also multiply by 1.5- poelma "correction"
                percinvalid(Cas,Lig,Bra,Agg)=sum(invalid(:))/numel(invalid)*100;  % record the percentage of invalid vectors
                  %% velocity profiles 
                uProfC=mean(ur,2);      % inner bit of profile
                sProfC=std(ur,[],2)*CI;          % get std along the profile           

                trim=3; % number of data points for extrapolation
                    if Lig==2         % for SM, set wall velocity to zero.
                    uprof=[0; uProfC; 0]';
                    elseif Lig==1     % for RBC 
                    uProfL=polyval(polyfit(X(2:trim+1)',uProfC(1:trim),1),X(1)); % extrapolate LHS
                    uProfR=polyval(polyfit(X(end-trim:end-1)',uProfC(end-trim+1:end),1),X(end)); % extrapolate LHS
                    uprof=[uProfL uProfC' uProfR] ;     %  note that this u, and the X array go  0 2 3 4 5....               
                    end
                sProfL=polyval(polyfit(X(2:trim+1)',sProfC(1:trim),1),X(1)); % extrapolate LHS
                sProfR=polyval(polyfit(X(end-trim:end-1)',sProfC(end-trim+1:end),1),X(end)); % extrapolate LHS

                sprof=[sProfL sProfC' sProfR] ;
                  
                Uprof=spline(X,uprof,XH);       % interpolate velocity onto haematcrit resolution
                Sprof=spline(X,sprof,XH);       % interpolate std onto haematcrit resolution

                UPROFfit(:,Cas,Lig,Bra,Agg)=Uprof;    % save the haematocrit resolution profile
                dUPROFfit(:,Cas,Lig,Bra,Agg)=Sprof;    % save the haematocrit resolution profile
                UPROFraw(:,Cas,Lig,Bra,Agg)=uprof;    % save the raw profile
                dUPROFraw(:,Cas,Lig,Bra,Agg)=sprof;    % save the raw std profile
                  
                if Bra==5                     
                    USym=(uProfC+flipud(uProfC))/2; % if it's the parent branch, build the symmetric profile
                    dUSym=sqrt(sProfC.^2+flipud(sProfC).^2);  % add the errors (95CI) in quadrature
                    if Lig==1;
                        [USymF,dUSymF,~,~,~,~,~,PIV(Cas,Lig,Agg)]=FitCoshU2(X(2:end-1)',USym,dUSym.^-2);
                        PIV(Cas,Lig,Agg)
                    elseif Lig==2;
                        [USymF,dUSymF]=FitCoshU(X(2:end-1)',USym,dUSym.^-2);
                    end
                    [Umax(Cas,Lig,Agg),im]=max(USymF);
                    dUmax(Cas,Lig,Agg)=dUSymF(im);
                    USYM(:,Cas,Lig,Agg)=USymF;      % normalise by max
                    dUSYM(:,Cas,Lig,Agg)=dUSymF;                     % normalise by max velocity
                    efitu(Cas,Lig,Agg)=RootMS(USym-spline(XH,USymF,X(2:end-1))')/trapz(XH,USymF)*100;
                end
                  
                %% z plane velocity distribution
                for l=1:length(XH)                                    % for each element in the profile    
                    u3(:,l)=Uprof(l)*USYM(:,Cas,Lig,Agg)/Umax(Cas,Lig,Agg);              % calculate the y/z profile
                    du3(:,l)=1./Umax(Cas,Lig,Agg).*sqrt((USYM(:,Cas,Lig,Agg).*Sprof(l)).^2+...
                                                    (Uprof(l).*dUSYM(:,Cas,Lig,Agg)).^2+...
                                                    (dUmax(Cas,Lig,Agg).*Uprof(l).*USYM(:,Cas,Lig,Agg)/Umax(Cas,Lig,Agg)).^2); % estimate the error on each value
%                     errorbar(XH,Uprof,(USYM(:,Cas,Lig,Agg).*Sprof(l)),'.');hold all
%                     errorbar(XH,Uprof,(Uprof(l).*dUSYM(:,Cas,Lig,Agg)),'.')
%                     errorbar(XH,Uprof,dUmax(Cas,Lig,Agg).*Uprof(l).*(USYM(:,Cas,Lig,Agg))/Umax(Cas,Lig,Agg),'.')     
                    
                    um3(l)=trapz(XH,u3(:,l)); 
                    dum3(l)=trapz(XH,du3(:,l)); 
                end
                U3(:,:,Cas,Lig,Bra,Agg)=u3;
                dU3(:,:,Cas,Lig,Bra,Agg)=du3;
                Umean(Cas,Lig,Bra,Agg)=trapz(XH,um3);                           % save average velocity - because 3D weighting has mean of
                dUmean (Cas,Lig,Bra,Agg)=trapz(XH,dum3);   
                %% z plane velocity distribution - exluding extrapolated region
                r=9:75;
                for l=r%1:length(XH)      % for each element in the profile    
                    u3r(:,l-r(1)+1)=Uprof(l)*USYM(r,Cas,Lig,Agg)/Umax(Cas,Lig,Agg);              % calculate the y/z profile
                    du3r(:,l-r(1)+1)=1./Umax(Cas,Lig,Agg).*sqrt((USYM(r,Cas,Lig,Agg).*Sprof(l)).^2+...
                                                    (Uprof(l).*dUSYM(r,Cas,Lig,Agg)).^2+...
                                                    dUmax(Cas,Lig,Agg).*Uprof(l).*(USYM(r,Cas,Lig,Agg))); % estimate the error on each value
                    um3r(l-r(1)+1)=trapz(XH(r),u3(r,l)); 
                    dum3r(l-r(1)+1)=trapz(XH(r),du3(r,l)); 
                end
 
                UmeanR(Cas,Lig,Bra,Agg)=trapz(XH(r),um3r);                           % save average velocity - because 3D weighting has mean of
                dUmeanR(Cas,Lig,Bra,Agg)=trapz(XH(r),dum3r); 
             end % velocity cases
        end % lighting
        %% Import haematocrit 
        cd([hctpath,'\',char(branches(Bra))])
        cases=struct2cell(dir);cases=cases(1,3:end);                        % list cases
        hcases=cases(1:2:end);
        dhcases=cases(2:2:end);
        
        cd([hctprofpath])
        cases=struct2cell(dir);cases=cases(1,3:end);                        % list cases
        hcasesprof=cases(1:2:end);
        dhcasesprof=cases(2:2:end);
        
        D=84;
        if Bra==5;Range=D:1.5*D-1; % 1:1.5D
        elseif Bra==3, Range=2.25*D:2.75*D-1; % 2.25:2.75D
        else Range=3*D:3.5*D-1;     %3:3.5D    
        end                             % this defines the averaging region. For parent, take start, for all others, take end
        for Cas=1:length(hcases)       % for each flow case
            H=double(imread([hctpath,'\',char(branches(Bra)),'\',char(hcases(Cas))]))*0.0001;H=H';
            dH=double(imread([hctpath,'\',char(branches(Bra)),'\',char(dhcases(Cas))]))*0.0001;dH=dH';
            %% get weighted mean
%             for k=1:size(H,1)
%                 Hprof(k)=
            [HProf,~,dHProf]=Wmean2(H(:,Range),(dH(:,Range)+0.00001).^-2,1,1);
            HPROF(:,Cas,Bra,Agg)=coerce1(HProf,0,1);
            dHPROF(:,Cas,Bra,Agg)=coerce1(dHProf,0,1);
            Hct(Cas,Bra,Agg)=trapz(XH,HProf);   % 2D haematocrit
            dHct(Cas,Bra,Agg)=trapz(XH,dHProf);
            
                 %%
            HSYM=dlmread([hctprofpath,'\',char(hcasesprof(Cas))]);
            dHSYM=dlmread([hctprofpath,'\',char(dhcasesprof(Cas))])/1.96;  % as confidence intervals in fit
             
            for l=1:length(XH)                                    % for each element in the profile    
                h3(:,l)=HProf(l)*HSYM./Hct(Cas,5,Agg);                     % calculate the y/z profile
                dh3(:,l)=1./Hct(Cas,Bra,Agg)*sqrt((HSYM.*dHProf(l)).^2+(HProf(l).*dHSYM).^2+(dHct(Cas,Bra,Agg)/Hct(Cas,Bra,Agg).*HProf(l).*HSYM).^2); % estimate the error on each value
                hm3(l)=trapz(XH,h3(:,l)); 
                dhm3(l)=trapz(XH,dh3(:,l)); 
            end
            H3(:,:,Cas,Bra,Agg)=h3;
            dH3(:,:,Cas,Bra,Agg)=dh3;
            Hmean(Cas,Bra,Agg)=trapz(XH,hm3);                           % save average velocity - because 3D weighting has mean of
            dHmean (Cas,Bra,Agg)=trapz(XH,dhm3);
               

        end
        %% Fluxes
        for Cas=1:length(hcases) 
            Haem(:,:)=H3(:,:,Cas,Bra,Agg);
            dHaem(:,:)=dH3(:,:,Cas,Bra,Agg);
            for Lig=1:2
                Vel(:,:)=U3(:,:,Cas,Lig,Bra,Agg);
                dVel(:,:)=dU3(:,:,Cas,Lig,Bra,Agg);
                if Lig==1;
                    Flux=Haem.*Vel;
                    dFlux=sqrt((dHaem.*Vel).^2+(dVel.*Haem).^2);
                else Flux=(1-Haem).*Vel;
                    dFlux=sqrt((dHaem.*Vel).^2+(dVel.*(1-Haem)).^2);
                end
                for l=1:length(XH)
                    f(l)=trapz(XH,Flux(:,l));
                    df(l)=trapz(XH,dFlux(:,l));
                end
                F(Cas,Lig,Bra,Agg)=trapz(XH,f);
                dF(Cas,Lig,Bra,Agg)=trapz(XH,df);
                
                r=9:75;
                for l=r%1:length(XH)                                    % for each element in the profile   
                    fr(l)=trapz(XH(r),Flux(r,l));
                    dfr(l)=trapz(XH(r),dFlux(r,l));
                end
                Fr(Cas,Lig,Bra,Agg)=trapz(XH(r),fr(r));
                dFr(Cas,Lig,Bra,Agg)=trapz(XH(r),dfr(r));
              
            end
            Ft(Cas,Bra,Agg)=F(Cas,1,Bra,Agg)+F(Cas,2,Bra,Agg); 
            dFt(Cas,Bra,Agg)=sqrt(dF(Cas,1,Bra,Agg)^2+dF(Cas,2,Bra,Agg)^2);
            
            Ftr(Cas,Bra,Agg)=Fr(Cas,1,Bra,Agg)+Fr(Cas,2,Bra,Agg); 
            dFtr(Cas,Bra,Agg)=sqrt(dFr(Cas,1,Bra,Agg)^2+dFr(Cas,2,Bra,Agg)^2);
            
            %% Relative and Total Velocity Profiles
            rbcprof=UPROFraw(2:end-1,Cas,1,Bra,Agg);%/Umean(Cas,1,Bra,Agg);     % ignore zeros at beginning and end
            drbcprof=dUPROFraw(2:end-1,Cas,1,Bra,Agg);%/Umean(Cas,1,Bra,Agg);
            smprof=UPROFraw(2:end-1,Cas,2,Bra,Agg);%/Umean(Cas,2,Bra,Agg);
            dsmprof=dUPROFraw(2:end-1,Cas,2,Bra,Agg);%/Umean(Cas,2,Bra,Agg);
            hctprof=HPROF(:,Cas,Bra,Agg);hctprof_r=spline(XH,hctprof,XR)';
            dhctprof=dHPROF(:,Cas,Bra,Agg);dhctprof_r=spline(XH,dhctprof,XR)';
            %
            RPROFraw(:,Cas,Bra,Agg)=rbcprof-smprof;
            dRPROFraw(:,Cas,Bra,Agg)=sqrt(drbcprof.^2+dsmprof.^2);
            TPROFraw(:,Cas,Bra,Agg)=rbcprof.*hctprof_r+smprof.*(1-hctprof_r);
            dTPROFraw(:,Cas,Bra,Agg)=sqrt((drbcprof.*hctprof_r).^2+(rbcprof.*dhctprof_r).^2+...
                        (dsmprof.*(1-hctprof_r)).^2+(smprof.*dhctprof_r).^2); 
                   
            rbcprof=UPROFfit(:,Cas,1,Bra,Agg);%/Umean(Cas,1,Bra,Agg);     % ignore zeros at beginning and end
            drbcprof=dUPROFfit(:,Cas,1,Bra,Agg);%/Umean(Cas,1,Bra,Agg);
            smprof=UPROFfit(:,Cas,2,Bra,Agg);%/Umean(Cas,2,Bra,Agg);
            dsmprof=dUPROFfit(:,Cas,2,Bra,Agg);%/Umean(Cas,2,Bra,Agg);
            hctprof=HPROF(:,Cas,Bra,Agg);
            dhctprof=dHPROF(:,Cas,Bra,Agg);

            RPROFfit(:,Cas,Bra,Agg)=rbcprof-smprof;
            dRPROFfit(:,Cas,Bra,Agg)=sqrt(drbcprof.^2+dsmprof.^2);
            TPROFfit(:,Cas,Bra,Agg)=rbcprof.*hctprof+smprof.*(1-hctprof);
            dTPROFfit(:,Cas,Bra,Agg)=sqrt((drbcprof.*hctprof).^2+(rbcprof.*dhctprof).^2+...
                        (dsmprof.*(1-hctprof)).^2+(smprof.*dhctprof).^2);
            
        end 
    end % branches
end % aggregating/non aggregating
%Umean(Cas,Lig,Bra,Agg)

[Q1dr(:,:),dQ1dr(:,:)]=DivideError(Umean(:,1,1,:),Umean(:,1,5,:),dUmean(:,1,1,:),dUmean(:,1,5,:));
[Q1or(:,:),dQ1or(:,:)]=DivideError(Umean(:,1,3,:),Umean(:,1,5,:),dUmean(:,1,3,:),dUmean(:,1,5,:));
[Q2dr(:,:),dQ2dr(:,:)]=DivideError(Umean(:,1,2,:),Umean(:,1,3,:),dUmean(:,1,2,:),dUmean(:,1,3,:));
[Q2or(:,:),dQ2or(:,:)]=DivideError(Umean(:,1,4,:),Umean(:,1,3,:),dUmean(:,1,4,:),dUmean(:,1,3,:));

[Q1d(:,:),dQ1d(:,:)]=DivideError(Ft(:,1,:),Ft(:,5,:),dFt(:,1,:),dFt(:,5,:));
[Q1o(:,:),dQ1o(:,:)]=DivideError(Ft(:,3,:),Ft(:,5,:),dFt(:,3,:),dFt(:,5,:));
[Q2d(:,:),dQ2d(:,:)]=DivideError(Ft(:,2,:),Ft(:,3,:),dFt(:,2,:),dFt(:,3,:));
[Q2o(:,:),dQ2o(:,:)]=DivideError(Ft(:,4,:),Ft(:,3,:),dFt(:,4,:),dFt(:,3,:));

[H1d(:,:),dH1d(:,:)]=DivideError(Hmean(:,1,:),Hmean(:,5,:),dHmean(:,1,:),dHmean(:,5,:));
[H1o(:,:),dH1o(:,:)]=DivideError(Hmean(:,3,:),Hmean(:,5,:),dHmean(:,3,:),dHmean(:,5,:));
[H2d(:,:),dH2d(:,:)]=DivideError(Hmean(:,2,:),Hmean(:,3,:),dHmean(:,2,:),dHmean(:,3,:));
[H2o(:,:),dH2o(:,:)]=DivideError(Hmean(:,4,:),Hmean(:,3,:),dHmean(:,4,:),dHmean(:,3,:));

Ur(:,:,:)=UmeanR(:,1,:,:)-UmeanR(:,2,:,:);      %% 
dUr(:,:,:)=sqrt(dUmeanR(:,1,:,:).^2+dUmeanR(:,2,:,:).^2);
[Fnr(:,:,:),dFnr(:,:,:)]=DivideError(Ur,Ftr,dUr,dFtr);

Ure(:,:,:)=Umean(:,1,:,:)-Umean(:,2,:,:);
dUre(:,:,:)=sqrt(dUmean(:,1,:,:).^2+dUmean(:,2,:,:).^2);
[Fne(:,:,:),dFne(:,:,:)]=DivideError(Ure,Ft,dUre,dFt);

for Agg=1:2
    Qt(:,:,Agg)=[Q1d(:,Agg) Q2d(:,Agg) Q1o(:,Agg) Q2o(:,Agg)];
    dQt(:,:,Agg)=[dQ1d(:,Agg) dQ2d(:,Agg) dQ1o(:,Agg) dQ2o(:,Agg)];
end
for Agg=1:2
    Ht(:,:,Agg)=[H1d(:,Agg) H2d(:,Agg) H1o(:,Agg) H2o(:,Agg)];
    dHt(:,:,Agg)=[dH1d(:,Agg) dH2d(:,Agg) dH1o(:,Agg) dH2o(:,Agg)];
end

clearvars -except ...
U3 dU3 Umean dUmean ...
UPROFfit dUPROFfit ...
UPROFraw dUPROFraw ...
R3 dR3 Rmean dRmean ...
RPROFfit dRPROFfit ...
RPROFraw dRPROFraw ...
TPROFfit dTPROFfit ...
TPROFraw dTPROFraw ...
H3 dH3 Hmean dHmean ...
HPROF dHPROF ...
F dF Ft dFt ...     % Fluxes of each component. Note Ft is total velocity
Fr dFr ...
Ur dUr Fnr dFnr ...     % Fluxes of each component. Note Ft is total velocity
Ure dUre Fne dFne ...     % Fluxes of each component. Note Ft is total velocity
XH XR ...
Q1dr Q1or Q2dr Q2or ...     % estimate of flow ratio from RBC data
dQ1dr dQ1or dQ2dr dQ2or ...     % estimate errors in flow ratio from RBC data
Q1d Q1o Q2d Q2o ...       % flow ratios from total velocity
dQ1d dQ1o dQ2d dQ2o ...
Qt dQt...
Ht dHt...
USYMN dUSYMN USYMM...
percinvalid efitu...
PIV

% Umean_px=Umean; dUmean_px=dUmean;
% clearvars -except Umean_px dUmean_px
% load('G:\New\Mats\PIV_2204');save('G:\New\Mats\PIV_2204');

load('G:\New\Mats\efith')
save('G:\New\Mats\PIV_1105')

%% use 6408_1608_0808 - 16 pixel pre-shift

PIA(:,:)=percinvalid(:,1,:,1);
PIN(:,:)=percinvalid(1:35,1,:,2);
'percent invalid'
[mean(PIA(:)) mean(PIN(:)) ]

'velocity profile fitting errors RBC'
[mean(efitu(:,1,1)) mean(efitu(1:35,1,2))]

'velocity profile fitting errors SM'
[mean(efitu(:,2,1)) mean(efitu(1:35,2,2))]

'haematocrit profile fitting errors'
[mean(efith(:,1)) mean(efith(1:35,2)) ]

'average and std haematocrits'
[HA,~,dHa]=Wmean(Hmean(:,5,1),dHmean(:,5,1).^-2,1);
[HN,~,dHn]=Wmean(Hmean(1:35,5,2),dHmean(1:35,5,2).^-2,1);







% [Q1dr(:,:),dQ1dr(:,:)]=DivideError(Umean(:,1,1,:),Umean(:,1,5,:),dUmean(:,1,1,:),dUmean(:,1,5,:));
% [Q1or(:,:),dQ1or(:,:)]=DivideError(Umean(:,1,3,:),Umean(:,1,5,:),dUmean(:,1,3,:),dUmean(:,1,5,:));
% [Q2dr(:,:),dQ2dr(:,:)]=DivideError(Umean(:,1,2,:),Umean(:,1,3,:),dUmean(:,1,2,:),dUmean(:,1,3,:));
% [Q2or(:,:),dQ2or(:,:)]=DivideError(Umean(:,1,4,:),Umean(:,1,3,:),dUmean(:,1,4,:),dUmean(:,1,3,:));
% 
% [Q1ds(:,:),dQ1ds(:,:)]=DivideError(Umean(:,2,1,:),Umean(:,2,5,:),dUmean(:,2,1,:),dUmean(:,2,5,:));
% [Q1os(:,:),dQ1os(:,:)]=DivideError(Umean(:,2,3,:),Umean(:,2,5,:),dUmean(:,2,3,:),dUmean(:,2,5,:));
% [Q2ds(:,:),dQ2ds(:,:)]=DivideError(Umean(:,2,2,:),Umean(:,2,3,:),dUmean(:,2,2,:),dUmean(:,2,3,:));
% [Q2os(:,:),dQ2os(:,:)]=DivideError(Umean(:,2,4,:),Umean(:,2,3,:),dUmean(:,2,4,:),dUmean(:,2,3,:));

