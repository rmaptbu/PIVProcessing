clear;clc;%close all
D=84;  
Agg=1;
% load('G:\New\Mats\ImageAdjust')
% image regions
DaughterX=4.5*D+1:5.5*D;
DaughterY=5*D+1:9.5*D;
Daughter2X=10.5*D+1:11.5*D;
Daughter2Y=5*D+1:9.5*D;
OutletX=5.5*D+1:10.5*D; % added extra half diameter compared to PIV
OutletY=9.5*D+1:10.5*D;
Outlet2X=11.5*D+1:16*D;
Outlet2Y=9.5*D+1:10.5*D;
ParentX=1:4.5*D;
ParentY=9.5*D+1:10.5*D;
YY={DaughterY Daughter2Y OutletY Outlet2Y ParentY};
XX={DaughterX Daughter2X OutletX Outlet2X ParentX};
savepaths={'Daughter';'Daughter2';'Outlet';'Outlet2';'Parent'};

% settings
D=84;
Trim=6;                                         % number of dodgy pixels to remove
X=(0.5:83.5)/84;XT=[X(Trim+1:end-Trim)];        % build new X array     
XH=(0.5:83.5)/84-0.5;
alpha=0.685;                                   % scaling for haematocrit
beta=9.244;   
dH=0.014;
L=6;                                            % sweep length - gives 50/84*13=7.74 micron
for Agg=1:2
    if Agg==1;
        basepath='G:\New\Images\A\PreProcess_Double\MeanCor\';
        imsavepath='G:\New\HctImages\A\';
        profsavepath='G:\New\HctProfiles\A\';
    else
        basepath='G:\New\Images\N\PreProcess_Double\MeanCor\';
        imsavepath='G:\New\HctImages\N\';
        profsavepath='G:\New\HctProfiles\N\';
    end
    cd(basepath)
    Cases=struct2cell(dir);Cases=Cases(1,3:end); % these are the corrected images
    for Cas=1:length(Cases)
        filename=char(Cases(Cas));
        Im=single(imread(filename));  %% import 
        for Bra=1:5
            ImBra=Im(cell2mat(YY(Bra)),cell2mat(XX(Bra)));          % select branch
            ImBra=single(ImBra);                                    % convert to single format
            MaxPx=4095;         
            ImStar=coerce2((MaxPx-ImBra)/MaxPx,0,alpha-0.0001);     % find I star            
         
            if Bra<3
                ImStar=imrotate(ImStar,-90);
            end
            ImSweep=zeros(size(ImStar));
            dImSweep=zeros(size(ImStar));
            for i=1:size(ImStar,2)
                if i<(L+1);a=1;else a=i-L;end                       % for first bits, just use forward points
                if i>size(ImStar,2)-L;b=size(ImStar,2);else b=i+L;end % for last bits, just use backward points
                ImSec=double(ImStar(:,a:b));
                IProfC=mean(ImSec(Trim+1:end-Trim,:),2);      % inner bit of profile
                HProfL=polyval(polyfit(XT(1:Trim)',IProfC(1:Trim),1),X(1:Trim)); % extrapolate LHS
                HProfR=polyval(polyfit(XT(end-Trim:end)',IProfC(end-Trim:end),1),X(end-Trim+1:end));
                IProf=[HProfL'; IProfC; HProfR'];
                dIProfC=std(ImSec(Trim+1:end-Trim,:),[],2);      % inner bit of profile - std
                dHProfL=polyval(polyfit(XT(1:Trim)',dIProfC(1:Trim),1),X(1:Trim)); % extrapolate LHS
                dHProfR=polyval(polyfit(XT(end-Trim:end)',dIProfC(end-Trim:end),1),X(end-Trim+1:end));
                dIProf=[dHProfL'; dIProfC; dHProfR'];
                dIProfS=dIProf/sqrt(size(ImSec,2));
                ImSweep(:,i)=coerce1(IProf,0,alpha-0.001);
                dImSweep(:,i)=coerce1(dIProfS,0,alpha-0.001);
            end
            ImHct=ApplyCalib([alpha,beta],ImSweep); 
            dImHct=ApplyCalib([alpha,beta],ImSweep+dImSweep)-ImHct; 
            dImHctC=sqrt(dImHct.^2+dH.^2); % SEM intensity over averaging is approx 2.5% of error due to calibration
            
            ImWrite=uint16(ImHct*10000); % each pixel value is worth 0.001 hct values
            if isdir([imsavepath,char(savepaths(Bra))]);else  mkdir([imsavepath,char(savepaths(Bra))]);end
            imwrite(ImWrite,[imsavepath,char(savepaths(Bra)),'\',filename(1:6),'_Hct.tif'])
            
            dImWrite=uint16(dImHct*10000); % each pixel value is worth 0.001 hct values
            if isdir([imsavepath,char(savepaths(Bra))]);else  mkdir([dimsavepath,char(savepaths(Bra))]);end
            imwrite(dImWrite,[imsavepath,char(savepaths(Bra)),'\',filename(1:6),'_dHct.tif'])
                                
            % select profiles
            if Bra==5
                ImPar=double(ImHct(:,1:D));
                ImSym=[ImPar; flipud(ImPar)];
                dImPar=double(dImHct(:,1:D));
                dImSym=[dImPar; flipud(dImPar)];
                
                [HProfC,~,dHProfC]=Wmean2(ImSym(:,Trim+1:end-Trim),dImSym(:,Trim+1:end-Trim).^-2,1,2);
                
                HProfL=polyval(polyfit(XT(1:Trim),HProfC(1:Trim),1),X(1:Trim)); % extrapolate LHS linearly
                HProfR=polyval(polyfit(XT(end-Trim:end),HProfC(end-Trim:end),1),X(end-Trim+1:end)); % extrapolate RHS linearly
                HProfSym=[HProfL HProfC HProfR]';    % compile profi
                
                dHProfL=polyval(polyfit(XT(1:Trim)',dIProfC(1:Trim),1),X(1:Trim)); % extrapolate LHS
                dHProfR=polyval(polyfit(XT(end-Trim:end)',dIProfC(end-Trim:end),1),X(end-Trim+1:end));
                dHProfSym=[dHProfL dHProfC dHProfR]';
                
                [SymProf,dSymProf,~,~,~,~,rmse]=FitCoshU2(XH',HProfSym,dHProfSym.^-2);

                Have=trapz(XH',HProfSym);
                dHave=trapz(XH',dHProfSym);
                
                efith(Cas,Agg)=RootMS(HProfSym-SymProf)/Have*100;

                
%                 [SymProfN,n]=normalise(SymProf,1);
%                 dSymProfN=dSymProf/n;

                dlmwrite([profsavepath,filename(1:6),'_HctSym.txt'],SymProf)
                dlmwrite([profsavepath,filename(1:6),'_dHctSym.txt'],dSymProf)
            end
        end % Branches
    end % Cases
end % agg/pbs
clearvars -except efith
save('G:\New\Mats\efith')
        