clear;clc;%close all
load('G:\New\Mats\Analytical');
XH=(0.5:83.5)/84-0.5;
UT=spline(Xt,VtN,XH);
UTN=UT'/trapz(XH,UT);
UTR=spline(Xt,VtN,XH);

clear Xt VtN

load('G:\New\Mats\PIV_2704')

npbs=1:35;
Lig=1;Bra=5;
%% Bluntness
   
low=32;
for Agg=1:2
    for Lig=1:3
        for Cas=1:44
                %% Velocity
            if Lig<3
                Uprof=UPROFfit(:,Cas,Lig,Bra,Agg);
                dUprof=dUPROFfit(:,Cas,Lig,Bra,Agg);
            else
                Uprof=TPROFfit(:,Cas,Bra,Agg);
                dUprof=dTPROFfit(:,Cas,Bra,Agg);
            end
%             errorbar(XH,Uprof,dUprof)
            
            un=trapz(XH,Uprof);
            UN=Uprof/un;
            dUN=dUprof/un;

            UB(Cas,Lig,Bra,Agg)=trapz(XH,abs(UTN-UN));
            dUB(Cas,Lig,Bra,Agg)=trapz(XH,dUN);

            %% Haematocrit
            if Lig==1
                Hprof=HPROF(:,Cas,Bra,Agg);
                dHprof=dHPROF(:,Cas,Bra,Agg);

                hn=trapz(XH,Hprof);
                dhn=trapz(XH,dHprof);

                hmid=trapz(XH(29:56),Hprof(29:56));
                dhmid=trapz(XH(29:56),dHprof(29:56));                
                HB(Cas,Bra,Agg)=1.5*(hn-hmid)./hn;
                [~,dhb]=DivideError(hmid,hn,dhmid,dhn);
                dHB(Cas,Bra,Agg)=1.5*dhb;
            end
        end
    end
end
            

for Agg=1:2
    for Lig=1:2
        for Bra=1:5
            for Cas=1:44
                Hprof=HPROF(:,Cas,Bra,Agg);
                dHprof=dHPROF(:,Cas,Bra,Agg);

                HprofL=Hprof(1:41);
                XHL=XH(1:41);
                dHprofL=dHprof(1:41);

                hnl=trapz(XHL,HprofL);
                dhnl=trapz(XHL,dHprofL);
                hn=trapz(XH,Hprof);
                dhn=trapz(XH,dHprof);

                HS(Cas,Bra,Agg)=abs(hnl/hn-0.5);
                dHS(Cas,Bra,Agg)=sqrt((dhnl/hn)^2+(hnl*dhn/hn^2)^2);
            end
        end
    end
end

clearvars -except Umean dUmean UB dUB HB dHB HS dHS ...
    Q1dr Q1or Q2dr Q2or dQ1dr dQ1or dQ2dr dQ2or
% 
save('G:\New\Mats\Parameters_2704');
            
            
            
  