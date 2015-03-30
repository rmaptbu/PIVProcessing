clear;clc;close all

load('G:\New\Mats\PIV_2704')
ca=1:44;
cn=[1:4,6:33,35];
%% continuity of spatially resolved fluxes - terms in haematocrit weighted volumetric flow
lig={'RBC','SM','Tot'};
col=[0.8 0 0; 0 0 0.8; 0 0.8 0];
for Agg=1:2
    if Agg==1;c=ca;else c=cn;end
    for Lig=1:3
        clear f df
        if Lig<3
            f(:,:)=F(c,Lig,:,Agg);
            df(:,:)=dF(c,Lig,:,Agg);
        else
            f(:,:)=Ft(c,:,Agg);
            df(:,:)=dFt(c,:,Agg);
        end
        c1=(f(:,5)-f(:,1)-f(:,3));
        c2=(f(:,3)-f(:,4)-f(:,2));
        ct=(f(:,5)-f(:,1)-f(:,2)-f(:,4));
        dc1=emag([df(:,5) df(:,1) df(:,3)]);
        dc2=emag([df(:,3) df(:,4) df(:,2)]);
        dct=emag([df(:,1) df(:,2) df(:,4) df(:,5)]);
        
        [C1(c,Lig,Agg),dC1(c,Lig,Agg)]=DivideError(c1,f(:,5),dc1,df(:,5));
        [C2(c,Lig,Agg),dC2(c,Lig,Agg)]=DivideError(c2,f(:,3),dc2,df(:,3));
        [CT(c,Lig,Agg),dCT(c,Lig,Agg)]=DivideError(ct,f(:,5),dct,df(:,5));
        
        errorbar(C1,dC1)
        
       [mC1(Lig,Agg),~,sC1(Lig,Agg)]=Wmean(C1(c,Lig,Agg),dC1(c,Lig,Agg).^-2,1);
       [mC2(Lig,Agg),~,sC2(Lig,Agg)]=Wmean(C2(c,Lig,Agg),dC2(c,Lig,Agg).^-2,1);  
       [mCT(Lig,Agg),~,sCT(Lig,Agg)]=Wmean(CT(c,Lig,Agg),dCT(c,Lig,Agg).^-2,1);  

    figure(Agg)
    if Lig <3; Up=Umean(c,Lig,5,Agg);else Up=Ft(c,5,Agg); end
    subplot(3,3,(Lig-1)*3+1);h=herrorbar(C1(c,Lig,Agg),Up,dC1(c,Lig,Agg),'x');set(h,'color',col(Lig,:));hold all;title([char(lig(Lig)),' Bif 1'])
    plot(mC1(Lig,Agg)*[1 1],[0 10],'k');axis([-0.15 0.15 0 10])
    plot((mC1(Lig,Agg)+1.96*sC1(Lig,Agg))*[1 1],[0 10],':k');
    plot((mC1(Lig,Agg)-1.96*sC1(Lig,Agg))*[1 1],[0 10],':k'); 
    
    subplot(3,3,(Lig-1)*3+2);h=herrorbar(C2(c,Lig,Agg),Up,dC2(c,Lig,Agg),'+');set(h,'color',col(Lig,:));hold all;title([char(lig(Lig)),' Bif 2'])
    plot(mC2(Lig,Agg)*[1 1],[0 10],'k');axis([-0.15 0.15 0 10])
    plot((mC2(Lig,Agg)+1.96*sC2(Lig,Agg))*[1 1],[0 10],':k');
    plot((mC2(Lig,Agg)-1.96*sC2(Lig,Agg))*[1 1],[0 10],':k'); 
    
    subplot(3,3,(Lig-1)*3+3);h=herrorbar(CT(c,Lig,Agg),Up,dCT(c,Lig,Agg),'.');set(h,'color',col(Lig,:));hold all
    plot(mCT(Lig,Agg)*[1 1],[0 10],'k');axis([-0.15 0.15 0 10]);title([char(lig(Lig)),' Total Bif'])
    plot((mCT(Lig,Agg)+1.96*sCT(Lig,Agg))*[1 1],[0 10],':k');
    plot((mCT(Lig,Agg)-1.96*sCT(Lig,Agg))*[1 1],[0 10],':k'); 
    
    end
end

Lig=1;Agg=1;
[mC1(Lig,Agg) sC1(Lig,Agg) mC2(Lig,Agg) sC2(Lig,Agg) mCT(Lig,Agg) sCT(Lig,Agg)]*100
Lig=1;Agg=2;
[mC1(Lig,Agg) sC1(Lig,Agg) mC2(Lig,Agg) sC2(Lig,Agg) mCT(Lig,Agg) sCT(Lig,Agg)]*100



clearvars -except mC1 mC2 mCT sC1 sC2 sCT ...
                    C1 C2 CT dC1 dC2 dCT
save('G:\New\Mats\MassCont_2004')