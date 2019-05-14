clear all
close all

drange=1:10:22000;
fGHz=0.850;
 
%%% Constants
rmaH=5; % Avg Building heights
hBS=35;
hUT=1.5; % UT height
C=3e8;
dBP=2*pi*hBS*hUT*fGHz*1e9/C;


C1 = min(0.03*(rmaH^ 1.72), 10);
C2 = min(0.044*(rmaH^ 1.72), 14.77);
C3 = 0.002 * log10(rmaH);
W=20;
%%% nlos constansts
C4=161.04-7.1*log10(W)+7.5*log10(rmaH);
C5=-(24.37-3.7*(rmaH/hBS)^2)*log10(hBS);
C6=(43.42-3.1*log10(hBS));
C7=20*log10(fGHz)-(3.2*(log10(11.75*hUT)^2)-4.97);

%%% Evalate P1
indx=1;
P1=[];
P2=[];
LOS=[];
NLOS=[];
NLOSeH=[];
NLOSeS=[];
NLOSeHS=[];
FS=[];
P1BP=20*log10(40*pi*dBP*fGHz/3)+C1*log10(dBP)-C2+C3*dBP;
for d=drange
    FS(indx)= 20*log10(d) + 20*log10(fGHz)+32.45;
    P1(indx)=20*log10(40*pi*d*fGHz/3)+C1*log10(d)-C2+C3*d;
    
    P2(indx)=P1BP+40*log10(d/dBP);
    
    P3(indx)=C4+C5+C6*(log10(d3D)-3)+C7;
    
    if d2D<dBP
        LOS(indx)=P1(indx);
    else
        LOS(indx)=P2(indx);
    end
    
    NLOS(indx)=max(LOS(indx),P3(indx));
    NLOSeH(indx)=max(LOS(indx),P3(indx)-12);
    NLOSeS(indx)=max(LOS(indx)+6,P3(indx)-12);
   if d2D<dBP
       NLOSeHS(indx)=NLOS(indx);
   else
       NLOSeHS(indx)=NLOSeS(indx);
    end
 
    indx=indx+1;
end
bpline=[dBP,-100;dBP,200];

figure;
olddrange=drange;
% drange=log10(drange/1000);

semilogx(drange ,P1,'LineStyle','--','LineWidth',2)  % d Normalized to 1km
hold all
semilogx(drange,P2,'LineStyle','--','LineWidth',2)  % d Normalized to 1km
semilogx(drange ,P3,'LineStyle','--','LineWidth',2)  % d Normalized to 1km

semilogx(drange ,LOS,  'LineWidth',1)  % d Normalized to 1km
semilogx(drange ,NLOS, 'LineWidth',2)  % d Normalized to 1km
semilogx(drange ,NLOSeH,'LineWidth',2)  % d Normalized to 1km
semilogx(drange ,NLOSeS,'LineWidth',2)  % d Normalized to 1km
% semilogx(drange ,NLOSeHS,'LineWidth',2)  % d Normalized to 1km
semilogx(drange ,FS,'k','LineWidth',2);
grid on;


legend('P1','P2','P3','LOS','NLOS','NLOSeH' ,'NLOSeS' , 'FSpace');
h=line(bpline(:,1),bpline(:,2));
set(h,'Color',[1,0,0])
set(h,'LineStyle','-.');
ylabel('PL [dB]')
xlabel('Distance log10(d)(m)')



figure;

bpline=[dBP,-100;dBP,200];
semilogx(drange/1000 ,LOS)  % d Normalized to 1km
hold all
grid on
semilogx(drange/1000 ,NLOS)  % d Normalized to 1km
% h=line(bpline(:,1),bpline(:,2));
% set(h,'Color',[1,0,0])
% set(h,'LineStyle','-.');
ylabel('PL [dB]')
xlabel('Distance log10(d)(m)')


