clear all
% close all
figure
drange=1:1:22000;
fGHz=3.5;



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
FreeSpace=[];


P1BP=20*log10(40*pi*dBP*fGHz/3)+C1*log10(dBP)-C2+C3*dBP;
for d=drange
    d3d=d;
    %% Freespace PL calculation
    FreeSpace(indx)= 20*log10(d) + 20*log10(fGHz)+32.45;
    
    
    %% LOS PL calculation
    P1(indx)=20*log10(40*pi*d3d*fGHz/3)+C1*log10(d3d)-C2+C3*d;    
    P2(indx)=P1BP+40*log10(d3d/dBP);
    if d<dBP
        LOS(indx)=P1(indx);
    else
        LOS(indx)=P2(indx);
    end

    %% NLOS PL calculation
    P3(indx)=C4+C5+C6*(log10(d3d)-3)+C7;
    

    NLOS(indx)=max(LOS(indx),P3(indx)-12); 
    indx=indx+1;
end
bpline=[dBP,-100;dBP,200];

figure;
olddrange=drange;
% drange=log10(drange/1000);

% semilogx(drange ,P1,'LineStyle','--','LineWidth',2)  % d Normalized to 1km
% hold all
% semilogx(drange,P2,'LineStyle','--','LineWidth',2)  % d Normalized to 1km
% semilogx(drange ,P3,'LineStyle','--','LineWidth',2)  % d Normalized to 1km

semilogx(drange ,LOS,'r','LineWidth',1)  % d Normalized to 1km
hold on;
semilogx(drange ,NLOS,'g','LineWidth',1)  % d Normalized to 1km
semilogx(drange ,FreeSpace,'k','LineWidth',2);
grid on;


ISD=6000;
R=ISD/sqrt(3);
celledge=[R,-100;R,200];
h=line(celledge(:,1),celledge(:,2));
set(h,'Color',[1,0,0],'LineStyle',':')
legend('LOS','NLOS','Free Space','Cellradius ISD=6km');
ylabel('PL [dB]')
xlabel('Distance d(m)')
title(sprintf('For freq=%f',fGHz))



deltad=[];
cnt=1;
rWindow=100:10:6000;
for r=rWindow
[v indx]=min(abs(drange-r));

pllos=LOS(indx)-3; % a 3dB gain
plnlos=NLOS(indx)-3;  % for a 3dB gain
plfspace=FreeSpace(indx)-3; % for a 3dB gain

[v dlosindx]=min(abs(LOS-pllos));
[v dnlosindx]=min(abs(NLOS-plnlos));
[v dFSindx]=min(abs(FreeSpace-plfspace));
deltad(cnt,:)=[drange(dlosindx)/r drange(dnlosindx)/r drange(dFSindx)/r];
cnt=cnt+1;
end



figure;
plot(rWindow,1-deltad);
legend 'LOS' 'NLOS' 'Freespace'
grid on;
ylim([0 1])
xlabel('Ref. Distance d(m)')
ylabel('Distance Gain (%) for 3dB')
title(sprintf('Distance Comparison freq=%f',fGHz))
