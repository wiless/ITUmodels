function d=BPdistance(fGHz) 
hBS=35;
hUT=1.5;
C=3e8;
d=2*pi*hBS*hUT*fGHz*1e9/C;
