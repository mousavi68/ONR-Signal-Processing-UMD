function [xDP,xSP] = IntTang(sig,st,nd)
% Function to calculate intersecting tangent 
dfSig = diff(sig(st:nd));
[m,loc] = max(dfSig);
loc = st+loc;
xDP = floor(1/m*(sig(st)-sig(loc))+loc);
xSP = floor(1/m*(sig(nd)-sig(loc))+loc);
end