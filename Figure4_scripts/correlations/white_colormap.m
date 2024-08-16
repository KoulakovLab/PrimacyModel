function m = white_colormap

m=zeros(64, 3);

dr = [8:16, 16*ones(1,24), 16:(-1):0 ]';
dr=dr/max(dr);
r = m(:,1);
r(1:length(dr))=dr;

b = r(end:(-1):1);

g = r*0;

dg = [zeros(1,15),  [1:16], 16*ones(1,3), [15:(-1):0] ]';
dg = dg/max(dg);

g(1:length(dg))=dg;

m = [b,g,r];


return
end