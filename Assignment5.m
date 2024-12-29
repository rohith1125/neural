xor.smat = [0 0;0 1;1 0;1 1];
xor.tmat = [0 1 1 1]';
xornet0 = initnet3(2,2,1,2,2,rs);
xoract0 = forw3(xornet0,xor);
xornet1k = bp3(xornet0,xor,1000,1,0,rs);
xornet1k.wih
xoract1k = forw3(xornet1k,xor);
kolor = [xor.tmat,zeros(size(xor.tmat)),1-xor.tmat];

figure
hold on
scatter(xoract0.hid(:,1),xoract0.hid(:,2),10,"filled")
scatter(xoract1k.hid(:,1),xoract1k.hid(:,2),10,kolor,"filled")
drawArrow = @(x1,y1,x2,y2) quiver(x1,y1,x2-x1,y2-y1,0);
drawArrow(xoract0.hid(:,1),xoract0.hid(:,2),xoract1k.hid(:,1),xoract1k.hid(:,2))
hold off

figure
x = 0:0.01:1;
y = 0:0.01:1;
stims=combvec(x,y)';
netwk=xornet1k;
hid=layersigpn(stims,netwk.wih,netwk.hbias) ;
out=layersig01(hid,netwk.whout,netwk.obias) ;
[X,Y]= meshgrid(x,y);
contour(X,Y,griddata(stims(:,1),stims(:,2),out,X,Y))

function netact=forw3(netwk,pats)
netact.stim=pats.smat;
netact.hid=layersigpn(netact.stim,netwk.wih,netwk.hbias) ;
netact.out=layersig01(netact.hid,netwk.whout,netwk.obias) ;
end

function squash=sigpn(x)
squash=(1-exp(-x))./(1+exp(-x)) ;
end

function finalnet=bp3(net0,patstr,niter,eta,nlev,rs)
rng(rs) ;
netk=net0;
for i=1:niter
    netk=cyc3(netk,patstr,eta,nlev) ;
end
finalnet=netk;
end

function lout=layersigpn(x,w,b)
lsz=size(x) ;
lll=b'*ones(lsz(1),1)'+w*x' ;
lout=sigpn(lll)';
end

function lout=nlayersigpn(x,w,b,noise)
lsz=size(x) ;
% lll=b'*ones(lsz(1),1)'+w*x'
bb=repmat(b',1,lsz(1)) ;
lll=bb+w*x' ;
lout=sigpn(lll)'+noise*(rand(size(lll'))-0.5) ;
end

function newstruct=cyc3(nstruct,pstruct,dt,noi)
newstruct=nstruct;
szs=size(pstruct.smat);
patk=ceil(szs(1)*rand());
activity=forw1p3(nstruct,pstruct,patk,noi);
odelt=(pstruct.tmat(patk,:)-activity.out); %output deltas
hdelt=0.5*(nstruct.whout'*odelt').*(1+activity.hid').*(1-activity.hid');%hid deltas
%adjust weights and biases
newstruct.whout=newstruct.whout+dt*odelt'*activity.hid;
newstruct.obias=newstruct.obias+dt*odelt ;
newstruct.wih=newstruct.wih+dt*hdelt*activity.stim;
newstruct.hbias=newstruct.hbias+dt*hdelt' ;
end

function netstruct=initnet3(n1,n2,n3,uamp,vamp,rs)
rng(rs) ;
netstruct.wih=uamp*(rand(n2,n1)-0.5) ;
netstruct.hbias=uamp*(rand(1,n2)-0.5) ;
netstruct.whout=vamp*(rand(n3,n2)-0.5) ;
netstruct.obias=vamp*(rand(1,n3)-0.5);
end

function squash=sig01(x)
squash=1./(1+exp(-x)) ;
end

function lout=layersig01(x,w,b)
lsz=size(x) ;
lll=b'*ones(lsz(1),1)'+w*x' ;
lout=sig01(lll)';
end

function lout=nlayersig01(x,w,b,noise)
lsz=size(x) ;
lll=b'*ones(lsz(1),1)'+w*x'
% lll=b'+w*x' ;
%lout=sig01(lll)'+noise*(rand(lsz(1),1)-0.5) ;
lout=sig01(lll)';
end

function netact=forw1p3(netwk,pats,patno,nois)
netact.stim=pats.smat(patno,:);
netact.hid=nlayersigpn(netact.stim,netwk.wih,netwk.hbias,nois) ;
netact.out=layersig01(netact.hid,netwk.whout,netwk.obias) ;
end