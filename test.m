rs = rng;
p8.smat = eye(8);
net0 = initnet3(8,2,8,2,2,rs);
p8.tmat = eye(8);
net_a = bp3(net0,p8,10000,1,0,rs);
act_a = forw3(net_a,p8)
i8 = eye(8) ;
p8.tmat = i8 + i8([2 3 4 5 6 7 8 1],:) + i8([8 1 2 3 4 5 6 7],:);
p8.tmat
net_b = bp3(net0,p8,10000,1,0,rs);
act_b = forw3(net_b,p8)
act_a.out
act_b.out
function netstruct=initnet3(n1,n2,n3,uamp,vamp,rs)
rng(rs) ;
netstruct.wih=uamp*(rand(n2,n1)-0.5) ;
netstruct.hbias=uamp*(rand(1,n2)-0.5) ;
netstruct.whout=vamp*(rand(n3,n2)-0.5) ;
netstruct.obias=vamp*(rand(1,n3)-0.5);
end

function netact=forw3(netwk,pats)
netact.stim=pats.smat;
netact.hid=layersigpn(netact.stim,netwk.wih,netwk.hbias) ;
netact.out=layersig01(netact.hid,netwk.whout,netwk.obias) ;
end

function finalnet=bp3(net0,patstr,niter,eta,nlev,rs)
rng(rs) ;
netk=net0;
for i=1:niter
    netk=cyc3(netk,patstr,eta,nlev) ;
end
finalnet=netk;
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

function netact=forw1p3(netwk,pats,patno,nois)
netact.stim=pats.smat(patno,:);
netact.hid=nlayersigpn(netact.stim,netwk.wih,netwk.hbias,nois) ;
netact.out=layersig01(netact.hid,netwk.whout,netwk.obias) ;
end

function lout=nlayersigpn(x,w,b,noise)
lsz=size(x) ;
% lll=b'*ones(lsz(1),1)'+w*x'
bb=repmat(b',1,lsz(1)) ;
lll=bb+w*x' ;
lout=sigpn(lll)'+noise*(rand(size(lll'))-0.5) ;
end

function squash=sigpn(x)
squash=(1-exp(-x))./(1+exp(-x)) ;
end

function lout=layersig01(x,w,b)
lsz=size(x) ;
lll=b'*ones(lsz(1),1)'+w*x' ;
lout=sig01(lll)';
end

function squash=sig01(x)
squash=1./(1+exp(-x)) ;
end

function lout=layersigpn(x,w,b)
lsz=size(x) ;
lll=b'*ones(lsz(1),1)'+w*x' ;
lout=sigpn(lll)';
end