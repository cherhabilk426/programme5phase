function tensionref=ondhyst(ias,ibs,ics,ifs,ies,iasref,ibsref,icsref,ifsref,iesref,int1,int2,int3,int4,int5)
pas=0.1;

e1=ias-iasref;
e2=ibs-ibsref;
e3=ics-icsref;
e4=ifs-ifsref;
e5=ies-iesref;


int1old=int1;

if e1>pas,
   int1=0;
elseif e1<-pas,
   int1=1;
else
   int1=int1old;
end;

int2old=int2;

if e2>pas,
   int2=0;
elseif e2<-pas,
   int2=1;
else
   int2=int2old;
end;

int3old=int3;

if e3>pas,
   int3=0;
elseif e3<-pas,
   int3=1;
else
   int3=int3old;
end;

int4old=int4;

if e4>pas,
   int4=0;
elseif e4<-pas,
   int4=1;
else
   int4=int4old;
end;

int5old=int5;

if e5>pas,
   int5=0;
elseif e5<-pas,
   int5=1;
else
   int5=int5old;
end;

vs=514;   
va=vs*( 2*int1-int2-int3-int4-int5)/5;
vb=vs*(-int1+2*int2-int3-int4-int5)/5;
vc=vs*(-int1-int2+2*int3-int4-int5)/5;
vf=vs*(-int1-int2-int3+2*int4-int5)/5;
ve=vs*(-int1-int2-int3-int4+2*int5)/5;

tensionref=[va;vb;vc;vf;ve];

