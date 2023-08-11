function [iqsref,diqsref]=Floue3_w(iqsref,ew,vew);

a1=0.5;

% fonction mu (NG,NM,NP,ZE,PP,PM,PG) pour l'erreur
muewNG=max(min(1,-ew/a1),0);
muewZE=max(min(1+ew/a1,1-ew/a1),0);
muewPG=max(min(ew/a1,1),0);

%fonction mu (NG,NM,NP,ZE,PP,PM,PG) pour la dérivée
muvewNG=max(min(1,-vew/a1),0);
muvewZE=max(min(1+vew/a1,1-vew/a1),0);
muvewPG=max(min(vew/a1,1),0);


% defuzzification : max-min

muqNG=max([min(muewNG,muvewNG);min(muewNG,muvewZE);min(muewZE,muvewNG)]);
muqZE=max([min(muewNG,muvewPG);min(muewZE,muvewZE);min(muewPG,muvewNG)]);
muqPG=max([min(muewZE,muvewPG);min(muewPG,muvewPG);min(muewPG,muvewZE)]);

som_muq=muqNG+muqZE+muqPG;


Iq=1;

s_PG=0.5*Iq;
s_ZE=0*Iq;
s_NG=-0.5*Iq;

diqsref=(s_PG*muqPG++s_ZE*muqZE+s_NG*muqNG)/som_muq;

Guq=6.5;
iqsref=iqsref+Guq*diqsref;
