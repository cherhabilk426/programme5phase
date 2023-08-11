function [idsref,didsref]=Floue3_phr(idsref,ephr,vephr);

a1=0.5;

% fonction mu (NG,NM,NP,ZE,PP,PM,PG) pour l'erreur

mueNG=max(min(1,-ephr/a1),0);
mueZE=max(min(1+ephr/a1,1-ephr/a1),0);
muePG=max(min(ephr/a1,1),0);

%fonction mu (NG,NM,NP,ZE,PP,PM,PG) pour la dérivée
muveNG=max(min(1,-vephr/a1),0);
muveZE=max(min(1+vephr/a1,1-vephr/a1),0);
muvePG=max(min(vephr/a1,1),0);


% defuzzification : max-min

musNG=max([min(mueNG,muveNG);min(mueNG,muveZE);min(mueZE,muveNG)]);
musZE=max([min(mueNG,muvePG);min(mueZE,muveZE);min(muePG,muveNG)]);
musPG=max([min(mueZE,muvePG);min(muePG,muvePG);min(muePG,muveZE)]);

som_mus=musNG+musZE+musPG;

Id=1;

s_PG=0.5*Id;
s_ZE=0*Id;
s_NG=-0.5*Id;

didsref=(s_PG*musPG+s_ZE*musZE+s_NG*musNG)/som_mus;

Gud=4.5;
idsref=idsref+Gud*didsref;
