default(parisizemax, 512000000);
units=[-3,-1,0,1,3];
emit(a)={my(e=ellinit(a));if(#e==0,return());my(m=ellminimalmodel(e));print(Str(m.a1,",",m.a2,",",m.a3,",",m.a4,",",m.a6,",",ellrootno(m,2),",",ellrootno(m,3)));};
for(p=2,3,for(i=0,8,for(j=0,9,for(u=1,#units,for(v=1,#units,emit([0,0,0,p^i*units[u],p^j*units[v]]))))));
quit;
