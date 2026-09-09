default(parisizemax, 512000000);
setrand(20260909);
emit(a)={my(e=ellinit(a)); if(#e==0,return()); my(m=ellminimalmodel(e),n=ellglobalred(e)[1],d2=elllocalred(m,2),d3=elllocalred(m,3)); print(Str(a[1],",",a[2],",",a[3],",",a[4],",",a[5],",",n,",",m.a1,",",m.a2,",",m.a3,",",m.a4,",",m.a6,",",d2[2],",",d3[2]));};
for(i=1,2000,emit(vector(5,j,random(81)-40)));
for(p=2,3,for(i=1,2000,emit([p*random(5),p*random(5),p^2*random(5),p^(2+random(3))*(random(21)-10),p^(3+random(4))*(random(21)-10)])));
quit;
