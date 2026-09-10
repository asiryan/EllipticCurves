\\ PARI/GP reference outputs only; no PARI implementation code is used by the library.
default(parisizemax,512000000);
default(realprecision,70);
setrand(20260910);
localrow(a)={my(e=ellinit(a));if(#e==0,return());my(m=ellminimalmodel(e),ps=factor(abs(m.disc))[,1]);for(i=1,#ps,my(p=ps[i],d=elllocalred(m,p));print(Str("L,",a[1],",",a[2],",",a[3],",",a[4],",",a[5],",",p,",",valuation(m.disc,p),",",d[1],",",d[2],",",d[4],",",ellrootno(m,p))));};
realrow(a)={my(e=ellinit(a));if(#e==0,return());my(m=ellminimalmodel(e),w=m.omega);print(Str("P,",a[1],",",a[2],",",a[3],",",a[4],",",a[5],",",real(w[1]),",",abs(imag(w[2])),",",real(w[1])*if(m.disc>0,2,1),",",abs(real(w[1])*imag(w[2]))));my(points=ellratpoints(m,15));for(i=1,min(2,#points),my(p=points[i]);print(Str("H,",m.a1,",",m.a2,",",m.a3,",",m.a4,",",m.a6,",",p[1],",",p[2],",",ellheight(m,p))));};
for(i=1,500,localrow(vector(5,j,random(61)-30)));
forprime(p=2,7,for(i=1,700,localrow([p*random(5),p*random(5),p^2*random(5),p^(2+random(3))*(random(21)-10),p^(3+random(4))*(random(21)-10)])));
for(i=1,100,realrow(vector(5,j,random(11)-5)));
for(a=-12,12,for(b=-8,8,if(random(12)==0,realrow([0,0,0,a,b]))));
realrow([0,0,0,-1,0]);realrow([0,0,1,0,0]);realrow([0,0,1,-1,0]);realrow([0,0,1,-7,6]);
quit;
