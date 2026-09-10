\\ Independent PARI reference outputs; no database or C# implementation is used.
default(realprecision,80);
default(parisizemax,512000000);
setrand(20260912);
heightrow(a)={my(e=ellinit(a));if(#e==0,return());e=ellminimalmodel(e);my(h=-log(e.area)/2,s=h+log(denominator(e.j)/abs(e.disc))/12);print(Str("H,",a[1],",",a[2],",",a[3],",",a[4],",",a[5],",",h,",",s))};
for(a=-8,8,for(b=-5,5,if(random(4)==0,heightrow([0,0,0,a,b]))));
heightrow([0,0,1,-1,0]);heightrow([0,-17,0,72,0]);heightrow([0,0,0,-1,0]);heightrow([0,0,0,0,1]);heightrow([0,0,1,-7,6]);
decode(n,g,p,k)=sum(i=0,k-1,(n\p^i%p)*g^i);
encode(v,g,p,k)={v=v*g^0;sum(i=0,k-1,lift(polcoef(v.pol,i))*p^i)};
coord(P,g,p,k)={if(#P==1,"-1,-1",Str(encode(P[1],g,p,k),",",encode(P[2],g,p,k)))};
fieldrows(p,m)={my(k=#m-1,f=Mod(1,p)*Polrev(m),g=ffgen(f,'t),q=p^k,tag=Str(p,",",strjoin(apply(x->Str(x),m),":")));if(!polisirreducible(f),error("Reducible fixture modulus"));for(i=1,12,my(a=random(q),b=random(q),u=decode(a,g,p,k),v=decode(b,g,p,k));print(Str("F,",tag,",",a,",",b,",",encode(u+v,g,p,k),",",encode(u*v,g,p,k),",",if(a,encode(1/u,g,p,k),-1))));for(s=1,5,my(c=vector(5,i,random(q)),e=ellinit(vector(5,i,decode(c[i],g,p,k)),g));if(#e==0,next());my(ps=List());for(i=0,q-1,for(j=0,q-1,my(P=[decode(i,g,p,k),decode(j,g,p,k)]);if(ellisoncurve(e,P),listput(ps,P);if(#ps==3,break(2)))));if(#ps<2,next());my(P=ps[1],Q=ps[2],n=s-3);print(Str("E,",tag,",",strjoin(apply(x->Str(x),c),","),",",ellcard(e),",",encode(e.disc,g,p,k),",",encode(e.c4,g,p,k),",",encode(e.c6,g,p,k),",",encode(e.j,g,p,k),",",coord(P,g,p,k),",",coord(Q,g,p,k),",",coord(elladd(e,P,Q),g,p,k),",",n,",",coord(ellmul(e,P,n),g,p,k))))};
fieldrows(2,[1,1,1]);fieldrows(2,[1,1,0,1]);fieldrows(2,[1,1,0,0,1]);fieldrows(3,[1,0,1]);fieldrows(3,[1,2,0,1]);fieldrows(5,[2,0,1]);fieldrows(7,[1,0,1]);fieldrows(5,[1,1,0,1]);
quit;
