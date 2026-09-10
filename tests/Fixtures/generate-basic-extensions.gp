\\ Independent PARI outputs only. The C# implementation does not use PARI at runtime.
default(parisizemax,512000000);
default(realprecision,70);
setrand(20260911);
logrow(a)={my(e=ellinit(a));if(#e==0,return());my(m=ellminimalmodel(e),ps=ellratpoints(m,12));for(i=1,min(6,#ps),my(p=ps[i],z=ellpointtoz(m,p),w=real(m.omega[1]),u=real(z));u-=floor(u/w)*w;print(Str("L,",m.a1,",",m.a2,",",m.a3,",",m.a4,",",m.a6,",",p[1],",",p[2],",",u,",",abs(imag(z)),",",w)))};
for(a=-6,6,for(b=-4,4,if(random(5)==0,logrow([0,0,0,a,b]))));
logrow([0,0,0,-25,0]);logrow([0,-17,0,72,0]);logrow([0,0,0,-1,0]);logrow([0,0,0,0,1]);logrow([0,0,1,-1,0]);logrow([0,0,1,-7,6]);
{my(e=ellinit([0,0,1,-1,0]));for(n=1,12,my(p=ellmul(e,[0,0],n),z=ellpointtoz(e,p));print(Str("L,0,0,1,-1,0,",p[1],",",p[2],",",real(z),",",abs(imag(z)),",",real(e.omega[1]))))};
isogrow(a,p)={my(e=ellinit(a),g=ellisogeny(e,p),q=ellinit(g[1]));print(Str("I,",a[1],",",a[2],",",a[3],",",a[4],",",a[5],",",p[1],",",p[2],",",ellorder(e,p),",",q.c4,",",q.c6))};
isogrow([0,0,0,0,1],[-1,0]);isogrow([0,0,0,0,1],[0,1]);isogrow([0,0,0,0,1],[2,3]);isogrow([0,-1,1,-10,-20],[5,5]);isogrow([0,-17,0,72,0],[6,6]);
quit;
