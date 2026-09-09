\\ Independent inputs and outputs only; no curve database is loaded.
default(parisizemax, 512000000);
setrand(20260910);
emit(a)={my(e=ellinit(a));if(#e==0,return());my(rk=ellrankinit(e),r=ellrank(rk,2),cov=ell2cover(rk));print(Str(a[1],",",a[2],",",a[3],",",a[4],",",a[5],",",r[1],",",r[2],",",#cov));};
emit([0,-1,1,-10,-20]);
emit([0,0,1,-1,0]);
emit([0,1,1,-2,0]);
emit([0,0,1,-7,6]);
emit([1,-1,0,-79,289]);
emit([0,-1,1,-929,-10595]);
emit([0,0,0,-113^2,0]);
emit([0,0,0,-34^2,0]);
for(i=1,120,emit([random(2),random(3)-1,random(2),random(31)-15,random(61)-30]));
\\ Reducible resolvents: full and partial rational 2-torsion, including J=0.
for(n=1,16,emit([0,0,0,-n^2,0]));
for(n=1,16,emit([0,random(11)-5,0,random(25)-12,0]));
quit;
