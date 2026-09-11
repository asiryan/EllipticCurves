\\ Independent full torsion sets, including every Mazur group and changed models.
default(parisizemax, 512000000);
setrand(20260911);
emit(a)={
  my(E=ellinit(a)); if(#E==0,return());
  my(T=elltors(E), P=vector(T[1],k,my(Q=[0],n=k-1);for(j=1,#T[2],Q=elladd(E,Q,ellmul(E,T[3][j],n%T[2][j]));n=n\T[2][j]);Q));
  my(group=if(#T[2],strjoin(vector(#T[2],i,Str("Z/",T[2][#T[2]+1-i],"Z"))," x "),"Z/1Z"));
  print(strjoin(vector(5,i,Str(a[i])),","),",",group,",",strjoin(vector(#P,i,if(#P[i]==1,"O",Str(P[i][1],":",P[i][2]))),";"));
};
curves=[ [11,30,30,0,0], [4,-2,-2,0,0], [6,6,6,0,0], [1,30,30,0,0], [11,10,10,0,0], [3,-2,-2,0,0], [-5,12,12,0,0], [1,4,4,0,0], [3,6,6,0,0], [-11,-12,-12,0,0], [1,8/3,8/3,0,0], [0,3,0,2,0], [1,-3,-3,0,0], [0,37,0,160,0], [0,337,0,20736,0] ];
for(i=1,#curves,E=ellinit(curves[i]);emit(curves[i]);emit(ellchangecurve(E,[2,0,0,0])[1..5]);emit(ellchangecurve(E,[-2/3,5,-2,7])[1..5]));
emit([0,-424/25,0,72,0]);
emit([0,-339/20,0,72,0]);
emit([0,-1697/100,0,72,0]);
emit([0,-17,0,72,0]);
for(i=1,32,emit(vector(5,j,random(11)-5)));
quit;
