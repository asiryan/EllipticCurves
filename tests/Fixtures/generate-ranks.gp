default(parisizemax, 512000000);
setrand(20260909);
for(i=1,300,a=random(41)-20;b=random(201)-100;if(b==0||a*a==4*b,next());e=ellinit([0,a,0,b,0]);r=ellrank(e);if(r[1]==r[2],print(Str(a,",",b,",",r[1]))));
quit;
