setrand(20260916);
print("p,q");
for(i=1,8,d=17+2*i;p=nextprime(10^(d-1)+random(9*10^(d-1)));q=nextprime(10^(d-1)+random(9*10^(d-1)));if(!isprime(p)||!isprime(q),error("Unproved fixture prime"));print(p,",",q));
quit
