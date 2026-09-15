default(parisizemax,536870912);
E=ellinit([1,0,0,-22442227122985365806140019829619195351264650303963380,1494792882769884194322140639447353804389891718364698900054413322736138280709776]);
print("VERSION=",version());
print("AINVS=",E[1..5]);
print("DISCRIMINANT=",E.disc);
print("C4=",E.c4);
print("GCD_C4_DISC=",gcd(E.c4,E.disc));
gettime();
R=ellglobalred(E);
print("CONDUCTOR=",R[1]);
print("CONDUCTOR_FACTORIZATION=",R[4]);
print("GLOBAL_REDUCTION=",R);
print("MILLISECONDS=",gettime());
print("PROVEN_PRIMES=",vector(matsize(R[4])[1],i,isprime(R[4][i,1])));
quit;

