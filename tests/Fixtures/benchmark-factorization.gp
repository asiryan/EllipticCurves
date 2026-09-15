\\ Run from the repository root. Use fresh ell structures to avoid result caching.
default(parisizemax,512000000);
default(nbthreads,1);
default(factor_proven,1);
print("case,digits,first_ms,median_ms");
{
  my(times=vector(4),e,result,started);
  for(j=1,4,
    setrand(20260916+j);
    started=getwalltime();
    e=ellinit([0,1,0,-221556180740323405132844117936,35386140191724122461245294467670188433973860]);
    result=ellglobalred(e)[1];
    times[j]=getwalltime()-started;
    if(result!=1103561624055499058867562340698878392772504928025988266715523317532246643920,error("Conductor mismatch")));
  print("conductor,90,",times[1],",",vecsort(times[2..4])[2]);
}
{
  my(rows=readstr("tests/Fixtures/factorization.csv"),fields,n,times,started,result);
  for(i=2,#rows,
    fields=strsplit(rows[i],",");
    n=eval(fields[1])*eval(fields[2]);
    times=vector(4);
    for(j=1,4,
      setrand(20260916+j);
      started=getwalltime();
      result=factor(n);
      times[j]=getwalltime()-started;
      if(factorback(result)!=n,error("Factorization mismatch")));
    print(i-1,",",#Str(n),",",times[1],",",vecsort(times[2..4])[2]));
}
quit
