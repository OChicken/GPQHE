q=1<<61;

a=0;
b=0;
for (i=0, 127, a+=(i+2)*x^i)
for (i=0, 127, b+=(i+3)*x^i)
r = liftall(Mod(Mod(a*b, x^128+1),q));
for (i=0, 127,\
  if (polcoef(r,i)>=q/2,\
    r-=q*x^i))
print("Mod(Mod(a*b, x^128+1),1<<61) = \n", r)

print("\n");

a=0;
b=0;
for (i=0, 127, a+=(576460752303434497-i-1)*x^i)
for (i=0, 127, b+=(576460752303436801-i-1)*x^i)
r = liftall(Mod(Mod(a*b, x^128+1),q));
for (i=0, 127,\
  if (polcoef(r,i)>=q/2,\
    r-=q*x^i))
print("Mod(Mod(a*b, x^128+1),1<<61) = \n", r)
