a=0;
b=0;
for (i=0, 127, a+=(i+2)*x^i)
print("a = ",a)
for (i=0, 127, b+=(i+3)*x^i)
print("b = ",b)
print("Mod(Mod(a*b, x^128+1),1048576) = \n", liftall(Mod(Mod(a*b, x^128+1),1048576)))

print("\n\n");

a=0;
b=0;
for (i=0, 127, a+=(576460752303434497-i-1)*x^i)
print("a = ",a)
for (i=0, 127, b+=(576460752303436801-i-1)*x^i)
print("b = ",b)
print("Mod(Mod(a*b, x^128+1),1048576) = \n", liftall(Mod(Mod(a*b, x^128+1),1048576)))


