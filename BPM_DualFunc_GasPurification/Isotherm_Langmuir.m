%% Langmuir isotherm
function g =Isotherm_Langmuir(c,b,qm)
c= max(c,0);
g= qm*(b*c)./(1+b*c);
end