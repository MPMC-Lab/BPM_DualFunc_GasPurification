%% Derivative of Langmuir isotherm
function gprime =Deriv_Langmuir(c,b,qm)
gprime= (b*qm)./(b*c + 1) - (b^2*c*qm)./(b*c + 1).^2;
end