%% IUQ Procedure Module of Two-Step Parameter Estimation in Dual-Functional Gas Purification: Step 1
%======================================================================================================================%==========================================================================================================
%> @details     This function, uq_DualFunc_GasPurification_Step1, is designed to execute the Inverse Uncertainty Quantification(IUQ) procedure within 2-step parameter estimation for dual-functional gas purification processes.
%>              It processes multiple data sets to estimate adsorption model parameters using Bayesian analysis.
%>              The function includes:

%>              (1) Parameter Transformation and Looping Over Data Sets:
%>                  - Transforms parameters from logarithmic to linear scale.
%>                  - Iterates over multiple experimental data curves.
%>
%>              (2) Parameter Calculation for Each Curve:
%>                  - Calculates adsorption-related parameters like mass transfer coefficients.
%>                  - Determines the breakthrough curve using a forward adsorption model.
%>
%>              (3) Data Interpolation and Differentiation:
%>                  - Interpolates the estimated concentration data to reference time points.
%>                  - Calculates the rate of change of concentration.
%>
%>              (4) Error Handling and Output Data Compilation:
%>                  - Implements checks for unrealistic model outputs.
%>                  - Compiles the differential rate and time-to-breakthrough data for each parameter set.
%======================================================================================================================%==========================================================================================================
function t_IUQ = uq_DualFunc_GasPurification_Step1(X,tref_multi,cref_multi,dt,Ngrid,dz,cFeed,epsilon,tf,v,num,rho,rho_b,rho_g,dp,mu,NumCurve,max_dcdt,max_tbr,Ndata,Max_C)

Npar = size(X, 1); % Number of parameter sets
for j = 1:Npar
    % Transform parameters from log scale to actual scale
    par_tmp =10.^X(j,1:6);
    
    for i = 1:NumCurve
        c_ref = cref_multi(1+Ndata*(i-1):Ndata+Ndata*(i-1));
        t_ref = tref_multi(1+Ndata*(i-1):Ndata+Ndata*(i-1));
        
        if num==1
            Dm = par_tmp(3);
            Sc = mu/(rho_g*Dm);
            Re = v*dp*rho_g/mu;
            Dz = Dm * (20+ 0.5*Sc*Re)/epsilon;
            kg = Dm/dp*(2.0+ 1.8 * Re^0.5 * Sc^(1/3));
            Rp = dp/2;
            De = par_tmp(4);
            Isotherm_param = [par_tmp(1),par_tmp(2)];
            qstar_phys= Isotherm_Langmuir(cFeed(i),par_tmp(1),par_tmp(2));
            kc= par_tmp(5);
            kd= par_tmp(6);
            q0star = qstar_phys;
            inv_K = (Rp*rho_b*q0star)./(3*kg*cFeed(i)*epsilon) + (Rp^2*rho_b*q0star)./(15*De*cFeed(i)*epsilon);
            K_G = 1/inv_K;
        end
        par = [K_G,Isotherm_param,Dz,kc,kd]; %disp(log10(par))

        [t_est,c_tmp,~]= Model_DualFunc_GasPurification(dt,Ngrid,dz,par,epsilon,cFeed(i),tf(i),num,v,rho);
        c_breakthrough=(c_tmp(Ngrid,:)+c_tmp(Ngrid+1,:))./2;
        c_est = c_breakthrough;
        
        idx1= find(c_est> 0,1,'first');
        idx2= length(c_est);
        idx3=find(c_est> 0 & c_est<= Max_C(i)* cFeed(i),1,'last');
        
        if isempty(idx1) || isempty(idx2) || isempty(idx3)
            idx1=1;
            idx2= idx1;
            idx3= idx1;
        end
        if length(t_est)< 3 ||max(c_est)>cFeed(i) ||min(c_est)< 0 || abs(idx2-idx1)<=2 || abs(idx3-idx1)<=2
            tmp_dcdt = 100*ones(1,length(c_ref)-1);
            tmp_tbr(i)= 100;
        else
            c_interp =  interp1(t_est,c_est,t_ref,'linear', 'extrap');
            dt = t_ref(2) - t_ref(1);
            dc_dt = diff(c_interp) / dt;
            tmp_dcdt=  dc_dt/max_dcdt(i);
            
            interp_func = @(c) interp1(c_est(idx1:idx3)/cFeed(i), t_est(idx1:idx3)/60, c, 'linear', 'extrap');
            specific_c = 0.01;
            specific_t = interp_func(specific_c);
            tmp_tbr(i)= specific_t;
        end
        
        if i==1
            tmp_stack =tmp_dcdt;
        else
            tmp_stack = [tmp_stack,tmp_dcdt];
        end
    end
    
    % Compile differential rate and time-to-breakthrough data
    t_IUQ(j, :) = [((Ndata - 1) / Ndata) * tmp_stack, (1 / Ndata) * tmp_tbr / max_tbr];
end
end