%% Forward Modeling for Dual-Functional Gas Purification Process
%=========================================================================================================================================================================================================================================================================
%> @details     This function, Model_DualFunc_GasPurification, performs forward modeling for dual-functional gas purification. It simulates concurrent mechanisms within fixed-bed, considering adsorption and catalytic decomposition, with various isotherm models 
%>              and solves the system using an implicit solver with adaptive time stepping. Key steps include:

%>              (1) Model Initialization:
%>                  - Transforms model parameters from logarithmic to linear scale.
%>                  - Sets up initial conditions and parameters based on the chosen isotherm model (In this case, Langmuir model).
%>
%>              (2) Time-stepping and Implicit Solving:
%>                  - Iterates through time steps, updating concentrations and removed quantities by concurrrent removal mechanisms (adsorption and catalytic decomposition) within dual-functional adsorbent.
%>                  - Applies a finite volume method with implicit solving for the nonlinear model.
%>
%>              (3) Boundary Conditions and Matrix Assembly:
%>                  - Implements boundary conditions and assembles the system matrix for solving.
%>
%>              (4) Newton Iteration and Convergence:
%>                  - Performs Newton iteration for nonlinear system convergence.
%>                  - Checks and records convergence metrics.
%>
%>              (5) Results Compilation and Breakthrough Detection:
%>                  - Compiles concentration and adsorption data at each time step.
%>                  - Detects breakthrough times based on concentration thresholds.
%=========================================================================================================================================================================================================================================================================
function [t,c,q]= Model_DualFunc_GasPurification(dt,ngrid,dz,par,epsilon,cFeed,tf,num,us,rho)
% Define model parameter 
if num==1 %num 1= Langmuir
    k_LDF= par(1);
    b_phys=par(2);
    qm_phys=par(3);
    D= par(4);
    kc= par(5);
    kd= par(6);    
end

% Time stepping and implicit solving
t_tmp= 0:dt:tf;
ntime = length(t_tmp);

idx= zeros(2,1);
t= zeros(ntime,1);
c = zeros(ngrid+1,ntime); 
c(1,1) = cFeed;
q = zeros(ngrid+1,ntime); 
tol_save = zeros(ntime,100);

alpha= rho*(1-epsilon)/epsilon;
ui= us/epsilon; % interstitial velocity
beta = ui*dt/(2*dz);
gamma = D*dt/(2*dz*dz); 
lambda = k_LDF*dt/2;

Aim = -(beta + gamma); 
Aip = - gamma;

for nt = 2:ntime
    t(nt)=t(nt-1)+dt;
    t_cn= (t(nt)+t(nt-1))/2;
    eta = (alpha*kc*dt)/(2*(1+kd*t_cn));
    
    % Preallocate  and initialize necessary vector
    n = ngrid-1;
    Amatrix= zeros(n);
    delta_old= zeros(n,1);
    F= zeros(n,1);
    
    c_old = c(:,nt-1);
    q_old= q(:,nt-1);
    
    if nt==2
        if num==1
            g_phys_old=Isotherm_Langmuir(c_old,b_phys,qm_phys);
            g_old = g_phys_old;
        end
    end
    
    % Newton iteration for solving the nonlinear governing equations 
    for iter =1:5
        if iter==1
            c_tmp =  c_old;
        end
        
        for nz=2:ngrid
            if num==1
                g_phys_tmp=Isotherm_Langmuir(c_tmp(nz),b_phys,qm_phys);
                gprime_phys_tmp=Deriv_Langmuir(c_tmp(nz),b_phys,qm_phys);
                g_tmp = g_phys_tmp;
                gprime_tmp = gprime_phys_tmp;
            end

            F_conv= beta*(c_tmp(nz)-c_tmp(nz-1)+c_old(nz)-c_old(nz-1));
            F_diff= gamma*(c_tmp(nz+1)-2*c_tmp(nz)+c_tmp(nz-1)+c_old(nz+1)-2*c_old(nz)+c_old(nz-1));
            q_tmp= (lambda*(g_tmp+ g_old(nz))+(1-lambda)*q_old(nz))/(1+lambda);
            F_kinetic= alpha*lambda*(g_tmp+g_old(nz)- q_tmp-q_old(nz));
            F_catalytic = eta*(c_tmp(nz)+c_old(nz));
            
            F(nz-1)= c_tmp(nz)-c_old(nz) + F_conv - F_diff + F_kinetic + F_catalytic;
            
            Aic = 1 + beta + 2*gamma + alpha*lambda*gprime_tmp/(1+lambda) + eta;
            Amatrix(nz-1,nz-1) = Aic;
            
            if(nz-1 ~= 1)
                Amatrix(nz-1,nz-2) = Aim;
                Amatrix(nz-2,nz-1) = Aip;
            end
            
            % Boundary treatment 
            Amatrix(1,1) = Amatrix(1,1) + gamma;
        end    
        
        RHS = -F;
        
        if iter==1
            max_Fold = max(F);
        end
        
        % Solving the system using Thomas algorithm
        delta = thomas(Amatrix,RHS);
        
        % Updating concentration and checking convergence
        c_tmp(2:ngrid)= c_tmp(2:ngrid)+delta;
        c_tmp(1) = 2*cFeed  - c_tmp(2);
        c_tmp(ngrid+1) = c_tmp(ngrid);
        tol= sqrt(sum((delta-delta_old).^2)/length(delta));
        tol_save(nt,iter)= tol;
        delta_old= delta;
        max_F= max(F);
        relF= abs(max_F/max_Fold);
        
       % Break from iteration upon convergence
        if (relF <= 1e-10 || iter>=5)
            c_fin= c_tmp;
            break;
        end
    end 
    
    % Update concentration and adsorption quantities for current time step    
    c(:,nt) = c_fin;
    
    if num==1
        g_phys_new=Isotherm_Langmuir(c(:,nt),b_phys,qm_phys);
        g_new = g_phys_new;
    end
    
    q(:,nt) = (lambda*(g_new+g_old) + (1-lambda)*q_old)/(1+lambda); % calculating new q
    g_old = g_new;
end 
end
%%
function x = thomas(A,d)
% Thomas algorithm for solving x  tridiagonal linear sys Ax=d
% determines n
n = length(d);

% preallocates all necessary vectors
a = zeros(n-1,1);
b = zeros(n,1);
c = zeros(n-1,1);
x = zeros(n,1);

% extracts first element of "b" from "A"
b(1) = A(1,1);
% forward loop
for i = 2:n
    % extract relevant elements of "a", "b", and "c" from "A"
    a(i-1) = A(i,i-1);
    b(i) = A(i,i);
    c(i-1) = A(i-1,i);
    % forward elimination
    w = a(i-1)/b(i-1);
    b(i) = b(i)-w*c(i-1);
    d(i) = d(i)-w*d(i-1);
end

% backward loop (backward substitution)
x(n) = d(n)/b(n);
for i = (n-1):(-1):1
    x(i) = (d(i)-c(i)*x(i+1))/b(i);
end
end