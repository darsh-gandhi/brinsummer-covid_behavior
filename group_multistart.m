function [psi_12, psi_21, a_1, a_2, beta] = group_multistart(data,initialvalues,tspan,lb,ub,w,nu)
 
% parameter order:
% z := [(1) pi, (2) alpha, (3) r, (4) sigma, (5) p, (6) gamma, (7) delta,
% (8) xi, (9) psi_12^max, (10) psi_21^max, (11) a_1, (12) a_2, (13) beta]

% How many runs for parameter estimation do you want:

NoStartPoints = 10;

%% Upper and Lower bounds for parameters you are fitting.

LowerBounds=lb;      %Lowerbounds for the parameters you are estimating
UpperBounds=ub;       %Upperbounds for the parameters you are estimating

xstart=.5*(LowerBounds+UpperBounds);                            %What initial parameter values you want to start with

%% MultiStart and fmincon - Fitting Part - Parallelization - not many comments ask Prof. Edholm for clarification if you want.

% Here we set-up the optimization problem, specifying we will use fmincon
% as the local solver, and the what model we want to minimize along with
% the specific measure down below in the SIR_RUN_ODE45 function. We give
% initial conditions and the bounds.

problem = createOptimProblem('fmincon','objective',@(z) group_functional(z,data,initialvalues, tspan,w,nu)...
    ,'x0',xstart,'lb',LowerBounds,'ub',UpperBounds);%,'Aineq',A,'bineq',b)%,'Aeq',Aeq,'beq',beq);

problem.options = optimoptions(problem.options,'MaxFunEvals',9999,'MaxIter',9999);%,'TolFun',0,'TolCon',0)
%problem.options = optimoptions(problem.options,'MaxFunEvals',inf,'MaxIter',inf,'TolFun',1e-10,'TolCon',0,'TolX',0,'MaxFunEvals',999999)

numstartpoints=NoStartPoints;                               % How many runs do you want?

% %  ms=MultiStart('Display','iter');                       %defines a multistart problem without parallel

ms=MultiStart('UseParallel',true,'Display','iter');         %defines a parallel multistart problem

%parpool %accesses the cores for parallel on your computer (laptop goes for 2-8, can be more specific)

[b,fval,exitflag,output,manymins]=run(ms,problem,numstartpoints);  %runs the multistart

% the following takes solutions from manymins and makes a matrix out of them


for i=1:length(manymins)
    Parameters(i,:)=manymins(i).X;       %what are the parameter values
end

for i=1:length(manymins)
    fvalues(i)=manymins(i).Fval;            %the minimization error
end

for i=1:length(manymins)
    ExitFlags(i)=manymins(i).Exitflag;      %how "good" is the solution, we want 1 or 2.
end


%delete(gcp('nocreate'))  %turns off the parallel feature


%% Plot the "best" solution

%%Outputs state variables for "best" fit
[t,y] = ode45(@(t,y) odesystem(t,y,Parameters(1,:)),tspan,initialvalues);

% S_1 = y(:,1);
% I_1 = y(:,2);
% R_1 = y(:,3);
% S_2 = y(:,4);
% I_2 = y(:,5);
% R_2 = y(:,6);
% H = y(:,7);
% D = y(:,8);

I_tot = y(:,2) + y(:,5);


figure
plot(tspan,I_tot,'LineWidth',1.5)
hold on
scatter(tspan,data, 'filled')
% title('Average Cases')
xlabel('Days')
hold off


psi_12 = Parameters(1,:);
psi_21 = Parameters(2,:);
a_1 = Parameters(3,:);
a_2 = Parameters(4,:);
beta = Parameters(5,:);

end



