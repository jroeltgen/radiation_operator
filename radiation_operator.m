%% Radiation operator function
% DESCRIPTIVE TEXT
function radfit=radiation_operator(element, dataTable, charge_state, weight_power, ...
    x0, xl, xu, ne, toplot, publish, changeV0, maxTe, minTe, te, lz)
if(~exist('toplot','var')||isempty(toplot))
    toplot = false;
end
if(~exist("ne",'var')||isempty(ne))
    ne=19;
end
if(~exist("x0",'var')||isempty(x0))
    x0=[.01,1000,2,3,-3];
end
if(~exist('xl','var')||isempty(xl))
    xl=[1e-15,.01,.01,.1,-Inf];
end
if(~exist('xu','var')||isempty(xu))
    xu=[Inf,Inf,Inf,80,80];
end
if(~exist('weight_power','var')||isempty(weight_power))
    weight_power=0.0;
end
if(~exist('charge_state','var')||isempty(charge_state))
    charge_state=1;
end
if(~exist('publish','var')||isempty(publish))
    publish=true;
end
if(~exist('maxTe','var')||isempty(maxTe))
    maxTe=Inf;
end
if(~exist('minTe','var')||isempty(maxTe))
    minTe=0;
end
global mass %#ok<GVMIS>
mass = 9.10938215e-31; % Same as Gkeyll

yfixed = -78;
if(exist('te','var')&&exist('lz','var')&&~isempty(te)&&~isempty(lz))
    y=lz;
    ind=1;
    ne=20;
else
    [te,myNe,myPlt]=getData(element, charge_state, dataTable);
    if(ne<0)
        ind=-ne;
        ne=myNe(ind);
    else
        [~,ind]=min(abs(myNe-ne));
    end

    y=myPlt(:,ind);%ne of 1e19
end

scale = 10^30;
if(~toplot)
    display='off';
else
    display='iter';
end
warning('off','MATLAB:nearlySingularMatrix')
problem2.solver = 'fmincon';
problem2.options = optimoptions('fmincon','Algorithm','interior-point','MaxIterations',...
    8000,'Diagnostics','off','Display',display,'FiniteDifferenceType','central','FunctionTolerance',...
    1e-6,'FunValCheck','on','MaxFunctionEvaluations',12000,'StepTolerance',1e-12,'FunctionTolerance',1e-12,...
    'OptimalityTolerance',1e-18,'ConstraintTolerance',1e-12,'FiniteDifferenceStepSize',sqrt(eps),...
    'DiffMaxChange',0.1);
[x,y2]=preprocess(te',y,scale,yfixed,maxTe,minTe);

%Constraints:
% alpha+gamma>0 -> -(x(2)+x(5))<0
% gamma -beta<0 -> (x(5)-x(3))<2
A=[0,-1,0,0,-1;0,0,-1,0,1];
b=[0;2];
problem2.A=A;
problem2.b=b;
problem2.objective=@(x0p) objective(x,y2, x0p, weight_power);
problem2.x0 = x0;
problem2.lb = xl;
problem2.ub = xu;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% fmincon %%%%%%%%%%%%%%%%%
successes=[false,false];
[~,max_ind]=max(y);
%while(x0(4)<10^te(max_ind) && all(~successes))
while(x0(4)<200 && all(~successes))
    [xn,fval,exitflag,output,lambda,grad,hessian]=fmincon(problem2);
    radiation=fun(x,xn)/scale;
    ratio = max(radiation./(y2/scale),y2/scale./radiation);    
    if(exist('changeV0','var')&&~isempty(changeV0)&&changeV0)        
        xu(4)=x0(4)+20;
        problem2.x0 = x0;
        [successes,~]=error_analysis(ratio,x,y2/scale);
        if(all(~successes))
            fprintf("Increasing x0 from %f to %f.  ",x0(4),x0(4)*1.25);
        end
        x0(4)=x0(4)*1.25;
    else
        x0(4)=1e10;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xn(1)=xn(1)/scale;
if(toplot)
    rad_plot(x,xn,te,y,radiation,element,charge_state,weight_power,y2,scale,publish);
end
N=3e2;
v = [eps,logspace(0,8,N)]; % V par & mu
radfit = create_rad_struct(fval,exitflag,output,lambda,grad,hessian,xn,problem2,ratio,weight_power, ...
    radiation,x,y2/scale,v,ind,ne,minTe,maxTe);
nu2=nuOfV(vpa(v(1:2)),vpa([xn(1),xn(2),xn(3),xn(4)*scaled_mass(),xn(5)]));
nuratio = nu2(1)/nu2(2);
if(nuratio>1)
    fprintf("Warning! Nu possibly going to infinity as v->0");
    fprintf("nu(eps)/nu(1): %s\n",nuratio);
end
warning('on','MATLAB:nearlySingularMatrix')
end

function [ residual ] = objective(x, y, params, weight_power)
%objective - The objective function to be minimized
    calcY = fun(x,params);
    residual = sum(((calcY-y).*y.^(weight_power-1)).^2);
end

function [x,y]= preprocess(te,Lz,scale,ymin,temax,minTe)

    if(~exist("temax","var")||isempty(temax))
        temax=Inf;
    end
    if(~exist("temin","var")||isempty(temin))
        temin=0;
    end
    te = 10.^te;
    % Rescale y and convert from log space
    mask = Lz>ymin+1 & te<temax & te>temin;
    x = te(mask);
    y = 10.^Lz(mask)*scale;
end

function F=fun(xdata,x)
% The function to be minimized, 0 lower bound since spherical coordinates
    p = x;
    F= lhs(p,xdata,0,0);
end

function integrated=lhs(p,T,vb,lower_bound)
   integrated = integral(@(v) integrand(v,p,T,vb),lower_bound,inf,'AbsTol',1e-8,'RelTol',...
    1e-8,'ArrayValued',true)./T.^(3/2);
end

function int = integrand(v,p,T,vb)
   int = v.^4.*nuOfV(v,p).*exp(-(v-vb).^2./T);
end


% Function to breakdown data from given table, charge state and element
% @param element: string with atomic symbol
% @param charge_state: charge state of desired fit
% @param dataTable: table of coefficients formatted in ADAS style
% @return myTe: temperature intervals
% @return myNe: density intervals
% @return myPlt: radiation coefficients
function [myTe,myNe,myPlt]=getData(element, charge_state, dataTable)
    arr = dataTable;
    %arr = table2array(dataTable);
    if(strcmp(element,"H"))
        te_intervals = 29;
        ne_intervals = 24;
    elseif(strcmp(element,"Li"))
        te_intervals = 25;
        ne_intervals = 16;
    elseif(strcmp(element,"B")||strcmp(element,"Ar")||strcmp(element,"F")||...
            strcmp(element,"S")||strcmp(element,"Cl"))
        te_intervals = 48;
        ne_intervals = 26;
    elseif(strcmp(element,"He")||strcmp(element,"C")||strcmp(element,"N")||...
            strcmp(element,"O")||strcmp(element,"Ne")||strcmp(element,"Si"))
        te_intervals = 30;
        ne_intervals = 24;
    elseif(strcmp(element,"Be"))
        te_intervals = 25;
        ne_intervals = 24;
    elseif(strcmp(element,"Al"))
        te_intervals = 40;
        ne_intervals = 24;
    else
        disp("Element not programmed.");
        exit;
    end
    count = 1;
    te = linspace(0,1,te_intervals);
    ne = linspace(0,1,ne_intervals);
    plt = zeros(charge_state,te_intervals,ne_intervals);
    for i=1:charge_state
        for j=1:te_intervals
            te(j)=arr(count,2);
            for k=1:ne_intervals
                ne(k)=arr(k,3);
                plt(i,j,k)=arr(count,4);
                count=count+1;
            end
        end
    end
    myTe=te;
    myNe=ne+6;
    myPlt=squeeze(plt(charge_state,:,:))-6;
end



