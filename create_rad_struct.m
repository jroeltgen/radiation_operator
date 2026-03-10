function radfit = create_rad_struct(fval,exitflag,output,lambda,grad,hessian,xn,problem2,ratio,weight_power, ...
    radiation,x,y2,v,ind,ne,min_te,max_te,note)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    % Assign output
    mass = 9.10938215e-31; % Same as Gkeyll
    electron_charge = 1.602176487e-19; 
    %c_const = 8*sqrt(pi)*electron_charge^(5/2)/mass;
    c_const = 8*electron_charge/sqrt(pi)*(2*electron_charge/mass)^(xn(5)/2);
    bmagmin = 2.2;
    bmag0 = 1;
    bmagmax = 7.4;
    radfit.fval=fval;
    radfit.exitflag=exitflag;
    radfit.output=output;
    radfit.lambda=lambda;
    radfit.grad=grad;
    radfit.hessian=hessian;
    radfit.coeffs=xn;
    radfit.problem_struct=problem2;
    radfit.ratio=ratio;
    radfit.weight_power=weight_power;
    radfit.radiation= radiation;
    radfit.te = x;
    radfit.Lz = y2;
    radfit.min_te=min_te;
    radfit.max_te=max_te;
    N = numel(v);
   
    radfit.vpar = v;
    radfit.vperp = v;
    radfit.mu = mass*radfit.vperp.^2/(2*bmag0);

   % radfit.mu0 = radfit.mu;
   % radfit.mumax = 2*mass*radfit.vperp.^2/bmagmax;
   % radfit.mumin = 2*mass*radfit.vperp.^2/bmagmin;
    radfit.vmag = zeros(N,N);
    for i=1:N
        radfit.vmag(i,:)=sqrt(radfit.vpar.^2+radfit.vperp(i).^2);
    end
    radfit.corresponding_nu = nuOfV(radfit.vmag,[xn(1),xn(2),xn(3),xn(4)*scaled_mass(),xn(5)])./(c_const);
    for i=1:N
        radfit.nuprime(i,:)=radfit.vpar.*radfit.corresponding_nu(i,:);
        radfit.nudoubleprime(:,i)=2*radfit.mu'.*radfit.corresponding_nu(:,i);
        %radfit.nudoubleprime0(:,i)=radfit.nudoubleprime(:,i);
        %radfit.nudoubleprimemax(:,i)=radfit.mumax'.*radfit.corresponding_nu(:,i);
        %radfit.nudoubleprimemin(:,i)=radfit.mumin'.*radfit.corresponding_nu(:,i);
    end
    [radfit.sucesses, radfit.max_error]=error_analysis(ratio,radfit.te,radfit.Lz);
    radfit.density_ind = ind;
    radfit.ne = ne;
    if(~exist("note","var")||isempty(note))
        radfit.note="";
    else
        radfit.note =note;
    end
end
