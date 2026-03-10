%function all_fits = radiation_wrapper()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

clearvars max_error dir;
%% USER INPUTS
element='Li'; %Character array
myplot=true;
weight_variation=.15;
maxTe=Inf;
minTe=0;
names=dir("C:\Users\Jonathan Roeltgen\Documents\MATLAB\radiation_operator\plt_data\plt*.txt");
start_over=true;
SIKE=false;
for charge_state=1:3
    fprintf("Now calculating charge state=%d, of %s\n",charge_state,element);
    if(charge_state==1||start_over)
         x0=[0.2e-1,8000,0.8,3,-4];
         changeV0=false;
         xu=[Inf,Inf,70,30,20];
    else
         x0=all_fits.(fit_name){charge_state-1,density_ind}.coeffs;
         x0(1)=x0(1)*1e30;
         changeV0=true;
         xu=[Inf,Inf,50,80,20];
    end
    changeV0=false;
    %x0(4)=72.759576;
    %x0(4)=x0(4)+5;
    %fit=all_fits_old.Al_fits{charge_state,end};
    x0=fit.coeffs; x0(1)=x0(1)*1e30;
    % x0=temp_fit(1).coeffs; x0(1)=x0(1)*1e30;   
    %x0(4)=22;
    %x0(1)=x0(1)*1e2;
    % x0(3)=.2;
    %xu(2)=80;
    xl=[1e-12,.01,.001,.1,-20]; 
    %xu(2)=1000;
    weight_powers=.0:.01:.3;%
    %weight_powers=linspace(weight_powers(18),weight_powers(20),20);
    ne=17;
    orig_ne=ne;
    toplot=false;
    publish=false;
    % End User Inputs
    %% Calculate fits
    fit_name=strcat(element,"_fits");
    max_error=1e30;
    if(SIKE)
        % elem=ArSIKELz;
        Te=log10(Ar_bimax_lj_unresolved.te);%log10(elem(:,1));
        Lz=log10((Ar_bimax_lj_unresolved.Lz(:,charge_state)))';%log10(elem(:,charge_state+1)*1e6)';
        % Te = log10(wpddatasets.Lithium2x);
        % Lz = log10(wpddatasets.Lithium2y)';
        % Te = log10(Arwpddatasets1.Argon11x);
        % Lz = log10(Arwpddatasets1.Argon11y)';
        %Te = log10(LiSIKELiLz.Te);
        %Lz = log10(LiSIKELiLz.Li1)';
        has_density_variation=false;
        dataTable=[];
    else
        for i=1:numel(names)
            if(numel(element)==1)
                if(names(i).name(7)==lower(element)&&names(i).name(8)=='_')
                    dataTable=readpltdata(strcat(names(i).folder,"\",names(i).name));
                    break;
                end
            else
                if(names(i).name(7:8)==lower(element))
                    dataTable=readpltdata(strcat(names(i).folder,"\",names(i).name));
                    break;
                end
            end
        end
        ne=dataTable(24,3)+6; % Specify a specific density index (first index in dataTable).
        ne_ints=get_ne_ints(element);
        has_density_variation=density_variation(dataTable,ne_ints);
        Te=[];
        Lz=[];
    end
    LogicalStr = {'false', 'true'};
    final=numel(weight_powers);
    first=1;
    best_fit=0;
    clearvars temp_fit all_error
    fprintf("alpha(0)=%e, beta(0)=%e, gamma(0)=%e, V0(0)=%e, A(0)=%e \n",x0(2),x0(3),x0(5),x0(4),x0(1));
    for i=first:1:final
        fprintf("Fit#=%d; Currently fitting with a weighting power=%f. \n",i,weight_powers(i));
        temp_fit(i)=radiation_operator(element, dataTable, charge_state, weight_powers(i), ...
            x0, xl, xu, ne, toplot, publish, changeV0, maxTe, minTe, Te, Lz); %#ok<SAGROW>
        fprintf("\t Pass test: ");
        for j=1:6
            if (temp_fit(i).sucesses(j))
                cprintf('green',"Pass ");
            else
                cprintf('red',"Fail ");
            end
        end
        fprintf("\n\t Fitting parameters found: alpha(0)=%e, beta(0)=%e, gamma(0)=%e, V0(0)=%e, A(0)=%e\t",temp_fit(i).coeffs(2),temp_fit(i).coeffs(3),temp_fit(i).coeffs(5),temp_fit(i).coeffs(4),temp_fit(i).coeffs(1));
        if(any(temp_fit(i).sucesses))
            if(changeV0)
                x0 = temp_fit(i).problem_struct.x0;
                x0(1)=x0(1)*1e30;
            end
        end
        fprintf("\n");
        if(all(temp_fit(i).sucesses))
            if(temp_fit(i).max_error(1)<max_error)
                best_fit=i;
                max_error=temp_fit(i).max_error(1);
            end
        end
        all_error(i,:)=temp_fit(i).max_error(:);
    end
    density_ind = temp_fit(1).density_ind;
    if(best_fit>0)
        temp_fit(best_fit).note="All pass";
        all_fits.(fit_name){charge_state,density_ind}=temp_fit(best_fit);
    else
        break;
    end
    if(myplot&&best_fit>0)
        figure;
        loglog(all_fits.(fit_name){charge_state,density_ind}.te,all_fits.(fit_name){charge_state,density_ind}.Lz,all_fits.(fit_name){charge_state,density_ind}.te,all_fits.(fit_name){charge_state,density_ind}.radiation);
        ylabel("Lz");
        xlabel("Te");
        yyaxis right;
        semilogx(all_fits.(fit_name){charge_state,density_ind}.te,all_fits.(fit_name){charge_state,density_ind}.ratio);
        ylabel("ratio");
        ylim([0 10])
        grid on;
    end
%%
    %return;
    if(best_fit>0&&has_density_variation)
        clearvars all_error;
        x0=all_fits.(fit_name){charge_state,density_ind}.coeffs;
        weight_powers=max((weight_powers(best_fit)-weight_variation),0):.01:(weight_powers(best_fit)+weight_variation);
        final=numel(weight_powers);
        for ne=-1:-1:-ne_ints
            best_fit=0;
            clearvars temp_fit;
            if(abs(ne)==orig_ne)
                continue;
            end
            max_error=1e30;
            fprintf("Currently fitting a density index=%f\n",abs(ne));
            for i=first:final
                fprintf("Fit#=%d; Currently fitting with a weighting power=%f\n",i, weight_powers(i));
                temp_fit(i)=radiation_operator(element, dataTable, charge_state, weight_powers(i), ...
                    x0, xl, xu, ne, toplot, publish, false, maxTe, minTe);
                fprintf("\tPass test: ");
                for j=1:6
                    fprintf("%s ",LogicalStr{temp_fit(i).sucesses(j)+1});
                end
                fprintf("Fitting parameters found: alpha(0)=%e, beta(0)=%e, gamma(0)=%e, V0(0)=%e, A(0)=%e\t",temp_fit(i).coeffs(2),temp_fit(i).coeffs(3),temp_fit(i).coeffs(5),temp_fit(i).coeffs(4),temp_fit(i).coeffs(1));
                fprintf("\n");
                if(all(temp_fit(i).sucesses))
                    if(temp_fit(i).max_error(1)<max_error)
                        best_fit=i;
                        max_error=temp_fit(i).max_error(1);
                    end
                end
                all_error(abs(ne),i,:)=temp_fit(i).max_error(:);
            end            
            if(best_fit>0)
                density_ind = temp_fit(1).density_ind;
                temp_fit(best_fit).note="All pass";
                all_fits.(fit_name){charge_state,density_ind}=temp_fit(best_fit);
            end
        end
    elseif(best_fit>0)
        all_fits.(fit_name){charge_state,density_ind}.note=strcat(all_fits.(fit_name){charge_state,density_ind}.note,"; No density dependence");
    end
end
if (SIKE)
    sike_fits=all_fits;
end

function has_density_variation=density_variation(dataTable,ne_ints)
    [~,ind]=min(abs(dataTable(:,1)-2));
    te_int = (ind-1)/ne_ints;
    has_density_variation=false;
    for i=1:te_int
        dev=std(dataTable((i-1)*ne_ints+1:i*ne_ints,4));
        if(dev>1e-8)
            has_density_variation=true;
        end
    end
end
function ne_ints = get_ne_ints(element)
    if(strcmp(element,"Li"))
        ne_ints = 16;
    elseif(strcmp(element,"B")||strcmp(element,"Ar")||strcmp(element,"F")||...
            strcmp(element,"S")||strcmp(element,"Cl"))
        ne_ints = 26;
    elseif(strcmp(element,"He")||strcmp(element,"C")||strcmp(element,"N")||...
            strcmp(element,"O")||strcmp(element,"Ne")||strcmp(element,"H")||...
            strcmp(element,"Be")||strcmp(element,"Al")||strcmp(element,"Si"))
        ne_ints = 24;
    else
        disp("Element not programmed.");
        exit;
    end
end