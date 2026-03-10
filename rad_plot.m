function ratio = rad_plot(x,xn,te,y,radiation,element,charge_state,weight_power,y2,scale,publish,ne)
    xns2(1) = string(num2str(xn(1),"%.1e"));
    for i=2:5
        xns2(i)=string(num2str(xn(i),'%.1f'));
    end
    legend(gca,'off')
    fstring2 = strcat("Fit=$\frac{",xns2(1),"(",xns2(2),"+",xns2(3),")}{",xns2(3),"(v/",xns2(4)...
        ,")^{-",xns2(2),"}+",xns2(2),"(v/",xns2(4),")^{",xns2(3),"}}v^{",xns2(5),"}$");
    exponent = floor(ne);
    coeff = 10^ne/10^exponent;
    data_legend=strcat("Data (n$_e=",num2str(coeff,"%.2f"),"\times 10^{",num2str(exponent,"%d"),"}m^{-3}$)");
    loglog(x,radiation,'DisplayName',fstring2);
    
    hold on;
    loglog(10.^te,10.^y,':x',"MarkerSize",8,'DisplayName',data_legend);
    legend(gca,'show',"Interpreter","latex",'FontSize',20);

    grid on;
    xlabel("T_e (eV)");
    ylabel("L_z (Wm^3)");
    set(gca,'FontSize',20);
    if(publish)
        title(strcat(element,"^{+",num2str(charge_state),"}"));
    else
        title(strcat(element,"^{+",num2str(charge_state-1),"} Weight power: ",num2str(weight_power))); %#ok<*UNRCH>
    end
    %
    toplot=true;
    if(~publish)
        ratio=error_plotting(x,y2/scale,radiation,weight_power,charge_state,element,toplot);
    else
        ratio = max(radiation./(y2/scale),y2/scale./radiation);
    end
end

function ratio=error_plotting(te,data,calculated,N,charge_state,element,plotting)
    ratio = max(calculated./data,data./calculated);
    if(plotting)
        figure;semilogx(te,ratio);
        yyaxis right;
        loglog(te,data./max(data));
        title(strcat(element,'+',num2str(charge_state-1),' Weights to power:',num2str(N)));
        xlabel("T_e");
        ylabel("Data/max(data)");
        yyaxis left;
        ylabel("Max(fit/data,data/fit), (error ratio)");
        ylim([0 10]);
        figure;
        loglog(data./max(data),ratio);
        xlabel("Data/max(data)");
        ylabel("ratio")
        title(strcat(element,'+',num2str(charge_state-1),' Weights to power:',num2str(N)));
        ylim([0 8]);
        grid on;
        xlim([1e-8 1]);
    end
    [perror] = abs((calculated-data)./data);
    [maxperror,ind]=max(perror); 
    [mratio,ind1] = max(ratio);
    [mTeratio,indTe] = max(ratio(te>1));
    [mLzratio,indLz] = max(ratio(data>1e-4*max(data)));
    success = true;
    if(any(ratio(and(te>1,data>1e-2*max(data)))>1.2))
        fprintf("Fails ratio <1.2 at L_z<0.01L_{z,max}\n");
        success =false;
    end
    if(any(ratio(and(te>1,and(data>1e-4*max(data),data<1e-2*max(data))))>1.4))
       % disp(te(ratio(and(te>1,data>1e-4*max(data)))>1.4));
       % disp(te(ratio(and(te>1,and(data>1e-8*max(data),data<1e-4*max(data))))>2));
        fprintf("TRY TO REDUCE ERROR! Fails for Lz*1e-4<data<Lz*1e-2 with Te>1\n");
        success = false;
    end
    if(any(ratio(and(te>1,and(data>1e-8*max(data),data<1e-4*max(data))))>2))
        fprintf("TRY TO REDUCE ERROR! Fails for Lz*1e-8<data<Lz*1e-4 with Te>1\n");
        success = false;
    end
    if(any(and(data<max(data)*1e-8,data>10.^(.5*log10(max(data)./data)))))
        fprintf("Signficant error at low L_z. Check plots.\n");
        success = false;
    end
    if(success)
        fprintf("Within Error tolerances.\n");
    end
    if(plotting)
        indLz = indLz + sum(data<1e-4*max(data));
        fprintf("Errors %e %f %f\n",mratio,mTeratio,mLzratio);
        fprintf("Lz at max ratio: %e %e %e\n",data(ind1)/max(data),...
            data(indTe+sum(te<=1))/max(data),data(indLz)/max(data));
        fprintf("Max ratio: %f Mean ratio: %f\n",mratio,mean(ratio));
    
        [perror] = max(calculated./data,data./calculated)-1;
        fprintf("Max of percent off from 1: %f\n",...
            max(perror));
        fprintf("Mean of percent off from 1: %f\n",...
            mean(perror));
    end
end