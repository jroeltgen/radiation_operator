function [sucesses, max_error]=error_analysis(ratio,te,Lz)
    Lz_ratio=Lz./max(Lz);
    Lz_ind_1e2=Lz_ratio>1e-2 & te>1;
    Lz_ind_1e4=Lz_ratio>1e-4 & te>1;
    Lz_ind_1e8=Lz_ratio>1e-8 & te>1;
    Lz_mag=abs(log10(Lz_ratio));
    sucesses=[true,true,true,true,true,true];
    if(any(ratio(Lz_ind_1e4)>2))
        sucesses(1)=false;
    end
    if(any(ratio(Lz_ind_1e8)>8))
        sucesses(2)=false;
    end
    if(any(and(log10(ratio)>1/2*Lz_mag,Lz_mag>log10(8))))
        sucesses(3)=false;
    end
    if(any(ratio(Lz_ind_1e2)>1.2))
        sucesses(4)=false;
    end
    if(any(ratio(Lz_ind_1e4)>1.4))
        sucesses(5)=false;
    end
    if(any(ratio(Lz_ind_1e8)>2))
        sucesses(6)=false;
    end
    max_error(1)=max(ratio(Lz_ind_1e2));
    max_error(2)=max(ratio(Lz_ind_1e4));
    max_error(3)=max(ratio(Lz_ind_1e8));
end