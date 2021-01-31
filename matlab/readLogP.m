function [logP]=readLogP(filename)

string_ = fileread(filename);

lenS=size(string_,2);

tempS_lo = []; 
tempS_up = []; 
lo_read =0;
up_read =0;
i=1;
while(string_(i)~='+'||(string_(i)=='+') && string_(i-1)=='e')
    if(string_(i)=='[')
         lo_read =1;
    elseif(string_(i)==',' && lo_read==1)
         lo_read=0;
         up_read=1;
         logI(1) = str2double( tempS_lo );
         tempS_lo = [];
    elseif(string_(i)==']' && up_read==1)
         up_read=0;
         logI(2) = str2double( tempS_up ); 
         tempS_up = []; 
    elseif(lo_read==1)
         tempS_lo = strcat( tempS_lo, string_(i));
    elseif(up_read==1)
         tempS_up = strcat( tempS_up, string_(i));     
    end
    i = i+1;
end
i = i+1;

logP = [];
tempS_lo = []; 
tempS_up = []; 
lo_read =0;
up_read =0;
while (i<=lenS)
    if(string_(i)=='[')
         lo_read =1;
    elseif(string_(i)==',' && lo_read==1)
         lo_read=0;
         up_read=1;
         temp_lo = str2double( tempS_lo );
         tempS_lo = [];
    elseif(string_(i)==']' && up_read==1)
         up_read=0;
         temp_up = str2double( tempS_up );
         temp_mid = (temp_lo+temp_up)/2;
         tempS_mid = num2str( temp_mid );
         tempS_up = []; 
         if(temp_mid<0)
             %tempS_mid= strcat('(',tempS_mid,')');
             if(strlength(logP)>0)
             logP = extractBefore(logP,strlength(logP));
             end
             logP = strcat(logP,tempS_mid);
         else
         logP = strcat(logP,tempS_mid);
         end
    elseif(lo_read==1)
         tempS_lo = strcat( tempS_lo, string_(i));
    elseif(up_read==1)
         tempS_up = strcat( tempS_up, string_(i));
    else
         logP = strcat(logP,string_(i));
    end

    i=i+1;
end

%%if logP is too long, split it into several subset. (long logP will cause issue in str2sym function)
    logP0 =logP;
    clear logP;
    k=1;
    while(strlength(logP0)>1000)
        i=1001;
        while(logP0(i)~='+')
            i=i+1;
        end
        logP{k}=extractBefore(logP0,i);
        logP0 = extractAfter(logP0,i);
        k=k+1;
    end
    logP{k}=logP0;   
end