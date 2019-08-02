function [Overindex ] = ica_fuse_falsemaxdetect(data, trendPara,LtrendPara)
% false maximun detect fucntion is to detect if a flase maximum occurs
% by checking entroy's trend along time.
% if there is  a decreaseing trend , then the false maximum occur;
% if osciilation for a long time without increasing occurs, Then the false
% maximum occur; 

if nargin<2
    LtrendPara= 1e-4;
    trendPara= -1e-3; % the parameter -1e-3 or -5e-4, or -1e-4; need to test on simulation for overfitting problem with  low correlation
elseif nargin<3
    LtrendPara= -1e-4;
end
if isempty(LtrendPara) LtrendPara= 0.0001 ;end
if isempty(trendPara) trendPara= -0.001; end
  
Overindex=0;

% if osciilation for a long time without increasing occurs, Then the false maximum occur; 
 n=length(data); 
% 
 if  n>60;
     datat=data(n-49:n)-mean(data(n-49:n));
     p=polyfit([1:50],datat,1);
     if abs(p(1))<LtrendPara 
        Overindex=1;
     end
     
 end
% 
if ~Overindex
    datat=data(n-4:n)-mean(data(n-4:n));
    p=polyfit([1:5],datat,1);
    r=datat-polyval(p,[1:5]);
    if p(1)<trendPara 
     %   if sum(r.^2)/10 < varPara
            Overindex=1;
     %   end
    end      
 
end
