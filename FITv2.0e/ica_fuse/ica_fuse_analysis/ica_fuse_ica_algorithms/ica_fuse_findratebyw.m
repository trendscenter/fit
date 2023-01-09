function [lrate2,A1,A2,change]=findratebyw(lrate1,A1,A2,oldweights,wchwindow,whiteM,indt);
% add constraint to A
sh=size(A1,1);
sp=size(A2,1);
chans=size(A1,2);
a1=mean(A1);
a2=mean(A2);
lrate2=lrate1;
A11=A1;
A22=A2;
    for i=indt
        a1=A1(:,i);
        a2=A2(:,i);
        % [h p]=ttest2(a1,a2);

        % fk1=-(mean(a1)-mean(a2))^2/(nanvar(a1)*(sh-1)+nanvar(a2)*(sp-1));
        commonterm=(mean(a1)-mean(a2))/(nanvar(a1)*(sh-1)+nanvar(a2)*(sp-1));
        deltafk1=2*(-commonterm/sh+commonterm.^2*(a1-mean(a1)));
        deltafk2=2*(commonterm/sp+commonterm.^2*(a2-mean(a2)));

        a1=a1-lrate2*deltafk1;%
        a2=a2-lrate2*deltafk2;%

        A11(:,i)=a1;
        A22(:,i)=a2;
        %  fk2=-(mean(a1)-mean(a2))^2/(nanvar(a1)*(sh-1)+nanvar(a2)*(sp-1));
    end

    AA=[A11;A22];
    weights=inv(whiteM*AA);
    oldwtchange = weights-oldweights;
    delta=reshape(oldwtchange,1,chans*chans);
    tempchange=delta*delta';

    wchwindow1=[wchwindow,tempchange];
    X=[ones(1,10);1:10]';
    a=X\wchwindow1';
    if a(2)>0.001
         lrate2=lrate2*0.1;
         AA=[A1;A2];
            weights=inv(whiteM*AA);
            oldwtchange = weights-oldweights;
            delta=reshape(oldwtchange,1,chans*chans);
            change=delta*delta';
    elseif a(2)>0.00001
           lrate2=lrate1*0.9;
             AA=[A1;A2];
            weights=inv(whiteM*AA);
            oldwtchange = weights-oldweights;
            delta=reshape(oldwtchange,1,chans*chans);
            change=delta*delta';
    else
        %lrate2=lrate1;
        change=tempchange;
        A1=A11;
        A2=A22;
    end


