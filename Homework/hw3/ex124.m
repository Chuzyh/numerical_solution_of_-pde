theta = linspace(0, 2*pi, 1000)*1i;
for p=1:4
    hold on;
    z=Adams_Bahsforth_r(p,theta);
    plot(real(z), imag(z));
end
legend('p=1','p=2','p=3','p=4') 
title('The RASs of Adams-Bashforth formulas') 
saveas(gcf, 'Adams-Bashforth.jpg');

cla
for p=3:5
    hold on;
    z=Adams_Moulton_r(p,theta);
    plot(real(z), imag(z));
end
legend('p=3','p=4','p=5') 
title('The RASs of Adams-Moulton formulas') 
saveas(gcf, 'Adams-Moulton.jpg');

cla
for p=1:4
    hold on;
    z=backforth_r(p,theta);
    plot(real(z), imag(z));
end
legend('p=1','p=2','p=3','p=4') 
title('The RASs of backforth formulas') 
saveas(gcf, 'backforth.jpg');
function [re]=Adams_Bahsforth_r(p,theta)
    if(p==1)
        re=exp(theta)-1;
    else if(p==2)
            re=2*(exp(2*theta)-exp(theta))./(3*exp(theta)-1);
        else
            if(p==3)
                re=12*(exp(3*theta)-exp(2*theta))./(23*exp(2*theta)-16*exp(theta)+5);
            else
                if(p==4)
                    re=24*(exp(4*theta)-exp(3*theta))./(55*exp(3*theta)-59*exp(2*theta)+37*exp(theta)-9);
                end
            end
        end
    end
end
function [re]=Adams_Moulton_r(p,theta)
    if(p==3)
        re=12*(exp(2*theta)-exp(theta))./(5*exp(2*theta)+8*exp(theta) - 1);
    else if(p==4)
            re=24*(exp(3*theta)-exp(2*theta))./(9*exp(3*theta)+19*exp(2*theta)-5*exp(theta)+1);
        else
            if(p==5)
                re=720*(exp(4*theta)-exp(3*theta))./(251*exp(4*theta)+646*exp(3*theta)-264*exp(2*theta)+106*exp(theta) - 19);
            end
        end
    end
end

function [re]=backforth_r(p,theta)
    if(p==1)
        re=(exp(theta)-1)./exp(theta);
    else if(p==2)
            re=(3*exp(2*theta)-4*exp(theta)+1)./(2*exp(2*theta));
        else
            if(p==3)
                re=(11*exp(3*theta)-18*exp(2*theta)+9*exp(theta)-2)./(6*exp(3*theta));
            else
                if(p==4)
                    re=(25*exp(4*theta)-48*exp(3*theta)+36*exp(2*theta)-16*exp(theta)+3)./(12*exp(4*theta));
                end
            end
        end
    end
end