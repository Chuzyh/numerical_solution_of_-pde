clc; clear;

format rat;
b=zeros(5,5);
for s=1:5
    a=zeros(s+1,1);
    co=zeros(s,s);
    rhs=zeros(s,1);
    a(s+1)=1;
    a(s)=-1;
    for q=1:s
        for j=0:s-1
            co(q,j+1)=j^(q-1)/factorial(q-1);
            rhs(q)=rhs(q)+j^(q)/factorial(q)*a(j+1);
        end
        rhs(q)=rhs(q)+s^(q)/factorial(q)*a(s+1);
    end
    b(s,1:s)=(co\rhs).';
end
AdamsBashforth_co=b

b=zeros(4,5);
for s=1:4
    a=zeros(s+1,1);
    co=zeros(s+1,s+1);
    rhs=zeros(s+1,1);
    a(s+1)=1;
    a(s)=-1;
    for q=1:s+1
        for j=0:s
            co(q,j+1)=j^(q-1)/factorial(q-1);
            rhs(q)=rhs(q)+j^(q)/factorial(q)*a(j+1);
        end
    end
    b(s,1:s+1)=(co\rhs).';
end
AdamsMoulton_co=b

for s=1:4
    co=zeros(s+2,s+2);
    rhs=zeros(s+2,1);
    for q=0:s
        for j=0:s
            co(q+1,j+2)=j^(q)/factorial(q);
        end
        if(q>=1)co(q+1)=-s^(q-1)/factorial(q-1);end
    end
    co(s+2,s+2)=1;
    rhs(s+2)=1;
    ab(s,1:s+2)=(co\rhs).';
end
Backwarddifferentiation_co=ab
    