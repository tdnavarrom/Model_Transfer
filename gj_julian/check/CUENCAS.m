clear all
clc

A=xlsread('CUENCAS.xlsx','TEPETITLAN');
for i=1:size(A,2)
    M=dmatriz(A(:,i),1,1);
    media(:,i)=nanmean(M)';
end

m=1;
for i=1:size(A,1)
    a=isnan(A(i,1));
    if a==1
        A(i,:)=media(m,:);
        m=m+1;
    else
        m=1;        
    end
end



dia=1;
mes=1;
anio=2005;

d1=datenum([anio,mes,dia]);
d2=datenum([2020,12,31]);
d=d1;
for i=1:d2-d1+1
    fecha(i,:)=datevec(d);
    d=d+1;    
end

%M1=[fecha(:,1:3) A(:,5)]

%M - your data matrix with 4 columns < Year Month Day data > - double type

% [a,~,c] = unique(M1(:,1:2),'rows');
% QM = ([a, accumarray(c,M1(:,4),[],@sum)]); %(m3/mes)


for s=1:size(A,2)
    SM=A(:,s);
    MM=dmatriz(SM,1,1);   
    
    for i=1:size(fecha,1)  
        SD(i,s)=MM(fecha(i,1)-anio+1,fecha(i,2))/(eomday(fecha(i,1),fecha(i,2)))*10^6;       
    end
end