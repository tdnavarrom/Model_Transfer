

clear all
clc

% convertir datos mensuales a diarios
A=xlsread('SOLIS.xlsx');
FINAL=xlsread('SOLIS.xlsx','ALEXA');

dia=1;
mes=1;
anio=2009;

d=datenum([anio,mes,dia]);
for i=1:size(A,1)
    fecha(i,:)=datevec(d);
    d=d+1;    
end

M1=[fecha(:,1:3) A(:,5)]

%M - your data matrix with 4 columns < Year Month Day data > - double type

[a,~,c] = unique(M1(:,1:2),'rows');
QM = ([a, accumarray(c,M1(:,4),[],@sum)]); %(m3/mes)


for s=1:size(FINAL,2)
    SM=FINAL(:,s);
    MM=dmatriz(SM,1,1);   
    
    for i=1:size(fecha,1)  
        SD(i,s)=MM(fecha(i,1)-anio+1,fecha(i,2))/(eomday(fecha(i,1),fecha(i,2)));       
    end
end






%% primero vamos a llenar con el ciclo anual
% 
% 
% for i=1:size(A,2)
%     DATOS=A(:,i);
%     M=dmatriz(DATOS,1,1);
%     p=nanmean(M);
%     
%     for l=1:size(M,1)
%         for m=1:size(M,2)
%             a=isnan(M(l,m));
%             if a==1
%                 M(l,m)=p(1,m);
%             end
%         end
%     end
%     d=M';
%     d=d(:);
%     FINAL(:,i)=d';    
% end
% 
% for s=1:size(FINAL,2)
%     SM=FINAL(:,s);
%     MM=dmatriz(SM,1,11);
%     
%     fecha=fecha1;
%     for i=1:(fecha2-fecha1+1)
%         f=datevec(fecha);
%         SD(i,s)=MM(f(1,1)-anio+1,f(1,2))*10^6/(eomday(f(1,1),f(1,2)));
%         fecha=fecha+1;
%     end
% end
%         
% % Convertir diario a Mensaul
