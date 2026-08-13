%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Crop Height and type     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hc,SAI,B,Ccrown,ZR95,Nreserve,Preserve,Kreserve,...
    Careserve,Mgreserve,Sireserve,...
    AgrHarNut,Rf] = CropHeightType(LAI,LAIdead,cc,ZR95,B,...
    Zs,CASE_ROOT,Ccrown,Crem,...
    Nreserve,Preserve,Kreserve,...
    Careserve,Mgreserve,Sireserve,...
    ManI,Mpar,Veg_Param_Dyn,OPT_SoilBiogeochemistry)
%%%% 
BRoot=B(:,cc,3); 
Btot = squeeze(sum(B, 3));
%Crop_type=Mpar.Crop_type; 
MHcrop =  Veg_Param_Dyn.MHcrop(cc); 
%%% INPUT
LAI = LAI + LAIdead/3; 
LAI_max=6.5; 
%%% CLM5.0 / AgroIBIS 
if LAI/(LAI_max-1)  < 1.0
    hc =  MHcrop*(LAI/(LAI_max-1))^2;
else
    hc = MHcrop; 
    %%%%%%
end
hc(hc<0.05)=0.05; 
SAI=max(0.15*LAI,0.001); 
%%%%%%%%%%%%%%%%%
AgrHarNut=[0 0 0 0 0 0]; 
if  ManI(cc)>0 %% Sowing
    Ccrown(cc)=Mpar(cc).Crop_crown(ManI(cc));
    
    if Ccrown(cc)==0
        B(:,cc,:)=zeros(1,8);
    end
    
    if sum(Ccrown) > 1 - Crem
        % free remaining land and redistribute it
        ix = setdiff(1:length(Ccrown), cc);
        
        for j = 1:length(ix)
           if Btot(ix(j)) == 0
               Ccrown(ix(j)) = 0;
           end
        end
        
    end
    
elseif ManI(cc) == -2 %%%%% Harvest
    
    if OPT_SoilBiogeochemistry==1
        AgrHarNut=[Nreserve(1,cc) Preserve(1,cc) Kreserve(1,cc) Careserve(1,cc) Mgreserve(1,cc) Sireserve(1,cc)];
        Nreserve(1,cc)=0;
        Preserve(1,cc)=0;
        Kreserve(1,cc)=0;

        Careserve(1,cc)=0;
        Mgreserve(1,cc)=0;
        Sireserve(1,cc)=0;
    end
    
    %%% Removing the crop
    % Ccrown(cc)=0;
    
end


if sum(Ccrown) < 1 - Crem
    
    ix = find(Btot == 0);
    Cfill = (1 - Crem) - sum(Ccrown);
    
    for j = 1:length(ix)
        Ccrown(ix(j)) = ( Cfill )/length(ix);
    end
    
end
        
%%%%%%%%%%%
%%%%% Potential update for rooting depth -- 
if  CASE_ROOT~= 1 
    disp('IN ORDER TO HAVE A VARIABLE ROOT DEPTH CASE_ROOT MUST BE 1')
    return
end
% %ZR=0.5*(2*BRoot)^r;
% ZR = 1283*BRoot^0.1713-1914; ZR(ZR<5)=5; %%[mm]
% ZR95(cc)=ZR;
% % %%%%
ZR50= NaN*ZR95;
ZRmax= NaN*ZR95;
% % %%%% Root depth Update
[Rf]=Root_Fraction_General(Zs,CASE_ROOT,ZR95,ZR50,0*ZR95,0*ZR50,ZRmax,0*ZRmax);
return
%%%%%%

