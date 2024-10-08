%% Processing Gilbert Metrics at Dragon's Back

savename = 'DragonsBackGilbertMetrics.mat';
DEM = GRIDobj('C:\Users\wtstr\OneDrive\Documents\Projects\DragonsBack\11_big.tif');
DEMf = fillsinks(DEM);
G=arcslope(DEM);

% load('DragonsBackFLOWobj.mat');
load('C:\Users\wtstr\OneDrive\Documents\Projects\DragonsBack\Coordinatesfor1m.mat');

x = xsnap;
y = ysnap;

criticalarea = 2000;  % What is the channel head critical area you want to use?
curvsmooth = 4;

FD = FLOWobj(DEM,'preprocess','carve');
A = flowacc(FD).*DEM.cellsize^2;
S = STREAMobj(FD,A>criticalarea);

% st = readtable('outlets.csv');
% x = st.x;
% y = st.y;

% Have discovered that catchments 56 and 61 aren't correctly snapped
% (HOW??), so we are going to go ahead and drop them. 
% x([56,61])=[];
% y([56,61])=[];

ix = coord2ind(DEM,x,y);

St = modify(S,'upstreamto',ix);

%% Gilbert metrics, if you only want the values per basin

% % First, recognize that for some of our "basins" of interests, our
% % reference drainage area might correspond to more than one point (b/c a
% % drainage network branches). We are going to be taking an average.
% L = drainagebasins(FD,ix);
% Relief_GO = GRIDobj(L);
% Grad_GO = GRIDobj(L);
% Elev_GO = GRIDobj(L);
% 
% for i = 1:length(ix) 
% 
%     LP = drainagebasins(FD,ix(i));
% 
%     Sm = modify(S,'upstreamto',ix(i));
%     Heads = streampoi(Sm,'channelheads','ix');
%     
%     RelfL=[];
%     GradL=[];
%     ElevL=[];
%     qq=1;
%     
%     for j = 1:length(Heads)
%         LL = drainagebasins(FD,Heads(j));
%         RelfL(j) = max(DEM.Z(LL.Z==1))-DEM.Z(Heads(j));
%         GradL(qq:qq+nnz(LL.Z)-1) = G.Z(LL.Z==1)';
%         ElevL(j) = DEM.Z(Heads(j))';
% %         GradL{j} = G.Z(LL.Z==1); 
% %         ElevL{j} = DEM.Z(LL.Z==1);
%         qq = qq+nnz(LL.Z);
%     end
% 
%     Relf(i) = mean(RelfL);
%     Relfmed(i) = median(RelfL);
%     Relfst(i) = std(RelfL);
%     Relf_num(i) = length(RelfL);
% 
%     Grad(i) = mean(GradL);
%     Gradmed(i) = median(GradL);
%     Gradst(i) = std(GradL);
%     Grad_num(i) = length(GradL);
% 
%     Elev(i) = mean(ElevL);
%     Elevmed(i) = median(ElevL);
%     Elevst(i) = std(ElevL);
%     Elev_num(i) = length(ElevL);
% 
%     Relief_GO.Z(LP.Z==1) = Relfmed(i);
%     Grad_GO.Z(LP.Z==1) = Gradmed(i);
%     Elev_GO.Z(LP.Z==1) = Elevmed(i);
% 
% end
%  
% DB=table;
% DB.x = x';
% DB.y = y';
% DB.Median_Relief = Relfmed';
% DB.Relief_st_err = Relfst'./sqrt(Relf_num)';
% DB.Median_Gradient = Gradmed';
% DB.Gradient_st_err = Gradst'./sqrt(Grad_num)';
% DB.Median_Elevation = Elevmed';
% DB.Elevation_st_err = Elevst'./sqrt(Elev_num)';
% DB.NumberOfHeadsAtRefArea = Relf_num';
% DB.NumberofPixelsForGradient = Grad_num';
% writetable(DB,'DB_Gilbert.csv');
% 
% Relief_GO.Z(Relief_GO.Z==0) = NaN;
% Grad_GO.Z(Grad_GO.Z==0) = NaN;
% Elev_GO.Z(Elev_GO.Z==0) = NaN;
% % 
% % figure
% % imageschs(DEM,Relief_GO);
% % 
% % figure
% % imageschs(DEM,Grad_GO);
% % 
% % figure
% % imageschs(DEM,Elev_GO,'caxis',[600,680]);


%% Gilbert metrics, if you want the values assigned to a ridgeline (True Gilbert Metrics)
load('C:\Users\wtstr\OneDrive\Documents\Projects\DragonsBack\SharedCatchments_1m.mat');

k=1;
nn = 1;
for i = 1:length(DivShare)-1
    if ~isempty(DivShare{i})
        tt = DivShare{i};
        
        for j = 1:length(tt)
            St1 = modify(S,'upstreamto',ix(i));   % Side 1
            St2 = modify(S,'upstreamto',ix(tt(j))); % Side 2
            LP1 = drainagebasins(FD,ix(i));   % Basin 1
            LP2 = drainagebasins(FD,ix(tt(j)));  % Basin 2

            Heads1 = streampoi(St1,'channelheads','ix');
            Heads2 = streampoi(St2,'channelheads','ix');
           
            RelfL1=[]; RelfL2=[];
            GradL1=[]; GradL2=[];
            ElevL1=[]; ElevL2=[];
            qq=1;
            
            % ksn
            ksnarray1=ksn(St1,DEMf,flowacc(FD),0.68,10);
            ksnarray2=ksn(St2,DEMf,flowacc(FD),0.68,10);
    
            d1 = getnal(St1,LP1);
            d2 = getnal(St2,LP2);

            if ~isempty(ksnarray1)
                ks1 = accumarray(d1,ksnarray1,[],@mean);
                ksmed1 = accumarray(d1,ksnarray1,[],@median);
                ks_st1 = accumarray(d1,ksnarray1,[],@std);
                ks_num1 = length(ksnarray1);
            else
                ks1=NaN;
                ksmed1=NaN;
                ks_st1=NaN;
                ks_num1 = NaN;
            end

            if ~isempty(ksnarray2)
                ks2 = accumarray(d2,ksnarray2,[],@mean);
                ksmed2 = accumarray(d2,ksnarray2,[],@median);
                ks_st2 = accumarray(d2,ksnarray2,[],@std);
                ks_num2 = length(ksnarray2);
            else
                ks2=NaN;
                ksmed2=NaN;
                ks_st2=NaN;
                ks_num2 = NaN;
            end

            ksnGilb_med{i,j} = ks1-ks2;
            ksnGilb_std = sqrt(ks_st1.^2 + ks_st2.^2);
            ksnGilb_st_err{i,j} = ksnGilb_std./sqrt(ks_num1+ks_num2);

            % Gilberts
            for pp = 1:length(Heads1)
                LL = drainagebasins(FD,Heads1(pp));
                RelfL1(pp) = max(DEM.Z(LL.Z==1))-DEM.Z(Heads1(pp));
                GradL1(qq:qq+nnz(LL.Z)-1) = G.Z(LL.Z==1)';
                ElevL1(pp) = DEM.Z(Heads1(pp))';
                qq = qq+nnz(LL.Z);
            end

            qq=1;
            for pp = 1:length(Heads2)
                LL = drainagebasins(FD,Heads2(pp));
                RelfL2(pp) = max(DEM.Z(LL.Z==1))-DEM.Z(Heads2(pp));
                GradL2(qq:qq+nnz(LL.Z)-1) = G.Z(LL.Z==1)';
                ElevL2(pp) = DEM.Z(Heads2(pp))';
                qq = qq+nnz(LL.Z);
            end

            RelfGilb_med{i,j} = median(RelfL1)-median(RelfL2);
            RelfGilb_std = sqrt(std(RelfL1).^2 + std(RelfL2).^2);
            RelfGilb_st_err{i,j} = RelfGilb_std./sqrt(length(RelfL1)+length(RelfL2));

            GradGilb_med{i,j} = median(GradL1)-median(GradL2);
            GradGilb_std = sqrt(std(GradL1).^2 + std(GradL2).^2);
            GradGilb_st_err{i,j} = GradGilb_std./sqrt(length(GradL1)+length(GradL2));

            ElevGilb_med{i,j} = median(ElevL1)-median(ElevL2);
            ElevGilb_std = sqrt(std(ElevL1).^2 + std(ElevL2).^2);
            ElevGilb_st_err{i,j} = ElevGilb_std./sqrt(length(ElevL1)+length(ElevL2));

            
            k=k+1;
        end        
    else
            RelfGilb_med{i,j} = [];
            GradGilb_med{i,j} = [];
            ElevGilb_med{i,j} = [];
            RelfGilb_st_err{i,j} = [];
            GradGilb_st_err{i,j} = [];
            ElevGilb_st_err{i,j} = [];

            ksnGilb_med{i,j} = [];
            ksnGilb_st_err{i,j} = [];
              
    end
    disp(['Done with catchment ',num2str(i),' of ',num2str(length(DivShare))]);
end

save(savename,'RelfGilb_med','RelfGilb_st_err','GradGilb_med','GradGilb_st_err',...
     'ElevGilb_med','ElevGilb_st_err','-v7.3');

