% Main Divide Processing Script

% This script outputs ALL ridgelines along DBPR that correspond to mapped 
% catchment outlets and their corresponding "mirrored" curvature values along
% DBPR. Divides used in the analysis can be seen in Fig. 1 in the paper where
% only divides that spanned from the outlets of catchments all the way to the 
% top of DBPR at the San Andreas Fault were used. 

% To run this script, TopoToolbox must be installed an on the file path
% (for large DEMs, processing may take a while. This script can be run on a
% server). 

% Files needed to run:
% A DEM (DBPR DEM is named '11_big.tif');
% A .mat file with the coordinates of catchment outlets (the provided file
%   includes outlet coordinates that are snapped to the stream network
%   produced here. ('Coordinatesfor1m.mat') 
% A csv file with the location of maximum ridgeline deviation (this is not
%   necessary for the analysis, but is a relict portion of the code).

% The output files are in .mat format. 

% Written by Will Struble, Univiersity of Houston, University of Arizona,
% 2024

addpath(genpath('topotoolbox-master')); % If running on a server.

savename = 'MirroredDividesOutput_1m_doubleHT_4smooth_untrimmed.mat';
maxdev = 1;  % 1 if getting Casym where ridgelines max offset, 0 if only whole ridgeline
    % NOTE! Whether maxdev==1 or 0, the code stil outputs results for the
    % whole ridgeline. 
savenamemax = 'MirroredDividesOutput_1m_doubleHT_4smooth_maxdev_untrimmed.mat';

% DEM = GRIDobj('C:\Users\wtstr\OneDrive\Documents\Projects\DragonsBack\11_big.tif');
DEM = GRIDobj('11_big.tif');

% load('DragonsBackFLOWobj.mat');
load('Coordinatesfor1m.mat');

x = xsnap;
y = ysnap;

criticalarea = 500;  % What is the channel head critical area you want to use?
curvsmooth = 4;

FD = FLOWobj(DEM,'preprocess','carve');
A = flowacc(FD).*DEM.cellsize^2;
S = STREAMobj(FD,A>criticalarea);

% st = readtable('outlets.csv');
% x = st.x;
% y = st.y;

ix = coord2ind(DEM,x,y);

if maxdev==1
    T = readtable('MaxDeviations_locations.csv');
    xneg = T.ridgexneg;
    yneg = T.ridgeyneg;
    xpos = T.ridgexpos;
    ypos = T.ridgeypos;
    ridgeid = T.ridge_id;

%     TrimT = [1;2;3;4;5;7;9;14;16;21;24;26;27;29;33;34;35;40;39;42;46;48;50;54;52;56;58;59;60;61;...
%     62;63;64;65;66;67;68;69;70;71;72;73]';

end

%% Identifying which catchments share a divide with each other

catches = 1:length(ix);

SE = strel('square',5);  % We'll dilate by 5 pixels. 
for i = 1:length(ix)-1
    L = drainagebasins(FD,ix(i));
    LDil = dilate(L,SE); % Dilate each basin by SE
    k = 1;
    Lshare=[];
    for j = i+1:length(ix)
        
        Lt = drainagebasins(FD,ix(j));
        LDilt = dilate(Lt,SE);
        
        Lsum = LDil+LDilt;

        if nnz(Lsum.Z==2)>0
            Lshare(k)=j;
            k = k+1;
        else
        end
        
    end

    DivShare{i} = Lshare;

    disp(['Through Catchment ',num2str(i),'.']);
end
DivShare{length(ix)}=[];

disp('Done checking for catchment neighbors.');

save('SharedCatchments_1m.mat','DivShare','-v7.3');    

load('C:\Users\wtstr\OneDrive\Documents\Projects\DragonsBack\SharedCatchments_1m.mat');
%%
S = STREAMobj(FD,A>criticalarea);
disp('Starting Mirroring')
k=1;
nn = 1;
for i = 1:length(DivShare)-1
    if ~isempty(DivShare{i})
        tt = DivShare{i};
        
        for j = 1:length(tt)
            St = modify(S,'upstreamto',ix([i,tt(j)]));
%             St = klargestconncomps(St,2); % You shouldn't need this...
            dx=1;
            paralleldx=0.43;  % Note that this value will depend on Lh below (it is a fraction of Lh)
            Lh=18;
            [OverlapDivix,OverlapDiv1ix,OverlapDiv2ix,z1t,z2t,zft,...
                IX_flowpath_1,IX_flowpath_2,Cht_1,Cht_2,Cht_full,DivTrim,DivTrim1,DivTrim2,rmv] = ...
            mirror_hilltops(DEM,FD,St,curvsmooth,dx,paralleldx,Lh);
                      
            z1{i,j} = z1t;
            z2{i,j} = z2t;
            zfull{i,j} = zft;
            flowpath1{i,j} = IX_flowpath_1;
            flowpath2{i,j} = IX_flowpath_2;
            Cht1{i,j} = Cht_1;
            Cht2{i,j} = Cht_2;
            Chtf{i,j} = Cht_full;
            OverlapDivides{i,j} = OverlapDivix;
            OverlapDivides1{i,j} = OverlapDiv1ix;
            OverlapDivides2{i,j} = OverlapDiv2ix;
            
            %if we want maxdev, & this loop corresponds to a ridge we
            %already want & that ridge is in our trimmed list, then... (This was not done in the paper)        
            if maxdev == 1 && ismember(k,ridgeid) 

                DivOff = GRIDobj(DEM);  % Define Max Offset GRIDobj
                div_ix_p = coord2ind(DEM,xpos(nn),ypos(nn)); % Positive offset
                div_ix_n = coord2ind(DEM,xneg(nn),yneg(nn)); % Negative offset
                nn=nn+1;
                
                DivOff.Z(div_ix_p) = 1;  % Coords of max offset are 1
                DivOff.Z(div_ix_n) = 1;

                SE = strel('square',30);  %Dilate them by 30 pixels (arbitrary for now)
                DivDilate = dilate(DivOff,SE);
                DivDilateT = find(DivDilate);  % Find which pixels are true.
 
                emptybefore(nn) = isempty(Cht_1);

                z1_maxdev{i,j} = z1t(:,ismember(OverlapDiv1ix,DivDilateT));
                z2_maxdev{i,j} = z2t(:,ismember(OverlapDiv1ix,DivDilateT));
                zfull_maxdev{i,j} = zft(:,ismember(OverlapDiv1ix,DivDilateT));
                flowpath1_maxdev{i,j} = IX_flowpath_1(ismember(OverlapDiv1ix,DivDilateT));
                flowpath2_maxdev{i,j} = IX_flowpath_2(ismember(OverlapDiv1ix,DivDilateT)); % defined by side 1
                Cht1_maxdev{i,j} = Cht_1(ismember(OverlapDiv1ix,DivDilateT));
                Cht2_maxdev{i,j} = Cht_2(ismember(OverlapDiv1ix,DivDilateT));
                Chtf_maxdev{i,j} = Cht_full(ismember(OverlapDiv1ix,DivDilateT));
                OverlapDivides_maxdev{i,j} = OverlapDivix(ismember(OverlapDivix,DivDilateT));
                OverlapDivides1_maxdev{i,j} = OverlapDiv1ix(ismember(OverlapDiv1ix,DivDilateT));
                OverlapDivides2_maxdev{i,j} = OverlapDiv2ix(ismember(OverlapDiv2ix,DivDilateT));
                
                emptyafter(nn) = isempty(Cht1_maxdev{i,j});

            elseif maxdev==1 && ~ismember(k,ridgeid)

%                 z1_maxdev{i,j} = [];
%                 z2_maxdev{i,j} = [];
%                 zfull_maxdev{i,j} = [];
%                 flowpath1_maxdev{i,j} = [];
%                 flowpath2_maxdev{i,j} = [];
%                 Cht1_maxdev{i,j} = [];
%                 Cht2_maxdev{i,j} = [];
%                 Chtf_maxdev{i,j} = [];
%                 OverlapDivides_maxdev{i,j} = [];
%                 OverlapDivides1_maxdev{i,j} = [];
%                 OverlapDivides2_maxdev{i,j} = [];

                z1_maxdev{i,j} = NaN;
                z2_maxdev{i,j} = NaN;
                zfull_maxdev{i,j} = NaN;
                flowpath1_maxdev{i,j} = NaN;
                flowpath2_maxdev{i,j} = NaN;
                Cht1_maxdev{i,j} = NaN;
                Cht2_maxdev{i,j} = NaN;
                Chtf_maxdev{i,j} = NaN;
                OverlapDivides_maxdev{i,j} = NaN;
                OverlapDivides1_maxdev{i,j} = NaN;
                OverlapDivides2_maxdev{i,j} = NaN;
       
            end 
            k=k+1;
        end        
    else
        z1{i,j} = [];
        z2{i,j} = [];
        zfull{i,j} =[];
        flowpath1{i,j} = [];
        flowpath2{i,j} = [];
        Cht1{i,j} = [];
        Cht2{i,j} = [];
        Chtf{i,j} = [];
        OverlapDivides{i,j} = [];
        OverlapDivides1{i,j} = [];
        OverlapDivides2{i,j} = [];
       
        if maxdev == 1
            z1_maxdev{i,j} = [];
            z2_maxdev{i,j} = [];
            zfull_maxdev{i,j} = [];
            flowpath1_maxdev{i,j} = [];
            flowpath2_maxdev{i,j} = [];
            Cht1_maxdev{i,j} = [];
            Cht2_maxdev{i,j} = [];
            Chtf_maxdev{i,j} = [];
            OverlapDivides_maxdev{i,j} = [];
            OverlapDivides1_maxdev{i,j} = [];
            OverlapDivides2_maxdev{i,j} = [];
        end
        
    end
    disp(['Done with catchment ',num2str(i),' of ',num2str(length(DivShare))]);
    save(savename,'DivShare','z1','z2','zfull','flowpath1','flowpath2',...
        'Cht1','Cht2','Chtf','OverlapDivides','OverlapDivides1','OverlapDivides2','-v7.3');
    
    if maxdev == 1
    save(savenamemax,'DivShare','z1_maxdev','z2_maxdev','zfull_maxdev',...
        'flowpath1_maxdev','flowpath2_maxdev','Cht1_maxdev','Cht2_maxdev',...
        'Chtf_maxdev','OverlapDivides_maxdev','OverlapDivides1_maxdev',...
        'OverlapDivides2_maxdev','emptybefore','emptyafter','-v7.3');
    end
    
end

save(savename,'DivShare','z1','z2','zfull','flowpath1','flowpath2','Cht1',...
    'Cht2','Chtf','OverlapDivides','OverlapDivides1','OverlapDivides2','-v7.3');

if maxdev==1
save(savenamemax,'DivShare','z1_maxdev','z2_maxdev','zfull_maxdev',...
    'flowpath1_maxdev','flowpath2_maxdev','Cht1_maxdev','Cht2_maxdev',...
    'Chtf_maxdev','OverlapDivides_maxdev','OverlapDivides1_maxdev',...
    'OverlapDivides2_maxdev','emptybefore','emptyafter','-v7.3');
end
