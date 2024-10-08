function [OverlapDivix,OverlapDiv1ix,OverlapDiv2ix,z1,z2,zfull,...
    IX_flowpath_1,IX_flowpath_2,Cht_1,Cht_2,Cht_full,DivTrim,DivTrim1,DivTrim2,rmv] = ...
    mirror_hilltops(DEM,FD,S,curvsmooth,dx,paralleldx,Lh)

%Extracts hillslope profiles on opposite sides of divides and calcualtes
%curvature asymmetry
%
% Input
%
%     DEM         digital elevation model (GRIDobj)
%     FD          flow direction raster (FLOWobj)
%     S           stream network (STREAMobj)
%                   THIS FUNCTION IS WRITTEN ASSUMING THERE ARE TWO (!!) DRAINAGE BASINS ONLY 
%     curvsmooth  Curvature smoothing scale (meters)
%     dx          How often do you want to sample the ridgeline pixels
%                   (dx=1 is sample every single pixel, dx=2 means skip every other)
%     paralleldx  To remove parallel flowpaths on mirrored hillslopes, what
%                   is the percentage that dictates how many pixels can be within the ridgeline window?
%     Lh          Length of hillslope profiles in PIXELS (not sure this is right, but: must be odd to ensure
%                   single pixel ridgelines)
% Output
%
%     OverlabDivix(B)     The shared divide (B has pixels labeled by each corresponding basin)
%          AS WRITTEN, THIS WILL OUTPUT THE LINEAR INDICES FOR THESE!! 
%     z1, z2            The mirrored hillslope profiles for each side
%     IX_flowpath_*     The flowpath for each side of the divide for each
%                           profile
%     Cht_*             The mirrored hilltop curvature value of each side.
%
%
% Author: Will Struble, University of Houston, University of Arizona

% First see if there are only two basins

SP = streampoi(S,'outlets','ix');
if length(SP)>2 || length(SP)<2
    error('Paired divide mapping is only written for two connected stream components.');
end

Lbig = drainagebasins(FD,S);  % The main drainage basins
[RL,DB] = FindRidgeLines(DEM,FD,S);  % Map the ridgelines and basins for S

RL_log = find(RL);       % Indices of ridgeline pixels
Lshare = GRIDobj(Lbig);  %L  
Lshareb = GRIDobj(Lbig); %L
Lnum = Lbig; %L

for i = 1:max(Lbig.Z,[],'all') %L
    
    Ltemp = find((Lbig.Z==i)); %L  % Find indices where each basin is
    
    if nnz(ismember(Ltemp,RL_log))~=0   % Does this catchment share pixels with the main divide?
        maindiv(i) = 1;         % Logical; not currently used later... may drop
        Lshare.Z(Ltemp) = 1;    % Logical basin map
        Lshareb.Z(Ltemp) = Lbig.Z(Ltemp);  % Label basin by basin number
    else
        maindiv(i) = 0;
        Lshare.Z(Ltemp) = 0;
        Lshareb.Z(Ltemp) = 0;
        Lnum.Z(Ltemp) = 0;
    end
    
end

% Divides currently clipped to stream network main boundaries. Now we want
% just the shared boundary    

SE = strel('square',5);  % We'll dilate by 5 pixels. 
Lshareb1 = Lshareb;
Lshareb1.Z(Lshareb1.Z==1) = 0; % Set basin 1 == 0 (counterintuitive, I know...)
Lshareb2 = Lshareb;
Lshareb2.Z(Lshareb2.Z==2) = 0; % Set basin 2 == 0 (... but should be arbitrary)

DivDil1 = dilate(Lshareb1,SE); % Dilate each basin by SE
DivDil2 = dilate(Lshareb2,SE); % 

Div = DivDil1 + DivDil2; % Add basins together
Div.Z = (Div.Z==3); % Where basins intersect is shared divide. 

Divix = find(Div); % Now, find the indices of this shared divide
ShareDivide = GRIDobj(Lbig); % Set up GRIDobj for shared divide
ShareDivideB = Lnum;      % Set up GRIDobj where basin numbers are saved

% Here we are just confirming that each basin does indeed
% share a boundary with the mapped divide 
for i = 1:max(Lnum.Z,[],'all')    
    Ltemp = find((Lnum.Z==i));
    if nnz(ismember(Ltemp,Divix))~=0
        ShareDivide.Z(Ltemp) = 1;       
    else
        ShareDivideB.Z(Ltemp) = 0;
    end
end

% Find where Div (the dilated divides) and RL (the actual ridgelines for
% the full catchments) overlap. (Basically removes extra pixels from "dilated"
% divide)
OverlapDivix = intersect(Divix,find((RL.Z==1)));
OverlapDiv = GRIDobj(RL);
OverlapDiv.Z(OverlapDivix) = 1;

% What basins are each divide pixel in?
OverlapDivB = OverlapDiv;
OverlapDivB.Z(OverlapDiv.Z==1 & DB.Z~=0) = DB.Z(OverlapDiv.Z==1 & DB.Z~=0); % Some might not actually correspond to DB basins
OverlapDivix = find((OverlapDivB.Z~=0));

% For some rare, weird basin pairs, they might touch when dilated, but don't
% actually share a divide (i.e., close, but not that close). 
if isempty(OverlapDivix)
    OverlapDivix=[];
    OverlapDiv1ix=[];
    OverlapDiv2ix=[];
    z1=[];
    z2=[];
    IX_flowpath_1=[];
    IX_flowpath_2=[];
    Cht_1=[];
    Cht_2=[];
    warning('Basins don''t actually share a divide. Returning empty vectors.');
    return
end

% Now, we only want divide pixels that share a pixel neighbor that drains
% to the other catchment (this should be most already mapped, but there could be some stragglers) 
for i = 1:length(OverlapDivix)
    % Put in "try" here in case pixels are on DEM edge?
    [row,col] = ind2sub(size(OverlapDivB.Z),OverlapDivix(i));
    ul=OverlapDivB.Z(row-1,col-1); us=OverlapDivB.Z(row-1,col);  ur=OverlapDivB.Z(row-1,col+1);
    sl=OverlapDivB.Z(row,col-1);                                 sr=OverlapDivB.Z(row,col+1);
    dl=OverlapDivB.Z(row+1,col-1); ds=OverlapDivB.Z(row+1,col);  dr=OverlapDivB.Z(row+1,col+1);
    
    neighbors = [ul,us,ur,sl,sr,dl,ds,dr];
    if OverlapDivB.Z(OverlapDivix(i)) == 1  % If the divide pixel is in catchment #1...
        if ~ismember(2,neighbors)   % And the neighboring pixel is NOT in catchment 2...
            remove(i) = 1;          % then remove
        else
            remove(i) = 0;
        end
            
    elseif OverlapDivB.Z(OverlapDivix(i)) == 2   % If the divide pixel is in catchment #2...
        if ~ismember(1,neighbors)       % And the neighboring pixel is NOT in catchment #1...
            remove(i) = 1;              % then remove
        else
            remove(i) = 0;
        end
    elseif OverlapDivB.Z(OverlapDivix(i)) == 0  % If the divide pixel is in no catchment
    end
end

remove=logical(remove); % Make logical
% Remove these isolated ridgelines from the map/list (should be very few, given the way we've defined the ridgelines)
OverlapDivB.Z(OverlapDivix(remove))=0;  
OverlapDivix(remove)=[];

% We seem to have undertrimmed OverlapDiv (even though we did OverlapDivB). So overwriting and defining. 
OverlapDiv = GRIDobj(OverlapDiv);
OverlapDiv.Z(OverlapDivix)=1;

% Here, we are denoting which basin each divide pixel resides in (after
% trimming)
OverlapDivB1 = GRIDobj(OverlapDivB); OverlapDivB2 = GRIDobj(OverlapDivB);
OverlapDivB1.Z(OverlapDivB.Z==1)=1;
OverlapDiv1ix = find(OverlapDivB1);
OverlapDivB2.Z(OverlapDivB.Z==2)=1;
OverlapDiv2ix = find(OverlapDivB2);             

% IX_flowpath_s = [];
% distance_s =[];
% figure
% imageschs(DEM)
% hold on
% for i = 1:length(OverlapDivix)
%     [IX_flowpath,distance] = flowpathextract(FD,OverlapDivix(i)); %,A,10000);
%     % trim to just headwaters, ~ 150 pixels ~~150 m (touch more, given
%     % diagonal movements 
%     IX_flowpath = IX_flowpath(1:150); distance = distance(1:150);
%     IX_flowpath_s = [IX_flowpath_s, IX_flowpath];
%     distance_s = [distance_s, distance];
%     [xpath,ypath] = ind2coord(DEM,IX_flowpath);
%     plot(xpath,ypath,'k-');
% end

% figure('units','normalized','outerposition',[0 0 1 1])
% imageschs(DEM,OverlapDivB,'colorbar',true);
% pause
% close
% return

%% Now, take paired divide and make paired profiles
% We will just move along one side of the divide. If it shares >1 pixel, we
% will just take a single value (the dataset will get too big otherwise!)

IX_flowpath_1 = []; IX_flowpath_2 = [];
distance_1 =[]; distance_2 =[];

%dx=8;
runlength = rem(nnz(OverlapDivB.Z==1),dx);   % Added this and also made multiple of dx below
k=1;

if dx>nnz(OverlapDivB.Z==1) % Just in case dx is bigger than the length of our divide, we still want 1 profile. 
    dxloop = nnz(OverlapDivB.Z==1);
    warning('Selected dx is greater than length of ridgeline. Using last ridgeline pixel.')
else
    dxloop = dx;
end

for i = 1:dxloop:nnz(OverlapDivB.Z==1) %-runlength  % For all pixels on the watershed #1 side of the divide...
    
    [row,col] = ind2sub(size(OverlapDivB.Z),OverlapDiv1ix(i));  %ind2sub(size(OverlapDivB.Z),OverlapDivix(i));
    ul=OverlapDivB.Z(row-1,col-1); us=OverlapDivB.Z(row-1,col);  ur=OverlapDivB.Z(row-1,col+1);
    sl=OverlapDivB.Z(row,col-1);                                 sr=OverlapDivB.Z(row,col+1);
    dl=OverlapDivB.Z(row+1,col-1); ds=OverlapDivB.Z(row+1,col);  dr=OverlapDivB.Z(row+1,col+1);
    
    neighbors = [ul,us,ur,sl,sr,dl,ds,dr];
    
    neighborshare=find((neighbors==2));   % Find a single neighbor pixel on the opposite side of the divide
    neighborshare = neighborshare(1);     % We used to just pick the first value...   
    neighborshare = neighborshare(randi(length(neighborshare))); % but I decided to randomize it.

    if neighborshare==1
        jr = row-1;
        jc = col-1;
    elseif neighborshare==2
        jr = row-1;
        jc = col;
    elseif neighborshare==3
        jr = row-1;
        jc = col+1;
    elseif neighborshare==4
        jr = row;
        jc = col-1;
    elseif neighborshare==5
        jr = row;
        jc = col+1;
    elseif neighborshare==6
        jr = row+1;
        jc = col-1;
    elseif neighborshare==7
        jr = row+1;
        jc = col;
    elseif neighborshare==8
        jr = row+1;
        jc = col+1;
    end    
           
    % Flowpath of side 1
    [IX_flowpath,distance] = flowpathextract(FD,OverlapDiv1ix(i));

        %    
    if length(IX_flowpath)<Lh
        IX_flowpath=ones(Lh,1)*NaN;
        distance = ones(Lh,1)*NaN;
    else
        IX_flowpath = IX_flowpath(1:Lh);
        distance = distance(1:Lh); 
    end

    HillslopeBuffer = Lh-0.9*curvsmooth; % We don't want to include short hillslope profiles that are close to the channel. 

    % Removing hillslope profiles that overlap the channels (i.e., concave) within curvature window
    if nnz(ismember(IX_flowpath,S.IXgrid)) >= HillslopeBuffer 
        HS_Chan_over(i) = nnz(ismember(IX_flowpath,S.IXgrid));
        test(i) = 1;
        IX_flowpath=ones(Lh,1)*NaN;
        distance = ones(Lh,1)*NaN;
    else
        HS_Chan_over(i) = nnz(ismember(IX_flowpath,S.IXgrid));
        test(i) = 0;
    end

    IX_flowpath_1 = [IX_flowpath_1, IX_flowpath];
    distance_1 = [distance_1, distance];
    
    % Flowpath of side 2
    side2ix = sub2ind(size(DEM.Z),jr,jc);
    [IX_flowpath,distance] = flowpathextract(FD,side2ix);
    if length(IX_flowpath)<Lh || nnz(isnan(IX_flowpath_1(:,i)))>0  % changed an i to k
        IX_flowpath=ones(Lh,1)*NaN;
        distance = ones(Lh,1)*NaN;
        IX_flowpath_1(:,i) = NaN; % In case side 1 wasn't NaN, but side 2 is.
        distance_1(:,i) = NaN;
    else
        IX_flowpath = IX_flowpath(1:Lh);
        distance = distance(1:Lh); % We'll check Lh pixels downstream
    end
 
 % Removing hillslope profiles that overlap the channels (i.e., concave)
 %      within curvature window on side 2
    if nnz(ismember(IX_flowpath,S.IXgrid)) >= HillslopeBuffer 
        HS_Chan_over(i) = nnz(ismember(IX_flowpath,S.IXgrid));
        test2(i) = 1;
        IX_flowpath=ones(Lh,1)*NaN;
        distance = ones(Lh,1)*NaN;
        
        % Then also need to set side 1 to NaN
        IX_flowpath_1(:,i) = ones(Lh,1)*NaN;
        distance_1(:,i) = ones(Lh,1)*NaN;

    else
        HS_Chan_over(i) = nnz(ismember(IX_flowpath,S.IXgrid));
        test2(i) = 0;
    end

    IX_flowpath_2 = [IX_flowpath_2, IX_flowpath];
    distance_2 = [distance_2, distance];
    
    k=k+1;
end

% test=test+test2;
% for i = 1:length(test)
%     if test(i) > 1
%         test(i) = 1;
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%
% figure
% imageschs(DEM)
% hold on
% plot(S,'w');
% for i = 1:length(test)
%     
%     if test(i) == 0
%         [x1,y1] = ind2coord(DEM,IX_flowpath_1(:,i));
%         [x2,y2] = ind2coord(DEM,IX_flowpath_2(:,i));
%     
%         plot(x1,y1,'k-');
%         plot(x2,y2,'k-');       
%     end
% end
% 
% for i = 1:length(test)
%     
%     if test(i) == 1
%         tt = IX_flowpath_1(:,i);
%         tt1=tt(~ismember(tt,S.IXgrid));
%         [xt1,yt1] = ind2coord(DEM,tt1);
% 
%         tt = IX_flowpath_2(:,i);
%         tt2=tt(~ismember(tt,S.IXgrid));
%         [xt2,yt2] = ind2coord(DEM,tt2);
% 
%         [x1,y1] = ind2coord(DEM,IX_flowpath_1(:,i));
%         [x2,y2] = ind2coord(DEM,IX_flowpath_2(:,i));
%         plot(x1,y1,'m-');
%         plot(x2,y2,'m-');
% 
%         plot(xt1,yt1,'y-');
%         plot(xt2,yt2,'y-');
% 
%     end
% end
%%%%%%%%%%%%%%%%%%%%

clear distance

%% Dilate main divide to make sure that we aren't including "parallel" hillslopes
SE = strel('square',5);
DivDilate = dilate(OverlapDiv,SE); 
DivDilateT = find(DivDilate);

DivTrim = OverlapDivix;
DivTrim1 = OverlapDiv1ix;
DivTrim2 = OverlapDiv2ix;

for i = 1:size(IX_flowpath_2,2)   %length(IX_flowpath_2)  
    if nnz(ismember(IX_flowpath_2(:,i), DivDilateT))/length(IX_flowpath_2(:,i))>paralleldx ||...
            nnz(ismember(IX_flowpath_1(:,i), DivDilateT))/length(IX_flowpath_1(:,i))>paralleldx ||...
            nnz(isnan(IX_flowpath_1(:,i)))>0 || nnz(isnan(IX_flowpath_2(:,i)))>0
        

% Commented all this out so that our curvature vectors are the same length as the flowpaths, etc.         
%         [~,tt]=ismember(IX_flowpath_1(1,i),DivTrim);
%         if tt~=0
%             DivTrim(tt) = [];
%         end
% 
%         [~,tt]=ismember(IX_flowpath_2(1,i),DivTrim);
%         if tt~=0
%             DivTrim(tt) = [];
%         end
%         
%         [~,tt]=ismember(IX_flowpath_1(1,i),DivTrim1);
%         if tt~=0
%             DivTrim1(tt) = [];
%         end
% 
%         [~,tt]=ismember(IX_flowpath_2(1,i),DivTrim2);
%         if tt~=0
%             DivTrim2(tt) = [];
%         end
        
        rmv(i) = 1;
    else
        rmv(i)=0;
    end
end

% rmvix1 = IX_flowpath_1(:,rmv==1); rmvix1 = rmvix1(1,:);
% rmvix2 = IX_flowpath_2(:,rmv==1); rmvix2 = rmvix2(1,:);
% IX_flowpath_1(:,rmv==1)=[];
% IX_flowpath_2(:,rmv==1)=[];

IX_flowpath_1(:,rmv==1)=NaN;
IX_flowpath_2(:,rmv==1)=NaN;

if isempty(IX_flowpath_1)
    z1=[];
    z2=[];
    IX_flowpath_1=[];
    IX_flowpath_2=[];
    Cht_1=[];
    Cht_2=[];
    warning(['Empty hillslope flowpaths, likely because dx is too large.',...
        ' Returning empty vectors.']);
    return
end

%% Get elevations, make profiles, split and mirror
z1=[]; z2=[]; xprofile1=[]; yprofile1=[]; xprofile2=[]; yprofile2=[];
% k=1; % b/c IX_flowpath vectors are now shorter than rmv (but will save by i)
for i = 1:length(rmv)  %size(IX_flowpath_2,2)   %length(IX_flowpath_2)
    
    if rmv(i)==1  % Discovered that if we want to trim out hilltop pixels later, helps to do this.
        xprofile1(:,i) = ones(length(IX_flowpath_1(:,1)),1)*NaN;
        yprofile1(:,i) = ones(length(IX_flowpath_1(:,1)),1)*NaN;
        xprofile2(:,i) = ones(length(IX_flowpath_1(:,1)),1)*NaN;
        yprofile2(:,i) = ones(length(IX_flowpath_1(:,1)),1)*NaN;
        z1(:,i) = ones(length(IX_flowpath_1(:,1)),1)*NaN;
        z2(:,i) = ones(length(IX_flowpath_1(:,1)),1)*NaN;
    else
        [xprofile1(:,i),yprofile1(:,i)] = ind2coord(DEM,IX_flowpath_1(:,i));
        [xprofile2(:,i),yprofile2(:,i)] = ind2coord(DEM,IX_flowpath_2(:,i));
    
        z1(:,i) = double(interp(DEM,xprofile1(:,i),yprofile1(:,i)));
        z2(:,i) = double(interp(DEM,xprofile2(:,i),yprofile2(:,i)));
        k=k+1;
    end
    
end

xprofilefull = [flip(xprofile1); xprofile2];
yprofilefull = [flip(yprofile1); yprofile2];

zfull = [flip(z1); z2];

% disp('Waiting!')
% pause 

z1f = flipud(z1);
% z1f(end,:) = [];  % We don't want the hilltop pixels repeated at the new mirrored hilltops
%                                 %, so delete one of them (still in z1).
z1f(1,:)=[];
z1 = [z1f; z1];  % This is now the MIRRORED profile for side 1.

z2f = flipud(z2);
% z2f(end,:) = [];  % Same as above, but z2 (still in z2).
z2f(1,:)=[];
z2 = [z2f; z2];   % This is now the MIRRORED profile for side 2.


distancefull = [(-1)*flipud(distance_1)-DEM.cellsize; distance_2];

distance_1 = [(-1)*flipud(distance_1)-DEM.cellsize; distance_1];
distance_2 = [(-1)*flipud(distance_2)-DEM.cellsize; distance_2];

%% Now, extend into 2D and calculate curvature!
z1_2D = {};
z2_2D = {};
zfull_2D = {};

a = curvsmooth/(sqrt(2)*pi*DEM.cellsize);
Cht_1 = []; Cht_2 = []; Cht_full = [];

for i = 1:size(z1,2)   
    if isnan(z1(:,i))
        Cht_1(i) = NaN;
        Cht_2(i) = NaN;
        Cht_full(i) = NaN;
    else
        z1_2D{i} = repmat(z1(:,i),1,101);
        z2_2D{i} = repmat(z2(:,i),1,101);
        zfull_2D{i} = repmat(zfull(:,i),1,101);
    
        [C1] = conv2_mexh_curv(z1_2D{i},a,DEM.cellsize);
        [C2] = conv2_mexh_curv(z2_2D{i},a,DEM.cellsize);
        [C3] = conv2_mexh_curv(zfull_2D{i},a,DEM.cellsize);
    
        % Since mirrored, values are (mostly, except near edges) constant along the hilltop,
        % so taking the middle value at the top of the hilltop
        Cht_1(i) = C1(round(size(C1,1)/2),50); % Since dimensions are odd, middle pixel is length/2 + 0.5
        Cht_2(i) = C2(round(size(C2,1)/2),50); 
        Cht_full(i) = C3(round(size(C3,1)/2),50); 
    end
    
end

%%


    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RidgeLines,DB]=FindRidgeLines(DEM,FD,S,varargin)
	% Helper function to find ridge line locations, defined as top of
	% a drainage basin (i.e. ridge lines are two grids thick, with true divide between them)
    % Modified version of Adam Forte's FindRidgeLines
    
	p=inputParser;
	p.FunctionName = 'FindRidgeLines';
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));

	addParamValue(p,'minimum_order',3,@(x) isscalar(x) && isnumeric(x));
    addParamValue(p,'river_mouths',[1 1],@(x) isnumeric(x) && size(x,2)==2);

	parse(p,FD,S,varargin{:});
	FD=p.Results.FD;
	S=p.Results.S;
    rm=p.Results.river_mouths;

	so=p.Results.minimum_order;

	
    COix=streampoi(S,'outlets','ix');
    DB=drainagebasins(FD,COix);
	
	[db,X,Y]=GRIDobj2mat(DB);
	db=double(db);

	rl_l=diff(db,1,2);
	rl_d=diff(db,1,1);

	rl_l(rl_l~=0)=1;
	rl_d(rl_d~=0)=1;

	size_left=size(rl_l);
	size_down=size(rl_d);

	rl_u=vertcat(zeros(1,size_down(:,2)),rl_d);
	rl_r=horzcat(zeros(size_left(:,1),1),rl_l);

	rl_d=vertcat(rl_d,zeros(1,size_down(:,2)));
	rl_l=horzcat(rl_l,zeros(size_left(:,1),1));

	rdgs=rl_d+rl_u+rl_r+rl_l;
	rdgs(rdgs~=0)=1;
% 	RidgeLines=GRIDobj(X,Y,rdgs);
    RidgeLines=GRIDobj(DEM);
    RidgeLines.Z=rdgs;
end

end