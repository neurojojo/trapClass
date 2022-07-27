classdef analyzeTrapsClass < handle
    % This class takes tracksTableClass and segsTableClass objects as inputs
    %
    % showImmobileConfined
    % locateTraps (creates immobile_confined_mask)
    % locateTrappedSegs
    % locateTrappedTracks
    %
    
    properties
        metadata
        tracksTable
        segsTable
        % Tables created within the class
        trap_rois_tbl
        traps_segs_join_tbl
        traps_segs_tracks_tbl
        % Images created within the class
        immobile_confined_mask
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % CONSTRUCTOR FUNCTION %
        %%%%%%%%%%%%%%%%%%%%%%%%
        function obj = analyzeTrapsClass( tracksTable, segsTable )
            
            obj.metadata = tracksTable.metadata;
            obj.tracksTable = tracksTable.tracksTable;
            
            % In a preliminary step, we will clean up any NaNs in the data
            % before we assign this to the instantiation of the class
            segstbl = segsTable.segsTable;
            segstbl = segstbl(~isnan(segstbl.segType),:);
            [segstbl.xSeg, segstbl.ySeg] = deal( ...
                cellfun(@(x) inpaint_nans(x), segstbl.xSeg , 'UniformOutput', false),...
                cellfun(@(x) inpaint_nans(x), segstbl.ySeg , 'UniformOutput', false) );
            segstbl( or( segstbl.segType==0,segstbl.segType==1 ), : ).segType = ... 
                repmat(0,sum(or( segstbl.segType==0,segstbl.segType==1 )),1);
            segstbl( or( segstbl.segType==2,segstbl.segType==3 ), : ).segType = ...
                repmat(1,sum(or( segstbl.segType==2,segstbl.segType==3 )),1);
            obj.segsTable = segstbl;
            
            
        end

        function showImmobileConfined(obj)
            mytbl = obj.segsTable;
            % Looking only at immobile or confined segments
            i_or_c_segs = mytbl( or(mytbl.segType==1,mytbl.segType==2),: );
            % A figure showing the overlap of x and y locations
            figure
            ax = axes('nextplot','add');
            arrayfun(@(x) plot( ax,i_or_c_segs(x,:).xSeg{1},i_or_c_segs(x,:).ySeg{1}),...
                [1:size(i_or_c_segs,1)] )
        end
        
        % Use the next function to ask where traps might possibly be located
        % using the locations of immobilized particles
        function obj = locateTraps(obj)

            % (1) Creates a table of segments with sequences in the column
            % (2) takes the immobile and confined segments for every track
            % (3) shows every single track

            segment_xy_by_track = table();
            tmp_tbl = obj.segsTable;
            % Only take 0 or 1 valued segments (which are immobile or confined)
            tmp_tbl = tmp_tbl( tmp_tbl.segType==0, : );
            for i = unique( tmp_tbl.trackIdx )'

                % Build a track from the matching segType values
                %wesSegs{expIdx}.segsTable( wesSegs{expIdx}.segsTable.trackIdx==i, : ).segType;
                segment_xy_by_track = [segment_xy_by_track;...
                    table( i,...
                           {cell2mat(tmp_tbl( tmp_tbl.trackIdx==i, : ).xSeg')},...
                           {cell2mat(tmp_tbl( tmp_tbl.trackIdx==i, : ).ySeg')}, ...
                    'VariableNames', {'trackIdx','x','y'} ) ];
            end

            % Take the immobile and confined segments for every track (Create trap rois table)
            [x_,y_] = deal( segment_xy_by_track.x, segment_xy_by_track.y );

            mask_ = zeros(256,256);
            for i = 1:numel(x_)
                [x_tmp,y_tmp] = deal( round(inpaint_nans(round(x_{i}'))), round(inpaint_nans(round(y_{i}'))) );
                % This step adds a one to every location where 
                tmp = accumarray( [x_tmp,y_tmp], ones(numel(x_{i}),1), [256,256] );
                tmp = tmp>0;
                mask_ = mask_+tmp;
            end
            
            trap_rois = bwconncomp( mask_>1 );
            Ntraps = numel( trap_rois.PixelIdxList );
            %
            trap_rois_xy = cellfun( @(x) obj.ind2sub_array(trap_rois.ImageSize,x), trap_rois.PixelIdxList,...
            'UniformOutput', false);

            trap_rois_tbl = table();
            for i = 1:numel( trap_rois_xy )
                trap_rois_tbl = [ trap_rois_tbl;...
                    array2table( mean(trap_rois_xy{i},1) ) ]
            end
            trap_rois_tbl.trapIdx = [1:size(trap_rois_tbl,1)]';
            trap_rois_tbl.Properties.VariableNames = {'CentroidX','CentroidY','trapIdx'};
            trap_rois_tbl.CentroidIndex = rowfun(@(x,y) sub2ind([256,256],x,y), trap_rois_tbl(:,{'CentroidX','CentroidY'}) );
            
            obj.trap_rois_tbl = trap_rois_tbl;
            obj.immobile_confined_mask = mask_;
        end
        
        % Locate the indices of all of the segments that become trapped in the identified traps
        function obj = locateTrappedSegs(obj)
            
            minDistance = 5;
            
            tracks_tbl = obj.tracksTable;
            segs_tbl = obj.segsTable;
            trap_rois_tbl = obj.trap_rois_tbl;
            
            traps_tbl = table();
            for i = 1:size(segs_tbl,1)
                distances_to_traps = ...
                    pdist2( [trap_rois_tbl.CentroidX,trap_rois_tbl.CentroidY],...
                            [segs_tbl(i,:).xSeg{1}',segs_tbl(i,:).ySeg{1}'] ); 
                % Distance between all traps and all points along this line segment
                % For a line segment that is immobile or confined, comparing to all points along the line
                % is a bit overkill -- but it doesn't take too long and might work for some edge cases
                [min_distance_to_trap,trap_approached] = min(distances_to_traps,[],1);
                if min_distance_to_trap < minDistance
                    unique_traps_approached = unique(trap_approached)';
                    N_traps = numel( unique_traps_approached );
                    traps_tbl = [traps_tbl; table( ...
                        repmat( segs_tbl(i,:).segIdx,N_traps,1),...
                        unique_traps_approached )];
                end
            end
            traps_tbl.Properties.VariableNames = {'segIdx','trapIdx'};

            traps_segs_join = innerjoin( traps_tbl, segs_tbl, 'Key', 'segIdx' );
            
            obj.traps_segs_join_tbl = traps_segs_join;
            
        end
        
        function trapbias = computeTrappedSegStats(obj)
            
            trapstats = varfun(@(x) histcounts(x,[0:4]), obj.traps_segs_join_tbl,...
                'InputVariable', 'segType', 'GroupingVariable', 'trapIdx' );
            trapbias = ( sum(trapstats.Fun_segType(:,[1:2]),2) - ... % Sum of IC per trap
                sum(trapstats.Fun_segType(:,[3:4]),2) )./ ... % Subtract Sum of DIFF per trap
                ( sum(trapstats.Fun_segType(:,[1:2]),2) + ... % Sum of IC per trap
                sum(trapstats.Fun_segType(:,[3:4]),2) ); ... % Add Sum of DIFF per trap
                
%             
%             traps_segs_join = obj.traps_segs_join;
%             
%             traps_segs_join.segLifetime = arrayfun(@(x) numel(x{1}), traps_segs_join.xSeg)
%             segment_lifetime_tbl = varfun(@mean, traps_segs_join, 'InputVariables','segLifetime', 'GroupingVariables', 'trapIdx' )
% 
%             segstbl.segLifetime = arrayfun(@(x) numel(x{1}), segstbl.xSeg);
%             
%             % Percentiles for the lifetimes of particles that go into traps versus all segments
%             stats = [prctile( traps_segs_join.segLifetime, 25 ),...
%                 prctile( segstbl.segLifetime, 25 ),...
%                 size(traps_segs_join,1)/size(segstbl,1) ]
%             
%         end
%         
%         function obj = locateTrappedTracks(obj)
%             
%             tracks_tbl = obj.tracksTable;
%             traps_segs_join = obj.traps_segs_join_tbl;
%             trap_rois_tbl = obj.trap_rois_tbl;
%             
%             tracks_tbl.Properties.VariableNames{1}='trackIdx';
%             traps_tracks_join = innerjoin(traps_segs_join,tracks_tbl,'Key','trackIdx');
%             traps_tracks_join = innerjoin(traps_segs_join,trap_rois_tbl(:,[1:3]),'Key','trapIdx');
% 
%             % Compute the furthest that a track gets away from the center of a trap
%             traps_tracks_join.track_to_trap_maxd = rowfun( @(centroidX,centroidY,trackx,tracky) ...
%                 max(pdist2( [centroidX,centroidY],[trackx{1}',tracky{1}'])),...
%                 traps_tracks_join(:,{'CentroidX','CentroidY','xSeg','ySeg'}),'OutputFormat','uniform' );
% 
%             output = varfun( @(x) [min(x),max(x),range(x),numel(x)], traps_segs_join, ...
%                 'InputVariables', 'segStart', 'GroupingVariables', 'trapIdx' )
% 
%             traps_segs_tracks_tbl = innerjoin( traps_segs_join,...
%                 array2table( [output.trapIdx, output.Fun_segStart],...
%                 'VariableNames', {'trapIdx','TrapStart','TrapEnd','TrapLifetime','NsegsTrapped'} ),...
%                 'Keys', 'trapIdx' )
% 
%             trap_durations = unique( traps_segs_tracks_tbl(:,{'trapIdx','TrapStart','TrapEnd','TrapLifetime','NsegsTrapped'}) );
%             
%             obj.traps_segs_tracks_tbl = traps_segs_tracks_tbl;

        end
        
        function makeTrapFigures(obj)
        %
        % Emphasize variability among dim objects
        factor = -0.008;
        cmap = jet(1000);
        Y_map1 = 999*rescale( 1./(1+exp(factor*[1:1000])) );
        cmap1 = cmap(1+fix(Y_map1),:);
        cmap1(1,:) = [0,0,0];
        figure('color','w'); imagesc( obj.immobile_confined_mask );
        colormap(cmap1);
        colorbar();
            %figure(); ax = axes('nextplot','add');
            %trap_rois_xy = cellfun( @(x) ind2sub_array(trap_rois.ImageSize,x), trap_rois.PixelIdxList,...
            %    'UniformOutput', false);
%             
%         %
% %         
%             figure('color','w'); imagesc( mask_>1 ); colormap( hot(256) )
%             figure('color','w'); imagesc( mask_ ); colormap( hot(256) )
%             figure('color','w'); imagesc( mask_ ); colormap( hot(256) ); colorbar;
%             colormap([1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;jet(256)])
        
            % Show every single track
%             figure('color','w'); ax = axes('nextplot','add');rowfun( @(x,y,c) ...
%                 plot(ax,x{1},y{1},'color',colors(c+1,:)),...
%                 segstbl(:,{'xSeg','ySeg','segType'}) )
%         
%             figure(); ax=axes('nextplot','add');
%             colors = lines(100);
%             rowfun( @(trackIdx) plot(ax,tracks_tbl(tracks_tbl.id==trackIdx,:).x{1},...
%                 tracks_tbl(tracks_tbl.id==trackIdx,:).y{1}), traps_segs_join(:,{'trackIdx'}) );
%             rowfun( @(x,y) plot(ax,x{1},y{1},'color',[0,0,0,.1]), tracks_tbl(:,{'x','y'}) );
%             title(...
%                 sprintf(['Colored Tracks = Tracks which become stuck in a trap,',...
%                 '\nGrey tracks = Tracks that do not']) )
% figure; scatter( segment_lifetime_tbl.mean_segLifetime,segment_lifetime_tbl.GroupCount );
%             ylabel('# of segments trapped');xlabel('Lifetime of segment in trap')
% 
% figure();ax=axes('nextplot','add');
% rowfun( @(trapIdx,start_,end_,lt,NsegsTrapped)...
%     plot( [start_,end_],[trapIdx,trapIdx]), trap_durations )
% rowfun( @(trapIdx,start_,end_,lt,NsegsTrapped)...
%     plot( end_,trapIdx,'.','markersize',NsegsTrapped), trap_durations )
% ylim([0,Ntraps]);
% title('Number of segments trapped versus lifetime of the trap')
% 
% 
% figure();ax=axes('nextplot','add');
% scatter( trap_durations.NsegsTrapped, trap_durations.TrapLifetime )
% title('Number of segments trapped versus lifetime of the trap')
% 
% trap_maxradius = varfun(@max, tracks_traps_join, 'InputVariable', 'track_to_trap_maxd', 'GroupingVariable', 'trapIdx' );
% 
% figure; h1 = histogram( trap_maxradius.max_track_to_trap_maxd, 'BinMethod', 'integers' )
% set(gca,'XTick',[ceil(h1.BinLimits(1)):1:floor(h1.BinLimits(2))])
% %% How far apart in time are trapped tracks
        end
        
    end
    
    methods (Static)
        
        function xy = ind2sub_array( size_, indx_ )

            [x,y] = ind2sub( size_, indx_ );
            xy = [x,y];

        end
        
    end
end

