classdef simulationClass < handle
    
    % Demo
%     for i = 1:10; molecules{i} = moleculeClass(); end
%     for i = 1:10; traps{i} = trapClass(); end
%     simulation = simulationClass();
%     simulation.addObjects( molecules );
%     simulation.addObjects( traps );
    
    properties
        molecules % Either one moleculeClass object or a cell array of moleculeClass objects
        Nmolecules = 0;
        
        traps % Either one trapClass object or a cell array of trapClass objects
        Ntraps = 0;
        
        trapRadius = 5; % Distance around a trap that will slow down a molecule 
        
        % Default frame size and number
        frame_height = 1000;
        frame_width = 2000;
        frames = 1000;
    end
    
    methods
        function obj = simulationClass(varargin)
            if nargin>0
               [obj.frame_height,obj.frame_width] = deal(varargin{1},varargin{2}); 
            end
        end
        
        function obj = addObjects(obj,input)
            if iscell( input ) % Add every molecule
                
                for i = 1:numel(input)
                    if strcmp( class( input{i} ) , 'moleculeClass' )
                        obj.Nmolecules = obj.Nmolecules + 1;
                        obj.molecules{ obj.Nmolecules } = input{i}; 
                    end
                    
                    if strcmp( class( input{i} ) , 'trapClass' )
                        obj.Ntraps = obj.Ntraps + 1;
                        obj.traps{ obj.Ntraps } = input{i}; 
                    end
                end
                
            else % Add a single molecule
                if strcmp( class( input ) , 'trapClass' )
                    obj.Ntraps = obj.Ntraps + 1;
                    obj.traps{ obj.Ntraps } = input;
                end
                
                if strcmp( class( input ) , 'moleculeClass' )
                    obj.Nmolecules = obj.Nmolecules + 1;
                    obj.traps{ obj.Nmolecules } = input;
                end
            end
        end
        
        function outputSimulation(obj)
           
            % Extract trap x,y locations
            [trap_x,trap_y] = deal( cellfun( @(x) x.x0 , obj.traps )',...
                                    cellfun( @(x) x.y0 , obj.traps )' );
            
            for i = 1 : obj.frames
                for j = 1: obj.Nmolecules
                    degree1 = 2*pi()*rand();
                    d = obj.molecules{j}.d;
                    dist_to_traps = pdist2( [obj.molecules{j}.x(i),obj.molecules{j}.y(i)],...
                                            [trap_x,trap_y] );
                    
                    % Check if the molecules current position is near enough to a trap
                    if dist_to_traps > obj.trapRadius % Trap radius
                        obj.molecules{j}.x = ...
                            [ obj.molecules{j}.x,...
                            obj.molecules{j}.x(i) + d*cos(degree1) ];
                        obj.molecules{j}.y = ...
                            [ obj.molecules{j}.y,...
                            obj.molecules{j}.y(i) + d*sin(degree1) ];
                    else
                        obj.molecules{j}.x = ...
                            [ obj.molecules{j}.x, ...
                            obj.molecules{j}.x(i) + (rand()>0.9)*5 + (d/10)*cos(degree1) ];
                        obj.molecules{j}.y = ...
                            [ obj.molecules{j}.y, ...
                            obj.molecules{j}.y(i) + (rand()>0.9)*5 + (d/10)*sin(degree1) ];
                    end
                    
                end % Loop over obj.molecules
            end % Loop over obj.frames
            
        end
        
        function watchSimulation(obj,varargin)
            
            1
            
            f = figure('color','k','position',[680   558   1120   420]);
            ax = axes('color','k','nextplot','add',...
                'xcolor','k','ycolor','k','xlim',[0,obj.frame_width],'ylim',[0,obj.frame_height]);
            
            for i = 1:obj.frames
                arrayfun(@(Nmolecules) plot( ax, obj.molecules{Nmolecules}.x(1:i),...
                    obj.molecules{Nmolecules}.y(1:i) ), [1:obj.Nmolecules] );
                pause(0.1);
                cla;
            end
            
        end
        
    end
end

