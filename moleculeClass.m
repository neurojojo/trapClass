classdef moleculeClass < handle
    
    properties
        d % Diffusion coefficient  (round)
        x0 % Initial x position (integer)
        y0 % Initial y position (integer)
        
        x % Current x position
        y % Current y position
        
        % Constant properties %
        d_rand_mu = 5;
        d_rand_sigma = 3;
        frame_height = 1000;
        frame_width = 2000;
    end
    
    methods
        function obj = moleculeClass(varargin)
            if nargin==3
               [obj.d,obj.x0,obj.y0] = deal(varargin{1},varargin{2},varargin{3});
            else
               [obj.d,obj.x0,obj.y0] =  ...
                   deal( round( normrnd( obj.d_rand_mu,obj.d_rand_sigma ) ), ... % Random diffusion coeff
                     randi( obj.frame_width ),...
                     randi( obj.frame_height ) );
            end
            [obj.x,obj.y]=deal(obj.x0,obj.y0);
            
        end
    end
end

