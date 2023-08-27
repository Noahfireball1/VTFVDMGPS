classdef PropagateGravity < handle
    %PROPAGATEENGINE Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = private)
        gravity;
    end

    methods
        function  [obj,gravityForces,gravityMoments] = PropagateGravity(states,mass,cg)

            obj.gravity = 9.8066;

            gravityForces = mass*[-obj.gravity*sin(states(11));obj.gravity*sin(states(10))*cos(states(11));obj.gravity*cos(states(10))*cos(states(11))];
            gravityMoments = cross(cg,gravityForces);
        end

    end
end

