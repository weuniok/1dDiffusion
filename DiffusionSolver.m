classdef DiffusionSolver < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    % x range - 0:1

    properties
        numVolumes
        numTimes 

        faceArea
        lambda
    end

    properties
        xVector
        dx
        U
        dU
        diffusionThroughFace
        dt
        currentTime
        numFaces

        V
        alpha
        A

        DirichletBCs = []

        isSCSet = false
    end

    methods
        function obj = DiffusionSolver(numVolumes, numTimes, lambda)
            obj.numVolumes = numVolumes;
            obj.numTimes = numTimes;
            obj.numFaces = obj.numVolumes+1;
            obj.lambda = lambda;

            obj.dx = 1/(obj.numVolumes);
            obj.xVector = 0+obj.dx/2:obj.dx:1-obj.dx/2;
            obj.dt  = obj.dx^2/2/max(obj.lambda);

            % TODO: volumes of different face areas 
            obj.faceArea = 1;
            obj.V = obj.faceArea*obj.dx;

            % Array size
            obj.U = zeros(obj.numVolumes, obj.numTimes);
            obj.dU = zeros(obj.numVolumes, obj.numTimes);
            obj.currentTime=zeros(1, obj.numTimes);
        end

        function solve(obj)
            if obj.isSCSet == false
                obj.setStartingCondition()
            end
            obj.calculateWallDiffusion()
            obj.calculateSolution()
        end

        function setStartingCondition(obj, SC)
            % Applies function SC form 0 to 1 as starting U
            % default SC - standard deviation centered on 0.5
            if nargin < 2
                mu = 0.5;
                sigma = 0.05;

                obj.U = zeros(obj.numVolumes, obj.numTimes);

                for i=1:obj.numVolumes
                    obj.U(i,1) = ...
                        exp(-(obj.xVector(i)-mu)^2/(2*sigma^2)) / sqrt(2*pi*sigma^2);
                end
            else
                for i=1:obj.numVolumes
                    obj.U(i,1) = SC(obj.xVector(i));
                end
            end
            obj.isSCSet = true;
        end

        function calculateWallDiffusion(obj)
            lambdaR = 2./([1./obj.lambda, 0] + [0, 1./obj.lambda]);
            lambdaR = lambdaR(2:end-1);

            obj.alpha = ones(1, obj.numVolumes-1) .* lambdaR * obj.faceArea / obj.dx;
            obj.A = diag(obj.alpha, -1) + diag(obj.alpha, 1) + ...
                - diag([obj.alpha, 0]) - diag([0, obj.alpha]);
        end

        function calculateSolution(obj)
            discretizationFactor = obj.dt / obj.V;

            for timeStep=1:obj.numTimes-1
                obj.currentTime(timeStep+1)=obj.currentTime(timeStep)+obj.dt;

                obj.applyBCs(timeStep);
                obj.dU(:, timeStep) = obj.dU(:, timeStep) + obj.A*obj.U(:,timeStep);
                obj.U(:,timeStep+1) = obj.dU(:, timeStep) * ...
                    discretizationFactor + obj.U(:, timeStep);
            end
        end

        function applyBCs(obj, timeStep)
            for i = 1:length(obj.DirichletBCs)
                index = obj.DirichletBCs(i).key;
                value = obj.DirichletBCs(i).value;
                obj.U(index, timeStep) = value;
            end
        end

        function setDirichletBC(obj, index, value)
            % Sets permanent U[index] value
            condition.key = index;
            condition.value = value;
            obj.DirichletBCs = [obj.DirichletBCs, condition];
        end

        function setSinglePointNeumannBC(obj, index, value, firstTime, lastTime)
            % TODO: make actual du/dt = value instead of = value + ...
            obj.dU(index, firstTime:lastTime) = value;
        end

        function plot(obj, timeStep)
            % Plot a single frame
            % or plot start, middle and end of simulation on default
            if nargin < 2
                hold on;
                plot(obj.xVector, obj.U(:, 1))
                plot(obj.xVector, obj.U(:, round(obj.numTimes/2)))
                plot(obj.xVector, obj.U(:, obj.numTimes))
                legend(strcat(...
                    num2str(...
                    [obj.currentTime(1);...
                    obj.currentTime(round(obj.numTimes/2)); ...
                    obj.currentTime(obj.numTimes)], 1) ...
                    , 's'))
                xlabel('x');
                ylabel('U(x,t)');
                hold off
            else
                plot(obj.xVector, obj.U(:, timeStep))
            end

        end

        function plotAnimated(obj, speed)
            % speed: number of timesteps between each frame
            clf
            if nargin < 2 || speed > obj.numTimes
                speed = 100;
            end
            for timestep = 1:speed:obj.numTimes
                temporaryU = obj.U(:, timestep);
                plot(obj.xVector, temporaryU)
                ylim([0 max(obj.U(:,1))])
                title('Timestep:', timestep)
                xlabel('x');
                ylabel('U(x,t)');
                drawnow;
            end
        end

        function plotIntegral(obj)
            S = sum(obj.U);
            plot(obj.currentTime, S)
        end

    end

end

