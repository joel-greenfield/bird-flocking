classdef Flock
    %FLOCK Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        boids
        lattice_size
        step_counter=1;
        
        center_x = zeros(300,1);
        center_y = zeros(300,1);
        spread = zeros(300,1);
    end
    
    methods
        
        function obj = Flock(boids,lattice_size)
            obj.boids=boids;
            obj.lattice_size=lattice_size;
        end
        
        function run(obj, plane)
            while true
                obj = flock(obj);
                obj = update_boids(obj);
                obj = borders(obj);
                
                % getting center and std of distances
                cent_x = 0;
                cent_y = 0;
                
                for i=1:length(obj.boids)
                    if ~isnan(obj.boids(i).position(1))
                        cent_x = cent_x + obj.boids(i).position(1);
                        cent_y = cent_y + obj.boids(i).position(2);
                    end
                end
                
                obj.center_x(obj.step_counter) = cent_x/length(obj.boids);
                obj.center_y(obj.step_counter) = cent_y/length(obj.boids);
                    
                positions = zeros(2,length(obj.boids));
                for i=1:1:length(obj.boids)
                    if ~isnan(obj.boids(i).position(1))
                        positions(:,i) = obj.boids(i).position;
                    end
                end
                d = pdist([[obj.center_x(obj.step_counter),obj.center_y(obj.step_counter)]; positions']);
                d = d(1:length(obj.boids));
                obj.spread(obj.step_counter) = std(d);

                [obj,plane] = render(obj,plane);
                
            end
        end
        
        function obj = update_boids(obj)
            for i=1:length(obj.boids)
                obj.boids(i)=obj.boids(i).update();
            end
        end
        
        function obj = flock(obj)
            for i=1:length(obj.boids)
                obj.boids(i)=obj.boids(i).flock(obj.boids);
            end
        end
        
        function [obj,plane] = render(obj,plane)
            fprintf('Rendering %s \n',num2str(obj.step_counter))
            csvwrite('myCenters.txt',[obj.center_x, obj.center_y])
            csvwrite('mySpread.txt',obj.spread)
            
            for i=1:length(obj.boids)
                delete(plane.boids_figure_handles(i));
                theta = atan2(norm(cross([obj.boids(i).velocity 0],[1 0 0])),dot([obj.boids(i).velocity 0],[1 0 0]));
                x = [obj.boids(i).position(1)-2.5 obj.boids(i).position(1)+2.5 obj.boids(i).position(1)-2.5 obj.boids(i).position(1)-2.5];
                y = [obj.boids(i).position(2)-1.5 obj.boids(i).position(2) obj.boids(i).position(2)+1.5 obj.boids(i).position(2)+1.5];
                plane.boids_figure_handles(i) =  patch(x,y,'k');
                rotate(plane.boids_figure_handles(i), [0 0 1], rad2deg(theta), [obj.boids(i).position(1) obj.boids(i).position(2) 0]);
            end
            delete(plane.boids_figure_handles(length(obj.boids)+1))
            pos = [obj.center_x(obj.step_counter)-1.5, obj.center_y(obj.step_counter)-1.5, 3, 3];
            plane.boids_figure_handles(length(obj.boids)+1) = rectangle('Position',pos,'Curvature',[1,1],'Facecolor','Blue');
            drawnow;
            
            obj.step_counter=obj.step_counter+1;
        end
        
        function obj = borders(obj)
            for i=1:length(obj.boids)
                obj.boids(i) = obj.boids(i).borders(obj.lattice_size);
            end
            
            if obj.center_x(obj.step_counter) < 0
                obj.center_x(obj.step_counter)=obj.lattice_size(1);
            end
            
            if obj.center_y(obj.step_counter) < 0
                obj.center_y(obj.step_counter)=obj.lattice_size(2);
            end
            
            if obj.center_x(obj.step_counter) > obj.lattice_size(1)
                obj.center_x(obj.step_counter)=0;
            end
            
            if obj.center_y(obj.step_counter) > obj.lattice_size(2)
                obj.center_y(obj.step_counter)=0;
            end
            
        end
        
    end
    
end

