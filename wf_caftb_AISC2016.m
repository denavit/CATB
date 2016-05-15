classdef wf_caftb_AISC2016 < wf_caftb_base
    
    properties

    end
    
    methods
        function obj = wf_caftb_AISC2016(shape_name,Fy)
            obj = obj.set_shape_properties(shape_name);
            obj.Fy = Fy;
        end
        
        function Fcr = Fcr(obj,Fe)
            if obj.Fy/Fe < 2.25
                Fcr = (0.658^(obj.Fy/Fe))*obj.Fy;
            else
                Fcr = 0.877*Fe;
            end
        end
        
        function Ae = Ae(obj,Fcr)
            % flange
            c1 = 0.22;
            c2 = (1-sqrt(1-4*c1))/(2*c1);
            lambda = obj.bf_over_2tf;
            lambda_r = 0.56*sqrt(obj.E/obj.Fy);
            Fel = (c2*lambda_r/lambda)^2*obj.Fy;
            if lambda <= lambda_r*sqrt(obj.Fy/Fcr)
                Ar_flange = 0;
            else
                b = obj.bf/2;
                be = b*(1-c1*sqrt(Fel/Fcr))*sqrt(Fel/Fcr);
                Ar_flange = 4*(b-be)*obj.tf;
            end
            
            % web
            c1 = 0.18;
            c2 = (1-sqrt(1-4*c1))/(2*c1);
            lambda = obj.h_over_tw;
            lambda_r = 1.49*sqrt(obj.E/obj.Fy);
            Fel = (c2*lambda_r/lambda)^2*obj.Fy;
            if lambda <= lambda_r*sqrt(obj.Fy/Fcr)
                Ar_web = 0;
            else
                b = obj.h_over_tw*obj.tw;
                be = b*(1-c1*sqrt(Fel/Fcr))*sqrt(Fel/Fcr);
                Ar_web = (b-be)*obj.tw;
            end
            
            Ae = obj.A - Ar_flange - Ar_web;
        end
        
        function Pnx = Pnx(obj,L,K)
            KLr = K*L/obj.rx;
            Fe  = pi^2*obj.E/KLr^2;
            Fcr = obj.Fcr(Fe);
            Pnx = Fcr*obj.Ae(Fcr);
        end
        
        function Pny = Pny(obj,L,K)
            KLr = K*L/obj.ry;
            Fe  = pi^2*obj.E/KLr^2;
            Fcr = obj.Fcr(Fe);
            Pny = Fcr*obj.Ae(Fcr);
        end
        
        function Pnz = Pnz(obj,L,K)
            Fe = (pi^2*obj.E*obj.Cw/(K*L)^2 + obj.G*obj.J)/(obj.Ix+obj.Iy);
            Fcr = obj.Fcr(Fe);
            Pnz = Fcr*obj.Ae(Fcr);
        end
        
        function Pnca = Pnca(obj,L,K,a)
            if nargin < 4
                a = obj.ho/2;
            end
            %Fe = 0.9*(((pi^2*obj.E*(obj.Cw+obj.Iy*(obj.d/2)^2))/(K*L)^2)+obj.G*obj.J)*(1/(obj.Ix+obj.Iy+(obj.d/2)^2*obj.A))
            ro2 = obj.rx^2+obj.ry^2+a^2;
            Fe = 0.9*((pi^2*obj.E*obj.Iy)/(K*L)^2*(obj.ho^2/4+a^2)+obj.G*obj.J)/(obj.A*ro2);
            Fcr = obj.Fcr(Fe);
            Pnca = Fcr*obj.Ae(Fcr);
        end
        
        function [Pnys,tau,x] = Pnys(obj,Pr,L)
            Py = obj.Py;

            % Compute x         
            options = struct;
            options.Display = 'off';
            try
                [Ae,~,exitflag] =  fsolve(...
                    @(Ae)Ae-obj.Ae(Pr/Ae),...
                    obj.A,options);
                if exitflag <= 0
                    error('fsolve could not find solution');
                end                        
            catch err
                rethrow(err);
            end            
            x = Ae/obj.A;            
            
            % Compute tau
            if Pr/(x*Py) <= 0.39
                tau = x;
            else
                tau = min(max(-2.724*Pr/Py*log(Pr/(x*Py)),0),1);
            end
            
            % Pny
            Pnys = 0.877*tau*(pi^2*obj.E*obj.Iy)/L^2;
            
            if nargout < 2
                clear tau;
            end
        end
 
        function Pno = Pno_lb(obj)
            Pno = obj.Fy*obj.Ae(obj.Fy);
        end
    end
end

