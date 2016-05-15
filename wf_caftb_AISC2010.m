classdef wf_caftb_AISC2010 < wf_caftb_base
    
    properties

    end
    
    methods
        function obj = wf_caftb_AISC2010(shape_name,Fy)
            obj = obj.set_shape_properties(shape_name);
            obj.Fy = Fy;
        end
        
        function Qs = Qs(obj)
            if obj.bf_over_2tf < 0.56*sqrt(obj.E/obj.Fy)
                Qs = 1;
            elseif obj.bf_over_2tf >= 1.03*sqrt(obj.E/obj.Fy)
                Qs = (0.69*obj.E)/(obj.Fy*obj.bf_over_2tf^2);
            else
                Qs = 1.415-0.74*obj.bf_over_2tf*sqrt(obj.Fy/obj.E);
            end
        end
        
        function Qa = Qa(obj,f)
            h = obj.h_over_tw*obj.tw;
            if obj.h_over_tw < 1.49*sqrt(obj.E/f)
                be = obj.h_over_tw*obj.tw;
            else
                be = 1.92*obj.tw*sqrt(obj.E/f)*(1-((0.34/obj.h_over_tw)*sqrt(obj.E/f)));
                be = min(be,h);
            end
            Ae = obj.A -(h-be)*obj.tw;
            Qa = Ae/obj.A;
        end
        
        function Fcr = Fcr(obj,Fe,Q)
            if Q*obj.Fy/Fe < 2.25
                Fcr = (0.658^(Q*obj.Fy/Fe))*Q*obj.Fy;
            else
                Fcr = 0.877*Fe;
            end
        end
        
        function Pnx = Pnx(obj,L,K)
            KLr = K*L/obj.rx;
            Fe = pi^2*obj.E/KLr^2;
            Q = obj.Qs*obj.Qa(obj.Fcr(Fe,1));
            Fcr = obj.Fcr(Fe,Q);
            Pnx = Fcr*obj.A;
        end
        
        function Pny = Pny(obj,L,K)
            KLr = K*L/obj.ry;
            Fe = pi^2*obj.E/KLr^2;
            Q = obj.Qs*obj.Qa(obj.Fcr(Fe,1));
            Fcr = obj.Fcr(Fe,Q);
            Pny = Fcr*obj.A;
        end
        
        function Pnz = Pnz(obj,L,K)
            Fe = (pi^2*obj.E*obj.Cw/(K*L)^2 + obj.G*obj.J)/(obj.Ix+obj.Iy);
            Q = obj.Qs*obj.Qa(obj.Fcr(Fe,1));
            Fcr = obj.Fcr(Fe,Q);
            Pnz = Fcr*obj.A;
        end
        
        function Pnca = Pnca(obj,L,K)           
            Fe = 0.9*(((pi^2*obj.E*(obj.Cw+obj.Iy*(obj.d/2)^2))/(K*L)^2)+obj.G*obj.J)*(1/(obj.Ix+obj.Iy+(obj.d/2)^2*obj.A));
            Q = obj.Qs*obj.Qa(obj.Fcr(Fe,1));
            Fcr = obj.Fcr(Fe,Q);
            Pnca = Fcr*obj.A;
        end
        
        function [Pnys,tau] = Pnys(obj,Pr,L)
            Py = obj.Py;
            
            % Compute Q (with f = P/A)
            Q = obj.Qs*obj.Qa(Pr/obj.A);
            
            % Compute tau
            if Pr/(Q*Py) <= 0.39
                tau = 1;
            else
                tau = min(max(-2.724*Pr/(Q*Py)*log(Pr/(Q*Py)),0),1);
            end
            
            % Pny
            Pnys = 0.877*tau*(pi^2*obj.E*obj.Iy)/L^2;
            
            if nargout < 2
                clear tau;
            end
        end
        
        function Pno = Pno_lb(obj,Qtype)
            if nargin < 2
                Qtype = 'standard';
            end
            switch Qtype
                case 'standard'
                    Pno = obj.Qs*obj.Qa(obj.Fy)*obj.Fy*obj.A;
                case 'iterative'
                    fguess = obj.Fy;
                    options = struct;
                    options.Display = 'off';
                    [f,~,exitflag] =  fsolve(...
                        @(f)f-obj.Qs*obj.Qa(f)*obj.Fy,...
                        fguess,options);
                    if exitflag <= 0
                        error('fsolve could not find solution');
                    end
                    Pno = f*obj.A;
                otherwise
                    error('Unknown Qtype: %s',Qtype);
            end
        end
    end
end

