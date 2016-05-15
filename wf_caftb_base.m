classdef wf_caftb_base
    
    properties
        % Section Properties
        label
        A
        d
        tw
        Ix
        Iy
        J
        rx
        ry
        ho
        Cw
        bf_over_2tf
        h_over_tw
        
        % Material Properties
        Fy
        E = 29000;
        G = 11200;
        
        % Stiffness Resistance Factor
        phi_s = 0.75;
        
        % Compression Resistance Factor
        phi_c = 0.90;
    end
    
    methods
        function obj = set_shape_properties(obj,shape_name)
            load('wShapeData.mat');
            labels = {wShapeData(:).label};
            i = find(strcmpi(shape_name,labels));
            assert(isscalar(i),'Could not find: %s',shape_name);
            obj.label       = wShapeData(i).label;
            obj.A           = wShapeData(i).A;
            obj.d           = wShapeData(i).d;
            obj.tw          = wShapeData(i).tw;
            obj.Ix          = wShapeData(i).Ix;
            obj.Iy          = wShapeData(i).Iy;
            obj.J           = wShapeData(i).J;
            obj.rx          = wShapeData(i).rx;
            obj.ry          = wShapeData(i).ry;
            obj.ho          = wShapeData(i).ho;
            obj.Cw          = wShapeData(i).Cw;
            obj.bf_over_2tf = wShapeData(i).bf_over_2tf;
            obj.h_over_tw   = wShapeData(i).h_over_tw;
        end
        
        function Py = Py(obj)
            Py = obj.A*obj.Fy;
        end

        function beta_T = beta_T(obj,Pr,L)
            rs2  = obj.rx^2 + obj.ry^2 + (obj.d/2)^2;
            [Pnys,tau] = obj.Pnys(Pr,L);
            
            % beta
            if (Pnys*obj.d^2)/2 >= Pr*rs2
                beta_T = 0;
            else
                beta_T = (1.5/obj.phi_s)*(((Pr*rs2)-((Pnys*obj.d^2)/2))^2/(tau*obj.E*obj.Iy*obj.d^2));
            end
        end         

        function beta_Tb = beta_Tb(obj,Pr,L)
            beta_T   = obj.beta_T(Pr,L);
            beta_sec = (3.3*obj.E*obj.tw^3)/(12*obj.ho);
            
            if beta_T >= beta_sec
                beta_Tb = nan;
            else
                beta_Tb = beta_T/(1-beta_T/beta_sec);
            end
        end 

        function Pr = Pr_given_betaTb(obj,L,beta_Tbg)
            if beta_Tbg == 0
                rs2  = obj.rx^2 + obj.ry^2 + (obj.d/2)^2;
                
                Pguess = obj.Pnca(L,1);
                options = struct;
                options.Display = 'off';
                try
                    [Pr,~,exitflag] =  fsolve(...
                        @(P)P*rs2-0.5*obj.Pnys(P,L)*obj.d^2,...
                        Pguess,options);
                    if exitflag <= 0
                        error('fsolve could not find solution');
                    end                        
                catch err
                    fprintf('%s - L = %g - beta = %g\n',obj.label,L,beta_Tbg);
                    rethrow(err);
                end
                
            else
                beta_sec = (3.3*obj.E*obj.tw^3)/(12*obj.ho);
                beta_Tg = 1/(1/beta_Tbg+1/beta_sec);

                Pguess = obj.Pr_given_betaTb(L,0)+1e-2;
                options = struct;
                options.Display = 'off';
                try
                    [Pr,~,exitflag] =  fsolve(...
                        @(P)beta_Tg-obj.beta_T(P,L),...
                        Pguess,options);
                    if exitflag <= 0
                        [Pr,~,exitflag] =  fsolve(...
                            @(P)beta_Tg-obj.beta_T(P,L),...
                            max([Pguess+1 1.1*Pguess]),options);
                        if exitflag <= 0
                            Pguess
                            error('fsolve could not find solution');
                        end
                    end                        
                catch err
                    fprintf('%s - L = %g - beta = %g\n',obj.label,L,beta_Tbg);
                    rethrow(err);
                end
            end
        end
        
        function Pr = Pr_given_betaTb_and_phiMn(obj,L,beta_Tbg,phiMn) 
            if L ~= 0
                beta_Tbg = min(beta_Tbg,2*phiMn/(L/(1000*obj.ho)+phiMn/beta_Tbg));
            end
            Pr = obj.Pr_given_betaTb(L,beta_Tbg);
        end
        
    end
end

