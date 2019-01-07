function [ vA1, vT1, vA2, vT2 ] = velocity_evolution_ConstantAngle( triangle_mid, r_mid, zero, first, second, r )
    
    %% side 0
    if zero
        % no meaning
    else
        
        %% side 1
        if first
            sinalpha1sq = sin(triangle_mid.alpha1)^2;
            tanalpha1 = tan(triangle_mid.alpha1);
            vA1mid = triangle_mid.v1A;
            K1 = vA1mid * r_mid ^ sinalpha1sq;
            
            vA1 = @(r) K1 ./ (r.^sinalpha1sq);
            vT1 = @(r) vA1(r) * tanalpha1;
            
            if(nargin == 6)
                vA1 = vA1(r);
                vT1 = vT1(r);
            end
            
            vA2 = [];
            vT2 = [];
        end
        
        %% side 2
        if second
            sinalpha2sq = sin(triangle_mid.alpha2)^2;
            tanalpha2 = tan(triangle_mid.alpha2);
            vA2mid = triangle_mid.v2A;
            K2 = vA2mid * r_mid ^ sinalpha2sq;
            
            vA2 = @(r) K2 ./ (r.^sinalpha2sq);
            vT2 = @(r) vA2(r) * tanalpha2;
            
            if(nargin == 6)
                vA2 = vA2(r);
                vT2 = vT2(r);
            end
            
            vA1 = [];
            vT1 = [];
        end
        
    end
end

