function [rocket] = get_Rocket(rocket)
% Fills new rocket parameters, based on initial OpenRocket.
% Final version will have rand ranges around initial values

%% Geo
geo = rocket.geo;
    % rocket.geo.nose = struct('D',[],'L',[],'t',[],'x',0,'A_p',[],'A_w',[]);
        geo.nose.D = 7.8/12;
        geo.nose.L = 30/12;
        geo.nose.t = 0.079;
        geo.nose.A_p = 0.5 * geo.nose.D * geo.nose.L;
    % rocket.geo.body = struct('D',[],'L',[],'t',[],'x',[],'A_p',[],'K',1.5,'A_w',[]);
        geo.body.D = geo.nose.D;
        geo.body.L = 100/12;
        geo.body.t = 0.079;
        geo.body.x = geo.nose.L;
        geo.body.A_p = geo.body.D * geo.body.L;
    % rocket.geo.tail = struct('D',[],'L',[],'t',[],'x',[],'A_p',[],'R1',[],'R2',[],'A_w',[]);
%         geo.tail.D = ;
%         geo.tail.L = ;
%         geo.tail.t = ;
%         geo.tail.x = ;
%         geo.tail.A_p = ;
%         geo.tail.R1 = ;
%         geo.tail.R2 = ;
    % rocket.geo.fins = struct('RC',[],'TC',[],'SS',[],'SL',[],'MC',[],'N',[],'t',[],'x',[]),'A_w',[];
        % ignore for now
    % rocket.geo.press_t = struct('D',[],'L',[],'t',[],'x',[]);
    % rocket.geo.ox_t    = struct('D',[],'L',[],'t',[],'x',[]);
    % rocket.geo.CC      = struct('D',[],'L',[],'t',[],'x',[]);
    % rocket.geo.fuel    = struct('D_o',[],'D_i',[],'L',[]);
        geo.fuel.D_i = 2/12;
rocket.geo = geo;

%% Prop
prop = rocket.prop;
%     rocket.prop = struct('T_avg',[],'Isp',[],'t_b',[],'I',[],'OF',[]);
        prop.T_avg = 4423;
        prop.Isp = 250; % will fix
        prop.I = 9200; % will fix
        prop.t_b = prop.I / prop.T_avg;
        prop.OF = 6; % not sure if this is varied
        prop.P_c = 500*144; % [psf]
rocket.prop = prop;

%% Weight
W = rocket.weight;
        rocket = get_PropSys(rocket);
%     rocket.weight.nose = struct('W',[],'I',[],'I_pt',[],'CG',[]);
        W.nose.W = 2.57 * (geo.nose.L*geo.nose.D*geo.nose.t)/(2.5*0.65*0.0066); % Assumed fiberglass, conical nosecone, no tip
        W.nose.CG = 2*geo.nose.L/3;
%     rocket.weight.body = struct('W',[],'I',[],'I_pt',[],'CG',[]);
        W.body.W = 111.24 * pi*geo.body.L*(geo.body.D*geo.body.t - geo.body.t^2);
        W.body.CG = geo.body.x + geo.body.L/2;
%     rocket.weight.tail = struct('W',[],'I',[],'I_pt',[],'CG',[]);
%     rocket.weight.fins = struct('W',[],'I',[],'I_pt',[],'CG',[]);
        W.fins.W = 4; % Guess from Endurance
%     rocket.weight.press_t = struct('W',[],'I',[],'I_pt',[],'CG',[]);
%     rocket.weight.ox_t = struct('W',[],'I',[],'I_pt',[],'CG',[]);
        
%     rocket.weight.CC = struct('W',[],'I',[],'I_pt',[],'CG',[]);
%     rocket.weight.misc = struct('W',[],'I',[],'I_pt',[],'CG',[]); 
rocket.weight = W;

end