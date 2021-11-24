%% 0. Initialize Parameters
n = 1250; % Number of locations to evaluate bridge failure
L = 1250; % Length of bridge
x = linspace(0, L, n); % Define x coordinate
P = zeros(1,n); % Initializes Loads 

%% 1. Point Loading Analysis (SFD, BMD)
load = -318;
[SFD_PL, BMD_PL, P] = ApplyPL(550, load, x, P, n); % Construct SFD, BMD
[SFD_PL, BMD_PL, P] = ApplyPL(L, load, x, P, n); % Construct SFD, BMD

%% 2. Define cross-sections
% There are many (more elegant ways) to construct cross-section objects
xc = [0 550 L]; % Location, x, of cross-section change
bft = [100 100 100]; % Top Flange Width
tft = [2.54 2.54 2.54]; % Top Flange Thickness
hw = [100 120 100]; % Web Height
tw = [1.27 1.27 1.27]; % Web Thickness (Assuming 2 separate webs)
bfb = [80 80 80]; % Bottom Flange Width
tfb = [1.27 1.27 1.27]; % Bottom Flange Thickness
a = [400 400 400]; % Diaphragm Spacing
ws = [77.46 77.46 77.46]; % Web Spacing
rtw = [10 10 10]; % little rectanges at top of webs for glue area width
rtt = [1.27 1.27 1.27]; % little rectanges at top of webs for glue area thickness
rbw = [0 0 0]; % little rectanges at bottom of webs for glue area width
rbt = [0 0 0]; % little rectanges at bottom of webs for glue area thickness

% Optional but you need to ensure that your geometric inputs are correctly implemented
VisualizeBridge(xc, bft, tft, hw, tw, ws, bfb, tfb, rtw, rtt, rbw, rbt); 

%% 3. Define Material Properties
SigT = 30;
SigC = 6;
E = 4000;
TauU = 4;
TauG = 2;
mu = 0.2;

[Y_bridge, I_bridge, Q_bridge] = SectionProperties(xc, bft, tft, hw, tw, ws, bfb, tfb, rtw, rtt, rbw, rbt, n );

y = [1:124];
Q_bridge(1, 64)
plot(Q_bridge(1, 1:end))
Y_bridge(1)

%%
function [ SFD, BMD, Loads ] = ApplyPL( xP, P, x, Loads, n )
% Constructs SFD and BMD from application of 1 Point Load. Assumes fixed location of supports
%   Input: location and magnitude of point load. The previous Loads can be entered as input to 
%       construct SFD of multiple point loads
%   Output: SFD, BMD, Loads, all 1-D arrays of length n
    Loads(xP) = P;
    moment = 0;
    force = 0;
    
    % Calculates sum of forces and moments for use in reaction force
    % calculations
    for i = 1:length(Loads)
        if and(i ~= 1060, i ~= 1)
            moment = Loads(i) * i + moment;
            force = Loads(i) + force;
        end
    end
    
    % calculates reaction forces
    b = -moment/1060;
    a = -b - force;
    
    % Adds reaction forces to Loads array
    Loads(1) = a;
    Loads(1060) = b;
    
    % Initializes new SFD
    SFD = zeros(1, n);
    
    % Calculates shear force at every point of bridge
    shear = 0;
    for i = 1:length(Loads)
        if Loads(i) ~= 0
            shear = shear + Loads(i);
        end
        SFD(i) = shear;
    end
    
    % Integrates SFD
    BMD = cumtrapz(x, SFD);
    
    plot(x, SFD)
    set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
    figure;
    plot(x, BMD)
    set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'YDir', 'reverse')
end

%%
function [ ] = VisualizeBridge( csc, tfw, tft, wh, wt, ws, bfw, bft, rw, rh, rbw, rbh)
    % Provides a graphical interpretation of user geometric inputs
    %    Input: cross section changes x cordinates, top flange width, top flange thinckness, web height, web thickness, web spacing,
    %        bottom flange width, bottom flange thickness
    %    Output: polygons of crossections
    
    for i = 1:length(csc)
        % Calculates basic coordinates of features
        y_top_flang_top_cord = bft(i) + wh(i) + tft(i);
        y_top_flang_bot_cord = y_top_flang_top_cord - tft(i);
        
        x_right_web_cord = tfw(i) - ((tfw(i) - ws(i)) / 2) + wt(i);
        x_left_web_cord = ((tfw(i) - ws(i)) / 2) - wt(i);
        
        x_left_bot_flang_cord = ((tfw(i) - bfw(i)) / 2);
        x_right_bot_flang_cord = tfw(i) - ((tfw(i) - bfw(i)) / 2);

        outside_x_cords = [x_left_bot_flang_cord x_left_bot_flang_cord x_left_web_cord x_left_web_cord 0 0 tfw(i) tfw(i) x_right_web_cord x_right_web_cord x_right_bot_flang_cord x_right_bot_flang_cord ];
        outside_y_cords = [0 bft(i) bft(i) y_top_flang_bot_cord y_top_flang_bot_cord y_top_flang_top_cord y_top_flang_top_cord y_top_flang_bot_cord y_top_flang_bot_cord bft(i) bft(i) 0];
        
        x_left_inside_web = x_left_web_cord + wt(i);
        x_right_inside_web = x_right_web_cord - wt(i);

        x_left_top_rectangle = x_left_inside_web + rw(i);
        x_right_top_rectangle = x_right_inside_web - rw(i);

        x_left_bot_rectangle = x_left_inside_web + rbw(i);
        x_right_bot_rectangle = x_right_inside_web - rbw(i);


        y_rectangle_top = y_top_flang_bot_cord - rh(i);
        y_rectangle_bot = bft(i) + rbh(i);

        inside_x_cords = [x_left_bot_rectangle x_left_bot_rectangle x_left_inside_web x_left_inside_web x_left_top_rectangle x_left_top_rectangle x_right_top_rectangle x_right_top_rectangle x_right_inside_web x_right_inside_web x_right_bot_rectangle x_right_bot_rectangle];
        inside_y_cords = [bft(i) y_rectangle_bot y_rectangle_bot y_rectangle_top y_rectangle_top y_top_flang_bot_cord y_top_flang_bot_cord y_rectangle_top y_rectangle_top y_rectangle_bot y_rectangle_bot bft(i)];

        % Draws rectangles for features
        figure;
        cross_section_shape = polyshape(outside_x_cords, outside_y_cords);
        cross_section_shape = addboundary(cross_section_shape, inside_x_cords, inside_y_cords);

        plot(cross_section_shape)
        [p, pp] = centroid(cross_section_shape)
        
    end
end

%%
function [ Y_bar, I, Q, b ] = SectionProperties( csc, tfw, tft, wh, wt, ws, bfw, bft, rtw, rtt, rbw, rbt, n )
    % Calculates important sectional properties. Including but not limited to ybar, I, Q, etc.
    %    Input: cross section changes x cordinates, top flange width, top flange thinckness, web height, web thickness, web spacing,
    %        bottom flange width, bottom flange thickness
    % Output: Sectional Properties at every value of x. Each property is a 1-D array of length n

    Y_bar = zeros(1,n);
    I = zeros(1, n);

    for i = 1:length(csc)
         heights(i) = (tft(i) + wh(i) + bft(i));
    end
    
    max_height = ceil(max(heights));
    Q = zeros(n, max_height);

    for i = 1:length(csc)
        % Areas for features
        A_top = tfw(i) * tft(i);
        A_webs = 2 * wt(i) * wh(i);
        A_bot = bft(i) * bfw(i);
        A_r_top = 2 * (rtw(i) * rtt(i));
        A_r_bot = 2 * (rbw(i) * rbt(i));
        
        % Local centriods for features
        Y_top = (tft(i)/2 + wh(i) + bft(i));
        Y_webs = (wh(i)/2 + bft(i));
        Y_bot = (bft(i)/2);
        Y_r_top = (wh(i) + bft(i) - rtt(i)/2);
        Y_r_bot = (bft(i) + rbt(i)/2)
        
        % I for features
        I_top = (tfw(i) * (tft(i) ^ 3)) /12;
        I_webs = 2 * (wt(i) * (wh(i) ^ 3)) /12;
        I_bot = (bfw(i) * (bft(i) ^ 3)) /12;
        I_r_top = 2 * (rtw(i) * (rtt(i) ^ 3)) /12;
        I_r_bot = 2 * (rbw(i) * (rbt(i) ^ 3)) /12;

        % Calculates I, Y and Q for cross-section
        Y_section = ((A_top * Y_top) + (A_webs * Y_webs) + (A_bot * Y_bot ) + (A_r_top * Y_r_top) + (A_r_bot * Y_r_bot)) / (A_bot + A_top + A_webs + A_r_bot + A_r_top);
        I_section = (I_top + (A_top * (Y_top ^ 2))) + (I_webs + (A_webs * (Y_webs ^ 2))) + (I_bot + (A_bot * (Y_bot ^ 2))) + (I_r_top + (A_r_top * (Y_r_top ^ 2))) + (I_r_bot + (A_r_bot * (Y_r_bot ^ 2)));
 
        % Adds I and Y_bar for every x value of cross-section
        if i ~= length(csc)
            for j = csc(i):(csc(i+1)-1)
                Y_bar(j+1) = Y_section;
                I(j+1) = I_section;

                for y = 1:(tft(i) + wh(i) + bft(i));
                    if y <= bft(i);
                        sub_area = y * bfw(i);
                        sub_centriod = (bft(i) - y)/2;
                        d = abs((Y_bar(csc(i) + 1) - sub_centriod));
                        Q(j+1, y) = sub_area * d;
                    elseif bft(i) < y && y <= (bft(i) + wh(i));
                        web_height = y - bft(i);
                        sub_area = A_bot + (2 * web_height  * wt(i));
                        sub_centriod = ((A_bot * Y_bot) + ((2 * web_height * wt(i)) * ((web_height / 2) + bft(i)))) / sub_area;
                        d = abs(Y_bar(csc(i) + 1) - sub_centriod);
                        Q(j+1, y) = sub_area * d;
                    else
                        top_height = y - wh(i) - bft(i);
                        sub_area = A_bot + A_webs + (top_height * tfw(i));
                        sub_centriod = ((A_bot * Y_bot) + (A_webs * Y_webs) + (top_height * tfw(i) * ((top_height / 2) + wh(i) + bft(i)))) / sub_area;
                        d = abs(Y_bar(csc(i) + 1) - sub_centriod);
                        Q(j+1, y) = sub_area * d;
                    end
                end
            end
        else
            for j = csc(i):(n)
                Y_bar(j) = Y_section;
                I(j) = I_section;

                for y = 1:(tft(i) + wh(i) + bft(i))
                    if y <= bft(i)
                        sub_area = y * bfw(i);
                        sub_centriod = (bft(i) - y)/2;
                        d = abs(Y_bar(csc(i)) - sub_centriod);
                        Q(j, y) = sub_area * d;
                    elseif bft(i) < y && y <= (bft(i) + wh(i))
                        web_height = y - bft(i);
                        sub_area = A_bot + (2 * web_height  * wt(i));
                        sub_centriod = ((A_bot * Y_bot) + ((2 * web_height * wt(i)) * ((web_height / 2) + bft(i)))) / sub_area;
                        d = abs(Y_bar(csc(i)) - sub_centriod);
                        Q(j, y) = sub_area * d;
                    else
                        top_height = y - wh(i) - bft(i);
                        sub_area = A_bot + A_webs + (top_height * tfw(i));
                        sub_centriod = ((A_bot * Y_bot) + (A_webs * Y_webs) + (top_height * tfw(i) * ((top_height / 2) + wh(i) + bft(i)))) / sub_area;
                        d = abs(Y_bar(csc(i)) - sub_centriod);
                        Q(j, y) = sub_area * d;
                    end
                end
                
            end
        end
    end
end

function [ V_fail ] = Vfail( I, b, TauU )
    % Calculates shear forces at every value of x that would cause a matboard shear failure
    % Input: Sectional Properties (list of 1-D arrays), TauU (scalar material property)
    % Output: V_fail a 1-D array of length n
    I = {Sectional Properties};
    b = {Sectional Properties};
    Qcent = {Sectional Properties};
    
    V_fail = TauU .* I .* b ./ Qcent;


    end    