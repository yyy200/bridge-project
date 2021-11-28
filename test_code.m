%% 0. Initialize Parameters
n = 1250; % Number of locations to evaluate bridge failure
L = 1250; % Length of bridge
x = linspace(0, L, n); % Define x coordinate
P = zeros(1,n); % Initializes Loads 

%% 1. Point Loading Analysis (SFD, BMD)
load = -205;
[SFD_PL, BMD_PL, P] = ApplyPL(550, load, x, P, n); % Construct SFD, BMD
[SFD_PL, BMD_PL, P] = ApplyPL(L, load, x, P, n); % Construct SFD, BMD

%% 2. Define cross-sections
% There are many (more elegant ways) to construct cross-section objects
% xc = [0 350 415 600 650 800 950 980 1100, L]; % Location, x, of cross-section change
% tfw = [100 100 100 100 100 100 100 100 100 100]; % Top Flange Width
% tft = [1.27 1.27*1.5 1.27*2 1.27*1.5 1.27  1.27 1.27 1.27 1.27 1.27]*2; % Top Flange Thickness
% hw = [72.46 72.46 72.46 72.46 72.46 72.46 72.46 72.46 72.46 72.46]; % Web Height
% tw = [1.27 1.27 1.27 1.27 1.27 2*1.27  2*1.27 2*1.27  2*1.27 2*1.27]*2; % Web Thickness (Assuming 2 separate webs)
% bfw = [80 80 80 80 80 90 90 90 90 85]; % Bottom Flange Width
% bft = [1.27 1.27 1.27 1.27 1.27 1.27*1.5 1.27*1.5 1.27*2.5 1.27*2 1.27]*2; % Bottom Flange Thickness
% a = [0 275 300 550 670 805 900 990 1155 L]; % Diaphragm x coords
% ws = [74.92 74.92 74.92 74.92 74.92 74.92 74.92 74.92 74.92 74.92]; % Web Spacing
% rtw = [10 10 10 10 10 10 10 10 10 10]; % little rectanges at top of webs for glue area width
% rtt = [1.27 1.27 1.27 1.27 1.27 1.27 1.27 1.27 1.27 1.27]; % little rectanges at top of webs for glue area thickness
% rbw = [10 10 10 10 10 10 10 10 10 15]; % little rectanges at bottom of webs for glue area width
% rbt = [1.27 1.27 1.27 1.27 1.27 1.27 1.27 1.27 1.27 1.27]; % little rectanges at bottom of webs for glue area thickness

xc = [0, 800, L]
tfw = [100 100 100]
tft = [0, 1.27*4 1.27]
hw = [72.46 72.46 72.46]
tw = [1.27 1.27*1.5 2*1.27]*2
bfw = [0 0 90]
bft = [0 0 1.27*2]
a = xc
ws = [74.92 74.92 74.92]
rtw = [10 10 10]
rtt = [1.27 1.27 1.27]
rbw = [0 0 10]
rbt = [0 0 1.27]

% Optional but you need to ensure that your geometric inputs are correctly implemented
VisualizeBridge(xc, tfw, tft, hw, tw, ws, bfw, bft, rtw, rtt, rbw, rbt); 

%% 3. Define Material Properties
SigT = 30;
SigC = 6;
E = 4000;
TauU = 4;
TauG = 2;
mu = 0.2;

[Y_bridge, I_bridge, Q_bridge, b_bridge, section_heights] = SectionProperties(xc, tfw, tft, hw, tw, ws, bfw, bft, rtw, rtt, rbw, rbt, n );
Y_bridge(1)
I_bridge(1)


[v_fail] = Vfail(I_bridge, b_bridge, Y_bridge, TauU, Q_bridge);
[ v_buck ] = VfailBuck(xc, tw, a, hw, E, mu, n, I_bridge, b_bridge, Y_bridge, Q_bridge);
[ m_mat_tension ] = MfailMatT( I_bridge, Y_bridge, section_heights, SigT, BMD_PL);
[ m_mat_compression ] = MfailMatC( I_bridge, Y_bridge, section_heights, SigC, BMD_PL );
[ M_Buck1, M_Buck2, M_Buck3 ] = MfailBuck( xc, bfw, bft, tfw, tft, ws, tw, section_heights, Y_bridge, I_bridge, E, mu, BMD_PL);
[ V_GlueTF V_GlueBF V_GlueTW V_GlueBW] = VglueFail(I_bridge, Q_bridge, b_bridge, TauG, tft, bft, section_heights, xc);
V_GlueTW(1)
% if NUM = 1: matt, 2: matc, 3: buck1, 4: buck2, 5: buck3, 6: vmat, 7: vbuck, 8: vgluetf, 9: vgluebf, 10: vgluetw, 11: vgluebw

% FAILURE TYPES:
% matt = moment matboard tension failure
% matc = moment matboard compression failure
% buck1 = case 1 buckling of top or bottom flange (center part between webs)
% buck2 = case 2 buckling of top or bottom flange (part that sticks out past the webs)
% buck3 = case 3 buckling of left and right webs
% vmat = matboard shear failure
% vbuck = matboard shear buckling failure
% vgluetf = shear glue failure of top flange
% vgluebf = shear glue failure of bottom flange
% vgluetw = shear glue failure of top flange/web connection
% vgluebw = shear glue failure of bottom flange/web connection

[ Pf, failure_mode ] = FailLoad( load, SFD_PL, BMD_PL, v_fail, v_buck, m_mat_tension, m_mat_compression, M_Buck1, M_Buck2, M_Buck3, V_GlueTF, V_GlueBF, V_GlueTW, V_GlueBW, true, "Design0");
VisulizePL(x, load, Pf, SFD_PL, BMD_PL, v_fail, v_buck, m_mat_tension, m_mat_compression, M_Buck1, M_Buck2, M_Buck3, V_GlueTF, V_GlueBF, V_GlueTW, V_GlueBW);
Pf
failure_mode

[ material_ok ] = MaterialCheck(xc, tfw, tft, hw, tw, ws, bfw, bft, rtw, rtt, rbw, rbt); 
"Material is: " +  material_ok

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
        cross_section_shape = polyshape(outside_x_cords, outside_y_cords);
        cross_section_shape = addboundary(cross_section_shape, inside_x_cords, inside_y_cords);

        subplot(ceil(length(csc)/ 2), ceil(length(csc)/ 2), i)
        plot(cross_section_shape)
        if i ~= 1
            title("Cross Section from x = " + csc(i-1) + " to " + csc(i))
        else
            title("Cross Section at x = 0")
        end
        axis equal
        
    end
end

function out = getVarName(var)
    out = inputname(1);
end

%%
function [ Y_bar, I, Q, b, heights ] = SectionProperties( csc, tfw, tft, wh, wt, ws, bfw, bft, rtw, rtt, rbw, rbt, n )
    % Calculates important sectional properties. Including but not limited to ybar, I, Q, etc.
    %    Input: cross section changes x cordinates, top flange width, top flange thinckness, web height, web thickness, web spacing,
    %        bottom flange width, bottom flange thickness
    % Output: Sectional Properties at every value of x. Each property is a 1-D array of length n

    Y_bar = zeros(1,n);
    I = zeros(1, n);

    for i = 1:length(csc)-1
         heights(csc(i)+1:csc(i+1)) = (tft(i) + wh(i) + bft(i));
    end
    max_height = ceil(max(heights));
    Q = zeros(n, max_height);
    b = zeros(n, max_height);

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
        Y_r_bot = (bft(i) + rbt(i)/2);
        
        % I for features
        I_top = (tfw(i) * (tft(i) ^ 3)) /12;
        I_webs = 2 * (wt(i) * (wh(i) ^ 3)) /12;
        I_bot = (bfw(i) * (bft(i) ^ 3)) /12;
        I_r_top = 2 * (rtw(i) * (rtt(i) ^ 3)) /12;
        I_r_bot = 2 * (rbw(i) * (rbt(i) ^ 3)) /12;

        % Calculates I and Y for cross-section
        Y_section = ((A_top * Y_top) + (A_webs * Y_webs) + (A_bot * Y_bot ) + (A_r_top * Y_r_top) + (A_r_bot * Y_r_bot)) / (A_bot + A_top + A_webs + A_r_bot + A_r_top);
        I_section = (I_top + (A_top * ((Y_section - Y_top) ^ 2))) + (I_webs + (A_webs * ((Y_section - Y_webs) ^ 2))) + (I_bot + (A_bot * ((Y_section - Y_bot) ^ 2))) + (I_r_top + (A_r_top * ((Y_section - Y_r_top) ^ 2))) + (I_r_bot + (A_r_bot * ((Y_section - Y_r_bot) ^ 2)));

        % Calculates Q, Adds I, Q and Y_bar for every x value of cross-section
        if i ~= length(csc)
            for j = csc(i):(csc(i+1)-1)
                Y_bar(j+1) = Y_section;
                I(j+1) = I_section;
                
                for y = 1:((tft(i) + wh(i) + bft(i)));

                    if y <= bft(i)
                        b(j+1, y) = bft(i);
                        sub_area = y * bfw(i);
                        sub_centriod = (bft(i) - y)/2;
                        d = abs(Y_bar(csc(i) + 1) - sub_centriod);
                        Q(j+1, y) = sub_area * d;

                    elseif bft(i) < y && y <= (bft(i) + rbt(i))
                        b(j+1, y) = 2 * (rbw(i) + wt(i));
                        web_height = y - bft(i);
                        sub_area = A_bot + (2 * web_height  * (wt(i) + rbw(i)));
                        sub_centriod = ((A_bot * Y_bot) + ((2 * web_height * wt(i)) * ((web_height / 2) + bft(i)))) / sub_area;
                        d = abs(Y_bar(csc(i) + 1) - sub_centriod);
                        Q(j+1, y) = sub_area * d;

                    elseif (bft(i) + rbt(i)) < y && y <= (bft(i) + wh(i) - rtt(i))
                        b(j+1, y) = 2 * wt(i);
                        web_height = y - bft(i) - rbt(i);
                        sub_area = A_bot + (2 * web_height  * wt(i)) + A_r_bot;
                        sub_centriod = ((A_bot * Y_bot) + ((2 * web_height * wt(i)) * ((web_height / 2) + bft(i) + rbt(i))) + (A_r_bot * Y_r_bot)) / sub_area;
                        d = abs(Y_bar(csc(i) + 1) - sub_centriod);
                        Q(j+1, y) = sub_area * d;

                    elseif  (bft(i) + wh(i) - rtt(i)) < y && y <= (bft(i) + wh(i))
                        b(j+1, y) = 2 * (rtw(i) + wt(i));
                        web_height = y - wh(i) + rtt(i) - bft(i);
                        main_web_area =  (2 * (wh(i) - rtt(i)) * wt(i));
                        sub_area = A_bot + (2 * web_height  * (wt(i) + rtw(i))) + A_r_bot + main_web_area;
                        sub_centriod = ((A_bot * Y_bot) + (2 * web_height * (wt(i) + rtw(i)) * ((web_height / 2) + bft(i) + (wh(i) - rtt(i)))) + (A_r_bot * Y_r_bot) + (main_web_area * (((wh(i) - rtt(i))/2) + bft(i)))) / sub_area;
                        d = abs(Y_bar(csc(i) + 1) - sub_centriod);
                        Q(j+1, y) = sub_area * d;

                    else
                        b(j+1, y) = tfw(i);
                        top_height = y - wh(i) - bft(i);
                        sub_area = A_bot + A_webs + (top_height * tfw(i)) + A_r_bot + A_r_top;
                        sub_centriod = ((A_bot * Y_bot) + (A_webs * Y_webs) + (A_r_top * Y_r_top) + (A_r_bot * Y_r_bot) + (top_height * tfw(i) * ((top_height / 2) + wh(i) + bft(i)))) / sub_area;
                        d = abs(Y_bar(csc(i) + 1) - sub_centriod);
                        Q(j+1, y) = sub_area * d;
                    end
                end
            end
        else
            for j = csc(i):(n)
                Y_bar(j) = Y_section;
                I(j) = I_section;
                for y = 1:((tft(i) + wh(i) + bft(i)))

                    if y <= bft(i)
                        b(j, y) = bft(i);
                        sub_area = y * bfw(i);
                        sub_centriod = (bft(i) - y)/2;
                        d = abs(Y_bar(csc(i)) - sub_centriod);
                        Q(j, y) = sub_area * d;

                    elseif bft(i) < y && y <= (bft(i) + rbt(i))
                        b(j, y) = 2 * (rbw(i) + wt(i));
                        web_height = y - bft(i);
                        sub_area = A_bot + (2 * web_height  * (wt(i) + rbw(i)));
                        sub_centriod = ((A_bot * Y_bot) + ((2 * web_height * wt(i)) * ((web_height / 2) + bft(i)))) / sub_area;
                        d = abs(Y_bar(csc(i)) - sub_centriod);
                        Q(j, y) = sub_area * d;

                    elseif (bft(i) + rbt(i)) < y && y <= (bft(i) + wh(i) - rtt(i))
                        b(j, y) = 2 * wt(i);
                        web_height = y - bft(i) - rbt(i);
                        sub_area = A_bot + (2 * web_height  * wt(i)) + A_r_bot;
                        sub_centriod = ((A_bot * Y_bot) + ((2 * web_height * wt(i)) * ((web_height / 2) + bft(i) + rbt(i))) + (A_r_bot * Y_r_bot)) / sub_area;
                        d = abs(Y_bar(csc(i)) - sub_centriod);
                        Q(j, y) = sub_area * d;

                    elseif  (bft(i) + wh(i) - rtt(i)) < y && y <= (bft(i) + wh(i))
                        b(j, y) = 2 * (rtw(i) + wt(i));
                        web_height = y - wh(i) + rtt(i) - bft(i);
                        main_web_area =  (2 * (wh(i) - rtt(i)) * wt(i));
                        sub_area = A_bot + (2 * web_height  * (wt(i) + rtw(i))) + A_r_bot + main_web_area;
                        sub_centriod = ((A_bot * Y_bot) + (2 * web_height * (wt(i) + rtw(i)) * ((web_height / 2) + bft(i) + (wh(i) - rtt(i)))) + (A_r_bot * Y_r_bot) + (main_web_area * (((wh(i) - rtt(i))/2) + bft(i)))) / sub_area;
                        Q(j, y) = sub_area * d;

                    else
                        b(j, y) = tfw(i);
                        top_height = y - wh(i) - bft(i);
                        sub_area = A_bot + A_webs + (top_height * tfw(i)) + A_r_bot + A_r_top;
                        sub_centriod = ((A_bot * Y_bot) + (A_webs * Y_webs) + (A_r_top * Y_r_top) + (A_r_bot * Y_r_bot) + (top_height * tfw(i) * ((top_height / 2) + wh(i) + bft(i)))) / sub_area;
                        d = abs(Y_bar(csc(i)) - sub_centriod);
                        Q(j, y) = sub_area * d;
                    end
                end
                
            end
        end
    end
end

function [ V_fail ] = Vfail( I, b, Y_bar, TauU, Q)
    % Calculates shear forces at every value of x that would cause a matboard shear failure
    % Input: Sectional Properties (list of 1-D arrays), TauU (scalar material property)
    % Output: V_fail a 1-D array of length n
    Qcent = zeros(1,length(Y_bar));
    bcent = zeros(1,length(Y_bar));
    for i = 1:length(Y_bar)
        Qcent(i) = Q(i, round(Y_bar(i)));
        bcent(i) = b(i, round(Y_bar(i)));
    end
    
    V_fail = TauU .* I .* bcent ./ Qcent;
    
end    

function [ V_Buck ] = VfailBuck(csc, wt, ds, wh, E, mu, n, I, b, Y_bar, Q)
    % Calculates shear forces at every value of x that would cause a shear buckling failure in the web
    %   Input: Sectional Properties (list of 1-D arrays), E, mu (material property)
    %   Output: V_Buck a 1-D array of length n

        
    factor = (5 * (pi ^ 2) * E ) / (12 * (1 - (mu ^ 2) ));

    a_values = zeros(1, n);
    web_thicknesses = zeros(1, n);
    web_heights = zeros(1, n);
    for i = 1:(length(csc)-1)
        web_thicknesses((csc(i)+1):(csc(i+1))) = wt(i);
        web_heights((csc(i)+1):(csc(i+1))) = wh(i);
        a_values((csc(i)+1):(csc(i+1))) = ds(i+1) - ds(i);
    end
       
    TauCrit = factor .* (((web_thicknesses ./ web_heights) .^ 2) + ((web_thicknesses ./ a_values) .^ 2));

    [ V_Buck ] = Vfail( I, b, Y_bar, TauCrit, Q);

end

function [] = VisulizePL(x, P, Pfail, SFD, BMD, V_fail, V_buck, M_MatT, M_MatC, M_Buck1, M_Buck2, M_Buck3, V_GlueTF, V_GlueBF, V_GlueTW, V_GlueBW)
    
    figure('WindowState', 'maximized')

    % Plots SFD for Pfail
    subplot(2,3,1);
    hold on;
    plot(x, SFD);

    set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'origin');
    title(("SFD from P = " + P + "N, Pfail = " +  Pfail + "N"))
    xlim([0,x(end)])
    grid on;
    grid minor
    box on;
    yline(0)
    xlabel('Position on Bridge (mm)')
    ylabel('Shear Force (N)')


    % Plots BMD vs V_buck
    subplot(2,3,2);
    hold on;
    plot(x, SFD);
    plot(x, V_fail(V_fail~=0), 'r-');
    plot(x, -V_fail(V_fail~=0), 'r-');
    
    plot(x, V_GlueTF(V_GlueTF~=0), 'g-')
    plot(x, -V_GlueTF(V_GlueTF~=0), 'g-')
    
    plot(x, V_GlueBF(V_GlueBF~=0), 'b-')
    plot(x, -V_GlueBF(V_GlueBF~=0), 'b-')
    
    plot(x, V_GlueTW(V_GlueTW~=0), 'm-')
    plot(x, -V_GlueTW(V_GlueTW~=0), 'm-')

    plot(x, V_GlueBW(V_GlueBW~=0), 'c-')
    plot(x, -V_GlueBW(V_GlueBW~=0), 'c-')



    set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'origin');
    title('SFD vs Material Shear Failure')
    xlim([0,x(end)])
    grid on;
    grid minor
    box on;
    ylim([-max(V_fail) - 100, max(V_fail) + 100])
    yline(0)
    legend('', '', 'Matboard Shear Failure', '', 'Glue Shear (TF) Failure', '',  'Glue Shear (BF) Failure', '', 'Glue Shear (TW) Failure', '', 'Glue Shear (BW) Failure', 'FontSize', 7)
    xlabel('Position on Bridge (mm)')
    ylabel('Shear Force (N)')

    % Plots SFD vs V_buck
    subplot(2,3,3);
    hold on;
    
    plot(x, SFD);
    plot(x, V_buck, 'r-');
    plot(x, -V_buck, 'r-');    
    
    set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'origin');
    title('SFD vs Shear Buckling Failure')
    
    xlim([0,x(end)]);
    ylim([-max(V_buck) - 100, max(V_buck) + 100]);
    yline(0);
    legend('', '', 'Web Shear Buckling Failure', '', 'FontSize', 7)
    xlabel('Position on Bridge (mm)')
    ylabel('Shear Force (N)')

    grid on;
    grid minor;
    box on;
    
    % Plots BMD for Pfail
    subplot(2,3,4);
    hold on;
    plot(x, BMD);

    set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'origin', 'YDir', 'reverse');
    title("BMD from P = " + P + "N, Pfail = " +  Pfail + "N")
    xlim([0,x(end)])
    grid on;
    grid minor
    box on;
    yline(0)
    xlabel('Position on Bridge (mm)')
    ylabel('Bending Moment (Nmm)')

    % Plots BMD vs Material Moment Failures
    subplot(2,3,5);
    hold on;
    plot(x, BMD);
    plot(x, M_MatT)
    plot(x, M_MatC)
    
    set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'origin', 'YDir', 'reverse');
    title('BMD vs Material Moment Failures ')
    xlim([0,x(end)])
    grid on;
    grid minor
    box on;
    yline(0)
    legend('', 'Matboard Tension Failure', 'Matboard Compression Failure', '', 'Location', 'Northwest', 'FontSize', 7)
    xlabel('Position on Bridge (mm)')
    ylabel('Bending Moment (Nmm)')

    subplot(2,3,6);
    hold on;
    plot(x, BMD);
    plot(x, M_Buck1);
    plot(x, M_Buck2);
    plot(x, M_Buck3);
    
    set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'origin', 'YDir', 'reverse');
    title('BMD vs Material Buckling Failures ')
    xlim([0,x(end)])
    grid on;
    grid minor
    box on;
    yline(0)
    legend({'', 'Mid Flange Buckling', 'Side Flange Buclking', 'Web Compression Buckling', ''}, 'Location', 'Northwest', 'FontSize', 7)
    xlabel('Position on Bridge (mm)')
    ylabel('Bending Moment (Nmm)')

end

function [ M_MatT ] = MfailMatT( I, Y_bar, heights, SigT, BMD )
    % Calculates bending moments at every value of x that would cause a matboard tension failure
    % Input: Sectional Properties (list of 1-D arrays), SigT (material property), BMD (1-D array)
    % Output: M_MatT a 1-D array of length n
    for i = 1 : length(BMD)
        if BMD(i) > 0 % If the moment is positive, the tension failure will be at the bottom
            M_MatT(i) = SigT * I(i) / (Y_bar(i));
        elseif BMD(i) < 0 % If the moment is negative, the tension failure will be at the top
            M_MatT(i) = -SigT * I(i) / (heights(i) - Y_bar(i));
        end
    end
end

function [ M_MatC ] = MfailMatC( I, Y_bar, heights, SigC, BMD ) % Similar to MfailMatT
    % Calculates bending moments at every value of x that would cause a matboard Compression failure
    %   Input: Sectional Properties (list of 1-D arrays), SigC (material property), BMD (1-D array)
    %   Output: M_MatC a 1-D array of length n
    for i = 1 : length(BMD)
        if BMD(i) > 0 % If the moment is positive, the compression failure will be at the top
            M_MatC(i) = SigC * I(i) / (heights(i) - Y_bar(i));
            
        elseif BMD(i) < 0 % If the moment is negative, the compression failure will be at the bottom
            M_MatC(i) = -SigC * I(i) / Y_bar(i);
        end
    end
end

function [ M_Buck1 M_Buck2 M_Buck3 ] = MfailBuck( csc, bfw, bft, tfw, tft, ws, wt, heights, Y_bar, I, E, mu, BMD)
    % Calculates bending moments at every value of x that would cause a buckling failure
    % Input: Sectional Properties (list of 1-D arrays), E, mu (material property), BMD (1-D array)
    % Output: M_MatBuck a 1-D array of length n
    for c = 1:3
        if c == 1
            factor = (4 * (pi ^ 2) * E) / (12 * (1 - (mu ^ 2)));
            for i = 1:length(BMD);
                z = find(csc <= i, 1, 'last');
                if BMD(i) < 0 % if moment negative, compression on bottom
                    t = bft(z);
                    y = - Y_bar(i);
                    b =  ws(z);
                    M_Buck1(i) = (factor * ((t / b) ^ 2)) * I(i) / y;
                elseif BMD(i) > 0 % if moment positive, compression on top
                    t = tft(z);
                    y = heights(i) - Y_bar(i);
                    b = ws(z);
                    M_Buck1(i) = (factor * ((t / b) ^ 2)) * I(i) / y;
                end
            end
        elseif c == 2
            factor = (0.425 * (pi ^ 2) * E) / (12 * (1 - (mu ^ 2)));
            for i = 1:length(BMD);
                z = find(csc <= i, 1, 'last');
                if BMD(i) < 0 % if moment negative, compression on bottom
                    t = bft(z);
                    y = - Y_bar(i);
                    b = (bfw(z) - (2 * wt(z)) - ws(z)) / 2;
                    if b ~= 0 && bfw(z) ~= 0
                        M_Buck2(i) = (factor * ((t / b) ^ 2)) * I(i) / y;
                    else
                        M_Buck2(i) = 0;
                    end
                elseif BMD(i) > 0 % if moment positive, compression on top
                    t = tft(z);
                    y = heights(i) - Y_bar(i);
                    b = (tfw(z) - ws(z) - (2 * wt(z))) / 2;
                    if b ~= 0 && tfw(z) ~= 0
                        M_Buck2(i) = (factor * ((t / b) ^ 2)) * I(i) / y;
                    else
                        M_Buck2(i) = 0;
                    end
                end
            end
        else
            factor = (6 * (pi ^ 2) * E) / (12 * (1 - (mu ^ 2)));
            for i = 1:length(BMD);
                z = find(csc <= i, 1,'last');
                if BMD(i) < 0 % if moment negative, compression on bottom
                    t = wt(z);
                    y = bft(z) - Y_bar(i) ;
                    b = abs(y);
                    if b ~= 0
                        M_Buck3(i) = (factor * ((t / b) ^ 2)) * I(i) / y;
                    end
                elseif BMD(i) > 0 % if moment positive, compression on top
                    t = wt(z);
                    y = heights(i) - Y_bar(i) - tft(z);
                    b = abs(y);
                    if b ~= 0
                        M_Buck3(i) = (factor * ((t / b) ^ 2)) * I(i) / y;
                    end
                end
            end
        end
    end
end

function [ Pf, failure_mode ] = FailLoad( P, SFD, BMD, V_Mat, V_Buck, M_MatT, M_MatC, M_Buck1, M_Buck2, M_Buck3, V_GlueTF, V_GlueBF, V_GlueTW, V_GlueBW, write_to_xls, design_name )
    % Calculates the magnitude of the load P that will cause one of the failure mechanisms to occur
    %   Input: SFD, BMD under the currently applied points loads (P) (each 1-D array of length n)
    %       {V_Mat, V_Glue, … M_MatT, M_MatC, … } (each 1-D array of length n)
    %   Output: Failure Load value Pf

    SFD_new = SFD ./ abs(P);
    BMD_new = BMD ./ abs(P);

    fail_v_mat = V_Mat ./ SFD_new;
    fail_v_buck = V_Buck ./ SFD_new;

    fail_m_matt = M_MatT ./ BMD_new;
    fail_m_matc = M_MatC ./ BMD_new;
    fail_m_buck1 = M_Buck1 ./ BMD_new;
    fail_m_buck2 = M_Buck2 ./ BMD_new;
    fail_m_buck3 = M_Buck3 ./ BMD_new;

    fail_v_gluetf = V_GlueTF ./ SFD_new;
    fail_v_gluebf = V_GlueBF ./ SFD_new;
    fail_v_gluetw = V_GlueTW ./ SFD_new;
    fail_v_gluebw = V_GlueBW ./ SFD_new;

    matt = min(fail_m_matt(fail_m_matt>0));
    matc = min(fail_m_matc(fail_m_matc>0));
    buck1 = min(fail_m_buck1(fail_m_buck1>0));
    buck2 = min(fail_m_buck2(fail_m_buck2>0));
    buck3 = min(fail_m_buck3(fail_m_buck3>0));

    vmat = min(fail_v_mat(fail_v_mat>0));
    vbuck = min(fail_v_buck(fail_v_buck>0));
    
    vgluetf = min(fail_v_gluetf(fail_v_gluetf>0));
    vgluebf = min(fail_v_gluebf(fail_v_gluebf>0));
    vgluetw = min(fail_v_gluetw(fail_v_gluetw>0));
    vgluebw = min(fail_v_gluebw(fail_v_gluebw>0));

    if write_to_xls == true
        filename = design_name + '.xls';
        writematrix(design_name, filename);

        writematrix("Moment Matboard Tension Failure", filename, 'range', 'A2')
        writematrix(matt, filename, 'range', 'B2:C2');

        writematrix("Moment Matboard Compression Failure", filename, 'range', 'A3')
        writematrix(matc, filename, 'range', 'B3:C3');

        writematrix("Case 1 Buckling of Top or Bottom Flange (center part between webs)", filename, 'range', 'A4')
        writematrix(buck1, filename, 'range', 'B4:C4');
       
        writematrix("Case 2 Buckling of Top or Bottom Flange (part that sticks out past the webs)", filename, 'range', 'A5')
        writematrix(buck2, filename, 'range', 'B5:C5');

        writematrix("Case 3 Buckling of Left and Right Webs", filename, 'range', 'A6')
        writematrix(buck3, filename, 'range', 'B6:C6');
        
        writematrix("Matboard Shear Failure", filename, 'range', 'A7')
        writematrix(vmat, filename, 'range', 'B7:C7');

        writematrix("Matboard Shear Buckling Failure", filename, 'range', 'A8')
        writematrix(vbuck, filename, 'range', 'B8:C8');

        writematrix("Shear Glue Failure of Top Flange", filename, 'range', 'A9')
        writematrix(vgluetf, filename, 'range', 'B9:C9');

        writematrix("Shear Glue Failure of Bottom Flange", filename, 'range', 'A10')
        writematrix(vgluebf, filename, 'range', 'B10:C10');

        writematrix("Shear Glue Failure of Top Flange/Web Connection", filename, 'range', 'A11')
        writematrix(vgluetw, filename, 'range', 'B11:C11');

        writematrix("Shear Glue Failure of Bottom Flange/Web Connection", filename, 'range', 'A12')
        writematrix(vgluebw, filename, 'range', 'B12:C12');
    end

    results = ["Tension Failure", "Compression Failure", "Center Flange Buckling", "Protruding Flange Buckling", "Web Buckling", "V material", "vbuck", "Glue Shear tf", "Glue Shear bottom", "glue shear top width", "glue shear bottom width"];
    [Pf pp] = min([abs(matt) abs(matc) abs(buck1) abs(buck2) abs(buck3) abs(vmat) abs(vbuck) abs(vgluetf) abs(vgluebf) abs(vgluetw) abs(vgluebw)]); 
    % if NUM = 1: matt, 2: matc, 3: buck1, 4: buck2, 5: buck3, 6: vmat, 7: vbuck, 8: vgluetf, 9: vgluebf, 10: vgluetw, 11: vgluebw
    
    failure_mode = results(pp)

end

function [material_ok] = MaterialCheck( xc, tfw, tft, wh, wt, ws, bfw, bft, rw, rh, rbw, rbh)
    matboard_thickness = 1.27;

    top = tfw .* tft;
    top_glues = 2 * rh .* rw;
    bottom_glues = 2 * rbh .* rbw;
    webs = 2 * wh .* wt;
    bottom  = bfw .* bft;

    total_cross_section_width = (top + top_glues + bottom_glues + webs + bottom)/matboard_thickness;

    safety_factor = .8;

    total_area  = 0;
    for i = 2:length(xc)
        delta_x = xc(i) - xc(i - 1);
        total_area  = total_area + delta_x * total_cross_section_width(i);
    end

    total_material_area = 813 * 1016;
    material_ok = total_area < total_material_area * safety_factor;
end

function [ V_GlueTF V_GlueBF V_GlueTW V_GlueBW ] = VglueFail(I, Q, b, TauG, tft, bft, heights, csc)
    V_GlueBF = NaN(1,length(heights));
    V_GlueTF = NaN(1,length(heights));
    V_GlueTW = NaN(1,length(heights));
    V_GlueBW = NaN(1,length(heights));
    for i = 1 : length(heights)
        z = find(csc <= i, 1, 'last');
        % Top Flange Glue
        if (tft(z) / 1.27) > 1
            glue_y = ceil(heights(i) - 1.27);
            V_GlueTF(i) = (TauG * b(i, glue_y) * I(i)) / ((Q(i, floor(heights(i) - 1.27)) + Q(i, ceil(heights(i) - 1.27))) / 2);
        end
        % Bottom Flange Glue
        if (bft(z) / 1.27) > 1
            glue_y = floor(1.27);
            V_GlueBF(i) = (TauG * b(i, glue_y) * I(i)) / ((Q(i, floor(1.27)) + Q(i, ceil(1.27)))/2);
        end

        % Top web/flange glue (there will always be a top flange therefore no if statement needed) 
        %   (note: rounds up heights(i) - tft(i) in order to get smaller b value and larger Q value in order to get lower V_GlueTW for saftey
        V_GlueTW(i) = (TauG * b(i, floor(heights(i) - tft(z))) * I(i)) / ((Q(i, floor(heights(i) - tft(z))) + Q(i, ceil(heights(i) - tft(z)))) / 2);
        
        % Bottom web/flange glue (note: rounds up bft(i) in order to get smaller b value and larger Q value in order to get lower V_GlueBW for saftey
        if bft(z) > 0  
            V_GlueBW(i) = (TauG * b(i, ceil(bft(z))) * I(i)) / ((Q(i, floor(bft(z))) + Q(i, ceil(bft(z))))/2);
        end
    end
end