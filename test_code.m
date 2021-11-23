%% 0. Initialize Parameters
n = 1250; % Number of locations to evaluate bridge failure
L = 1250; % Length of bridge
x = linspace(0, L, n); % Define x coordinate
P = zeros(1,n); % Initializes Loads 

%% 1. Point Loading Analysis (SFD, BMD)
load = -318;
[SFD_PL, BMD_PL, P] = ApplyPL(550, load, x, P); % Construct SFD, BMD
[SFD_PL, BMD_PL, P] = ApplyPL(L, load, x, P); % Construct SFD, BMD

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
ws = [77.46 77.46 77.46] % Web Spacing

% Optional but you need to ensure that your geometric inputs are correctly implemented
VisualizeBridge(xc, bft, tft, hw, tw, ws, bfb, tfb ); 
%% 3. Define Material Properties
SigT = 30;
SigC = 6;
E = 4000;
TauU = 4;
TauG = 2;
mu = 0.2;

function [ SFD, BMD, Loads ] = ApplyPL( xP, P, x, Loads )
% Constructs SFD and BMD from application of 1 Point Load. Assumes fixed location of supports
%   Input: location and magnitude of point load. The previous Loads can be entered as input to 
%       construct SFD of multiple point loads
%   Output: SFD, BMD, Loads, all 1-D arrays of length n
    global n
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
    set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
end

function [ ] = VisualizeBridge( csc, tfw, tft, wh, wt, ws, bfw, bft )
    % Provides a graphical interpretation of user geometric inputs
    %    Input: cross section changes x cordinates, top flange width, top flange thinckness, web height, web thickness, web spacing,
    %        bottom flange width, bottom flange thickness
    %    Output: none
    
    % initial cooridates of cross section
    start_cord = [5 5];
    for i = 1:length(csc)
        % Calculates basic coordinates of features
        top_coordinate = start_cord(2) + bft(i) + wh(i) + tft(i);
        left_web_cord = ((tfw(i) - ws(i)) / 2) - wt(i);
        bot_flang_cord = start_cord(1) + ((tfw(i) - bfw(i)) / 2);
        
        % Draws rectangles for features
        figure;
        rectangle('position', [ start_cord(1) (top_coordinate - tft(i)) tfw(i) tft(i) ]);
        rectangle('position', [ (start_cord(1) + left_web_cord) (start_cord(2) + bft(i)) wt(i) wh(i)]);
        rectangle('Position', [ (start_cord(1) + tfw(i) - (left_web_cord) - wt(i)) (start_cord(2) + bft(i)) wt(i) wh(i)]);
        rectangle('position', [ bot_flang_cord start_cord(2) bfw(i) bft(i)]);

        % Defines axis
        axis([0 (2 * start_cord(1) + tfw(i)) 0 (2 * start_cord(2) + top_coordinate)]);
    end
end