%% 0. Initialize Parameters
n = 1250;               % Number of locations to evaluate bridge failure
L = 1250;               % Length of bridge
x = linspace(0, L, n)   % Define x coordinate
SFD_PL = zeros(1, n);   % Initialize SFD(x)

%% 1. Point Loading Analysis (SFD, BMD)
P = -318;
[SFD_PL, BMD_PL] = ApplyPL(550, P, x, SFD_PL); % Construct SFD, BMD
[SFD_PL, BMD_PL] = ApplyPL(L, P, x, SFD_PL); % Construct SFD, BMD
SFD_PL = Make_SFD(SFD_PL)
plot(x, SFD_PL)


% %% 2. Define cross-sections
%     % There are many (more elegant ways) to construct cross-section objects
% xc = [0 550 L];         % Location, x, of cross-section change
% bft = [100 100 100];    % Top Flange Width
% tft = [2.54 2.54 2.54]; % Top Flange Thickness
% hw = [100 120 100];     % Web Height
% tw = [1.27 1.27 1.27];  % Web Thickness (Assuming 2 separate webs)
% bfb = [80 80 80];       % Bottom Flange Width
% tfb = [1.27 1.27 1.27]; % Bottom Flange Thickness
% a = [400 400 400];      % Diaphragm Spacing
%     % Optional but you need to ensure that your geometric inputs are correctly implemented
% VisualizeBridge( {CrossSectionInputs} ); 

% %% 3. Define Material Properties
% SigT = 30;  % Matboard: Tensile strength
% SigC = 6;   % Matboard: Compressive strenght
% E = 4000;   % Matboard: Young's Modulus`
% TauU = 4;   % Matboard: Shear strength
% TauG = 2;   % Glue: Shear strength 
% mu = 0.2;   % Poisson's Ratio

% %% 4. Calculate Failure Moments and Shear Forces
% V_Mat = Vfail({CrossSectionInputs}, TauU);
% V_Glue = VfailGlue({CrossSectionInputs}, TauU);
% V_Buck = VfailBuck({CrossSectionInputs}, E, mu );
% M_MatT = MfailMatT({CrossSectionInputs}, SigT);
% M_MatC = MfailMatC({CrossSectionInputs}, SigC);
% M_Buck1 = MfailBuck({CrossSectionInputs}, E, mu, 1 );
% M_Buck2 = MfailBuck({CrossSectionInputs}, E, mu, 2 );
% M_Buck3 = MfailBuck({CrossSectionInputs}, E, mu, 3 );

%% 4.7 Calculate Failure Load
% Pf = FailLoad(P, SFD_PL, BMD_PL, V_Mat, V_Glue, V_Buck, M_MatT, M_MatC, M_Buck1, M_Buck2, M_Buck3);

%% Visualization
% VisualizePL(x, P, SFD_PL, BMD_PL, V_Mat, V_Glue, V_Buck, M_MatT, M_MatC, M_Buck1, M_Buck2, M_Buck3, Pf);

%% 5. Curvature, Slope, Deflections
% Defls = Deflections(x, BMD_PL, I, E);
%%
function [ SFD, BMD ] = ApplyPL( xP, P, x, SFD )
% Constructs SFD and BMD from application of 1 Point Load. Assumes fixed location of supports
%   Input: location and magnitude of point load. The previous SFD can be entered as input to construct SFD of multiple point loads
%   Output: SFD, BMD both 1-D arrays of length n
    % 1060 is the distance between supports
    
    xP_B = 1060;
    P_B = (-P * xP / xP_B) + SFD(xP_B) - SFD(xP_B-1);
    P_A = -P - P_B;
    SFD(1) = P_A;
    SFD(xP) = P;
    SFD(xP_B) = P_B;
    % for i = 2:length(SFD)
    %     if SFD(i) == SFD(i-1):
    %         SFD(i) = 0 
    %     % modifier = 0
    %     % if i == xP
    %     %     modifer = P
    %     % end
    %     % if i == xP_B 
    %     %     modifer = modifer + P_B
    %     % end
    %     SFD(i) = SFD(i) + SFD(i - 1);
    % end
    % SFD = @x 
    % for i = x
    %     BMD = integral(SFD(i), x(1), x(i))
    % end
    BMD = 0
end

function [SFD] = Make_SFD(SFD)
    for i = 2:length(SFD)
        SFD(i) = SFD(i) + SFD(i - 1);
    end
end 

%%

% for i == 1:length(SFD)
%     ex_SFD= SFD(i) 
%     if ex_SFD == SFD(i+1):
%         ex_SFD = SFD(i+1)
%         SFD(i+1) = 0

function [ ] = VisualizeBridge( {Geometric Inputs} )
% Optional. Provides a graphical interpretation of user geometric inputs
end
%%
function [ {Sectional Properties} ] = SectionProperties( {Geometric Inputs} )
% Calculates important sectional properties. Including but not limited to ybar, I, Q, etc.
%   Input: Geometric Inputs. Format will depend on user
%   Output: Sectional Properties at every value of x. Each property is a 1-D array of length n


end
%%

function [ V_fail ] = Vfail( {Sectional Properties}, TauU )
% Calculates shear forces at every value of x that would cause a matboard shear failure
%   Input: Sectional Properties (list of 1-D arrays), TauU (scalar material property)
%   Output: V_fail a 1-D array of length n
    I = {Sectional Properties};
    b = {Sectional Properties};
    Qcent = {Sectional Properties};
    V_fail = TauU .* I .* b ./ Qcent;
end
%%

function [ V_Buck ] = VfailBuck( {Sectional Properties}, E, mu )
% Calculates shear forces at every value of x that would cause a shear buckling failure in the web
%   Input: Sectional Properties (list of 1-D arrays), E, mu (material property)
%   Output: V_Buck a 1-D array of length n
end
%%

function [ M_MatT ] = MfailMatT( {Sectional Properties}, SigT, BMD )
% Calculates bending moments at every value of x that would cause a matboard tension failure
%   Input: Sectional Properties (list of 1-D arrays), SigT (material property), BMD (1-D array)
%   Output: M_MatT a 1-D array of length n
    [I, ybot, ytop] = {Sectional Properties};
    for i = 1 : length(x)
        if BMD(i) > 0 % If the moment is positive, the tension failure will be at the bottom
            M_MatT(i) = SigT * I(i) / ybot(i);
        elseif BMD(i) < 0 % If the moment is negative, the tension failure will be at the top
            M_MatT(i) = -SigT * I(i) / ytop(i);
        end
    end
end
%%
function [ M_MatT ] = MfailMatC( {Sectional Properties}, SigC, BMD ) % Similar to MfailMatT

end
%%

function [ M_Buck ] = MfailBuck( {Sectional Properties}, E, mu, BMD )
% Calculates bending moments at every value of x that would cause a buckling failure
%   Input: Sectional Properties (list of 1-D arrays), E, mu (material property), BMD (1-D array)
%   Output: M_MatBuck a 1-D array of length n
end
%%

function [ Pf ] = FailLoad( P, SFD, BMD, V_Mat, V_Buck, M_MatT, M_MatC, M_Buck1, M_Buck2, M_Buck3 )
% Calculates the magnitude of the load P that will cause one of the failure mechanisms to occur
%   Input: SFD, BMD under the currently applied points loads (P) (each 1-D array of length n)
%       {V_Mat, V_Glue, … M_MatT, M_MatC, … } (each 1-D array of length n)
%   Output: Failure Load value Pf
end
%%
function [] = VisualizePL(x, SFD, BMD, V_Mat, V_Buck, M_MatT, M_MatC, M_Buck1, M_Buck2,…, Pf)
% Plots all outputs of design process
end
%%

function [ Defls ] = Deflections( x, BMD, I, E )
% Calculates deflections
%   Input: I(1-D arrays), E (material property), BMD (1-D array)
%   Output: Deflection for every value of x (1-D array) or for the midspan only
end