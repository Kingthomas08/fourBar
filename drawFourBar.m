% drawFourBar(l12, l23, l34, l14, l35, plotTitle)
%
%   Takes input of four-bar linkage link lengths and plot title, then draws
%   linkage's ranges of motion. Can receive dimensions of arbitrary number
%   of sets of dimensions, but all vectors must have length 1 or n, with
%   'n' being length of longest vector provided.
%
%   'l12' is crank, 'l23' is coupler, 'l34' is output, 'l14' is ground, and
%   'l35' is extension to coupler, as shown here:
%
%       2---3--5
%      /   /
%   ==1===4==
%
%   'l35' and 'plotTitle' are optional.
%
%   Written by Travis Llado (travis@travisllado.com). He is aware that it
%   has bugs related to inversions. Last edited 2018-03-18.

function drawFourBar(l12, l23, l34, l14, l35, plotTitle)
    % Check which argument were received
    if ~exist('l12', 'var') % If no arguments provided, display demo
        l12 = 100*(3-sqrt(7));
        l23 = 100;
        l34 = 100;
        l14 = 100*(5-sqrt(7))/3;
        l35 = 100;
        plotTitle = "Chebyshev's Lambda Mechanism";
    end
    
    if exist('l35', 'var')
        ll = max([length(l12) length(l23) length(l34) length(l14) length(l35)]);
        l25(1:ll) = l23 + l35;
    else
        ll = max([length(l12) length(l23) length(l34) length(l14)]);
    end
    
    if ~exist('plotTitle', 'var')
        plotTitle = "Four Bar Mechanism Range";
    end

    l12(1:ll) = l12;
    l23(1:ll) = l23;
    l34(1:ll) = l34;
    l14(1:ll) = l14;

    % Specify range to rotate input/th2
    res = 0.01; % rad
    th2 = 0:res:2*pi;   % full circle
    rr = length(th2);

    % Initialize intermediate variables
    x = zeros(5, rr, ll);
    y = x;

    % Calculate forward kinematics
    for ii = 1:ll   % For each mechanism ...
        for jj = 1:rr   % For each position of this mechanism ...
            x(2,jj,ii) = l12(ii)*cos(th2(jj));
            y(2,jj,ii) = l12(ii)*sin(th2(jj));

            l24 = sqrt(l12(ii)^2 + l14(ii)^2 - 2*l12(ii)*l14(ii)*cos(th2(jj)));
            th142 = acos((l14(ii)^2 + l24^2 - l12(ii)^2)/(2*l14(ii)*l24)) ...
                *sign(-y(2, jj,ii));
            th243 = acos((l24^2 + l34(ii)^2 - l23(ii)^2)/(2*l24*l34(ii)));
            th143 = th142 + th243;

            x(3,jj,ii) = l14(ii) + l34(ii)*cos(th143 + pi);
            y(3,jj,ii) = l34(ii)*sin(th143 + pi);
            
            x(4,jj,ii) = l14(ii);
            y(4,jj,ii) = 0;
            
            if exist('l35', 'var')
                x(5,jj,ii) = x(2,jj,ii) + (x(3,jj,ii)-x(2,jj,ii))*l25(ii)/l23(ii);
                y(5,jj,ii) = y(2,jj,ii) + (y(3,jj,ii)-y(2,jj,ii))*l25(ii)/l23(ii);
            end
        end
    end

    % Plot results
    figure
    hold on
    nn = 30;    % random number so our mechanism is at an odd angle, not flat
    
    for ii = 1:ll   % Draw ground links first so they're background
        plot([x(1,nn,ii)-x(4,nn,ii)*0.2 x(4,nn,ii)*1.2], [y(1,nn,ii) y(4, ...
            nn, ii)], 'color', [.6 .25 .25], 'LineWidth', 5);
    end
    
    for ii = 1:ll   % For each mechanism ...
        % Draw links
        plot(x(1:4,nn,ii), y(1:4,nn,ii), 'k-o', 'LineWidth', 5);
        if exist('l35', 'var')
            plot([x(2,nn,ii) x(5,nn,ii)], [y(2,nn,ii) y(5,nn,ii)], 'k-o', 'LineWidth', 5);
        end
        
        % Draw joints
        plot(x(2,nn,ii), y(2,nn,ii), 'ro', 'LineWidth', 5);
        plot(x(3,nn,ii), y(3,nn,ii), 'go', 'LineWidth', 5);
        if exist('l35', 'var')
            plot(x(5,nn,ii), y(5,nn,ii), 'bo', 'LineWidth', 5);
        end
        
        % Draw range curves
        plot(x(2,:,ii)', y(2,:,ii)', 'r',  'LineWidth', 2);
        plot(x(3,:,ii)', y(3,:,ii)', 'g',  'LineWidth', 2);
        if exist('l35', 'var')
            plot(x(5,:,ii)', y(5,:,ii)', 'b',  'LineWidth', 2);
        end
        
        % Format plot
        daspect([1 1 1]);
        xlabel('x');
        ylabel('y');
        title(plotTitle);
    end
    hold off
end

% End of file