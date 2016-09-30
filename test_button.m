figh = figure;

global IS_ABORTED;

IS_ABORTED = false;


btn = uicontrol('style', 'pushb', 'string', 'Abort', ...
                'callback', @doAbort);

% simulate a long-running computation with PAUSE
%while 1
%    checkForAbort();
%    disp(i); pause(1);
%end

%% =====


