%mus28-Task: trial definition
function [trl, event] = mus28_trialDefinition(cfg)   
    hdr              = ft_read_header(cfg.dataset);
    event            = ft_read_event(cfg.dataset, 'detectflank', 'down', 'chanindx', find(ismember(hdr.label, 'DC3')));
    %event            = ft_read_event(cfg.dataset, 'detectflank', 'any', 'chanindx', find(ismember(hdr.label, 'DC3')));
    idx              = [];
    for e = 1:numel(event) %numel = number of array elements
      if ~isequal(event(e).type, 'DC3')
        idx = [idx e]; % events to be tossed
      end
    end
    event(idx) = [];
    %trigs = [event.sample]';

    % determine the number of samples before and after the trigger
    pretrig  =  round(cfg.trialdef.pre  * hdr.Fs);
    postrig = round(cfg.trialdef.post * hdr.Fs);
    trl = [];
    eventStart = event([event().value] == 0); %I think its down
    eventEnd = event([event().value] == 1); %I think its up
    eventSize = size(eventStart);
   
    for i = 1:eventSize(2) %ensure physiologically plausible stuff
        trlbegin = eventStart(i).sample - pretrig; % start of segment (MUST HAVE -pretrig)
        trlend = eventEnd(i).sample + postrig; % end of segment
        offset = -pretrig; % how many samples prestimulus 
        rt = (eventEnd(i).sample-eventStart(i).sample)/hdr.Fs; %Trial Duration
        newtrl   = [trlbegin trlend offset rt];      
        trl      = [trl; newtrl]; %#ok<AGROW>
    end
    %If you match the offset and the pre setting, t=0 will be the marker itself

    %Adding real time (i.e., seconds) 
    trl(:,5:6) = trl(:,1:2)/hdr.Fs;

    x = ones(1, eventSize(2));
    figure('Name','Stimulation duration during each trial','position',[0 0 600 500],'color','w','GraphicsSmoothing','on');
    ylabel('Time (s) between DC-up and DC-down flanking')
    swarmchart(x, trl(:,4),'^','filled', 'ColorVariable',[0.4660 0.6740 0.1880])
    text_annotation = sprintf('Mean = %f, std = %f', mean(trl(:,4)), std(trl(:,4)));
    text(0.6,1,text_annotation,'Units','normalized', 'FontSize',15)
    set(gca,'FontSize',18)
end

