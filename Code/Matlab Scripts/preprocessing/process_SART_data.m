function process_SART_data(filename)

%% READ DATA
SART_data = readtable(filename);

%% PRELIMINARIES
EV = unique(SART_data.EV);
subjid = extractBetween(filename, '/mus', '/');
music_session = unique(SART_data.MusicSession); 
digit_displayed = unique(SART_data.DigitDisplayed);
white_idx = strcmp(SART_data.DigitDisplayed, 'white');
black_idx = strcmp(SART_data.DigitDisplayed, 'black');

%Flag: How should we define reaction time(?)
%% COMPUTE METRICS

%Logic: Per song we can estimate the mean response time and accuracy in two
%ways: overall and per digit displayed. Therefore, we need a nx4 table (song, rt, acc, split)

SART_summary = cell(numel(EV)*3, 4);
counter = 1;
for song = EV' %must be horizontal
    current_data = strcmp(SART_data.EV, song);
    %White
    SART_summary(counter,1) = song;
    SART_summary(counter,2) = {mean(str2double(SART_data.Accuracy(current_data & white_idx)), "omitnan")};
    SART_summary(counter,3) = {mean(str2double(SART_data.ResponseTime(current_data & white_idx)), "omitnan")};
    SART_summary(counter,4) = {'white'};
    counter = counter + 1;

    %Black
    SART_summary(counter,1) = song;
    SART_summary(counter,2) = {mean(str2double(SART_data.Accuracy(current_data & black_idx)), "omitnan")};
    SART_summary(counter,3) = {mean(str2double(SART_data.ResponseTime(current_data & black_idx)), "omitnan")};
    SART_summary(counter,4) = {'black'};
    counter = counter + 1;

    %All
    SART_summary(counter,1) = song;
    SART_summary(counter,2) = {mean(str2double(SART_data.Accuracy(current_data)), "omitnan")};
    SART_summary(counter,3) = {mean(str2double(SART_data.ResponseTime(current_data)), "omitnan")};
    SART_summary(counter,4) = {'all'};
    counter = counter + 1;
end
SART_summary = cell2table(SART_summary, "VariableNames",{'EV', 'Accuracy', 'Reaction_Time', 'Digit'});
writetable(SART_summary, ['/Users/camilocastelblanco/Dropbox (Dartmouth College)/ieeg-music/data/mus' cell2mat(subjid) '/ieeg/mus' cell2mat(subjid) '_SART_summary.xlsx'])

end

%% NOTES:
%Reaction Time
%  We are still recording the response time past when the circle
%  circle was onscreen (during the crosshair stimulus after). Users
%  will be asked to press space when the circle is white and refrain
%  from pressing space when the circle is black.

%responseTime = time.time() - appearanceTime

%Accuracy: