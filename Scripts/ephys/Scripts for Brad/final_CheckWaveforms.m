% Check the spike unit waveforms (modified to save individual plots)

clear; clc; close all;
working_dir = 'C:\Users\Lab User\Desktop\Ephys Data';

Mice = {'061519_A','061519_B','061519_C','061619_C','061719_A','061719_B','061719_C','061819_A','061819_B','061819_C','061819_D','062019_A','062019_B','062019_C','062019_D','062119_A','062119_B','063019_A','063019_B','063019_C','070119_A','070119_B','070119_C','070619_A','070619_B','070619_C','070719_A','070719_B','070719_C','070719_D','070819_A','070819_B','070819_C','072419_A'};
Regions = {'iCT','cCT','Som','ZI'};
Frequencies = {'6','10','20'};
Recordings = {'A','B','C','D'};
lineStyles = {'-r','-r','-b'};
check_unsorted = 0; % Whether you want to look at the unsorted units as well
plotCountLimit = 50; % Limit the number of spikes to plot for performance reasons, set to a really high number (e.g. 100,000) to disable
i_png = 1;
% % bad_waveforms = [37,40,41,42,43,44,46,48,50,51,80,83,84,89,90,92,93,95,96,98,121,135,136,139,141,142,143,144,146,147,148,149,150,151,187,188,189,190,191,192,193,196,198,199,200,201,203,204,205,208,209,210,277,288,338,346,347,349,353,354,355,356,358,359,361,367,368,369,370,384,402,430,436,441,444,457,458,463,465,530,531,533,537,539,566,567,568,570,573,575,576,577,580,582,583,584,586,587,590,591,592,597,598,603,604,605,606,608,609,610,611,614,615,616,631,641,644,648,655,656,657,658,659,660,661,663,664,688,697,700,703,706,708,731,733,735,741,742,746,749,752,755,761,765,787,833,837,847,848,863,864,872,875,892,896,899,930,931,937,939,941,956,969,973,974,978,979,990,994,1012,1015,1019,1020,1030,1031,1035,1040,1041,1054,1056,1057,1113,1148,1235,1248,1253,1273,1278,1279,1281,1678,1766,1836,1928,1931,1960,1965,1971,1978,1999,2005,2006,2017,2044,2047,2055,2060,2084,2085,2090,2091,2094,2141,2160,2175,2207,2229,2232,2251,2264,2265,2275,2294,2295,2304,2308,2338,2343,2373,2385,2386,2432,2469];
bad_waveforms = [18,23,24,25,26,29,30,32,38,39,45,46,47,48,52,69,70,72,76,77,78,79,81,82,84,89,90,94,95,96,98,109,118,119,141,142,143,144,145,146,149,169,171,173,174,182,183,187,191,196,198,199,200,201,204,205,206,207,208,209,210,265,266,267,268,269,270,273,277,282,288,289,290,346,348,353,354,355,363,365,366,384,392,393,394,395,396,402,429,430,431,432,433,434,435,436,437,438,442,443,455,456,457,458,459,460,461,462,463,497,512,513,514,515,516,535,537,538,539,549,577,584,585,587,588,589,590,594,595,603,604,606,607,609,617,635,636,637,638,639,640,642,645,655,656,657,658,659,660,661,665,681,682,683,684,685,690,691,697,699,701,702,705,706,726,729,733,745,759,764,766,767,768,775,780,781,782,784,788,795,800,801,810,811,812,813,827,828,830,832,833,837,858,861,870,871,879,880,881,882,886,887,899,901,902,908,925,929,931,933,934,946,949,956,970,971,972,979,988,989,991,1004,1005,1006,1013,1014,1015,1019,1020,1021,1063,1113,1114,1115,1122,1123,1124,1126,1138,1143,1199,1232,1233,1234,1236,1237,1248,1249,1250,1251,1253,1256,1275,1280,1281,1282,1283,1309,1310,1366,1417,1418,1419,1431,1432,1459,1466,1485,1486,1531,1539,1632,1634,1658,1659,1660,1677,1678,1679,1709,1737,1738,1740,1741,1742,1761,1762,1764,1765,1766,1767,1769,1772,1780,1781,1782,1810,1816,1817,1849,1852,1868,1869,1870,1871,1878,1885,1887,1921,1929,1931,1939,1940,1959,1961,1968,1969,1971,1972,1973,1974,1977,1978,1979,1980,1981,1982,1985,1989,1995,1996,1997,1998,2001,2002,2011,2044,2046,2047,2051,2053,2060,2063,2066,2077,2078,2084,2085,2086,2089,2091,2096,2097,2123,2128,2129,2130,2131,2142,2143,2144,2146,2148,2150,2154,2159,2162,2171,2172,2177,2221,2222,2223,2224,2225,2227,2252,2254,2256,2263,2265,2267,2273,2288,2294,2296,2313,2338,2340,2373,2382,2384,2416,2417,2432,2436,2467,2482,2484,2496,2542,2543,2544,2605,2607,2612,2613,2624,2643,2649,2652,2668,2689];

%%
Fig1 = figure('Position', get(0, 'Screensize'),'visible','off');
hold on;
for i_mouse = 1:length(Mice)
    mouse = Mice{i_mouse};
    for i_region = 1:length(Regions)
        region = sprintf('BF_%s',Regions{i_region});
        
        msg = sprintf('\n%s %s ',mouse,region(4:end));
        fprintf(msg);
        
        Waveforms = cell(3,32,26);
        clear wf_index wf_idx;      
        for i_freq = 1:length(Frequencies)
            frequency = Frequencies{i_freq};
            for i_recording = 1:length(Recordings)
                recording = Recordings{i_recording};
                Filename = fullfile(working_dir,mouse,sprintf('%s_%sHz_%s.mat',region,frequency,recording));
                if exist(Filename,'file')
                    % Concatenate the spike waveforms from across the 4 recordings
%                     fprintf('Loading: %s\n',Filename);
                    load(Filename);
                    for i_channel = 1:size(Waveform_mV,2)
                        for i_unit = 1:size(Waveform_mV,3)
                            temp = Waveform_mV{1,i_channel,i_unit};
                            Waveforms{i_freq,i_channel,i_unit} = [Waveforms{i_freq,i_channel,i_unit};temp];
                        end
                    end
                end
            end
        end
        
        % Plot the spike waveforms
        wf_index = ~cellfun('isempty',Waveforms);
        [wf_idx(1,:,:), wf_idx(2,:,:), wf_idx(3,:,:)] = ind2sub(size(wf_index),find(wf_index));
        for i_wf_idx = 1:length(wf_idx)
            i_unit = wf_idx(3,i_wf_idx);
            if i_unit > 1 || check_unsorted == 1
                i_freq = wf_idx(1,i_wf_idx);
                i_channel = wf_idx(2,i_wf_idx);
                frequency = Frequencies{i_freq};
                
                if ismember(i_png,bad_waveforms) && ~isequal(frequency,'20')
                    fprintf('[%d,%d];',i_channel,i_unit-1);
                end

% % %                 fprintf('Plotting/Saving: %s %s_%sHz - Ch.%d Unit %d\n',mouse,region,frequency,i_channel,i_unit-1);
% % %                 data = Waveforms{i_freq,i_channel,i_unit};
% % %                 for i=1:min(plotCountLimit,size(data,1))
% % %                     plot1 = plot(data(i,:)',lineStyles{i_freq});
% % %                     plot1.Color(4) = 0.3;
% % %                 end
% % %                 clear data;
                title({sprintf('%s',mouse),sprintf('%s_%sHz',region,frequency),sprintf('Channel %d Unit %d',i_channel,i_unit-1)},'Interpreter','none');
                xlim([1 32]);
                xticks([1:4:32]);
                xlabel('25us steps');
                ylabel('mV');
                if i_wf_idx == length(wf_idx) || ~isequal([i_channel,i_unit],[wf_idx(2,i_wf_idx+1),wf_idx(3,i_wf_idx+1)])
% % %                     saveas(Fig1,fullfile(working_dir,'Figures_Indiv','CheckWaveforms',sprintf('%d.png',i_png)),'png');
                    i_png = i_png+1;
%                     pause;
                    cla;
                end
            end
        end
    end
end
close all;