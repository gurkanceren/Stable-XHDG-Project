% function readyforprint([W H],fntsize,fgc,bgc,lwidth)
%
% this script resizes and reformats images for inclusion in documents,
% presentations or web pages.
%
% What this function modifies;
% - paper size (size of output document when output usin 'print')
% - font sizes
% - colours of axes, legend boxes and colour bar boxes
% - widths of lines in the plot, including the axes, legend and colour bar
% boxes and 
%
% What this function does not modify; 
% - marker sizes
% 
% Inputs are;
% [W H];    width and height of the image (inches)
% fntsize;	Font size (pts). All text will be set to this size.
% fgc;      foreground colour, either as 'r' or [0.5 0.5 0.5] style. This
%			colour is for lines in the foreground (data)
% bgc;      background colour, either as 'r' or [0.5 0.5 0.5] style. This
%			colour is for axes, legend boxes, etc.
% lwidth;	Width of lines in figures. use 0.5 for journal plots, 2 for
%			overhead presentations, etc.
%
% example;
% readyforprint([6 3],10,'k','w',0.5)
%
% Any or all inputs can be left empty, in which case a default set is used.
% This is equivalent to
% > readyforprint([6 3],10,'k','w',0.5)
% These defaults are approximately suitable for a full-width journal image.
%
% Suggested image sizes are;
% [7 4];        for a picture spanning a whole A4 or letter page
% [3.33 3.33];  for a picture half the width of a page (subfigure)
%
% Note that typically used paper sizes are;
% [8.5 11];     letter paper
% [8.27 11.69]; A4 paper
% [5.83 8.27];  A5 paper
%
% and in fact the paper size should be fitted to the text width, not the
% paper width. This is the user's concern, and may be specified by a
% Journal or University, etc, rather than the user having any say in the
% matter... A consistent use of the same widths or heights should help any
% document look much smarter.
%
% follow this command with something like;
% print('-depsc2','-loose','GLR_fig2a')     % for a colour eps
% print('-djpeg','GRL_fig2a')               % for a colour jpg
%
% This version written by Andy, January 2010 for submission to MatlabCentral

function readyforprint(wh,fs,fgcin,bgcin,bs,lw,ms)

%% input check

% first deal with empty ones
if isempty(wh)
	wh = [6 3];
elseif(prod(size(wh))) == 1
	disp('')
	disp('Please set the required width and height (in inches)!')
	disp('')
	error('Width and height not defined in inputs to size_for_pub!')
	return
end
if isempty(fs)
	fs = 10;
end
if isempty(fgcin)
	fgc = 'k';
else
	fgc = fgcin;
end
if isempty(bgcin)
	bgc = 'w';
else
	bgc = bgcin;
end
if isempty(bs)
	bs = 0.5;
end
if isempty(lw)
	lw = 0.5;
end
if isempty(ms)
	ms = 0.5;
end

%specify some figure properties
figure_props.InvertHardcopy = 'off';

% specify the paper properties
paper_props.Color = bgc;
paper_props.PaperUnits ='inches';
paper_props.PaperOrientation ='Portrait';
paper_props.PaperPosition = [0 0 wh(1) wh(2)];
paper_props.PaperPositionMode = 'manual';
paper_props.PaperSize = [wh(1) wh(2)];

% specify the axis properties
axes_props.FontSize = fs;
axes_props.LineWidth = bs;
axes_props.Units = 'normalized';
axes_props.Color = bgc;
axes_props.Xcolor = fgc;					% foreground colours
axes_props.Ycolor = fgc;
axes_props.Zcolor = fgc;

% specify the line properties
line_props.LineWidth = lw;
line_props.MarkerSize = ms;
line_props.Color = fgc;

% specify the text properties
% text_props.Color = fgc;					% text colour
text_props.BackgroundColor = 'none';		% fill colour
text_props.EdgeColor = bgc;				% box colour
text_props.FontSize = fs;

% specify the legend properties
legend_props.FontSize = fs;
legend_props.Xcolor = fgc;
legend_props.Ycolor = fgc;
legend_props.TextColor = fgc;
legend_props.Color = bgc;

% specify the title properties
title_props.FontSize = fs;
title_props.Color = fgc;
title_props.BackgroundColor = bgc;

% specify the label properties
label_props.Color =  fgc;
label_props.FontSize = fs;

%%%%%%%%%%%%%%% END OF INPUTS %%%%%%%%%%%%%%%

% tell the user what's happening
disp('modifying plot...');

% now need to concatenate these inputs into a useful format; prefered
% option for setting properties is a cell array of strings for the
% properties, and a cell array of properties

%%%%%%%%%%%%%%% START MODIFICATION %%%%%%%%%%%%%%%


%% PAPER PROPERTIES %%%%%%%%%%%%%%%
set(gcf,fieldnames(figure_props)',struct2cell(figure_props)')
set(gcf,fieldnames(paper_props)',struct2cell(paper_props)')


%% AXES
myaxes = findobj(gcf,'Type','Axes');

for myaxe = 1:length(myaxes)
	% set the focus onto the current axis
	axes(myaxes(myaxe));
	
	% get the current foreground and background colours
	old_fgc = get(gca,'Xcolor');
	old_bgc = get(gca,'Color');
	
	% get the current legend data
	if isempty(findobj(legend)) == 0
		legend_data = get(legend,'UserData');
		legend_plot_handles = legend_data.handles;
		legend_string = get(legend,'String');
		legend_location = get(legend,'Location');
		legend_box_state = get(legend,'Box');
		legend_interpreter = get(legend,'Interpreter');
	end
	
	% get the current title...
	if ~isempty(get(get(gca,'Title'),'string'))
		title_string =get(get(gca,'Title'),'String');
		title_interpreter = get(get(gca,'Title'),'Interpreter');
	end
	
	%%%%%%%%%%%%%%% UPDATE AXES PROPERTIES %%%%%%%%%%%%%%%
	% now need to identify the axes objects and try to change these
	set(gca,fieldnames(axes_props)',struct2cell(axes_props)');
	% note that the grid inherits the colours either of the axes
	% or the line properties, and we don't worry about setting this
	% explicitly.
	
	%%%%%%%%%%%%%%% UPDATE LABEL PROPERTIES %%%%%%%%%%%%%%%
	try
		% get the properties of the labels that we want to change
		mylabels = fieldnames(label_props);
		for mylabel = 1:length(mylabels)
			set(get(gca,'xlabel'), mylabels{mylabel}, label_props.(char(mylabels(mylabel,:))) );
			set(get(gca,'ylabel'), mylabels{mylabel}, label_props.(char(mylabels(mylabel,:))) );
		end
	end
	
	%%%%%%%%%%%%%%% LINE PROPERTIES %%%%%%%%%%%%%%%
	
	% find all line items on this axis
	mylines = findobj(gca,'Type','line','-or','Type','hggroup','-or','Type','Patch');
	for myline = 1:length(mylines)
		try
            set(mylines(myline),'LineWidth',line_props.LineWidth);
            set(mylines(myline),'MarkerSize',line_props.MarkerSize);
        end
	end
	
	% now finding text items on this axis; note that we use gca not gcf
	% as this explicitly ignores the legend. The legend is set and
	% controlled seperately
	mytexts = findobj(gca,'Type','text');
	for mytext = 1:length(mytexts)
		try
			if strcmp(get(mytexts(mytext),'BackgroundColor'),'none')
				text_bgc = 'none';
			end
			set(mytexts(mytext),fieldnames(text_props)',struct2cell(text_props)');
			if strcmp(text_bgc,'none')
				set(mytexts(mytext),'BackGroundColor','none')
			end
			if strcmp(get(mytexts(mytext),'EdgeColor'),'none')
				text_bgc = 'none';
			end
			set(mytexts(mytext),fieldnames(text_props)',struct2cell(text_props)');
			if strcmp(text_bgc,'none')
				set(mytexts(mytext),'EdgeColor','none')
			end
		end
	end
	
	%%%%%%%%%%%%%%% TITLE PROPERTIES %%%%%%%%%%%%%%%
	if ~isempty(get(get(gca,'Title'),'string'))
		% and do the mods
		set(get(gca,'Title'),fieldnames(title_props)',struct2cell(title_props)');
	end
	
	%%%%%%%%%%%%%%% LEGEND PROPERTIES %%%%%%%%%%%%%%%
	
	% now finding legend items on this axes as well
	try
		if isempty(findobj(legend)) == 0
			%refresh the legend
			legend('off')
			legend(legend_plot_handles,legend_string)
			set(legend,'Location',legend_location,...
				'Interpreter',legend_interpreter,...
				'Box',legend_box_state)
			% and do the mods
			set(legend,fieldnames(legend_props)',struct2cell(legend_props)');
			
			legend_text = findall(legend,'Type','text');
			legend_text_lim = get(legend_text(1),'Extent');
			% reduce the size of the lines in the legend
			legend_lines = findall(legend,'Type','Line');
			legend_line_length = get(findall(legend,'type','line'),'XData');
			for l_i = 1:length(legend_line_length)
				if length(legend_line_length{l_i}) == 2
					set(legend_lines(l_i),'XData',[0.05 min(0.3,0.95*legend_text_lim(1))]);
				end
				if length(legend_line_length{l_i}) == 1
					set(legend_lines(l_i),'XData',mean([0.05 min(0.3,0.95*legend_text_lim(1))]));
				end
			end
			
			% reposition the legend
			set(legend,'Location',legend_location,'Box',legend_box_state)
			
			% bring the legend to the front; irritating otherwise.
			axes(legend);
			disp('...legend modified');
		else
			axes(findobj(gcf,'Tag','legend'));
		end
	catch
		disp('...no legend objects found to modify');
	end
end

%% and before we finish, update the figure.
drawnow


%% code for testing; uncomment before use!
% figure
% plot(rand(10),rand(10),'ko','MarkerFaceColor','r')
% legend('test')
% set(legend,'location','best')
% saveas(gcf,'this_is_my_raw_figure','fig')
% readyforprint([4 3],10,[1 1 1],[0 0.1 0.8],3)
% print('-djpeg','this_is_my_printed_jpg')  % for a jpg
% print('-depsc2','-loose','this_is_my_printed_eps') % for an eps