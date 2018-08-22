function record(file,fig)
% Call:
% record(file,fig)
%
% Description:
% Generate a gif file or include an image to the end of the specified gif file
%
% Inputs:
%     file
%     fig
% Outputs:
%     Generates a new gif file or includes an image to the end of the gif file
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

frame = getframe(fig);
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
if exist(file,'file') == 2
      imwrite(imind,cm,file,'gif','WriteMode','append'); 
else
      imwrite(imind,cm,file,'gif','Loopcount',inf);     
end
end

