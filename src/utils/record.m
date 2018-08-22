function record(fig,k)
frame = getframe(fig);
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
if k == 1 
      imwrite(imind,cm,'example.gif','gif', 'Loopcount',inf); 
else 
      imwrite(imind,cm,'example.gif','gif','WriteMode','append'); 
end
end

