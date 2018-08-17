%% fdirichlet_fft2_test

n=[111,155]; % gridpoints in each direction excluding boundaries
h=rand(1,2); % meshstep

F=rand(n);

disp('test p=1')
U=fdirichlet_fft2(F,h,@(x)x);
FF=mlap(F,h);
R=U-FF;
err_PS=norm(R,inf)

disp('test p=-1')
U=fdirichlet_fft2(F,h,@(x)1./x);
FF=mlap(U,h);
R=F-FF;
err_PS=norm(R,inf)