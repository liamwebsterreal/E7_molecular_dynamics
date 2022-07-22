function animation(n)
% generate .gif or .avi
nFrames = n;
mov(1:nFrames) = struct('cdata',[],'colormap',[]);
for i=1:nFrames
    Img = imread(sprintf('images/%d.bmp',i*100));
    imshow(Img,[]);
    frame = getframe(gcf);
    im = frame2im(frame);
    [I,map]=rgb2ind(im,256);
    mov(i)= getframe(gcf);
    if(i==1)
        imwrite(I,map,'movefig.gif','DelayTime',2.0,'LoopCount',Inf)
    else
        imwrite(I,map,'movefig.gif','WriteMode','append','DelayTime',2.0)    
    end
end
v=VideoWriter('collide','MPEG-4');
v.FrameRate = 10;
open(v);
writeVideo(v,mov);
close(v)

v=VideoWriter('collide.avi','Uncompressed AVI');
v.FrameRate = 10;
open(v);
writeVideo(v,mov);
close(v)
end