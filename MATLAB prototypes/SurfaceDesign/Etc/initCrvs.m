function [curves,kn] = initCrvs(opt)

knots = linspace(0,1,opt.nBlocks+1); x=knots';
r=repmat(x,1,2)';r=r(:)'; 
kn = [0 0 r 1 1 ];

if mod(opt.nBlocks,2) ~= 0
    initlinsp = linspace(-((opt.nBlocks-1)/2)*opt.bW,opt.bW*((opt.nBlocks+1)/2),4+3*(opt.nBlocks-1));
    if opt.nBlocks > 1
        indx2rem = [];
        for i = 2:opt.nBlocks
            indx2rem = [indx2rem 4+3*(i-2)];
        end
        initlinsp(indx2rem) = [];
    end
    zerolist = repmat(0.0,1,numel(initlinsp));
    zerolist(1:2:end) = 0.0001;
    init_curve = [initlinsp; zerolist];
else
    initlinsp = linspace(-((opt.nBlocks)/2)*opt.bW,opt.bW*((opt.nBlocks)/2),4+3*(opt.nBlocks-1));
    if opt.nBlocks > 1
        indx2rem = [];
        for i = 2:opt.nBlocks
            indx2rem = [indx2rem 4+3*(i-2)];
        end
        initlinsp(indx2rem) = [];
    end
    zerolist = repmat(0.0,1,numel(initlinsp));
    zerolist(1:2:end) = 0.0001;
    init_curve = [initlinsp; zerolist];
    
end
for j=1:opt.nICurves
     curves(j,:) = init_curve(:)'; 
end

end