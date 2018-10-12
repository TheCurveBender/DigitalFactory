function redraw_fig2(n_blocks,block_size,block_size_h,C2)
global csize active coefs Mold boxcon; 
     if nargin <4
         C2 = false;
     end
     actv = coefs(:,active);
     hold on
     plot(actv(1),actv(2),'.','color',[0 0 0],'MarkerSize',9*csize)
     plot(coefs(1,:),coefs(2,:),'--ok')
     axis equal
     set_axes2
     if boxcon
         if C2
             for j = 1:n_blocks-1
             %Two lines
             xval = 0;
             hold on
             plot([xval xval],[-block_size_h/2,block_size_h/2],'--k')
             end
             xval = Mold(1,1);
             hold on
             plot([xval xval],[-block_size_h/2,block_size_h/2],'--k')
             xval = Mold(1,end);
             hold on
              plot([xval xval],[-block_size_h/2,block_size_h/2],'--k')

         end
         if ~C2
         for j = 1:n_blocks-1
             %Two lines
             xval = (Mold(1,4+(j-1)*3));
             hold on
             plot([xval xval],[-block_size_h/2,block_size_h/2],'--k')
         end
        xval = Mold(1,1);
        hold on
        plot([xval xval],[-block_size_h/2,block_size_h/2],'--k')
        xval = Mold(1,end);
        hold on
        plot([xval xval],[-block_size_h/2,block_size_h/2],'--k')
         end
end
end