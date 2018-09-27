function compute_curve(poly,n_curves)
        x0=poly(1,1); y0=poly(1,2);
        xL=poly(2,1); yL=poly(2,2);
        xR=poly(3,1); yR=poly(3,2);
        xe=poly(4,1); ye=poly(4,2);
        FirstGuessD([xL,yL],[xR,yR],[x0,y0],[xe,ye],1,n_curves);   %1 means plot!
end