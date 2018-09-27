function plot_data(n_blocks)

global active polys scplot
global G11 G12 angerrorsG12 angerrors

if active == 1 || active == 4+(n_blocks-2)+(n_blocks-1)+1
    Categ = 1;
elseif active-1 == 1 || active == 4+(n_blocks-2)+(n_blocks-1)
    Categ = 2;
else
    if mod(active,2) == 0
        Categ = 4;
    else
        Categ = 3;
    end
end

%%  From Categ plot curves in the database

if Categ == 1
if active == 1
    crv = spmak([0 0 0 0 1 1 1 1],polys(:,:,1)');
    fnplt(crv,'r');
else
    crv = spmak([0 0 0 0 1 1 1 1],polys(:,:,end)');
    fnplt(crv,'r');
end
end
if Categ == 2
if active == 2
    crv = spmak([0 0 0 0 1 1 1 1],polys(:,:,1)');
    fnplt(crv,'r');
else
    crv = spmak([0 0 0 0 1 1 1 1],polys(:,:,end)');
    fnplt(crv,'r');
end
end
if Categ == 3
    crv = spmak([0 0 0 0 1 1 1 1],polys(:,:,floor(active/2))');
    fnplt(crv,'r')
    crv = spmak([0 0 0 0 1 1 1 1],polys(:,:,floor(active/2)+1)');
    fnplt(crv,'g')
end
if Categ == 4
    crv = spmak([0 0 0 0 1 1 1 1],polys(:,:,floor(active/2))');
    fnplt(crv,'r')
    crv = spmak([0 0 0 0 1 1 1 1],polys(:,:,floor(active/2)-1)');
    fnplt(crv,'g')
end

%Figure out which database to use

[xG,yG] = compute_points(Categ,active, polys,[0.4 0. 0]);
warning('off','MATLAB:alphaShape:DupPointsBasicWarnId');
if Categ ==3 || Categ == 4

    totalx = [xG{1,1}(:,1)'];
    totaly = [yG{1,1}(:,1)'];
    


   k1 = alphaShape(totalx',totaly'); 
    k1_x = k1.Points(:,1);
    k1_y = k1.Points(:,2);

    totalx1 = [xG{2,1}(:,1)'];
    totaly1 = [yG{2,1}(:,1)'];


    k3 = alphaShape(totalx1',totaly1');
    k3_x = k3.Points(:,1);
    k3_y = k3.Points(:,2);

    %compute intersection
    t = find(inShape(k3,k1_x,k1_y));
    t2 = find(inShape(k1,k3_x,k3_y));
   % [xb,yb] = polybool('intersection',k1_x,k1_y,k3_x,k3_y); %totalx(k1),totaly(k1),totalx1(k3),totaly1(k3));
   xb = k1_x(t);
   yb = k1_y(t);
   xb = [xb' k3_x(t2)']';
   yb = [yb' k3_y(t2)']';

   shp = alphaShape(xb,yb,0.4); 
   hold on
   h = plot(shp,'FaceColor', [0.58 0.2 0.23],'EdgeColor','none');
   set(h,'facealpha',.5);
h.HitTest = 'off';

end
if Categ == 1 || Categ ==2 
    totalx = [xG(:,1)'];
    totaly = [yG(:,1)'];

   
   shp = alphaShape(totalx',totaly',0.4); 
   hold on
   h = plot(shp,'FaceColor',  [0.58 0.2 0.23],'EdgeColor','none');
   set(h,'facealpha',.5);

h.HitTest = 'off';

end

end