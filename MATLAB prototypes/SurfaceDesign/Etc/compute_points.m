function [xpts_out,ypts_out] = compute_points(Categ,active, polys,color,radius)
if nargin < 5
    radius = 0.26;
end
global scplot searchrad
global G11 G12 angerrorsG12 angerrorsG
angtol = 0.1;
angtol = 100;
scplot = false;

if Categ == 1
    %only 1 curve is affected
    if active == 1
          cp = polys(:,:,1)';
          cx = cp(1,:);
          cy = cp(2,:);
          p0 = [cx(1) cy(1)];
          p1 = [cx(2) cy(2)];
          p2 = [cx(3) cy(3)];
          p3 = [cx(4) cy(4)];
          point = 1;%1 or 2          

          [p0,p1,p2,p3,Database,th,t0,t1,Sc,m,angError] =TransformCrv(p0,p1,p2,p3,point);%find class and database

         %ang_indices = find(angError < angtol);
         x2m = Database(:,3); %Database(ang_indices,3);
         y2m = Database(:,4);%Database(ang_indices,4);
         x1m = Database(:,1);%Database(ang_indices,1);
         y1m = Database(:,2);%Database(ang_indices,2);
         
         
         x2 = p1(1);
         y2 = p1(2);

         indices=find((x2m-x2).^2+(y2m-y2).^2<(radius).^2);

         xpts = (x1m(indices));
         ypts = y1m(indices);
         
         [xpts,ypts] = invTransformCrv_only1pt(xpts,ypts,m,Sc,th,t0,t1);

         if scplot
         p1 =scatter(xpts,ypts,20,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor',color);
         uistack(p1,'bottom');set(p1,'HitTest','off');
         end
         xpts_out = xpts;
         ypts_out = ypts;
         

    else %last point
          cp = fliplr(polys(:,:,end)');
          cx = cp(1,:);
          cy = cp(2,:);
          p0 = [cx(1) cy(1)];
          p1 = [cx(2) cy(2)];
          p2 = [cx(3) cy(3)];
          p3 = [cx(4) cy(4)];

          point = 1;%1 or 2
          [p0,p1,p2,p3,Database,th,t0,t1,Sc,m,angError] =TransformCrv(p0,p1,p2,p3,point);%find class and database

         x2m = Database(:,3); %Database(ang_indices,3);
         y2m = Database(:,4);%Database(ang_indices,4);
         x1m = Database(:,1);%Database(ang_indices,1);
         y1m = Database(:,2);%Database(ang_indices,2);
        % ang_indices = find(angError < angtol);
%          x2m = Database(ang_indices,3);
%          y2m = Database(ang_indices,4);
%          x1m = Database(ang_indices,1);
%          y1m = Database(ang_indices,2);

         x2 = p1(1);
         y2 = p1(2);
        
        
         indices=find((x2m-x2).^2+(y2m-y2).^2<(radius).^2);
         
         xpts = (x1m(indices));
         ypts = y1m(indices);
         [xpts,ypts] = invTransformCrv_only1pt(xpts,ypts,m,Sc,th,t0,t1);

        if scplot
         p1 = scatter(xpts,ypts,20,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor',color); 
         uistack(p1,'bottom');set(p1,'HitTest','off');
        end
         xpts_out = xpts;
         ypts_out = ypts;

    end
end 
if Categ == 2
    %only 1 curve is affected
    if active == 2 %hertil
          cp = polys(:,:,1)';
          cx = cp(1,:);
          cy = cp(2,:);
          p0 = [cx(1) cy(1)];
          p1 = [cx(2) cy(2)];
          p2 = [cx(3) cy(3)];
          p3 = [cx(4) cy(4)];

          point = 2;
          [p0,p1,p2,p3,Database,th,t0,t1,Sc,m,angError] =TransformCrv(p0,p1,p2,p3,point);
        
         x1 = p0(1);
         y1 = p0(2);

%          ang_indices = find(angError < angtol);
%          x2m = Database(ang_indices,3);
%          y2m = Database(ang_indices,4);
%          x1m = Database(ang_indices,1);
%          y1m = Database(ang_indices,2);
         x2m = Database(:,3); %Database(ang_indices,3);
         y2m = Database(:,4);%Database(ang_indices,4);
         x1m = Database(:,1);%Database(ang_indices,1);
         y1m = Database(:,2);%Database(ang_indices,2);
          
 
         indices=find((x1m-x1).^2+(y1m-y1).^2<(radius).^2);
         xpts = (x2m(indices));
         ypts = y2m(indices);
         [xpts,ypts] = invTransformCrv_only1pt(xpts,ypts,m,Sc,th,t0,t1);

         if scplot
         p1 =scatter(xpts,ypts,20,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor',color);
         uistack(p1,'bottom');set(p1,'HitTest','off');
         end
         xpts_out = xpts;
         ypts_out = ypts;
    else 
          cp = fliplr(polys(:,:,end)');
          cx = cp(1,:);
          cy = cp(2,:);
          p0 = [cx(1) cy(1)];
          p1 = [cx(2) cy(2)];
          p2 = [cx(3) cy(3)];
          p3 = [cx(4) cy(4)];
          point = 2;%1 or 2
          [p0,p1,p2,p3,Database,th,t0,t1,Sc,m,angError] =TransformCrv(p0,p1,p2,p3,point);
%          ang_indices = find(angError < angtol);
%          x2m = Database(ang_indices,3);
%          y2m = Database(ang_indices,4);
%          x1m = Database(ang_indices,1);
%          y1m = Database(ang_indices,2);
         x2m = Database(:,3); %Database(ang_indices,3);
         y2m = Database(:,4);%Database(ang_indices,4);
         x1m = Database(:,1);%Database(ang_indices,1);
         y1m = Database(:,2);%Database(ang_indices,2);
         
         x1 = p0(1);
         y1 = p0(2);
         
%          radius=0.07+0.0005*((x1-0.5)^2+y1^2);
%          radius=radius*(0.5+min(0.5*sqrt(x1^2+y1^2),2));
%          radius = max(searchrad,radius);

         indices=find((x1m-x1).^2+(y1m-y1).^2<(radius).^2);
         xpts = (x2m(indices));
         ypts = y2m(indices);
         [xpts,ypts] = invTransformCrv_only1pt(xpts,ypts,m,Sc,th,t0,t1);

         xpts_out = xpts;
         ypts_out = ypts;
         if scplot
         p1=scatter(xpts,ypts,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor',color);  
         uistack(p1,'bottom');set(p1,'HitTest','off');
         end
    end
end
if Categ == 3
         xpts_out = cell(2,1);
         ypts_out = cell(2,1);
         
         %1st curve segment, we need to find the mirror solutions in this
         %case
         
         cp = fliplr(polys(:,:,floor(active/2))');
         cx = cp(1,:);
         cy = cp(2,:);
         p0 = [cx(1) cy(1)];
         p1 = [cx(2) cy(2)];
         p2 = [cx(3) cy(3)];
         p3 = [cx(4) cy(4)];
         point = 2;%1 or 2
         [p0,p1,p2,p3,Database,th,t0,t1,Sc,m,angError] =TransformCrv(p0,p1,p2,p3,point);

         x4 = p0(1); %poly1(1,4);
         y4 = p0(2); %poly1(2,4);
         
         Boorx4 = x4+(x4-p1(1));
         Boory4 = y4+(y4-p1(2));
%          ang_indices = find(angError < angtol);
%          x2m = Database(ang_indices,3);
%          y2m = Database(ang_indices,4);
%          x1m = Database(ang_indices,1);
%          y1m = Database(ang_indices,2);
         x2m = Database(:,3); %Database(ang_indices,3);
         y2m = Database(:,4);%Database(ang_indices,4);
         x1m = Database(:,1);%Database(ang_indices,1);
         y1m = Database(:,2);%Database(ang_indices,2);
         
         %Global control points in the database
         Boorx1m = x1m+(x1m-x2m);
         Boory1m = y1m+(y1m-y2m);
  
%          radius=0.07+0.0005*((Boorx4-0.5)^2+Boory4^2);
%          radius=radius*(0.5+min(0.5*sqrt(Boorx4^2+Boory4^2),2));
%          radius = max(searchrad,radius);
%          radius = 0.56;
         indices=find((Boorx1m-Boorx4).^2+(Boory1m-Boory4).^2<(radius).^2);

         xpts = (x2m(indices));
         ypts = y2m(indices);
        
        [xpts,ypts] = invTransformCrv_only1pt(xpts,ypts,m,Sc,th,t0,t1);
         
         %Mirror solutions

         
         x2m = -1*x2m; %Database(ang_indices,3);
         %y2m = Database(ang_indices,4);
         x1m = -1*x1m; %Database(ang_indices,1);
         %y1m = Database(ang_indices,2);
         
         %Global control points in the database
         Boorx1m = x1m+(x1m-x2m);
         Boory1m = y1m+(y1m-y2m);
         
         indices=find((Boorx1m-Boorx4).^2+(Boory1m-Boory4).^2<(radius).^2);

         xpts2 = (x2m(indices));
         ypts2 = y2m(indices);
        
         [xpts2,ypts2] = invTransformCrv_only1pt(xpts2,ypts2,m,Sc,th,t0,t1);

         xpts = [xpts; xpts2];
         ypts = [ypts; ypts2];
         
         if scplot
         p1 = scatter(xpts,ypts,20,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor',color);  
         uistack(p1,'bottom');set(p1,'HitTest','off');   
         end

         xpts_out{1,1} =[xpts];
         ypts_out{1,1} = [ypts];

%%Plot with respect to second curve segment

         cp = polys(:,:,floor(active/2)+1)';
         cx = cp(1,:);
         cy = cp(2,:);
         p0 = [cx(1) cy(1)];
         p1 = [cx(2) cy(2)];
         p2 = [cx(3) cy(3)];
         p3 = [cx(4) cy(4)];
         point = 1;%1 or 2
         [p0,p1,p2,p3,Database,th,t0,t1,Sc,m,angError] =TransformCrv(p0,p1,p2,p3,point);
%           ang_indices = find(angError < angtol);
%          x2m = Database(ang_indices,3);
%          y2m = Database(ang_indices,4);
%          x1m = Database(ang_indices,1);
%          y1m = Database(ang_indices,2);
         x2m = Database(:,3); %Database(ang_indices,3);
         y2m = Database(:,4);%Database(ang_indices,4);
         x1m = Database(:,1);%Database(ang_indices,1);
         y1m = Database(:,2);%Database(ang_indices,2);
         x2 = p1(1); y2 = p1(2);
         


%          radius=0.07+0.0005*((x2-0.5)^2+y2^2);
%          radius=radius*(0.5+min(0.5*sqrt(x2^2+y2^2),2));
%          radius = max(searchrad,radius);
%         radius = 0.56;
         indices=find((x2m-x2).^2+(y2m-y2).^2<(radius).^2);

         xpts = x1m(indices)+x1m(indices)-x2m(indices);
         ypts = y1m(indices)+y1m(indices)-y2m(indices);
         [xpts,ypts] = invTransformCrv_only1pt(xpts,ypts,m,Sc,th,t0,t1);


         if scplot
         p1 = scatter(xpts,ypts,20,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor',circshift(color,[0 1]));  
         uistack(p1,'bottom');set(p1,'HitTest','off');
         end
         
         xpts_out{2,1} = [xpts];
         ypts_out{2,1} = [ypts];
end

if Categ == 4 
    
    %% Segm1
    
    xpts_out = cell(2,1);
    ypts_out = cell(2,1);
    cp = fliplr(polys(:,:,floor(active/2)-1)');
    cx = cp(1,:);
    cy = cp(2,:);
    p0 = [cx(1) cy(1)];
    p1 = [cx(2) cy(2)];
    p2 = [cx(3) cy(3)];
    p3 = [cx(4) cy(4)];
    
     point = 1;%1 or 2
     [p0,p1,p2,p3,Database,th,t0,t1,Sc,m,angError] =TransformCrv(p0,p1,p2,p3,point);%find class and database

     x2 = p1(1); y2 = p1(2);

     x2m = Database(:,3); %Database(ang_indices,3);
     y2m = Database(:,4);%Database(ang_indices,4);
     x1m = Database(:,1);%Database(ang_indices,1);
     y1m = Database(:,2);%Database(ang_indices,2);
%      ang_indices = find(angError < angtol);
%      x2m = Database(ang_indices,3);
%      y2m = Database(ang_indices,4);
%      x1m = Database(ang_indices,1);
%      y1m = Database(ang_indices,2); 
%      radius=0.07+0.0005*((x2-0.5)^2+y2^2);
%     radius=radius*(0.5+min(0.5*sqrt(x2^2+y2^2),2));
%     radius = max(0.1,radius);
%     radius = 0.56;
%radius = 0.26;
    indices=find((x2m-x2).^2+(y2m-y2).^2<(radius).^2);


    xpts = x1m(indices)+x1m(indices)-x2m(indices);
    ypts = y1m(indices)+y1m(indices)-y2m(indices);
    [xpts,ypts] = invTransformCrv_only1pt(xpts,ypts,m,Sc,th,t0,t1);

     if scplot
     p1 = scatter(xpts,ypts,20,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor',color);
     uistack(p1,'bottom');set(p1,'HitTest','off');
     end

    xpts_out{1,1} =[xpts];
    ypts_out{1,1} = [ypts]; 
    
    %% Plot with respect to second curve segment
    
    
     cp = polys(:,:,floor(active/2))';
     
     cx = cp(1,:);
     cy = cp(2,:);
     p0 = [cx(1) cy(1)];
     p1 = [cx(2) cy(2)];
     p2 = [cx(3) cy(3)];
     p3 = [cx(4) cy(4)];
     point = 2;

     [p0,p1,p2,p3,Database,th,t0,t1,Sc,m,angError] =TransformCrv(p0,p1,p2,p3,point);%find class and databasex
     
         x2m = Database(:,3); %Database(ang_indices,3);
         y2m = Database(:,4);%Database(ang_indices,4);
         x1m = Database(:,1);%Database(ang_indices,1);
         y1m = Database(:,2);%Database(ang_indices,2);
%       ang_indices = find(angError < angtol);
%          x2m = Database(ang_indices,3);
%          y2m = Database(ang_indices,4);
%          x1m = Database(ang_indices,1);
%          y1m = Database(ang_indices,2);
     x1 = p0(1); y1 = p0(2);

     Boorx1 =x1+(x1-p1(1));
     Boory1 = y1+(y1-p1(2));

     
     Boorx1m = x1m+(x1m-x2m);
     Boory1m = y1m+(y1m-y2m);
%      radius=0.07+0.0005*((Boorx1-0.5)^2+Boory1^2);
%      radius=radius*(0.5+min(0.5*sqrt(Boorx1^2+Boory1^2),2));
%      radius = max(0.1,radius);
%      radius = 0.56;
%radius = 0.26;

     indices=find((Boorx1m-Boorx1).^2+(Boory1m-Boory1).^2<(radius).^2);
     
     xpts = x2m(indices);
     ypts = y2m(indices);
     [xpts,ypts] = invTransformCrv_only1pt(xpts,ypts,m,Sc,th,t0,t1);

     %Mirror solutions

              x2m = -1*x2m;
          x1m = -1*x1m;
          
%          y1m = Database(ang_indices,2);     
%          x2m = -1*Database(ang_indices,3);
%          y2m = Database(ang_indices,4);
%          x1m = -1*Database(ang_indices,1);
%          y1m = Database(ang_indices,2);
         
         %Global control points in the database
         Boorx1m = x1m+(x1m-x2m);
         Boory1m = y1m+(y1m-y2m);
         
         indices=find((Boorx1m-Boorx1).^2+(Boory1m-Boory1).^2<(radius).^2);

       
         xpts2 = (x2m(indices));
         ypts2 = y2m(indices);
        
         [xpts2,ypts2] = invTransformCrv_only1pt(xpts2,ypts2,m,Sc,th,t0,t1);

         xpts = [xpts; xpts2];
         ypts = [ypts; ypts2];

    if scplot
    p1 = scatter(xpts,ypts,20,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor',circshift(color,[0 1])); 
    uistack(p1,'bottom');set(p1,'HitTest','off');
    end
    
    xpts_out{2,1} =[xpts];
    ypts_out{2,1} = [ypts];
    
  %  xpts_out{1,1} =[xpts];
  %  ypts_out{1,1} = [ypts]; 
  
end

end