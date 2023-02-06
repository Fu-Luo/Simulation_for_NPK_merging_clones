clc
clear
rng(726);
testcount=0;
outputmatrix=zeros(100,4);
while testcount<=1
    testcount=testcount+1;

%simulation for npk merging events 

endpointgridsize=600; %final amount space that can be occupied 
initialseeding=28872; %initial observed cell amount e.g 3 weeks  mean(totalcell.np3weeks$totalcell)= 29877.8 median(totalcell.np3weeks$totalcell)=28872
% negative clones will grow faster
grid=zeros(endpointgridsize,endpointgridsize);
%seeding clones and cells
seedingclones=650;%median total # clones at 3 weeks=643.5 mean =  713.05

sizechart=[2 0.17;3 0.075;4 0.035;5 0.021;8 0.046; 13 0.016 ;18 0.008; 23 0.0048; 28 0.0030 ]; %percentage of clone sizes based on npk week3, size 8 (6-10),13(11-15) .etc 
sizechart=flip(sizechart,1);
sizechart(:,3)=ceil(sizechart(:,2)*seedingclones);

%initiaze one seeding coordinate for each clone
[X,Y]=meshgrid(2:endpointgridsize-1,2:endpointgridsize-1);
for x = 1:(endpointgridsize-2)*(endpointgridsize-2)
matrix(x,1)=X(x);
matrix(x,2)=Y(x);
end
shuffledmatrix = matrix(randperm(size(matrix,1)),:);
i = shuffledmatrix(:,1);
j = shuffledmatrix(:,2);
for a =1: seedingclones
    seedingclone_i(a)=i(a);
    seedingclone_j(a)=j(a);
end


%seed clone based on clone sizes
clones=cell(1,seedingclones);
for a=1:size(clones,2)
     clones{a}=[seedingclone_i(a),seedingclone_j(a)];
end
sizechart(:,2)=[];
sizechart=[sizechart;1,seedingclones-sum(sizechart(:,2))];
sizes = repelem(sizechart(:,1), sizechart(:,2));
for a=1:length(sizes)
    maxtry=0;
     if sizes(a)>1 & maxtry<100
         sizecounter=1;
        while sizecounter<=sizes(a)-1   & maxtry<100
            trycenter=randi(size(clones{a},1));
            maxtry=maxtry+1;
            i=clones{a}(trycenter,1);
            j=clones{a}(trycenter,2);
            if i>1&i<endpointgridsize&j>1&j<endpointgridsize
            neighborslocation=[i-1,j-1;i-1,j;i-1,j+1;i,j-1;i,j+1;i+1,j-1;i+1,j;i+1,j+1];%record the 8 locations
            %checkborder 
             for c=randperm(8)
                if grid(neighborslocation(c,1),neighborslocation(c,2))==0
                    sizecounter=sizecounter+1;
                    grid(neighborslocation(c,1),neighborslocation(c,2))=1;
                    clones{a}=[clones{a};neighborslocation(c,1),neighborslocation(c,2)];
                  if sizecounter==sizes(a)
                      break
                  end
                end
             end
            end
        end
     end
 end

%seeding rest non labeled cells 
cellcount=(initialseeding-sum(sum(grid)))/(endpointgridsize*endpointgridsize-sum(sum(grid)));
for i=1:endpointgridsize
    for j=1:endpointgridsize 
         if grid(i,j)==0
            if rand()<cellcount
                 grid(i,j)=1;
            end
        end
    end
end
sum(sum(grid))
            

%assign color for seeding clones
colorspec=zeros(endpointgridsize,endpointgridsize);%create a matrix of the same size for color of the cells
red=ceil(seedingclones*0.419); %color distribution based on 3 weeks NP
yellow=ceil(seedingclones*0.246);
blue=ceil(seedingclones*0.310);
green=seedingclones-(red+yellow+blue);
colors = repelem([1 2 3 4], [red blue yellow green]);

for i =1:seedingclones
    for a=1:size(clones{i},1)
        colorspec(clones{i}(a,1),clones{i}(a,2))=colors(i);
        clones{i}(a,3)=colors(i);
    end
end


%20% clones have slow dying cells, out of these clones 0.3% cell die each round
% rest clones (80%) joins non labeled cells and 14% of all die each round
% (from np 3 weeks tunel data and conf+ cells loss)
% dying signal clones colomn 4, slow dying clone=1, fast dying clones =0
for a= 1:length(clones)
    if rand()<=0.20
        clones{a}(:,4)=1;
    else
        clones{a}(:,4)=0;
    end
end

% 3% clones have proliferating cells (np 3weeks, 8/269), out of
%these 3% clones: 
% 50% clones(0.5) contain 75-100% ki67+; 25% clones (1) will contain 50-75% ki67+ cells; %12.5(2)will contain
%25%-50% cells; 12.5% (3)will contain 0-25% ki67+ cells 
% rest clones 1% of all proliferate each round. 
% proliferating signal clones colomn 5, proliferating clone=0.5/1/2/3, regular =0
for a= 1:length(clones)
    if rand()<=0.03
         test=rand()
         if test<=0.5
            clones{a}(:,5)=0.5;
         elseif  test<=0.75 & test>0.5
             clones{a}(:,5)=1;
         elseif test>0.75 & test<=0.875
             clones{a}(:,5)=2;
         else
             clones{a}(:,5)=3;   
         end
    else
        clones{a}(:,5)=0;
    end
end

%cell death and proliferation begins here - start timepoint count 

timepoint=0;
clonesizes=0;
%maxclonesize=1245 %min largest clones among 12 weeks section
maxclonesize=90 %mean 115 median 66 largest clones among 4 weeks section 
while max(clonesizes)<maxclonesize 
      timepoint=timepoint+1
      sum(sum(grid))
for a=1:length(clones)
    if sum(clones{a}(:,4))>0
        b=0;
       while b<size(clones{a},1)
           b=b+1;
           if rand()<=0.003
               grid(clones{a}(b,1),clones{a}(b,2))=0;
               colorspec(clones{a}(b,1),clones{a}(b,2))=0;
               clones{a}(b,:)=[];
           end
       end
    else
        c=0;
        while c<size(clones{a},1)
            c=c+1;
           if rand()<0.14
               grid(clones{a}(c,1),clones{a}(c,2))=0;
               colorspec(clones{a}(c,1),clones{a}(c,2))=0;
               clones{a}(c,:)=[];
           end
        end
    end
end

for i=1:endpointgridsize
    for j=1:endpointgridsize 
         if grid(i,j)==1
            if colorspec(i,j)==0
                if rand()<0.08
                 grid(i,j)=0;
                end
            end
        end
    end
end
% cell proliferation one round 
%20% clones have proliferating cells (npk 3weeks, avg=19.44%), out of these
%clones 27% cell proliferate each round (average ki67% in ki67+clones
% rest clones  1% of all proliferate each round. Non labeling 30%
% proliferating signal clones colomn 5, proliferating clone=1, reguloar =0

for a=1:length(clones)
    if sum(clones{a}(:,5))>0
        if clones{a}(1,5)==0.5
            proliferationrate=0.25*rand() + 0.75;
        elseif clones{a}(1,5)==1
            proliferationrate=0.25*rand() + 0.5;
        elseif clones{a}(1,5)==2
            proliferationrate=0.25*rand() + 0.25;
        elseif clones{a}(1,5)==3
            proliferationrate=0.25*rand();
        end
        proliferatingcount=min(ceil(size(clones{a},1)*proliferationrate),50);%most ki67 cells observed in one clone (upper limit=15 for 4weeks; upper limit=50 for 12weeks)
        upperlimit=0;
        while proliferatingcount>0 && upperlimit<1000 %limit trys 
            b=ceil(rand()*size(clones{a},1));
            i=clones{a}(b,1);
            j=clones{a}(b,2);
            upperlimit=upperlimit+1;
            if i>1&i<endpointgridsize&j>1&j<endpointgridsize
            neighborslocation=[i-1,j-1;i-1,j;i-1,j+1;i,j-1;i,j+1;i+1,j-1;i+1,j;i+1,j+1];%record the 8 locations
             for c=randperm(8)
                 if proliferatingcount>0
                  if grid(neighborslocation(c,1),neighborslocation(c,2))==0
                    proliferatingcount=proliferatingcount-1;
                    grid(neighborslocation(c,1),neighborslocation(c,2))=1;
                    clones{a}=[clones{a};neighborslocation(c,1),neighborslocation(c,2),clones{a}(b,3:5)];
                    colorspec(neighborslocation(c,1),neighborslocation(c,2))=clones{a}(b,3);
                  end
                 end
             end
            end
        end
    else
         proliferatingcount=min(ceil(size(clones{a},1)*0.01),15);
         upperlimit=0;
         while proliferatingcount>0 && upperlimit<1000
            b=ceil(rand()*size(clones{a},1));
            i=clones{a}(b,1);
            j=clones{a}(b,2);
            upperlimit=upperlimit+1;
            if i>1&i<endpointgridsize&j>1&j<endpointgridsize
            neighborslocation=[i-1,j-1;i-1,j;i-1,j+1;i,j-1;i,j+1;i+1,j-1;i+1,j;i+1,j+1];%record the 8 locations
             for c=randperm(8)
                if proliferatingcount>0
                 if grid(neighborslocation(c,1),neighborslocation(c,2))==0
                    proliferatingcount=proliferatingcount-1;                  
                    grid(neighborslocation(c,1),neighborslocation(c,2))=1;
                    clones{a}=[clones{a};neighborslocation(c,1),neighborslocation(c,2),clones{a}(b,3:5)];
                    colorspec(neighborslocation(c,1),neighborslocation(c,2))=clones{a}(b,3);
                 end
                end
             end
            end
         end 
    end
end
for i=2:endpointgridsize-1
    for j=2:endpointgridsize-1 
        if i>1&i<endpointgridsize&j>1&j<endpointgridsize
         if grid(i,j)==1
            if colorspec(i,j)==0
                if rand()<0.2
                  neighborslocation=[i-1,j-1;i-1,j;i-1,j+1;i,j-1;i,j+1;i+1,j-1;i+1,j;i+1,j+1];%record the 8 locations
                      for c=randperm(8)
                      if grid(neighborslocation(c,1),neighborslocation(c,2))==0
                         grid(neighborslocation(c,1),neighborslocation(c,2))=1;
                      end
                      break
                      end                      
                end
            end
         end
        end
    end
end



 for a= 1:length(clones)
    
  if size(clones{a},1)<mean(clonesizes)
        if rand()<0.05
            grid(clones{a}(:,1),clones{a}(:,2))=0;
               colorspec(clones{a}(:,1),clones{a}(:,2))=0;
               clones{a}=[];
        end 
  end
 end

clones = clones(~cellfun(@isempty, clones));%remove diminished clones from original clone list




colormaps=zeros(endpointgridsize,endpointgridsize,3);% plot colors based on color labels 
for i=1:endpointgridsize
    for j=1:endpointgridsize

      if  colorspec(i,j)==1;
            colormaps(i,j,:)=[1 0 0];%1 for red,2 for blue, 3 for yellow, 4 for green         
      elseif colorspec(i,j)==2;
        colormaps(i,j,:)=[0 0 1];
      elseif  colorspec(i,j)==3;
        colormaps(i,j,:)=[1 1 0]; 
      elseif colorspec(i,j)==4;
        colormaps(i,j,:)=[0 1 0];   
      elseif grid(i,j)==1 && colorspec(i,j)==0;
        colormaps(i,j,:)=[0.5 0.5 0.5]; %grey for unlabled cells 
      else
        colormaps(i,j,:)=[0 0 0];    % black for unoccupied space
      end
    end
end
   imshow(colormaps,'InitialMagnification', 400)%plot colors
   %if timepoint==1;
      % imshow(colormaps,'InitialMagnification', 400)
  % end

  h1=colorspec;
  h1(h1>1)=0;
  h2=colorspec;
  h2(h2<2 | h2>=3)=0;
  h3=colorspec;
  h3(h3<3 |h3==4)=0;
  h4=colorspec;
  h4(h4<4)=0; % seperate 4 colors into individual matrices 
        clones1=bwconncomp(h1);% use build in function of matlab to find connected components 
        clones2=bwconncomp(h2);  
        clones3=bwconncomp(h3);
        clones4=bwconncomp(h4);
        clonesizes1=cellfun(@numel,clones1.PixelIdxList); % extract sizes of connected components 
        clonesizes2=cellfun(@numel,clones2.PixelIdxList);
        clonesizes3=cellfun(@numel,clones3.PixelIdxList);      
        clonesizes4=cellfun(@numel,clones4.PixelIdxList);
        clonesizes=[clonesizes1,clonesizes2,clonesizes3,clonesizes4];
%monitor for largest clone appeared after each round 

end 


%detect merge events after ending condition is met: 
minsize=2; %look at large clones only 
mergingevent=zeros(1,3);% create a matrix to record merging event
% in each connected components, find the origin clone and record them 
        stats_1 = regionprops(clones1,'pixellist');
        stats_1=addclonenumber(stats_1,clones,minsize);
        
        stats_2 = regionprops(clones2,'pixellist');
        stats_2=addclonenumber(stats_2,clones,minsize);
        
        stats_3 = regionprops(clones3,'pixellist');
        stats_3=addclonenumber(stats_3,clones,minsize);
        
        stats_4 = regionprops(clones4,'pixellist');
        stats_4=addclonenumber(stats_4,clones,minsize);
        n=1;
 
        for i =1:size(stats_1,1)
            if size(stats_1(i).PixelList,1) >= minsize
            clonedistribution=tabulate(stats_1(i).PixelList(:,3));
            if max(clonedistribution(:,3))<80 %detect clones contributing with less than 80 % of the ending clones
                mergingevent(n,1)=size(stats_1(i).PixelList,1); %record this merging event 
                mergingevent(n,2)=length(unique(stats_1(i).PixelList(:,3)));
                mergingevent(n,3)=1;
                n=n+1;
            end
            end
        end
        for i =1:size(stats_2,1)
            if size(stats_2(i).PixelList,1) >= minsize
            clonedistribution=tabulate(stats_2(i).PixelList(:,3));
            if max(clonedistribution(:,3))<80
               mergingevent(n,1)=size(stats_2(i).PixelList,1);
                mergingevent(n,2)=length(unique(stats_2(i).PixelList(:,3)));
                mergingevent(n,3)=2;
                n=n+1;
            end
            end
        end
        for i =1:size(stats_3,1)
            if size(stats_3(i).PixelList,1) >= minsize
            clonedistribution=tabulate(stats_3(i).PixelList(:,3));
            if max(clonedistribution(:,3))<80
                mergingevent(n,1)=size(stats_3(i).PixelList,1);
                mergingevent(n,2)=length(unique(stats_3(i).PixelList(:,3)));
                mergingevent(n,3)=3;
                n=n+1;
            end
            end
        end
        for i =1:size(stats_4,1)
           if size(stats_4(i).PixelList,1) >= minsize
            clonedistribution=tabulate(stats_4(i).PixelList(:,3));
            if max(clonedistribution(:,3))<80
                mergingevent(n,1)=size(stats_4(i).PixelList,1);
                mergingevent(n,2)=length(unique(stats_4(i).PixelList(:,3)));
                mergingevent(n,3)=4;
                n=n+1;
            end
           end
        end

     
if mergingevent(1,1)>0
            %check all clones
            %check large clones clonesize > 15 (top 5% clone sizes in np
            %4weeks)
       outputmatrix(testcount,1)=size(mergingevent,1)/sum(clonesizes>15);
       outputmatrix(testcount,2)=size(mergingevent,1);
       outputmatrix(testcount,3)=sum(clonesizes>15);
       outputmatrix(testcount,4)=max(clonesizes);
        else
       outputmatrix(testcount,1)=0;
       outputmatrix(testcount,2)=0;
       outputmatrix(testcount,3)=sum(clonesizes>14)
       ;
       outputmatrix(testcount,4)=max(clonesizes); 
        end

end

writematrix(outputmatrix,"NP3to4largeclone15.csv")