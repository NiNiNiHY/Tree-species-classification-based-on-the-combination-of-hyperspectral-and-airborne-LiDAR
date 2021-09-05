

class = 9;%样本种类
% lidar = 23;  %雷达特征
 lidar = 42;  %雷达特征
lidar= lidar+426; 


fid=fopen('C:\Users\25297\Desktop\LiDAR\data\train.txt');
fid1=fopen('C:\Users\25297\Desktop\LiDAR\data\test.txt');
c1 = '%d ';
c2 = '%*d';
c3 = ':%f';

formatSpec(:,1)=c1;
for i=1:lidar
 formatSpec(:,i*2) =c2;
 formatSpec(:,i*2+1) =c3;
end


A=fscanf(fid,formatSpec,[lidar+1,inf]); %指定’。’为分隔符，如果不指定分隔符的话，就需要把formatSpec写成'%d。%d。%d。%d' 。
fclose(fid);
A=A.';
lidar_row = size(A,1);
lidar_col = size(A,2);

for i = 1:lidar_row
for j = 2:lidar_col
    if A(i,j) == inf  
        A(i,j) = 0 ;
    end
    if A(i,j) == -inf  
        A(i,j) = 0 ;
    end
end
end
A2=fscanf(fid1,formatSpec,[lidar+1,inf]); %指定’。’为分隔符，如果不指定分隔符的话，就需要把formatSpec写成'%d。%d。%d。%d' 。
fclose(fid1);
A2=A2.';
lidar_row1 = size(A2,1);
lidar_col1 = size(A2,2);

for i = 1:lidar_row1
for j = 2:lidar_col1
    if A2(i,j) == inf  
        A2(i,j) = 0 ;
    end
    if A2(i,j) == -inf  
        A2(i,j) = 0 ;
    end
end
end
% %将两组数据进行合并，并进行归一化处理。
% A2 = [A;A1] ;
% for i=2:lidar+1
%    max_lidar = max(A2(:,i)) ;
%    min_lidar = min(A2(:,i)) ;
%    A2(:,i)=(A2(:,i)-min_lidar)/( max_lidar -min_lidar);
%    if isnan(A2(:,i))
%        A2(:,i) = 0;
%    end
% end
% 
% A = A2(1:lidar_row,:) ;
% A1 = A2(lidar_row+1:lidar_row+lidar_row1,:) ;

A1 = A;
for i=2:lidar+1
   max_lidar = max(A(:,i)) ;
   min_lidar = min(A(:,i)) ;
   A1(:,i)=(A(:,i)-min_lidar)/( max_lidar -min_lidar);
   if isnan(A1(:,i))
       A1(:,i) = 0;
   end
end

for i = 1:lidar_row
for j = 2:lidar_col
    if isnan(A1(i,j))
        A1(i,j) = 0 ;
    end
end
end

%将两组数据进行合并，并进行归一化处理。

for i=2:lidar+1
   max_lidar1 = max(A(:,i)) ;
   min_lidar1 = min(A(:,i)) ;
   A2(:,i)=(A2(:,i)-min_lidar1)/( max_lidar1 -min_lidar1);
   if isnan(A2(:,i))
       A2(:,i) = 0;
   end
end

for i = 1:lidar_row1
for j = 2:lidar_col1
    if isnan(A2(i,j))
        A2(i,j) = 0 ;
    end
end
end

fid = fopen('train.txt','wt');   % data.txt为写入文件名
matrix = A1;                     % M为要存储的矩阵
[m,n]=size(matrix);                      
 for i=1:1:m
   for j=1:1:n
      if j==n
%         fprintf(fid,'%d:%f\n',j-1,matrix(i,j));
        fprintf(fid,'\n');
      else
        if j==1
             fprintf(fid,'%d\t',matrix(i,j));
%         elseif  (j==6)||(j==15)||(j==196)||(j==206)||(j==427)||(j==465)||(j==466)||(j==429)||(j==430)||(j==436)||(j==454)||(j==437)
        elseif  (j==6)||(j==15)||(j==196)||(j==206)||(j==427)||(j==465)||(j==466)||(j==429)||(j==430)||(j==436)||(j==454)||(j==437)
        fprintf(fid,'%d:%f\t',j-1,matrix(i,j));
        end
      end
   end
end
fclose(fid);

fid = fopen('test.txt','wt');   % data.txt为写入文件名
matrix = A2;                     % M为要存储的矩阵
[m,n]=size(matrix);                      
 for i=1:1:m
   for j=1:1:n
      if j==n
%       fprintf(fid,'%d:%f\n',j-1,matrix(i,j));
        fprintf(fid,'\n');
      else
        if j==1
             fprintf(fid,'%d\t',matrix(i,j));
        elseif (j==6)||(j==15)||(j==196)||(j==206)||(j==427)||(j==465)||(j==466)||(j==429)||(j==430)||(j==436)||(j==454)||(j==437)
        fprintf(fid,'%d:%f\t',j-1,matrix(i,j));
        end
      end
   end
end
fclose(fid);



fid = fopen('trainhsi.txt','wt');   % data.txt为写入文件名
matrix = A1;                     % M为要存储的矩阵
[m,n]=size(matrix); 
%n =469;
n =427;
 for i=1:1:m
   for j=1:1:n
      if j==n
        fprintf(fid,'%d:%f\n',j-1,matrix(i,j));
      else
        if j==1
             fprintf(fid,'%d\t',matrix(i,j));
        else
        fprintf(fid,'%d:%f\t',j-1,matrix(i,j));
        end
      end
   end
end
fclose(fid);

fid = fopen('testhsi.txt','wt');   % data.txt为写入文件名
matrix = A2;                     % M为要存储的矩阵
[m,n]=size(matrix); 
%n =469;
n =427;
 for i=1:1:m
   for j=1:1:n
      if j==n
        fprintf(fid,'%d:%f\n',j-1,matrix(i,j));
      else
        if j==1
             fprintf(fid,'%d\t',matrix(i,j));
        else
        fprintf(fid,'%d:%f\t',j-1,matrix(i,j));
        end
      end
   end
end
fclose(fid);





%将A中数据转到xlsx表格中，其中第一行需进行手动标特征名称
% X = A2();
% X(:,1) = [];
% Y =  A2(:,1);
% xlswrite('X.xlsx', X);
% xlswrite('Y.xlsx', Y);%添加第一行，值为0

%/********************************将A中数据进行按比例划分，用于训练和和测试*************************************/
% B= ones(class,lidarhsi)*nan; %记录每类树的标号
% class_number = ones(1,14)*0;  %记录每类树的标签数
% k=1;
% 
%  for j =1:class  
%     for i = 1:lidarhsi
%      if A(i,1)==j
%          B(j,k)=i;
%          k=k+1;
%          class_number(1,j)=class_number(1,j)+1;
%      end
%     end
%     k=1;
%  end
% 
%   lidar_train = ones(1,lidar+1)*nan;
%   lidar_test = ones(1,lidar+1)*nan;
%   sum_train=0;
%   sum_test =0;
% for j=1:class
%   train = round(class_number(1,j)*0.8);
%   test =  round(class_number(1,j)*0.2);
% random = randperm(class_number(1,j));
% lidar_train1 = ones(train,lidar+1)*nan;
% lidar_test1 = ones(test,lidar+1)*nan;
% for i = 1:class_number(1,j)
%     if i <=train
%      lidar_train1(i,:) = A(B(j,random(i)),:);   
%     end
%      if i >train
%      lidar_test1(i-train,:) = A(B(j,random(i)),:);   
%     end
% end
% lidar_train = [lidar_train;lidar_train1(:,:)];
% lidar_test = [lidar_test;lidar_test1(:,:)];
% end
% lidar_train(1,:)=[];
% lidar_test (1,:)=[];



lidar_train=A1;
lidar_test=A2;

class_number = ones(1,class)*0;  %记录每类树的标签数
user_number = ones(1,class)*0;  %记录每类树的分类总数
lable_number = ones(1,class)*0;  %记录每类树的正确预测数
user_numberhsi = ones(1,14)*0;  %记录每类树的分类总数
lable_numberhsi = ones(1,14)*0;  %记录每类树的正确预测数
PA = ones(1,14)*0;  %
UA = ones(1,14)*0;  %
PAhsi = ones(1,14)*0;  %
UAhsi = ones(1,14)*0;  %

 for j =1:class  
    for i = 1:size(A2,1)
     if A2(i,1)==j
         class_number(1,j)=class_number(1,j)+1;
     end
    end
 end
%lidar+hsi精度
n =0;
for i=1:size(lidar_test,1)
    mindistance = inf;
for j=1:size(lidar_train,1)
     sumdistance=0;
%      for k = 2:469
%       sumdistance = sumdistance+(lidar_test(i,k)-lidar_train(j,k))^2;
%     end

%光谱和特征都按照排序改
  
    for k = 6
      sumdistance = sumdistance+(lidar_test(i,k)-lidar_train(j,k))^2;
    end
    for k = 15
      sumdistance = sumdistance+(lidar_test(i,k)-lidar_train(j,k))^2;
    end
    for k = 196
      sumdistance = sumdistance+(lidar_test(i,k)-lidar_train(j,k))^2;
    end
    for k = 206
      sumdistance = sumdistance+(lidar_test(i,k)-lidar_train(j,k))^2;
    end
    for k = 427
      sumdistance = sumdistance+(lidar_test(i,k)-lidar_train(j,k))^2;
    end
    
    for k = 465:466
      sumdistance = sumdistance+(lidar_test(i,k)-lidar_train(j,k))^2;
    end
    for k = 429:430
      sumdistance = sumdistance+(lidar_test(i,k)-lidar_train(j,k))^2;
    end
    for k = 436
      sumdistance = sumdistance+(lidar_test(i,k)-lidar_train(j,k))^2;
    end
    for k = 454
      sumdistance = sumdistance+(lidar_test(i,k)-lidar_train(j,k))^2;
    end
    for k = 437
      sumdistance = sumdistance+(lidar_test(i,k)-lidar_train(j,k))^2;
    end

  if  sumdistance<mindistance
      mindistance = sumdistance;
      lable = j;
  end 
end
 user_number(1,lidar_train(lable,1)) = user_number(1,lidar_train(lable,1)) +1;
 if lidar_test(i,1)==lidar_train(lable,1)
     n = n+1;
     lable_number(1,lidar_test(i,1))= lable_number(1,lidar_test(i,1))+1;
 end
end

OA = n / size(lidar_test,1);
for i=1:class
PA(1,i) = lable_number(1,i)/(class_number(1,i));
UA(1,i) = lable_number(1,i)/(user_number(1,i));
end
xii =0;
xxii =0;
for i=1:class
    xii = xii+lable_number(1,i);
    xxii =xxii+class_number(1,i)*user_number(1,i);
end
Kappa = (size(lidar_test,1)*xii-xxii)/(size(lidar_test,1)*size(lidar_test,1)-xxii);

%高光谱精度
n =0;
for i=1:size(lidar_test,1)
    mindistance = inf;
for j=1:size(lidar_train,1)
     sumdistance=0;
   % for k = 2:427
    for k = 2:469
      sumdistance = sumdistance+(lidar_test(i,k)-lidar_train(j,k))^2;
    end
  if  sumdistance<mindistance
      mindistance = sumdistance;
      lable = j;
  end 
end
user_numberhsi(1,lidar_train(lable,1)) = user_numberhsi(1,lidar_train(lable,1)) +1;
 if lidar_test(i,1)==lidar_train(lable,1)
     n = n+1;
     lable_numberhsi(1,lidar_test(i,1))= lable_numberhsi(1,lidar_test(i,1))+1;
 end
end

OAhsi = n / size(lidar_test,1);
for i=1:class
PAhsi(1,i) = lable_numberhsi(1,i)/class_number(1,i);
UAhsi(1,i) = lable_numberhsi(1,i)/(user_numberhsi(1,i));
end
xii =0;
xxii =0;
for i=1:class
    xii = xii+lable_numberhsi(1,i);
    xxii = xxii+class_number(1,i)*user_numberhsi(1,i);
end
Kappahsi = (size(lidar_test,1)*xii-xxii)/(size(lidar_test,1)*size(lidar_test,1)-xxii);