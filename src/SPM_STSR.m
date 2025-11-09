function Xhat=SPM_STSR(NoisyM,I,option,par, parBSM)

%% parameter setting
Xhat=zeros(size(NoisyM));
border_ext=option.border_ext;
interval=option.interval;
dim=option.dim;
dims = dim/interval; 
if mod(dim,interval)~=0
    err('interval must be divided by patch size');
end
stride=option.stride;

nChannel=option.nChannel;
buffersize=option.buffersize; %adjacent frame
p_frame=option.p_frame;
Ko=option.p_frame_o;
% min_ob_num=option.min_ob_num;

param.search_window_size=ceil(option.search_window_size/interval);
param.K=p_frame;

% para_TR.max_iter = option.TRmaxiter;
% para_TR.stopc  = option.TRstopc;
% para_TR.disp = option.TRdisp;
% para_TR.adaptiverank=option.adaptiverank;


% initOpt.printitn = 1;
% initOpt.maxiters = 100;
% initOpt.tol = 1e-8;

ndim=length(size(NoisyM));
sizeD=[dim,dim,nChannel,p_frame];
R           = ones(1,ndim)*6;
R(1)        = sizeD(1)-26;
R(2)        = sizeD(2)-26;
% ttm1D          = @(X,U,k,sizeX,Um)   shiftdim( reshape(U*reshape(  shiftdim(X,k-1) , sizeX(k),[] ), [Um,sizeX(k+1:ndim),sizeX(1:k-1)]), ndim+1-k);
RankChange  = find(R>sizeD);
for i = RankChange(1:end)
    R(i) = min(R(i),sizeD(i));
end





[m, n, ~]=size(NoisyM);
border_ext_row_left=border_ext;
border_ext_row_right=border_ext-mod(n+border_ext*2,interval);
border_ext_col_up=border_ext;
border_ext_col_down=border_ext-mod(m+border_ext*2,interval);

if nChannel==1
    NoisyM=reshape(NoisyM,m,n,nChannel,[]);

    total_frame_num=size(NoisyM,ndims(NoisyM));
%     para_TR.r = option.TRrank*ones(ndims(NoisyM),1);
else
    total_frame_num=size(NoisyM,4);
%     para_TR.r = option.TRrank*ones(4,1);
end
NoisyM=squeeze(NoisyM);


%% initial

m=m+border_ext_col_up+border_ext_col_down;
n=n+border_ext_row_left+border_ext_row_right;
sm=kron(1:interval,ones(1,interval));
sn=kron(ones(1,interval),1:interval);
patchnew=cell(1,4);


%% video patch index
total_patch_size = [m-dim+interval, n-dim+interval];
total_buffer_size = [m-dim+interval, n-dim+interval, buffersize];
patch_ind_adjust = (m-dim+interval) * (n-dim+interval);
total_patch_size_s = [(m-dim+interval)/interval, (n-dim+interval)/interval];

%% initial patch index
m2=m-dim+interval;
n2=n-dim+interval;
border_ext2=floor(border_ext/2);

ind_ref_m =  border_ext2:stride:m2-border_ext2; 
if ind_ref_m(end)<m2-border_ext2
    ind_ref_m=[ind_ref_m m2-border_ext2];
end
ind_ref_n  =  border_ext2:stride:n2-border_ext2;
if ind_ref_n(end)<n2-border_ext2
    ind_ref_n=[ind_ref_n n2-border_ext2];
end

ind_initial=[];
for  i  =  1 : length(ind_ref_m)
    for  j  =  1 : length(ind_ref_n)
        row    =   ind_ref_m(i);
        col     =   ind_ref_n(j);
        ind_initial = [ind_initial sub2ind(total_patch_size, row, col)];
    end
end

num_patch  =  (m-dim+interval) * (n-dim+interval); 
idx_patch  =  1 : num_patch;    
idx_patch   =  reshape(idx_patch, m-dim+interval, n-dim+interval);
for i=1:interval^2   
    si=sm(i);sj=sn(i);     
    cc=mod(si-1,interval)+mod(sj-1,interval)*interval+1;
    blk_arr_temp=idx_patch(si:interval:end,sj:interval:end,:);
    blk_arr_part{cc}=blk_arr_temp(:);
end

%% initial variables
C=zeros(m,n);
blk_arr=zeros(buffersize,length(ind_initial),p_frame);
% Utr_TR=cell(1,length(ind_initial));

AonlineAs=cell(1,length(ind_initial));
AonlinePs=cell(1,length(ind_initial));
Acore=cell(1,length(ind_initial));
APcore=cell(1,length(ind_initial));
AQcore=cell(1,length(ind_initial));
% AonlineAs_N=cell(1,length(ind_initial));
Atmp=cell(1,length(ind_initial));

for t=1:1:total_frame_num
    
    tic
    disp(t)
    disp('**********')
%     if t==4
%         disp('**********')
%     end
    %% geneate ind
    
    buffernum=min(buffersize,t);
    
    if nChannel==1
        NoisyM_t=NoisyM(:,:,t);

    else
        NoisyM_t=NoisyM(:,:,:,t);

    end 
    NoisyM_t = padarray(NoisyM_t,[border_ext_col_up border_ext_row_left],'symmetric','pre');
    NoisyM_t = padarray(NoisyM_t,[border_ext_col_down border_ext_row_right],'symmetric','post');

%     NoisyMse=imdilate(NoisyM_t,se);
    [m, n, ~]=size(NoisyM_t); 
    NoisyMdouble=double(NoisyM_t)/255;

    
%     blk_arr_pre=blk_arr;
    
    if t>buffersize
        % remove first entry
        blk_arr(1,:,:)=[];
        % remove zero columns of blk_arr
        ind_zero=sum(squeeze(blk_arr(:,1:length(AonlineAs),1)),1)~=0;
%         Utr_TR=Utr_TR(ind_zero);
        
        AonlineAs=AonlineAs(ind_zero);
        AonlinePs=AonlinePs(ind_zero);
   
        Acore=Acore(ind_zero);
        APcore=APcore(ind_zero);
        AQcore=AQcore(ind_zero);
%         AonlineAs_N=AonlineAs_N(ind_zero);
        Atmp=Atmp(ind_zero);
        
        
        ind_zero=sum(squeeze(blk_arr(:,:,1)),1)~=0;
        ind_new=ind_new(ind_zero);
%         blk_arr_pre=blk_arr_pre(:,ind_zero,:);
        blk_arr=blk_arr(:,ind_zero,:); 
        % adjust patch index and 
        blk_arr=(blk_arr-patch_ind_adjust).*(blk_arr~=0);
        % add new entry
        blk_arr(buffersize,:,:)=0;
    end
       
    if t==1
        patch_ind = ind_initial;
        NoisyMinterval=NoisyM_t(1:interval:end,1:interval:end,:);
        NoisyMstack=uint8(zeros(size(NoisyMinterval,1),size(NoisyMinterval,2),size(NoisyMinterval,3),interval^2));

        for i=1:interval^2   
            si=sm(i);sj=sn(i);
            cc=mod(si-1,interval)+mod(sj-1,interval)*interval+1;
            NoisyMstack(:,:,:,cc)=NoisyM_t(si:interval:end,sj:interval:end,:);    

        end
        patchnew=gpuArray(image2patch(NoisyMstack, dims));

        else
        patch_ind = ind_new;
    end

    patchold = patchnew;

    ind_new=zeros(1,length(patch_ind));
    patch_ind_s=zeros(1,length(patch_ind));
    patch_ind_c=zeros(1,length(patch_ind));
       
    %% compute new patch locations    
    if t>1
        NoisyMinterval=NoisyM_t(1:interval:end,1:interval:end,:);
        NoisyMstack=uint8(zeros(size(NoisyMinterval,1),size(NoisyMinterval,2),size(NoisyMinterval,3),interval^2));

        for i=1:interval^2   
            si=sm(i);sj=sn(i);
            cc=mod(si-1,interval)+mod(sj-1,interval)*interval+1;
            NoisyMstack(:,:,:,cc)=NoisyM_t(si:interval:end,sj:interval:end,:);    

        end
        patchnew=gpuArray(image2patch(NoisyMstack, dims));

    end

    
    for i =1:length(patch_ind)        
        if patch_ind(i)~=0
            [xx,yy]=ind2sub(total_patch_size,patch_ind(i));
            xxn=floor((xx-1)/interval)+1;
            yyn=floor((yy-1)/interval)+1;
            patch_ind_c(i)=mod(xx-1,interval)+mod(yy-1,interval)*interval+1;
            patch_ind_s(i)=sub2ind(total_patch_size_s, xxn, yyn);           
        end    
    end
    
    blk_patch=zeros(interval^2*p_frame,length(patch_ind_s));
    err_patch=blk_patch;
   
    for i =1:length(patch_ind_s)        
       if patch_ind_s(i)~=0
           refpatch = patchold(:,patch_ind_s(i),patch_ind_c(i));


           [blk_patch_s,error_arr_s] = patch_BM_adjacent_frame_stack(patchnew,  refpatch,  total_patch_size_s, param, patch_ind_s(i)); 
           for s=1:1:interval^2     
                blk_patch((s-1)*p_frame+1:s*p_frame,i)=blk_arr_part{s}(blk_patch_s(:,s));
                err_patch((s-1)*p_frame+1:s*p_frame,i)=error_arr_s(:,s);
           end

       end
    end
    
    for i=1:length(patch_ind)        
        if patch_ind(i)~=0
            [error_sort, ind] = sort(err_patch(:,i));
            if t==1&&blk_patch(ind(1),i)~=patch_ind(i)
                blk_arr(buffernum,i,1)=patch_ind(i)+(buffernum-1)*patch_ind_adjust;
                blk_arr(buffernum,i,2:p_frame)=blk_patch(ind(1:p_frame-1),i)+(buffernum-1)*patch_ind_adjust;
                ind_new(i)=patch_ind(i);
            else
                blk_arr(buffernum,i,:)=blk_patch(ind(1:p_frame),i)+(buffernum-1)*patch_ind_adjust;
                ind_new(i)=blk_patch(ind(1),i);
                if error_sort(1)>option.mistrack||sum(error_sort)==0
                    ind_new(i)=-1;
                end
            end           
       end
    end

    blk_arr(:,ind_new==-1,:)=[];
%     Utr_TR(ind_new==-1)=[]; 
    
    AonlineAs(ind_new==-1)=[]; 
    AonlinePs(ind_new==-1)=[]; 

    Acore(ind_new==-1)=[]; 
    APcore(ind_new==-1)=[]; 
    AQcore(ind_new==-1)=[]; 
%     AonlineAs_N(ind_new==-1)=[]; 
    Atmp(ind_new==-1)=[];
    
    ind_new(ind_new==-1)=[];
    
    
    %% compute overlaps
    C=zeros(m,n);
    for i =1:length(ind_new)      
        if ind_new(i)~=0      
            [xx,yy]=ind2sub(total_patch_size,ind_new(i));
            xxn=xx:min(xx+dim-1,m);
            yyn=yy:min(yy+dim-1,n);        
            C(xxn,yyn)=C(xxn,yyn)+ones(length(xxn),length(yyn));              
        end
    end
          
    %% detect holes and add patches
    BW=imbinarize(C);
    [L,num]= bwlabel(~BW,8);   
    stats = regionprops(L,'BoundingBox');
    
    ind_add=[];
    
    for j=1:1:num
        bbx=struct2array(stats(j));
        ind_ref_m = max(floor(bbx(2)-border_ext2),border_ext2):stride:min(floor(bbx(2)+bbx(4)+border_ext2),m2-border_ext2); 
        ind_ref_n = max(floor(bbx(1)-border_ext2),border_ext2):stride:min(floor(bbx(1)+bbx(3)+border_ext2),n2-border_ext2);
        ind_ref_m=[ind_ref_m m2-border_ext2];
        ind_ref_n=[ind_ref_n n2-border_ext2];
        for ii = 1:length(ind_ref_m)
            for jj = 1:length(ind_ref_n)
            row = ind_ref_m(ii);
            col = ind_ref_n(jj);
            patchbw=BW(row:min(row+dim-1,m),col:min(col+dim-1,n));
            if sum(patchbw(:))<dim^2
                ind_add = [ind_add sub2ind(total_patch_size, row, col)];
            end
            end
        end
    end
        
    n_new=length(ind_new);
    ind_new2=[ind_new ind_add];
    
    
    %% recompute overlaps
    C=zeros(m,n);
    for i =1:length(ind_new2)      
        if ind_new2(i)~=0      
            [xx,yy]=ind2sub(total_patch_size,ind_new2(i));
            xxn=xx:min(xx+dim-1,m);
            yyn=yy:min(yy+dim-1,n);        
            C(xxn,yyn)=C(xxn,yyn)+ones(length(xxn),length(yyn));              
        end
    end
    
    
    over_ind=[];
    
    for i=1:length(ind_new2)        
        if ind_new2(i)~=0
            [xx,yy]=ind2sub(total_patch_size,ind_new2(i));
            xxn=xx:min(xx+dim-1,m);
            yyn=yy:min(yy+dim-1,n);
            if min(min(C(xxn,yyn)))>option.overlap
               C(xxn,yyn)=C(xxn,yyn)-ones(length(xxn),length(yyn)); 
               over_ind=[over_ind i];
            end
       end
    end
    
    over_ind_new=over_ind(over_ind<=n_new);
    over_ind_add=over_ind(over_ind>n_new)-n_new;
    ind_new(over_ind_new)=[];
    blk_arr(:,over_ind_new,:)=[];    
%     Utr_TR(over_ind_new)=[];
    
    AonlineAs(over_ind_new)=[];  
    AonlinePs(over_ind_new)=[];  

    Acore(over_ind_new)=[];  
    APcore(over_ind_new)=[];  
    AQcore(over_ind_new)=[];  
%     AonlineAs_N(over_ind_new)=[];  
    Atmp(over_ind_new)=[]; 
    
    ind_add(over_ind_add)=[];
    
    ind_add_s=zeros(1,length(ind_add));
    ind_add_c=zeros(1,length(ind_add));
    
    for i =1:length(ind_add)        
        if ind_add(i)~=0
            [xx,yy]=ind2sub(total_patch_size,ind_add(i));
            xxn=floor((xx-1)/interval)+1;
            yyn=floor((yy-1)/interval)+1;
            ind_add_c(i)=mod(xx-1,interval)+mod(yy-1,interval)*interval+1;
            ind_add_s(i)=sub2ind(total_patch_size_s, xxn, yyn); 
        end
    end
    
    blk_patch=zeros(interval^2*p_frame,length(ind_add_s));
    err_patch=blk_patch;
    
    for i =1:length(ind_add_s)        
       refpatch = patchnew(:,ind_add_s(i),ind_add_c(i));

       [blk_patch_s,error_arr_s] = patch_BM_adjacent_frame_stack(patchnew, refpatch, total_patch_size_s, param, ind_add_s(i)); 
       for s=1:1:interval^2     
            blk_patch((s-1)*p_frame+1:s*p_frame,i)=blk_arr_part{s}(blk_patch_s(:,s));
            err_patch((s-1)*p_frame+1:s*p_frame,i)=error_arr_s(:,s);
       end
    end
    
    for i=1:length(ind_add)        
        [xx,yy]=ind2sub(total_patch_size,ind_add(i));
        xxn=xx:min(xx+dim-1,m);
        yyn=yy:min(yy+dim-1,n);
        C(xxn,yyn)=C(xxn,yyn)+ones(length(xxn),length(yyn)); 
        [~, ind] = sort(err_patch(:,i));
        
        if blk_patch(ind(1),i)~=ind_add(i)
            blk_arr(buffernum,i+length(ind_new),1)=ind_add(i)+(buffernum-1)*patch_ind_adjust;
            blk_arr(buffernum,i+length(ind_new),2:p_frame)=blk_patch(ind(1:p_frame-1),i)+(buffernum-1)*patch_ind_adjust;
        else
            blk_arr(buffernum,i+length(ind_new),:)=blk_patch(ind(1:p_frame),i)+(buffernum-1)*patch_ind_adjust;  
        end    
        
    end
          
    ind_new=[ind_new ind_add];
    %% 
     
%     Utr_TR_pre=Utr_TR;
%     Utr_TR=cell(1,length(ind_new));
%     Utr_TR(1:length(Utr_TR_pre))=Utr_TR_pre;
    
    AonlineAs_pre=AonlineAs;
    AonlineAs=cell(1,length(ind_new));
    AonlineAs(1:length(AonlineAs_pre))=AonlineAs_pre;

    AonlinePs_pre=AonlinePs;
    AonlinePs=cell(1,length(ind_new));
    AonlinePs(1:length(AonlinePs_pre))=AonlinePs_pre;


    Acore_pre=Acore;
    Acore=cell(1,length(ind_new));
    Acore(1:length(Acore_pre))=Acore_pre;

    APcore_pre=APcore;
    APcore=cell(1,length(ind_new));
    APcore(1:length(APcore_pre))=APcore_pre;

    AQcore_pre=AQcore;
    AQcore=cell(1,length(ind_new));
    AQcore(1:length(AQcore_pre))=AQcore_pre;

%     AonlineAs_N_pre=AonlineAs_N;
%     AonlineAs_N=cell(1,length(ind_new));
%     AonlineAs_N(1:length(AonlineAs_N_pre))=AonlineAs_N_pre;
    
    Atmp_pre=Atmp;
    Atmp=cell(1,length(ind_new));
    Atmp(1:length(Atmp_pre))=Atmp_pre;
    
    
    if t==1
        blk_arr(buffersize,:,:)=blk_arr(1,:,:);
    end
   
    %% update all subtensors 
    EMdouble= NoisyMdouble;
    MaxIter=parBSM.maxIter;
    for  iter = 1 : 1
%     for  iter = 1 : par.deNoisingIter
        
        EMdouble          = EMdouble + par.delta*(NoisyMdouble - EMdouble);
        EMdouble_old=EMdouble;
        blk_arr_out=squeeze(blk_arr(buffersize,:,1)); 
        blk_arr_out_ind=find(blk_arr_out~=0);    


        
        
        if(iter==1)
            Sigma_arr   = par.SigLam*par.nSig * ones(1,length(ind_new)); % First Iteration use the input noise parameter
        else
            
            for kk=1:1:length(blk_arr_out_ind)
                
                [xx,yy,~]=ind2sub(total_buffer_size,blk_arr(buffersize,blk_arr_out_ind(kk),1));
                xxn=xx:min(xx+dim-1,m);
                yyn=yy:min(yy+dim-1,n);
                sliceDifference = (EMdouble(xxn,yyn,:) - NoisyMdouble(xxn,yyn,:)).^2;
                temp=sum(sliceDifference, 'all')/(dim*dim*size(EMdouble,3));
                Sigma_arr(kk)= par.SigLam*sqrt(abs(par.nSig^2- temp)); % estimate local noise variance

            end
            
            parBSM.maxIter = max(parBSM.maxIter-10,15);
    %         clear Npatch Nmsi
    
           %% update blk_arr

            for i=1:interval^2   
                si=sm(i);sj=sn(i);
                cc=mod(si-1,interval)+mod(sj-1,interval)*interval+1;
                NoisyMstack(:,:,:,cc)=EMdouble(si:interval:end,sj:interval:end,:);    

            end
            patchnew=gpuArray(image2patch(NoisyMstack, dims));



            for i =1:length(ind_new)        
                if ind_new(i)~=0
                    [xx,yy]=ind2sub(total_patch_size,ind_new(i));
                    xxn=floor((xx-1)/interval)+1;
                    yyn=floor((yy-1)/interval)+1;
                    patch_ind_c(i)=mod(xx-1,interval)+mod(yy-1,interval)*interval+1;
                    patch_ind_s(i)=sub2ind(total_patch_size_s, xxn, yyn);           
                end    
            end

            blk_patch=zeros(interval^2*p_frame,length(patch_ind_s));
            err_patch=blk_patch;

            for i =1:length(patch_ind_s)        

               refpatch = patchnew(:,patch_ind_s(i),patch_ind_c(i));


               [blk_patch_s,error_arr_s] = patch_BM_adjacent_frame_stack(patchnew,  refpatch,  total_patch_size_s, param, patch_ind_s(i)); 
               for s=1:1:interval^2     
                    blk_patch((s-1)*p_frame+1:s*p_frame,i)=blk_arr_part{s}(blk_patch_s(:,s));
                    err_patch((s-1)*p_frame+1:s*p_frame,i)=error_arr_s(:,s);
               end

            end

            for i=1:length(ind_new)        
                if ind_new(i)~=0
                    [~, ind] = sort(err_patch(:,i));
                    if t==1&&blk_patch(ind(1),i)~=ind_new(i)
                        blk_arr(buffernum,i,1)=ind_new(i)+(buffernum-1)*patch_ind_adjust;
                        blk_arr(buffernum,i,2:p_frame)=blk_patch(ind(1:p_frame-1),i)+(buffernum-1)*patch_ind_adjust;

                    else
                        blk_arr(buffernum,i,:)=blk_patch(ind(1:p_frame),i)+(buffernum-1)*patch_ind_adjust;

                    end           
                end
            end

            blk_arr(:,ind_new==-1,:)=[];

%             if t==1
%                 blk_arr(buffersize,:,:)=blk_arr(1,:,:);
%             end
    
        end
        
       if t==1
            blk_arr(buffersize,:,:)=blk_arr(1,:,:);
        end
 
        parfor kk=1:length(ind_new)
%          for kk=1:length(ind_new)
%             if (mod(kk,30)==0)
%                 disp(kk)
%             end
            tempSigma=Sigma_arr(kk);
            if ind_new(kk)~=0
                arr=squeeze(blk_arr(buffersize,kk,:))';
                arr=arr(:);               
                if nChannel==1
                    patch_noisy=zeros(dim,dim,p_frame);  
                else
                    patch_noisy=zeros(dim,dim,nChannel,p_frame); 
                end
                p_count=1;
                for i=1:1:length(arr)
                    if arr(i)==0
                        continue;
                    end
                    [xx,yy,~]=ind2sub(total_buffer_size,arr(i));
                    xxn=xx:min(xx+dim-1,m);
                    yyn=yy:min(yy+dim-1,n);
                    
                    if nChannel==1
                        patch_noisy(1:length(xxn),1:length(yyn),p_count)=EMdouble(xxn,yyn);
                    else
                        patch_noisy(1:length(xxn),1:length(yyn),:,p_count)=EMdouble(xxn,yyn,:);                           
                    end
                    
                    
%                     patch_noisy(1:length(xxn)*length(yyn),:,p_count)=reshape(EMdouble(xxn,yyn,:),[],nChannel);                         

                    p_count=p_count+1;
                end
                if ~isempty(AonlineAs{kk})&&iter==1
                    if nChannel==1
                        patch_noisy=patch_noisy(:,:,1:Ko);
                    else
                        patch_noisy=patch_noisy(:,:,:,1:Ko);                        
                    end
    %                 Utr_TR{kk} = TRSSD_update(patch_noisy, double(patch_noisy~=-1), Utr_TR{kk}, Ko); 
                    [AonlineAs{kk}, AonlinePs{kk}, Acore{kk}, APcore{kk}, AQcore{kk}, AonlineAlpha] = streaming_tensor_sparsity_regularization(patch_noisy, AonlineAs{kk}, AonlinePs{kk},  Acore{kk}, APcore{kk}, AQcore{kk},1/tempSigma,parBSM,MaxIter);
%                     AonlineAs_N{kk}(end+1:end+Ko,:) = AonlineAlpha;
                    Atmp{kk} = [AonlineAs{kk}, {AonlineAlpha}];
                else 
    %                 Utr_TR{kk} = TRSSD_initial(patch_noisy, double(patch_noisy~=-1), para_TR);  
%                     [core, initAs] = DTucker(tensor(patch_noisy), R, 1e-8, 100, 1);
                    [core, initAs,patch_noisy,R] = tensor_sparsity_regularization(patch_noisy,1/tempSigma,parBSM); % Perform ITS-based tensor recovery on each FBP goup
                    [AonlinePs{kk}, Acore{kk}, APcore{kk}, AQcore{kk}] = TuckerOline_initial(patch_noisy, initAs, tensor(core), R);
                    AonlineAs{kk} = initAs(1:end-1);
%                     AonlineAs_N{kk} = initAs{end};
                    Atmp{kk}=initAs;

                end       
            end
        end 
        
        %% output the first image in the buffer    
%     AonlineAs_N(end+1:end+minibatchSize,:) = AonlineAlpha;
%     Atmp = [AonlineAs, {AonlineAs_N}];
%     time(2, k) = runtime;
%     fitness(2, k) = 1-(norm(tensor(Yt)-full(ttensor(Acore, Atmp)))/norm(tensor(Yt)));
        I_hat=zeros(m,n,nChannel);
        C_hat=zeros(m,n);

%         blk_arr_out=squeeze(blk_arr(buffersize,:,1)); 
%         blk_arr_out_ind=find(blk_arr_out~=0);    

        for kk=1:1:length(blk_arr_out_ind)

    %         Utr_n=Utr_TR{blk_arr_out_ind(kk)};
            patch_hat=double(full(ttensor(Acore{kk}, Atmp{kk})));
    %         TCP2=TCP(Utr_n(2:end));
    %         patch_hat=reshape(T2M_n(Utr_n{1})*T2M_r(TCP2,2)',dim,dim,nChannel,size(Utr_n{end},2));
            for ii=1:1:size(patch_hat,ndims(patch_hat))
                [xx,yy,~]=ind2sub(total_buffer_size,blk_arr(buffersize,blk_arr_out_ind(kk),ii));
                xxn=xx:min(xx+dim-1,m);
                yyn=yy:min(yy+dim-1,n);
                I_hat(xxn,yyn,:)=I_hat(xxn,yyn,:)+patch_hat(1:length(xxn),1:length(yyn),:,ii);
                C_hat(xxn,yyn)=C_hat(xxn,yyn)+ones(length(xxn),length(yyn));     
            end
        end

        timecost=toc;

        if t==1
           fprintf('Output frame, time cost, number of patches: ');
        elseif timecost>10
           fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
        else
           fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
        end   
        fprintf('%3i, %.4fs, %5i', t, timecost, length(blk_arr_out_ind));     
        EMdouble=I_hat./C_hat(:,:);
        EMdouble(isnan(EMdouble))=EMdouble_old(isnan(EMdouble));
        if option.debug
            disp(psnr(EMdouble(border_ext_col_up+1:end-border_ext_col_down,border_ext_row_left+1:end-border_ext_row_right,:),I(:,:,:,t)));
        end
     
        

        
                
    end
    
    %% output the first image in the buffer    
%     AonlineAs_N(end+1:end+minibatchSize,:) = AonlineAlpha;
%     Atmp = [AonlineAs, {AonlineAs_N}];
%     time(2, k) = runtime;
%     fitness(2, k) = 1-(norm(tensor(Yt)-full(ttensor(Acore, Atmp)))/norm(tensor(Yt)));
  
end

fprintf('\n');

Xhat(isnan(Xhat))=0;

Xhat=squeeze(Xhat);

end

