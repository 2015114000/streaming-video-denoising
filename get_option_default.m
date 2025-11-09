function option=get_option_default(sigma)

option=[];

% framenum=opt.framenum;

% if ndims(MissM)==3&&framenum~=1
%     option.nChannel=1;
% else
%     option.nChannel=3;
% end

option.stopc=1e-3;
option.maxitr  = 200;
option.debug=1;
% option.dim=18;
option.interval=3;
option.search_window_size  =  60;
if sigma <= 25
    option.p_frame=30;
    option.dim=36;
elseif sigma <= 38
    option.p_frame=30;
    option.dim=36;
elseif sigma <= 51
    option.p_frame=20;
    option.dim=24;
elseif sigma <= 64
    option.p_frame=30;
    option.dim=9;
elseif sigma <= 76
    option.p_frame=30;
    option.dim=36;
else
    option.p_frame=30;
    option.dim=36;
end
option.stride=floor((option.dim/3*2)); 
option.p_frame_o=10;
option.se=strel('square',3);
option.buffersize=2;
option.overlap=4;
option.mistrack=(sigma*2.5)^2;
option.TRdisp =  0;

   
option.border_ext=20;


   



% option.TRrank = 8;
% option.TRmaxiter = 10;
% option.TRstopc  = 1e-2;

option.adaptiverank = 1;
option.min_ob_num=5;

        









end