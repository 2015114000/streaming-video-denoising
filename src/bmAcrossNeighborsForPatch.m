function [blk_patch_all, err_patch_all, frame_idx_all] = bmAcrossNeighborsForPatch( ...
        refpatch_vect, ref_global_topleft, patchnew_frames, blk_arr_part, total_patch_size, param, p_frame, ...
        us_local, vs_local, curFrameRange, frame_o_eff, interval, center_padded_size)
% bmAcrossNeighborsForPatch - 在多个邻帧上基于光流做局部 BM，返回合并后的候选（全局索引）
%
% 输入（与主流程约定）:
%   refpatch_vect        - 参考 patch（应该和你传给 patch_BM_adjacent_frame_stack 的 refpatch 形式一致）
%   ref_global_topleft   - 参考patch在 patch-grid 上的 top-left 全局索引 (sub2ind(total_patch_size,row,col))
%   patchnew_frames      - cell array，每个元素为对应邻帧的 extractPatchnew（即你原先传给 patch_BM_adjacent_frame_stack 的第一个参数）
%   blk_arr_part         - cell array, blk_arr_part{s} 映射（你原代码的那个）
%   total_patch_size     - [rows_possible, cols_possible] （patch-grid 的可选 top-left 范围）
%   param                - 结构体（包含 search_window_size, K = p_frame, ...）
%   p_frame              - 每帧选择的候选数（K）
%   us_local, vs_local   - cell array, 与 patchnew_frames 顺序对应，光流场（大小匹配 padded image）
%   curFrameRange        - 数组，当前考虑的帧索引（client 代码的 curFrameRange）
%   frame_o_eff          - 实际使用的邻帧数量（<= numel(patchnew_frames)），通常靠近中心优先
%   interval             - 子采样数（interval^2 个 sub-sample）
%   center_padded_size   - [Hpad, Wpad] padded image 的像素尺寸（用于 flow 索引安全性）
%
% 输出:
%   blk_patch_all   - 向量 (interval^2 * p_frame * frame_o_eff) x 1，包含全局 top-left 索引（0 表示空）
%   err_patch_all   - 相同长度的误差值（SSD）
%   frame_idx_all   - 相同长度，指示每个 entry 来自哪一个邻帧（1..frame_o_eff）
%
% 备注: 函数内部会在每个邻帧上调用 patch_BM_adjacent_frame_stack，并将返回的 pos_arr (K x interval^2)
%       用 blk_arr_part{s}(...) 映射为实际的全局 patch 索引，再按 (s, frame) 顺序拼成输出向量。

% 参数检查 & 初始化
if nargin < 12, error('bmAcrossNeighborsForPatch: too few args'); end
if isempty(frame_o_eff)
    frame_o_eff = numel(patchnew_frames);
else
    frame_o_eff = min(frame_o_eff, numel(patchnew_frames));
end

p_total = p_frame * frame_o_eff * (interval^2);
blk_patch_all = zeros(p_total,1,'like',ref_global_topleft);
err_patch_all = inf(p_total,1);
frame_idx_all = zeros(p_total,1);

% center top-left -> center pixel coords (in patch-grid -> pixel mapping)
% total_patch_size is number of valid top-left positions (rows,cols)
[ref_row, ref_col] = ind2sub(total_patch_size, ref_global_topleft);
% compute reference patch center pixel coords in padded image:
half = floor(( (interval>0) * 0) + ( ( (numel(refpatch_vect) > 0) * 0) )); % dummy, not used directly
% To map to flow, we require pixel center: assume your ref_row/ref_col are top-left in pixel units of padded image grid
% So center pixel:
ref_center_r = ref_row + floor((dim_from_refpatch(refpatch_vect)-1)/2);
ref_center_c = ref_col + floor((dim_from_refpatch(refpatch_vect)-1)/2);

% We'll fill in output by iterating over sub-samples s then frame bi
outPtr = 0;

% local copy of param so we can shrink search window if desired
base_search_window = param.search_window_size;

% iterate neighbor frames (按 curFrameRange 顺序)，但只处理前 frame_o_eff 个
for bi_local = 1:frame_o_eff
    % patchnew for this neighbor
    extractPatchnew_b = patchnew_frames{bi_local};
    % local flow fields (note: us_local/vs_local are cell arrays same size as patchnew_frames)
    u_field = [];
    v_field = [];
    if ~isempty(us_local) && numel(us_local) >= bi_local
        u_field = us_local{bi_local};
    end
    if ~isempty(vs_local) && numel(vs_local) >= bi_local
        v_field = vs_local{bi_local};
    end

    % compute mapped center (use flow at reference center). 如果流不存在，用 0 位移
    if ~isempty(u_field) && ~isempty(v_field) ...
            && size(u_field,1) >= ref_center_r && size(u_field,2) >= ref_center_c
        u = round(u_field(ref_center_r, ref_center_c));
        v = round(v_field(ref_center_r, ref_center_c));
    else
        u = 0; v = 0;
    end
    % mapped center pixel coords in padded image
    mapped_center_r = ref_center_r + v;
    mapped_center_c = ref_center_c + u;
    % clamp mapped center to within padded image bounds
    mapped_center_r = max(1, min(mapped_center_r, center_padded_size(1)));
    mapped_center_c = max(1, min(mapped_center_c, center_padded_size(2)));

    % compute approximate mapped top-left (in patch-grid coordinates)
    % We need top-left row/col for the patch-grid: subtract half patch size
    p_dim = dim_from_refpatch(refpatch_vect); % helper below, from refpatch size
    mapped_tl_r = mapped_center_r - floor((p_dim-1)/2);
    mapped_tl_c = mapped_center_c - floor((p_dim-1)/2);
    % clamp to valid top-left range
    mapped_tl_r = max(1, min(mapped_tl_r, total_patch_size(1)));
    mapped_tl_c = max(1, min(mapped_tl_c, total_patch_size(2)));
    mapped_patchind = sub2ind(total_patch_size, mapped_tl_r, mapped_tl_c);

    % Optionally shrink search window for this neighbor to local area around mapped center:
    param_local = param;
    % choose a smaller pixel radius (in patch-grid unit) to speed up, but keep >=1
    if isfield(param,'search_window_size') && ~isempty(param.search_window_size)
        % shrink to e.g. half of base, but at least 1
        param_local.search_window_size = max(1, min(param.search_window_size, base_search_window));
    end
    param_local.K = p_frame;

    % call YOUR existing block-matching for this neighbor (returns K x interval^2)
    try
        [pos_arr_b, error_arr_b] = patch_BM_adjacent_frame_stack(extractPatchnew_b, refpatch_vect, total_patch_size, param_local, mapped_patchind);
        % pos_arr_b: either K x interval^2 or K x 1 depending on implementation; handle both
    catch ME
        % 如果出错，回退为把 mapped_patchind 自身作为唯一候选（填充 p_frame 个）
        warning('bmAcrossNeighborsForPatch: patch_BM_adjacent_frame_stack failed for neighbor %d: %s', bi_local, ME.message);
        pos_arr_b = repmat(mapped_patchind, p_frame, interval^2);
        error_arr_b = zeros(p_frame, interval^2);
    end

    % now pos_arr_b is size (K x interval^2) or (K x 1) -- normalize to K x interval^2
    [rpos, cpos] = size(pos_arr_b);
    if cpos == 1 && interval^2 > 1
        % replicate across sub-samples
        pos_arr_b = repmat(pos_arr_b, 1, interval^2);
        error_arr_b = repmat(error_arr_b, 1, interval^2);
    end

    % For each sub-sample s, map the local patch indices (returned by patch_BM...) through blk_arr_part{s}
    K = size(pos_arr_b,1); % should equal p_frame
    for s = 1:interval^2
        idxs_local = pos_arr_b(:, s);      % K x 1
        errs_local = error_arr_b(:, s);    % K x 1

        % blk_arr_part{s} maps local-grid indices to global patch indices: apply it (guard empty / out-of-range)
        if s <= numel(blk_arr_part) && ~isempty(blk_arr_part{s})
            try
                mapped_global = blk_arr_part{s}(idxs_local);
            catch
                % fallback if indexing fails (e.g., out of range)
                mapped_global = idxs_local; % best-effort fallback
            end
        else
            mapped_global = idxs_local;
        end

        % write into output vector at correct offset:
        % ordering: for s =1:interval^2, for bi_local=1:frame_o_eff, we place p_frame entries each
        base = (s-1)* (p_frame * frame_o_eff) + (bi_local-1)*p_frame;
        for k = 1:K
            outPtr = outPtr + 1;
            if outPtr > p_total
                break;
            end
            blk_patch_all(outPtr) = mapped_global(k);
            err_patch_all(outPtr) = errs_local(k);
            frame_idx_all(outPtr) = bi_local;
        end
        if outPtr >= p_total, break; end
    end
    if outPtr >= p_total, break; end
end

% if fewer entries filled, remaining keep zero/inf as initialized
end

%% helper: infer patch dimension from refpatch vector shape
function d = dim_from_refpatch(refpatch)
    % refpatch may be (patchdim x 1 x interval^2) or (patchdim x 1)
    sz = size(refpatch);
    % patchdim should be e.g. dim*dim*nChannel
    patchdim = sz(1);
    % Try to infer square patch: find integer dim such that dim^2 divides patchdim (assume nChannel small)
    % We prefer to solve for dim assuming nChannel known (but we don't have nChannel here).
    % A reasonable heuristic: find largest integer d where d^2 <= patchdim and divides patchdim
    d = round(sqrt(patchdim));
    if d^2 ~= patchdim
        % fallback: search divisors up to sqrt
        found = false;
        for cand = floor(sqrt(patchdim)):-1:1
            if mod(patchdim, cand) == 0
                d = cand;
                found = true; break;
            end
        end
        if ~found, d = floor(sqrt(patchdim)); end
    end
end
