function FR = FUN_get_FRSTM_TH_PARAM()
% Cached builder for freestream/throat parameter data.
% Returning the named struct keeps the call sites readable and avoids
% repeated reconstruction of the same fixed case parameters.

persistent FR_cached

if isempty(FR_cached)
    FR_cached = FUN_set_freestream_throat_params();
end

FR = FR_cached;
end
