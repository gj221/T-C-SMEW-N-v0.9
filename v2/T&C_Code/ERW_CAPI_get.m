function s = ERW_CAPI_get()
%% Return the init-seed capture struct, or empty if nothing was recorded.
global ERW_CAPI %#ok<GVMIS>
s = ERW_CAPI;
if isempty(s); s = struct(); end
end
