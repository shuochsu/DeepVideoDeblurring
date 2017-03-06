function [v_i0, reg_nan_map] = warpToRef(v0, vi, flo_i0, params)

[v_i0, reg_nan_map] = backwardsWarp(vi, flo_i0, 0);                
v_i0(isnan(v_i0))=0;
