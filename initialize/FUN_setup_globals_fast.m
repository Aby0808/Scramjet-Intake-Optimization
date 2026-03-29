function FUN_setup_globals_fast(BL, FR)
% setter function for setting glaobal variables
% called by all parellel workers
    global BL_SHAPE_PARAM FRSTM_TH_PARAM
    BL_SHAPE_PARAM = BL;
    FRSTM_TH_PARAM = FR;
end
