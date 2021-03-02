function [TS_T_ann, extra_e] = PerturbT(TS_T_ann_interim, deltaT, info)
    
    % created 08/10/2019 by Keirnan Fowler, University of Melbourne
    
    % add the temperature
    TS_T_ann = table2array(TS_T_ann_interim)+deltaT;
    
    % convert TS_T_ann to a table
    TS_T_ann = array2table(TS_T_ann, 'VariableNames', info.SubareaList);
    
    extra_e = [];
    
end