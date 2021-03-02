function [TS_P_ann, extra_d] = PerturbPmean(TS_P_ann_interim, deltaP, info)

    % created 08/10/2019 by Keirnan Fowler, University of Melbourne
    
    % apply scaling factor
    TS_P_ann = table2array(TS_P_ann_interim)*(1+deltaP);
    
    % convert TS_P_ann to a table
    TS_P_ann = array2table(TS_P_ann, 'VariableNames', info.SubareaList);
    
    extra_d = [];
    
end