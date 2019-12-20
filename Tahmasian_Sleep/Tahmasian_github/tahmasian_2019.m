%  Tahmasian 2019
%
%  programs: 
%  surfstat ()
%  nifti
%  spearman
%  pls
%  colorbrewer

% Project directory
savepathsleep = '/Users/sofievalk/Dropbox/Juelich/2019.Sleep/Tahmasian_github/'
mkdir(savepathsleep)

load([savepathsleep, 'Schaefer_parcels200_17.mat'])
load([savepathsleep, 'fsaverage5_surf.mat'])
load([savepathsleep, 'yeo200.mat'])
load([savepathsleep, 'dlabel_HCP200_17.mat'])
load([savepathsleep, 'template.mat'])

% Data that should be downloaded (HCP): 
% unrestricted HCP data
% restricted HCP data
% cortical thickness data; following 10 % trimmed mean thickness in each
% parcel is computed using the dlabel Schaefer 200 atlas (dlabel_HCP200_17)
% in 64k dlabel HCP template

% Data that should be downloaded (eNKI): 
% cortical thickness data; following 10 % trimmed mean thickness in each
% parcel is computed using the dlabel Schaefer 200 atlas
% (Schaefer_parcels200_17) in fsaverage5
% behavioral data (age, gender, PSQI, WASI, BMI, BDI)

load([hiddenpathsleep, 'seed_hcp_200.mat'])
load([hiddenpathsleep, 'unres_hcp.mat'])
load([hiddenpathsleep, 'res_hcp.mat'])
load([hiddenpathsleep, 'seed_hcp_mean_ctx.mat'])
load([hiddenpathsleep, 'seed_enki_200.mat'])
load([hiddenpathsleep, 'behavior_enki.mat'])


for supplementary_table1 = 1
      %% posthoc scores
     for i = 1          
         var1       = unres.PSQI_AmtSleep;
         ctrl_var1  = res.DSM_Depr_Pct;
         ctrl_var2  = res.BMI;
         ctrl_var3  = unres.CogTotalComp_Unadj;
         ctrl_var4  = res.SSAGA_Educ;
         keep       = find(mean(HCP200_CT,2)>-1)
         agek       = res.Age_in_Yrs(keep);
         sexk       = cellstr(unres.Gender(keep));
         gb         = global_thick_hcp(keep);
         sex_num    = grp2idx(unres.Gender(keep));
 
         keep_vars = [88:96,116,118,120,121,125,127,129,143,144,156,158,160,164,166, 175:191,508:512]
         
         rx = zeros(length(keep_vars),1);
         px = rx;
         for i =1:length(keep_vars)
             try
                 [r p] = partialcorr(unres.PSQI_Score,unres{:,keep_vars(i)},[D1.Age_in_Yrs, agek.*agek, sex_num, sex_num.*agek],'type','spearman','rows','complete')
                 rx(i,:) = r;
                 px(i,:) = p;
             catch
             end
         end
         
         rx(isnan(rx))=0;
         %find top 20
         absCorrs = abs(rx);
         [sortedCorrs,sortIdx] = sort(absCorrs,'descend');
         top_absCorrs = rx(sortIdx(1:45));
         top_absCorrsp = px(sortIdx(1:45));
         top_vars = unres.Properties.VariableNames(keep_vars(sortIdx(1:45)))';
         
         % Sort again within the top values
         [sortedCorrs2,sortIdx2] = sort(top_absCorrs,'descend');
         mySortedTopCorrs = top_absCorrs(sortIdx2);
         mySortedTopCorrsp = top_absCorrsp(sortIdx2);
         
         mySortedTopVars = top_vars(sortIdx2);
                  
         for iter_corr = 1:45
             disp([mySortedTopVars{iter_corr} ' - r=' num2str(mySortedTopCorrs(iter_corr,1),'%0.2f') ' - p=' num2str(mySortedTopCorrsp(iter_corr,1),'%0.2d')]);
         end
         
         mySortedTopCorrsp_unres = mySortedTopCorrsp;
         
         keep_vars = [13:22, 91,93,95,97,99:101,102,104:108,179:201]
         
         rx = zeros(length(keep_vars),1);
         px = rx;
         for i =1:length(keep_vars)
             try
                 [r p] = partialcorr(unres.PSQI_AmtSleep(keep),res{keep,keep_vars(i)},[agek, agek.*agek, sex_num, sex_num.*agek],'type','spearman','rows','complete')
                 rx(i,:) = r;
                 px(i,:) = p;
             catch
             end
         end
         
         
         rx(isnan(rx))=0;
         
         absCorrs = abs(rx);
         [sortedCorrs,sortIdx] = sort(absCorrs,'descend');
         top_absCorrs = rx(sortIdx(1:46));
         top_absCorrsp = px(sortIdx(1:46));
         top_vars = res.Properties.VariableNames(keep_vars(sortIdx(1:46)))';
         
         % Sort again within the top values
         [sortedCorrs2,sortIdx2] = sort(top_absCorrs,'descend');
         mySortedTopCorrs= top_absCorrs(sortIdx2);
         mySortedTopCorrsp= top_absCorrsp(sortIdx2);
         
         mySortedTopVars = top_vars(sortIdx2);
         
         for iter_corr = 1:46
             disp([mySortedTopVars{iter_corr} ' - r=' num2str(mySortedTopCorrs(iter_corr,1),'%0.2f') ' - p=' num2str(mySortedTopCorrsp(iter_corr,1),'%0.2d')]);
         end
         
         mySortedTopCorrsp_res = mySortedTopCorrsp;
         
         [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh([mySortedTopCorrsp_unres; mySortedTopCorrsp_res],0.05)

     end
end

% behavioral analysis of table 1, as well as post-hoc analysis on clusters
% from figure 1 and supplemmentary figure 1
for table1 = 1
    % hcp unrelated sample
for HCP_unrelated = 1   
    % depression and sleep amount
    for HCP = 1
        [tst,idx] = unique(res.Family_ID);
        tmp3      = zeros(1206,1);
        tmp3(idx) = 1;    
        variabl   = unres.PSQI_AmtSleep;
        z_variabl = zscore(variabl);
        keep      = mintersect(find(mean(HCP200_CT,2)>0), find(tmp3>0), find(~isnan(res.DSM_Depr_Pct)))
        vark      = (variabl(keep));
        agek      = res.Age_in_Yrs(keep);
        sexk      = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = unres.PSQI_Score(keep);
        depre     = unres.PSQI_SleepMeds(keep);
         
        [r p] = partialcorri(res.DSM_Depr_Pct(keep),vark, [sex_num, agek.*agek, sex_num.*agek],'type','spearman','rows','complete')
          
        [a b c d e] = Spearman(res.DSM_Depr_Pct(keep),vark,0,0.05)
        
        [r p] = partialcorri(res.DSM_Depr_Pct(keep),res.BMI(keep), [sex_num, agek.*agek, sex_num.*agek],'type','spearman','rows','complete')
        [r p] = partialcorri(res.DSM_Depr_Pct(keep),unres.CogTotalComp_Unadj(keep), [sex_num, agek.*agek, sex_num.*agek],'type','spearman','rows','complete')
        [r p] = partialcorri(res.BMI(keep),unres.CogTotalComp_Unadj(keep), [sex_num, agek.*agek, sex_num.*agek],'type','spearman','rows','complete')      
    end
    
    %DEPRESSION - score
    for HCP = 1
        [tst,idx] = unique(res.Family_ID);
        tmp3      = zeros(1206,1);
        tmp3(idx) = 1;    
        variabl   = unres.PSQI_Score;
        z_variabl = zscore(variabl);
        keep      = mintersect(find(mean(HCP200_CT,2)>0), find(tmp3>0), find(~isnan(res.DSM_Depr_Pct)))
        vark      = (variabl(keep));
        agek      = res.Age_in_Yrs(keep);
        sexk      = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = unres.PSQI_Score(keep);
        depre     = unres.PSQI_SleepMeds(keep);
         
        [r p] = partialcorri(res.DSM_Depr_Pct(keep),vark, [sex_num, agek.*agek, sex_num.*agek],'type','spearman','rows','complete')      
        [a b c d e] = Spearman(res.DSM_Depr_Pct(keep),vark,0,0.05)
        
        [r p] = partialcorri(HCP200_CT(keep,88),res.DSM_Depr_Pct(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,88), res.DSM_Depr_Pct(keep),0,0.05)
        
        [r p] = partialcorri(HCP200_CT(keep,91), res.DSM_Depr_Pct(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,91), res.DSM_Depr_Pct(keep),0,0.05)
        
        [r p] = partialcorri(HCP200_CT(keep,192), res.DSM_Depr_Pct(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,192), res.DSM_Depr_Pct(keep),0,0.05)

         [r p] = partialcorri(HCP200_CT(keep,30), res.DSM_Depr_Pct(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,30), res.DSM_Depr_Pct(keep),0,0.05)

         [r p] = partialcorri(HCP200_CT(keep,106), res.DSM_Depr_Pct(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,106), res.DSM_Depr_Pct(keep),0,0.05)

         [r p] = partialcorri(HCP200_CT(keep,131), res.DSM_Depr_Pct(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,131), res.DSM_Depr_Pct(keep),0,0.05)

    end
    
    %BMI - amount
    for HCP = 1
        [tst,idx] = unique(res.Family_ID);
        tmp3 = zeros(1206,1);
        tmp3(idx) = 1;    
        variabl  = unres.PSQI_AmtSleep;
        z_variabl = zscore(variabl);
        keep    = mintersect(find(mean(HCP200_CT,2)>0), find(tmp3>0), find(res.BMI))
        vark   = (variabl(keep));
        agek   = res.Age_in_Yrs(keep);
        sexk   = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = (unres.PSQI_Score(keep));
        depre     = unres.PSQI_SleepMeds(keep);
         
        [r p] = partialcorri(res.BMI(keep),vark, [sex_num, agek.*agek, sex_num.*agek],'type','spearman','rows','complete')       
        [a b c d e] = Spearman(res.BMI(keep),vark,0,0.05)
    end
    %BMI - score
    for HCP = 1
        [tst,idx] = unique(res.Family_ID);
        tmp3 = zeros(1206,1);
        tmp3(idx) = 1;    
        variabl  = unres.PSQI_Score;
        z_variabl = zscore(variabl);
        keep    = mintersect(find(mean(HCP200_CT,2)>0), find(tmp3>0), find(res.BMI))
        vark   = (variabl(keep));
        agek   = res.Age_in_Yrs(keep);
        sexk   = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = (unres.PSQI_Score(keep));
        depre     = unres.PSQI_SleepMeds(keep);
         
        [r p] = partialcorri(res.BMI(keep),vark, [sex_num, agek.*agek, sex_num.*agek],'type','spearman','rows','complete')       
        [a b c d e] = Spearman(res.BMI(keep),vark,0,0.05)
        
        [r p] = partialcorri(HCP200_CT(keep,88), res.BMI(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,88), res.BMI(keep),0,0.05)
        
        [r p] = partialcorri(HCP200_CT(keep,91), res.BMI(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,91), res.BMI(keep),0,0.05)
        
        [r p] = partialcorri(HCP200_CT(keep,192), res.BMI(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,192), res.BMI(keep),0,0.05)
        
           [r p] = partialcorri(HCP200_CT(keep,30), res.BMI(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,30), res.BMI(keep),0,0.05)
  
           [r p] = partialcorri(HCP200_CT(keep,106), res.BMI(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,106), res.BMI(keep),0,0.05)
  
           [r p] = partialcorri(HCP200_CT(keep,131), res.BMI(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,131), res.BMI(keep),0,0.05)
  
    end
    
    %IQ
    for HCP = 1
        [tst,idx] = unique(res.Family_ID);
        tmp3 = zeros(1206,1);
        tmp3(idx) = 1;    
        variabl  = unres.PSQI_AmtSleep;
        z_variabl = zscore(variabl);
        keep    = mintersect(find(mean(HCP200_CT,2)>0), find(tmp3>0), find(~isnan(unres.CogTotalComp_Unadj)))
        vark   = (variabl(keep));
        agek   = res.Age_in_Yrs(keep);
        sexk   = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = (unres.PSQI_Score(keep));
        depre     = unres.PSQI_SleepMeds(keep);
         
        [r p] = partialcorri(unres.CogTotalComp_Unadj(keep),vark, [sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(unres.CogTotalComp_Unadj(keep),vark,0,0.05)
    end    
    for HCP = 1
        [tst,idx] = unique(res.Family_ID);
        tmp3 = zeros(1206,1);
        tmp3(idx) = 1;    
        variabl  = unres.PSQI_Score;
        z_variabl = zscore(variabl);
        keep    = mintersect(find(mean(HCP200_CT,2)>0), find(tmp3>0), find(~isnan(unres.CogTotalComp_Unadj)))
        vark   = (variabl(keep));
        agek   = res.Age_in_Yrs(keep);
        sexk   = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = (unres.PSQI_Score(keep));
        depre     = unres.PSQI_SleepMeds(keep);
         
        [r p] = partialcorri(unres.CogTotalComp_Unadj(keep),vark, [sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(unres.CogTotalComp_Unadj(keep),vark,0,0.05)
        
           
        [r p] = partialcorri(HCP200_CT(keep,88), unres.CogTotalComp_Unadj(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,88), unres.CogTotalComp_Unadj(keep),0,0.05)
        
         [r p] = partialcorri(HCP200_CT(keep,91), unres.CogTotalComp_Unadj(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,91), unres.CogTotalComp_Unadj(keep),0,0.05)
        
         [r p] = partialcorri(HCP200_CT(keep,192), unres.CogTotalComp_Unadj(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,192), unres.CogTotalComp_Unadj(keep),0,0.05)
        
           [r p] = partialcorri(HCP200_CT(keep,30), unres.CogTotalComp_Unadj(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,30), unres.CogTotalComp_Unadj(keep),0,0.05)
      
              [r p] = partialcorri(HCP200_CT(keep,106), unres.CogTotalComp_Unadj(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,106), unres.CogTotalComp_Unadj(keep),0,0.05)
    
              [r p] = partialcorri(HCP200_CT(keep,131), unres.CogTotalComp_Unadj(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,131), unres.CogTotalComp_Unadj(keep),0,0.05)
    
        [r p] = partialcorr([res.DSM_Depr_Pct(keep),res.BMI(keep), unres.CogTotalComp_Unadj(keep)],[sex_num, agek.*agek, sex_num.*agek],'rows','complete','type','spearman')

    end    
end
    % enki sample
for eNKI = 1
    % depression score
    for k = 1
        variabl  = eNKI_beh.PSQI_04;
        variabl2 = eNKI_beh.depression;
        keep     = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl2)),  find(~isnan(variabl))); %,find(eNKI_beh.age>20),find(eNKI_beh.age<40)
        vark   = (variabl(keep));
        agek   = eNKI_beh.age(keep);
        sexk   = eNKI_beh.gender(keep);
        gb     = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));
        ctrl    = eNKI_beh.DAY_LAG(keep);
       
        [r p] = partialcorri(vark, variabl2(keep),[sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(vark, variabl2(keep),0,0.05)     
    end
    for k = 1
        variabl  = eNKI_beh.PSQI_18;
        variabl2 = eNKI_beh.depression;
        keep     = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl2)),  find(~isnan(variabl))); %,find(eNKI_beh.age>20),find(eNKI_beh.age<40)
        vark   = (variabl(keep));
        agek   = eNKI_beh.age(keep);
        sexk   = eNKI_beh.gender(keep);
        gb     = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));
       
        [r p] = partialcorri(vark, variabl2(keep),[sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(vark, variabl2(keep),0,0.05) 
        
        M = 1+term(gb) + term(agek) + term(sex_num) + term(agek)*term(agek) + term(agek)*term(sex_num);
        slm = SurfStatLinMod(eNKI200_CT(:,keep)',M);
        eNKI200_CTgb = eNKI200_CT(:,keep)' - slm.X*slm.coef;

        [r p] = partialcorri(variabl2(keep),eNKI200_CT(88,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CTgb(:,88),0,0.05)   

        [r p] = partialcorri(variabl2(keep),eNKI200_CT(91,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CT(91,keep)',0,0.05)   

        [r p] = partialcorri(variabl2(keep),eNKI200_CT(192,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CTgb(:,192),0,0.05)   
        
        [r p] = partialcorri(variabl2(keep),eNKI200_CT(30,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CT(30,keep)',0,0.05)   
        
        [r p] = partialcorri(variabl2(keep),eNKI200_CT(106,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CT(106,keep)',0,0.05)   
        
        [r p] = partialcorri(variabl2(keep),eNKI200_CT(131,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CT(131,keep)',0,0.05)   

    end
    % bmi 
    for k = 1
        variabl  = eNKI_beh.PSQI_04;
        variabl2 = eNKI_beh.bmi;
        keep     = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl2)), find((variabl2>0)), find(~isnan(variabl))); %,find(eNKI_beh.age>20),find(eNKI_beh.age<40)
        vark   = (variabl(keep));
        agek   = eNKI_beh.age(keep);
        sexk   = eNKI_beh.gender(keep);
        gb     = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));
       
        [r p] = partialcorri(vark, variabl2(keep),[sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(vark, variabl2(keep),0,0.05)     
    end
    for k = 1
        variabl  = eNKI_beh.PSQI_18;
        variabl2 = eNKI_beh.bmi;
        keep     = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl2)), find((variabl2>5)), find(~isnan(variabl))); %,find(eNKI_beh.age>20),find(eNKI_beh.age<40)
        vark   = (variabl(keep));
        agek   = eNKI_beh.age(keep);
        sexk   = eNKI_beh.gender(keep);
        gb     = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));
       
        [r p] = partialcorri(vark, variabl2(keep),[sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(vark, variabl2(keep),0,0.05)
        
        M = 1+term(gb) + term(agek) + term(sex_num) + term(agek)*term(agek) + term(agek)*term(sex_num);
        slm = SurfStatLinMod(eNKI200_CT(:,keep)',M);
        eNKI200_CTgb = eNKI200_CT(:,keep)' - slm.X*slm.coef;

        [r p] = partialcorri(variabl2(keep),eNKI200_CT(88,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CTgb(:,88),0,0.05)   

        [r p] = partialcorri(variabl2(keep),eNKI200_CT(91,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CTgb(:,91),0,0.05)   

        [r p] = partialcorri(variabl2(keep),eNKI200_CT(192,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CTgb(:,192),0,0.05)   
        
        [r p] = partialcorri(variabl2(keep),eNKI200_CT(30,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CTgb(:,30),0,0.05)   
        
        [r p] = partialcorri(variabl2(keep),eNKI200_CT(106,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CTgb(:,106),0,0.05)   
        
        [r p] = partialcorri(variabl2(keep),eNKI200_CT(131,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CTgb(:,131),0,0.05)   

    end
    % iq
    for k = 1
        variabl  = eNKI_beh.sleep_amount;
        variabl2 = eNKI_beh.iq;
        keep     = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl2)), find((variabl2>0)), find(~isnan(variabl))); %,find(eNKI_beh.age>20),find(eNKI_beh.age<40)
        vark   = (variabl(keep));
        agek   = eNKI_beh.age(keep);
        sexk   = eNKI_beh.gender(keep);
        gb     = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));
       
        eNKI_iq = variabl2(keep);
        [r p] = partialcorri(vark, variabl2(keep),[sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(vark, variabl2(keep),0,0.05)     
    end
    for k = 1
        variabl  = eNKI_beh.sleep_total;
        variabl2 = eNKI_beh.iq;
        keep     = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(eNKI_beh.depression)), find(~isnan(eNKI_beh.bmi)), find(~isnan(eNKI_beh.iq)),find((variabl2>0)), find(~isnan(variabl))); %,find(eNKI_beh.age>20),find(eNKI_beh.age<40)
        vark   = (variabl(keep));
        agek   = eNKI_beh.age(keep);
        sexk   = eNKI_beh.gender(keep);
        gb     = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));
       
        [r p] = partialcorri(vark, variabl2(keep),[sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(vark, variabl2(keep),0,0.05)     
        
        M = 1+term(gb) + term(agek) + term(sex_num) + term(agek)*term(agek) + term(agek)*term(sex_num);
        slm = SurfStatLinMod(eNKI200_CT(:,keep)',M);
        eNKI200_CTgb = eNKI200_CT(:,keep)' - slm.X*slm.coef;

      
        [r p] = partialcorri(variabl2(keep),eNKI200_CT(88,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CTgb(:,88),0,0.05)   

        [r p] = partialcorri(variabl2(keep),eNKI200_CT(91,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CTgb(:,91),0,0.05)   

        [r p] = partialcorri(variabl2(keep),eNKI200_CT(192,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CTgb(:,192),0,0.05)   
        
        [r p] = partialcorri(variabl2(keep),eNKI200_CT(30,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CTgb(:,30),0,0.05)   
        
        [r p] = partialcorri(variabl2(keep),eNKI200_CT(106,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CTgb(:,106),0,0.05)   
        
        [r p] = partialcorri(variabl2(keep),eNKI200_CT(131,keep)', [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(variabl2(keep),eNKI200_CTgb(:,131),0,0.05)   

        [r p] = partialcorr([eNKI_beh.depression(keep),eNKI_beh.bmi(keep), eNKI_beh.iq(keep)],[sex_num, agek.*agek, sex_num.*agek],'rows','complete','type','spearman')
    end
end
    % hcp related sample
for HCPl = 1   
    %DEPRESSIOn - amount
    for HCP = 1
        [tst,idx] = unique(res.Family_ID);
        tmp3 = zeros(1206,1);
        tmp3(idx) = 1;    
        variabl   = unres.PSQI_AmtSleep;
        z_variabl = zscore(variabl);
        keep      = mintersect(find(mean(HCP200_CT,2)>0),  find(~isnan(res.DSM_Depr_Pct)))
        vark      = (variabl(keep));
        agek      = res.Age_in_Yrs(keep);
        sexk      = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = (unres.PSQI_Score(keep));
        depre     = unres.PSQI_SleepMeds(keep);
         
        [r p] = partialcorri(res.DSM_Depr_Pct(keep),vark, [sex_num, agek.*agek, sex_num.*agek],'type','spearman','rows','complete')
          
        [a b c d e] = Spearman(res.DSM_Depr_Pct(keep),vark,0,0.05)
        
    end
    %DEPRESSION - score
    for HCP = 1
        [tst,idx] = unique(res.Family_ID);
        tmp3 = zeros(1206,1);
        tmp3(idx) = 1;    
        variabl   = unres.PSQI_Score;
        z_variabl = zscore(variabl);
        keep      = mintersect(find(mean(HCP200_CT,2)>0),  find(~isnan(res.DSM_Depr_Pct)))
        vark      = (variabl(keep));
        agek      = res.Age_in_Yrs(keep);
        sexk      = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = (unres.PSQI_Score(keep));
        depre     = unres.PSQI_SleepMeds(keep);
        depre_HCP = res.DSM_Depr_Pct(keep);
         
        [r p] = partialcorri(res.DSM_Depr_Pct(keep),vark, [sex_num, agek.*agek, sex_num.*agek],'type','spearman','rows','complete')      
        [a b c d e] = Spearman(res.DSM_Depr_Pct(keep),vark,0,0.05)
        
        [r p] = partialcorri(HCP200_CT(keep,88),res.DSM_Depr_Pct(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,88), res.DSM_Depr_Pct(keep),0,0.05)
        
        [r p] = partialcorri(HCP200_CT(keep,91), res.DSM_Depr_Pct(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,91), res.DSM_Depr_Pct(keep),0,0.05)
        
        [r p] = partialcorri(HCP200_CT(keep,192), res.DSM_Depr_Pct(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,192), res.DSM_Depr_Pct(keep),0,0.05)

         [r p] = partialcorri(HCP200_CT(keep,30), res.DSM_Depr_Pct(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,30), res.DSM_Depr_Pct(keep),0,0.05)

         [r p] = partialcorri(HCP200_CT(keep,106), res.DSM_Depr_Pct(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,106), res.DSM_Depr_Pct(keep),0,0.05)

         [r p] = partialcorri(HCP200_CT(keep,131), res.DSM_Depr_Pct(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CT(keep,131), res.DSM_Depr_Pct(keep),0,0.05)

        
    keep0    = keep;
    for i = 1
        cd('/Volumes/BnBfaSSD/USER/Valk/Genetics/Sleep/')
        filename = ['CT_200nodes_DSMdepr.csv']
        id       = res.Subject((keep0))
        age      = res.Age_in_Yrs((keep0));
        sex      = unres.Gender((keep0));
        node     = HCP200_CT(keep0,:);
        total    = unres.PSQI_Score(keep0);
        amount   = unres.PSQI_AmtSleep(keep0);
        sadness  = res.DSM_Depr_Pct(keep0);
        gb       = mean(CTX_fs32k(keep0,:),2);

        T = table(id, age, sex, amount, sadness, total, node, gb)
        writetable(T, filename);         
    end
    end
    
    %BMI - amount
    for HCP = 1
        [tst,idx] = unique(res.Family_ID);
        tmp3 = zeros(1206,1);
        tmp3(idx) = 1;    
        variabl  = unres.PSQI_AmtSleep;
        z_variabl = zscore(variabl);
        keep    = mintersect(find(mean(HCP200_CT,2)>0), find(~isnan(res.BMI)))
        vark   = (variabl(keep));
        agek   = res.Age_in_Yrs(keep);
        sexk   = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = (unres.PSQI_Score(keep));
        depre     = unres.PSQI_SleepMeds(keep);
         BMI_HCP =  res.BMI(keep)
        [r p] = partialcorri(res.BMI(keep),vark, [sex_num, agek.*agek, sex_num.*agek],'type','spearman','rows','complete')       
        [a b c d e] = Spearman(res.BMI(keep),vark,0,0.05)
    end
    %BMI - score
    for HCP = 1
        [tst,idx] = unique(res.Family_ID);
        tmp3 = zeros(1206,1);
        tmp3(idx) = 1;    
        variabl  = unres.PSQI_Score;
        z_variabl = zscore(variabl);
        keep    = mintersect(find(mean(HCP200_CT,2)>0), find(~isnan(res.BMI)))
        vark   = (variabl(keep));
        agek   = res.Age_in_Yrs(keep);
        sexk   = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = (unres.PSQI_Score(keep));
        depre     = unres.PSQI_SleepMeds(keep);
         
        [r p] = partialcorri(gb,vark, [sex_num, agek.*agek, sex_num.*agek],'type','spearman','rows','complete')       

         M = 1+term(gb) + term(agek) + term(sex_num) + term(agek)*term(agek) + term(agek)*term(sex_num);
        slm = SurfStatLinMod(HCP200_CT(keep,:),M);
        HCP200_CTgb = HCP200_CT(keep,:) - slm.X*slm.coef;

        [r p] = partialcorri(res.BMI(keep),vark, [sex_num, agek.*agek, sex_num.*agek],'type','spearman','rows','complete')       
        [a b c d e] = Spearman(res.BMI(keep),vark,0,0.05)
        
        [r p] = partialcorri(HCP200_CT(keep,88), res.BMI(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CTgb(:,88), res.BMI(keep),0,0.05)
        
        [r p] = partialcorri(HCP200_CT(keep,91), res.BMI(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CTgb(:,91), res.BMI(keep),0,0.05)
        
        [r p] = partialcorri(HCP200_CT(keep,192), res.BMI(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CTgb(:,192), res.BMI(keep),0,0.05)
        
          [r p] = partialcorri(HCP200_CT(keep,30), res.BMI(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CTgb(:,30), res.BMI(keep),0,0.05)
      
          [r p] = partialcorri(HCP200_CT(keep,106), res.BMI(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CTgb(:,106), res.BMI(keep),0,0.05)
      
          [r p] = partialcorri(HCP200_CT(keep,131), res.BMI(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CTgb(:,131), res.BMI(keep),0,0.05)
      
            keep0    = keep;
    for i = 1
        cd('/Volumes/BnBfaSSD/USER/Valk/Genetics/Sleep/')
        filename = ['CT_200nodes_BMI.csv']
        id       = res.Subject((keep0))
        age      = res.Age_in_Yrs((keep0));
        sex      = unres.Gender((keep0));
        node     = HCP200_CT(keep0,:);
        total    = unres.PSQI_Score(keep0);
        amount   = unres.PSQI_AmtSleep(keep0);
        bmi      = res.BMI(keep0);
        gb       = mean(CTX_fs32k(keep0,:),2);

        T = table(id, age, sex, amount, bmi, total, node, gb)
        writetable(T, filename);         
    end
    end
    
    
    %IQ
    for HCP = 1
        [tst,idx] = unique(res.Family_ID);
        tmp3 = zeros(1206,1);
        tmp3(idx) = 1;    
        variabl  = unres.PSQI_AmtSleep;
        z_variabl = zscore(variabl);
        keep    = mintersect(find(mean(HCP200_CT,2)>0), find(~isnan(unres.CogTotalComp_Unadj)))
        vark   = (variabl(keep));
        agek   = res.Age_in_Yrs(keep);
        sexk   = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = (unres.PSQI_Score(keep));
        depre     = unres.PSQI_SleepMeds(keep);
         
        [r p] = partialcorri(unres.CogTotalComp_Unadj(keep),vark, [sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(unres.CogTotalComp_Unadj(keep),vark,0,0.05)
    end    
    for HCP = 1
        [tst,idx] = unique(res.Family_ID);
        tmp3 = zeros(1206,1);
        tmp3(idx) = 1;    
        variabl  = unres.PSQI_Score;
        z_variabl = zscore(variabl);
        keep    = mintersect(find(mean(HCP200_CT,2)>0),  find(~isnan(unres.CogTotalComp_Unadj)))
        vark   = (variabl(keep));
        agek   = res.Age_in_Yrs(keep);
        sexk   = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = (unres.PSQI_Score(keep));
        depre     = unres.PSQI_SleepMeds(keep);
        
%         H =
%         
%         1
%         
%         
%         P =
%         
%         4.5099e-165
%         
%         
%         CI =
%         
%         18.6800
%         21.2598
%         
%         
%         STATS =
%         
%         struct with fields:
%         
%         tstat: 30.3628
%         df: 1877
%         sd: 14.0558
         
        [r p] = partialcorri(unres.CogTotalComp_Unadj(keep),vark, [sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(unres.CogTotalComp_Unadj(keep),vark,0,0.05)
        
        M = 1+term(gb) + term(agek) + term(sex_num) + term(agek)*term(agek) + term(agek)*term(sex_num);
        slm = SurfStatLinMod(HCP200_CT(keep,:),M);
        HCP200_CTgb = HCP200_CT(keep,:) - slm.X*slm.coef;

      
           
        [r p] = partialcorri(HCP200_CT(keep,88), unres.CogTotalComp_Unadj(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CTgb(:,88), unres.CogTotalComp_Unadj(keep),0,0.05)
        
         [r p] = partialcorri(HCP200_CT(keep,91), unres.CogTotalComp_Unadj(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CTgb(:,91), unres.CogTotalComp_Unadj(keep),0,0.05)
        
         [r p] = partialcorri(HCP200_CT(keep,192), unres.CogTotalComp_Unadj(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CTgb(:,192), unres.CogTotalComp_Unadj(keep),0,0.05)
        
        [r p] = partialcorri(HCP200_CT(keep,30), unres.CogTotalComp_Unadj(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CTgb(:,30), unres.CogTotalComp_Unadj(keep),0,0.05)
        
       
        [r p] = partialcorri(HCP200_CT(keep,106), unres.CogTotalComp_Unadj(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CTgb(:,106), unres.CogTotalComp_Unadj(keep),0,0.05)
        
       
        [r p] = partialcorri(HCP200_CT(keep,131), unres.CogTotalComp_Unadj(keep), [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman')
        [a b c d e] = Spearman(HCP200_CTgb(:,131), unres.CogTotalComp_Unadj(keep),0,0.05)
        
       
        [r p] = partialcorr([res.DSM_Depr_Pct(keep),res.BMI(keep), unres.CogTotalComp_Unadj(keep)],[sex_num, agek.*agek, sex_num.*agek],'rows','complete','type','spearman')

    for i = 1
        cd('/Volumes/BnBfaSSD/USER/Valk/Genetics/Sleep/')
        filename = ['CT_200nodes_IQ.csv']
        id       = res.Subject((keep0))
        age      = res.Age_in_Yrs((keep0));
        sex      = unres.Gender((keep0));
        node     = HCP200_CT(keep0,:);
        total    = unres.PSQI_Score(keep0);
        amount   = unres.PSQI_AmtSleep(keep0);
        iq       = unres.CogTotalComp_Unadj(keep0);
        gb       = mean(CTX_fs32k(keep0,:),2);

        T = table(id, age, sex, amount, iq, total, node, gb)
        writetable(T, filename);         
    end
    
    keep0    = keep;
    for i = 1
        cd('/Volumes/BnBfaSSD/USER/Valk/Genetics/Sleep/')
        filename = ['CT_200nodes_Depre_BMI_IQ.csv']
        id       = res.Subject((keep0))
        age      = res.Age_in_Yrs((keep0));
        sex      = unres.Gender((keep0));
        iq      = unres.CogTotalComp_Unadj(keep0);
        depr    = res.DSM_Depr_Pct(keep);
        bmi     = res.BMI(keep0); 
        T = table(id, age, sex, amount, sadness, total, node, gb)
        writetable(T, filename);         
    end
 
    end    
    
    % IQ _ BMI _ DPR
    
   for HCP = 1
        [tst,idx] = unique(res.Family_ID);
        tmp3 = zeros(1206,1);
        tmp3(idx) = 1;    
        variabl  = unres.PSQI_Score;
        z_variabl = zscore(variabl);
        keep    = mintersect(find(mean(HCP200_CT,2)>0),  find(~isnan(unres.CogTotalComp_Unadj)), find(~isnan(res.BMI)), find(~isnan(res.DSM_Depr_Pct)))
        vark   = (variabl(keep));
        agek   = res.Age_in_Yrs(keep);
        sexk   = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = (unres.PSQI_Score(keep));
        depre     = unres.PSQI_SleepMeds(keep);
         

    keep0    = keep;
    for i = 1
        cd('/Volumes/BnBfaSSD/USER/Valk/Genetics/Sleep/')
        filename = ['Depre_BMI_IQ.csv']
        id       = res.Subject((keep0))
        age      = res.Age_in_Yrs((keep0));
        sex      = unres.Gender((keep0));
        iq      = unres.CogTotalComp_Unadj(keep0);
        depr    = res.DSM_Depr_Pct(keep0);
        bmi     = res.BMI(keep0); 
        T = table(id, age, sex, iq, depr, bmi)
        writetable(T, filename);         
    end
 
    end    

end

end

for figure1 = 1
    
    % HCP: amount unrelated
    for k = 1
        [tst,idx]   = unique(res.Family_ID);
        tmp3        = zeros(1206,1);
        tmp3(idx)   = 1;    
        variabl     = unres.PSQI_AmtSleep;
        z_variabl   = zscore(variabl);
        keep        = mintersect(find(mean(HCP200_CT,2)>0), find(tmp3>0));
        vark        = (variabl(keep));
        agek        = res.Age_in_Yrs(keep);
        sexk        = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb          = mean(global_thick_hcp(keep,:),2);
        sex_num     = grp2idx(unres.Gender(keep));
        ctrl        = unres.PSQI_Score(keep);
        depre       = unres.PSQI_SleepMeds(keep);
        
        f = figure, 
        scatter(vark, unres.PSQI_Score(keep), 10, 'k', 'filled'),lsline 
        exportfigbo(f,[RPATH 'Amount_vs_Score_hist.png'],'png', 10)
        close(f)
      
        f = figure, 
        hist(vark), 
        exportfigbo(f,[RPATH 'Amount_hist.png'],'png', 10)
        close(f)
 
        [r p] = partialcorri(HCP200_CT(keep,:),(vark), [gb, sex_num, agek, agek.*agek, sex_num.*agek],'type','spearman','rows','complete');
        min(p)
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.02);
        max(h)
      
        sleep_amountRu = r; 
        sleep_amountPu = p;
        sleep_amountFDRu= h;
  
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN;
            p_valN(p<(0.05/200)) = 1;
            p_val = h;
            
            ct_thick = zeros(1,20484);
              
            sleep = zeros(1,64984);
            for i = 1:200
                sleep(:,find(dlabel_HCP200_17.cdata==i)) = (h(i).*r(i));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(sleep,m.template, '')
            SurfStatColLim([-0.1 0.1])
            colormap(flipud(ctx2))
            exportfigbo(f,[RPATH 'Amount_FDR.png'],'png', 10)
            close(f)
        end
        
        for plot_pvalue_trends = 1
            p_valN = zeros(size(p));
            p_valP = p_valN;
            p_valN(p<=(0.01)) = 1;
            p_val = p_valN + p_valP;
            
            sleep = zeros(1,64984);
            for i = 1:200
                sleep(:,find(dlabel_HCP200_17.cdata==i)) = (p_val(i)).*(r(i));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(sleep,m.template, '')
            SurfStatColLim([-0.1 0.1])
            colormap(flipud(ctx2))
            exportfigbo(f,[RPATH 'Amount_01.png'],'png', 10)
            close(f)
        end
        
    end   
    
    % HCP: score unrelated
    for k = 1
        [tst,idx]   = unique(res.Family_ID);
        tmp3        = zeros(1206,1);
        tmp3(idx)   = 1;    
        variabl     = unres.PSQI_Score;
        z_variabl   = zscore(variabl);
        keep        = mintersect(find(mean(HCP200_CT,2)>0), find(tmp3>0))
        vark        = (variabl(keep));
        agek        = res.Age_in_Yrs(keep);
        sexk        = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = unres.PSQI_Score(keep);
        depre     = unres.PSQI_SleepMeds(keep);
        
        f = figure, 
        hist(vark), 
        exportfigbo(f,[RPATH 'F1.Score_HCP_hist.png'],'png', 10)
        close(f)
        
        [r p] = partialcorri(HCP200_CT(keep,:),vark, [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman','rows','complete');
        min(p)
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05);
        max(h)
      
        sleep_scoreRu = r; 
        sleep_scorePu = p;
        sleep_scoreFDRu = h;
  
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN;
            p_valN(p<(0.05/200)) = 1;
            p_val = h;
            
            ct_thick = zeros(1,20484);
              
            sleep = zeros(1,64984);
            for i = 1:200
                sleep(:,find(dlabel_HCP200_17.cdata==i)) = (h(i).*r(i));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(sleep,m.template, '')
            SurfStatColLim([-0.1 0.1])
            colormap(flipud(ctx2))
            exportfigbo(f,[RPATH 'F1.Score_FDR.png'],'png', 10)
            close(f)
        end
        
        for plot_pvalue_trends = 1
            p_valN = zeros(size(p));
            p_valP = p_valN;
            p_valN(p<=(0.01)) = 1;
            p_val = p_valN + p_valP;
            
            sleep = zeros(1,64984);
            for i = 1:200
                sleep(:,find(dlabel_HCP200_17.cdata==i)) = (p_val(i)).*(r(i));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(sleep,m.template, '')
            SurfStatColLim([-0.1 0.1])
            colormap(flipud(ctx2))
            exportfigbo(f,[RPATH 'F1.Score_01.png'],'png', 10)
            close(f)
        end
        
          % Compute correlation HCP 
         [a, b, c,d,e] = Spearman(unres.PSQI_AmtSleep(keep),unres.PSQI_Score(keep),1,0.05)
 
    end  
    
    % eNKI sleep amount
    for k = 1
        variabl = eNKI_beh.sleep_amount;
        keep    = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl)));%,find(eNKI_beh.age>20),find(eNKI_beh.age<40))
        vark    = (variabl(keep));
        agek    = eNKI_beh.age(keep);
        sexk    = eNKI_beh.gender(keep);
        gb      = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));

        [r p]   = partialcorri(eNKI200_CT(:,keep)',vark,[gb, eNKI_beh.depression(keep), sex_num, agek, agek.*agek, sex_num.*agek],'type','spearman', 'rows','complete');
        min(p)
        
        f = figure, 
        scatter(vark, eNKI_beh.sleep_total(keep), 10, 'k', 'filled'),lsline 
        exportfigbo(f,[RPATH 'Amount_vs_Score_hist_eNKI.png'],'png', 10)
        close(f)
      
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.02);
        max(h)
       
        f = figure, 
        hist(vark), 
        exportfigbo(f,[RPATH 'F1.eNKI.Amount_hist.png'],'png', 10)
        close(f)
        
        eNKI_amount_phenoCT_FDR(:,k) = h;
        eNKI_amount_phenoCT_P(:,k) = p;
        eNKI_amount_phenoCT_R(:,k) = r;
            
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN;
            p_valN(p<(0.05/200)) = 1;
            p_val = h;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'F2pheno', num2str(k), 'neo.FDR.eNKIs.png'],'png', 10)
            close(f)
        end
        
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN
            p_valN(p<(0.01)) = 1;
            p_val = p_valN + p_valP;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'F2pheno', num2str(k), 'neo.p01.eNKIs.png'],'png', 10)
            close(f)
        end
    end
  
    % eNKI score
    for k = 1
        variabl  = eNKI_beh.sleep_total;
        keep     = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl))); %,find(eNKI_beh.age>20),find(eNKI_beh.age<40)
        vark    = (variabl(keep));
        agek    = eNKI_beh.age(keep);
        sexk    = eNKI_beh.gender(keep);
        gb      = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));

        [r p]  = partialcorri(eNKI200_CT(:,keep)',vark,[gb, sex_num, agek, agek.*agek, sex_num.*agek],'type','spearman');
        min(p)

        f = figure, 
        hist(vark), 
        xlim([0 20])
        exportfigbo(f,[RPATH 'F1.eNKI.Score_hist.png'],'png', 10)
        close(f)
      
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05);
        max(h)
        
        eNKI_score_phenoCT_FDR(:,k) = h;
        eNKI_score_phenoCT_P(:,k) = p;
        eNKI_score_phenoCT_R(:,k) = r;
        
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN;
            p_valN(p<(0.05/200)) = 1;
            p_val = h;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'F1.score.FDR.eNKIs.png'],'png', 10)
            close(f)
        end
        
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN
            p_valN(p<(0.01)) = 1;
            p_val = p_valN + p_valP;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'F1.score.p01.eNKIs.png'],'png', 10)
            close(f)
        end
        
        [a, b, c,d,e] = Spearman(eNKI_beh.sleep_amount(keep),eNKI_beh.sleep_total(keep),1,0.05)
    end

end

for figure2 = 1
    for heri200CTX = 1
        tmp = csvread('/Volumes/BnB_TEMP/Sofie/Genetics/Sleep/FN_CTX_her/CTX_200_sleep.csv')
        reodgr = zeros(200,1);
        for i = 1:200
            reodgr(i,:) = find(tmp(:,1)==i);
        end
        
        ctx200_her = tmp(reodgr,:)
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(ctx200_her(:,4),0.05)
        max(h)
    end
    
    % CTX amount
    for k = 1
        sleepy = csvread(['/Volumes/BnB_TEMP/Sofie/Genetics/Sleep/FN_CTX_NODEPR/amount_no_depr.csv']);
        
        skeepy = zeros(size(sleepy));
        for i = 1:200
            tmpi = find(i==sleepy(:,1));
            skeepy(i,:) = sleepy(tmpi,:);
        end
        
        skeepy_amount = skeepy;
        
        ct_thick = zeros(1,20484);
        
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(skeepy(:,5),0.025);
        max(h) 
        
        sleep = zeros(1,64984);
        for i = 1:200
            sleep(:,find(dlabel_HCP200_17.cdata==i)) = (h(i).*skeepy(i,4));
        end
        
        ctx2 = cbrewer('div','RdBu',11)
        f = figure,
        BoSurfStatViewData(sleep,m.template, '')
        SurfStatColLim([-0.5 0.5])
        colormap(flipud(ctx2))
        exportfigbo(f,[RPATH 'F2.FDR_CT_COHER_AMOUNT.png'],'png', 10);
        
        p_mask = skeepy(:,5)<(0.01);
        sleep = zeros(1,64984);
        for i = 1:200
            sleep(:,find(dlabel_HCP200_17.cdata==i)) = skeepy(i,4).*p_mask(i);
        end
        
        f = figure,
        BoSurfStatViewData(sleep,m.template, '')
        SurfStatColLim([-0.5 0.5])
        colormap(flipud(ctx2))
        exportfigbo(f,[RPATH 'F2.CT_COHER_AMOUNT.01.png'],'png', 10);
        
        %% depr
        sleepy = csvread(['/Volumes/BnB_TEMP/Sofie/Genetics/Sleep/FN_CTX/amount_depr.csv']);
        
        skeepy = zeros(size(sleepy));
        for i = 1:200
            tmpi = find(i==sleepy(:,1));
            skeepy(i,:) = sleepy(tmpi,:);
        end
        
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(skeepy(:,5),0.05);
        max(h) 
        
        skeepy_amount = skeepy;
        
        %%BMI
        sleepy = csvread(['/Volumes/BnB_TEMP/Sofie/Genetics/Sleep/FN_BMI/amount_bmi.csv']);
        
        skeepy = zeros(size(sleepy));
        for i = 1:200
            tmpi = find(i==sleepy(:,1));
            skeepy(i,:) = sleepy(tmpi,:);
        end
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(skeepy(:,5),0.05);
        max(h) % 88 %192
        
        %%IQ
        sleepy = csvread(['/Volumes/BnB_TEMP/Sofie/Genetics/Sleep/FN_CTX_IQ/amount_iq.csv']);
        
        skeepy = zeros(size(sleepy));
        for i = 1:200
            tmpi = find(i==sleepy(:,1));
            skeepy(i,:) = sleepy(tmpi,:);
        end
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(skeepy(:,5),0.05);
        max(h) % 88 %192

    end
    
    % CTX total
    for k = 1
        sleepy = csvread(['/Volumes/BnB_TEMP/Sofie/Genetics//Sleep/FN_CTX/total_FN.csv']);
        
        skeepy = zeros(size(sleepy));
        for i = 1:200
            tmpi = find(i==sleepy(:,1));
            skeepy(i,:) = sleepy(tmpi,:);
        end
        
        skeepy_score = skeepy;
        
        ct_thick = zeros(1,20484);
        
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(skeepy(:,5),0.05);
        max(h)
        
        sleep = zeros(1,64984);
        for i = 1:200
            sleep(:,find(dlabel_HCP200_17.cdata==i)) = (h(i).*skeepy(i,2));
        end
        
        ctx2 = cbrewer('div','RdBu',11)
        f = figure,
        BoSurfStatViewData(sleep,m.template, '')
        SurfStatColLim([-0.5 0.5])
        colormap(flipud(ctx2))
        exportfigbo(f,[RPATH 'F2.FDR_CT_COHER_TOTAL01.png'],'png', 10);
        
        
        p_mask = skeepy(:,5)<(0.01);
        sleep = zeros(1,64984);
        for i = 1:200
            sleep(:,find(dlabel_HCP200_17.cdata==i)) = skeepy(i,4).*p_mask(i);
        end
        
        f = figure,
        BoSurfStatViewData(sleep,m.template, '')
        SurfStatColLim([-0.5 0.5])
        colormap(flipud(ctx2))
        exportfigbo(f,[RPATH 'F2.CT_COHER_TOTAL01.png'],'png', 10);
    end
    
end

% supplementary figure 7
for interrelated = 1
    
    %phenotypical
    
    [a,b,c,d,e] = Spearman(sleep_amountRu,eNKI_amount_phenoCT_R,1,0.05)
    [a,b,c,d,e] = Spearman(sleep_amountRall,eNKI_amount_phenoCT_R,1,0.05)
    [a,b,c,d,e] = Spearman(sleep_amountRu,sleep_amountRall,1,0.05)

    [a,b,c,d,e] = Spearman(sleep_scoreRu,eNKI_score_phenoCT_R,1,0.05)
    [a,b,c,d,e] = Spearman(sleep_scoreRall,eNKI_score_phenoCT_R,1,0.05)
    [a,b,c,d,e] = Spearman(sleep_scoreRall,sleep_scoreRu,1,0.05)

    %genetic
    
    [a,b,c,d,e] = Spearman(sleep_amountRu,skeepy_amount(:,4),1,0.05)
    
    [a,b,c,d,e] = Spearman(eNKI_amount_phenoCT_R,skeepy_amount(:,4),1,0.05)
  
    [a,b,c,d,e] = Spearman(sleep_amountRall,skeepy_amount(:,4),1,0.05)
    
    [a,b,c,d,e] = Spearman(sleep_scoreRall,skeepy_score(:,4),1,0.05)
    [a,b,c,d,e] = Spearman(sleep_scoreRu,skeepy_score(:,4),1,0.05)

    [a,b,c,d,e] = Spearman(eNKI_score_phenoCT_R,skeepy_score(:,4),1,0.05)
 
end

% pre-pls
for k = 1
    variabl  = eNKI_beh.sleep_amount;
    ctrl_var1 = eNKI_beh.depression;
    ctrl_var2 = eNKI_beh.bmi;
    ctrl_var3 = eNKI_beh.iq;
    
    keep     = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl)), find(~isnan(ctrl_var1)), find(~isnan(ctrl_var2)), find(~isnan(ctrl_var3)));%,find(eNKI_beh.age>20),find(eNKI_beh.age<40))
    vark    = (variabl(keep));
    agek    = eNKI_beh.age(keep);
    sexk    = eNKI_beh.gender(keep);
    gb      = mean(eNKI200_CT(:,keep))';
    sex_num = grp2idx(eNKI_beh.gender(keep));
    
    sleep_vars_enkik = [variabl(keep), ctrl_var1(keep), ctrl_var2(keep),ctrl_var3(keep),eNKI_beh.sleep_total(keep)];
    
    M   = 1+term(gb) + term(agek) + term(sex_num) + term(agek)*term(agek) + term(agek)*term(sex_num);
    slm = SurfStatLinMod(eNKI200_CT(:,keep)',M);
    eNKI200_CTgb = eNKI200_CT(:,keep)' - slm.X*slm.coef;
end

 % run PLS
run('PLS_code/myPLS_PR4NI.m')    
    
for post_pls_enki = 1
    Lx_enkis = eNKI200_CTgb * V; % Brain scores : original imaging data projected on brain saliences
    Ly_enkis = sleep_vars_enkik * U; % Behavior scores : original behavior data projected on behavior saliences
    
    [r p] = corr(Lx_enkis,Ly_enkis,'type','spearman')
    
    [r p] = corr(Ly_enkis(:,1),sleep_vars_enkik,'type','spearman')
    
    sleep = zeros(1,64984);
    for i = 1:200
        sleep(:,find(dlabel_HCP200_17.cdata==i)) =  Vres(i,1);
    end
    
    ctx2 = cbrewer('div','Spectral',10);
    f = figure,
    BoSurfStatViewData(sleep,m.template, '')
    SurfStatColLim([-3 -2])
    colormap(flipud(ctx2))
    exportfigbo(f,[RPATH, 'Sleep_comp_1n.png'],'png', 10)
    close(f)
    
    sleep = zeros(1,64984);
    for i = 1:200
        sleep(:,find(dlabel_HCP200_17.cdata==i)) =  Vres(i,2);
    end
    
    ctx2 = cbrewer('div','Spectral',10);
    f = figure,
    BoSurfStatViewData(sleep,m.template, '')
    SurfStatColLim([-3 -2])
    colormap(flipud(ctx2))
    exportfigbo(f,[RPATH, 'Sleep_comp_2n.png'],'png', 10)
    close(f)
    
    for y = 1:2
        heri_ct = zeros(1,20484);
        for i = 1:100
            heri_ct(:,find(parcels200_17==i+1)) = Vres(i,y);
        end
        for i = 1:100
            heri_ct(:,find(parcels200_17==i+1001)) = Vres(i+100,y);
        end
        
        heri_ctp = heri_ct<-2;
        for i = 1:100
            xy(i) = max(heri_ctp(find(parcels200==i+1)));
        end        
        for i = 1:100
            xy(i+100) = max(heri_ctp(find(parcels200==i+1001)));
        end
        %
        length(find(yeo200(find(xy>0))==1))
        
        f= figure,
        bar([    (length(find(yeo200(find(xy>0))==1))./29)...
            ,    (length(find(yeo200(find(xy>0))==2))./35)...
            ,    (length(find(yeo200(find(xy>0))==3))./26)...
            ,    (length(find(yeo200(find(xy>0))==4))./22)...
            ,    (length(find(yeo200(find(xy>0))==5))./12)...
            ,    (length(find(yeo200(find(xy>0))==6))./30)...
            ,    (length(find(yeo200(find(xy>0))==7))./46);nan(1,7)],'stacked')
        exportfigbo(f,[RPATH, 'Sleep_comp_', num2str(i),'_stacked_n.png'],'png', 10)
        close(f)
        
        f =figure,
        boxplot(heri_ct(ct_yeo7_200>0)',ct_yeo7_200(ct_yeo7_200>0)')
        exportfigbo(f,[RPATH, 'Sleep_comp_', num2str(i),'.png'],'png', 10)
        close(f)
        
    end
end

% pls is run with PSL_code/myPLS_PR4NI.m

%pls of eNKI in HCP
for k = 1
    [tst,idx] = unique(res.Family_ID);
    tmp3 = zeros(1206,1);
    tmp3(idx) = 1;
    var1 = unres.PSQI_AmtSleep;
    ctrl_var1 = res.DSM_Depr_Pct;
    ctrl_var2 = res.BMI;
    ctrl_var3 = unres.CogTotalComp_Unadj;
    ctrl_var4 = res.SSAGA_Educ;
    keep        = find(mean(HCP200_CT,2)>-1)
    agek   = res.Age_in_Yrs(keep);
    sexk   = cellstr(unres.Gender(keep));
    gb        = mean(global_thick_hcp(keep,:),2);
    sex_num   = grp2idx(unres.Gender(keep));
    
    sleep_vars = [var1, ctrl_var1, ctrl_var2, ctrl_var3, unres.PSQI_Score];
    sleep_varsk = sleep_vars(keep,:)
    M = 1+term(gb) + term(agek) + term(sex_num) + term(agek)*term(agek) + term(agek)*term(sex_num);
    slm = SurfStatLinMod(HCP200_CT(keep,:),M);
    HCP200_CTgb = HCP200_CT(keep,:) - slm.X*slm.coef;
    
     Lx_hcp = HCP200_CTgb * V; % Brain scores : original imaging data projected on brain saliences
     Ly_hcp = sleep_varsk * U; % Behavior scores : original behavior data projected on behavior saliences

     [ r p] = corr(Lx_hcp,Ly_hcp,'type','spearman')
    
     f = figure, 
     scatter(Lx_hcp(:,1),Ly_hcp(:,1),'k','filled'),lsline
     exportfigbo(f,[RPATH, 'Sleep_comp_1_HCP.png'],'png', 10)
     close(f)
   
     f = figure, 
     scatter(Lx_hcp(:,2),Ly_hcp(:,2),'k','filled'),lsline
     exportfigbo(f,[RPATH, 'Sleep_comp_2_HCP.png'],'png', 10)
     close(f)
  
     keep_hcp   = keep;
     for i = 1
         filename = ['CTX200_PLS.csv']
         id       = res.Subject(keep_hcp)
         age      = res.Age_in_Yrs(keep_hcp);
         sex      = unres.Gender(keep_hcp);
         node     = HCP200_CT(keep_hcp,:);
         iq       = unres.CogTotalComp_Unadj(keep_hcp);
         F1 = Lx_hcp(:,1); %brain
         F2 = Lx_hcp(:,2);
         F3 = Ly_hcp(:,1); %behavior
         F4 = Ly_hcp(:,2);
         
         gb     = mean(global_thick_hcp(keep,:),2);
         
         T = table(id, age, sex, node, F1, F2, F3, F4, gb, iq)
         writetable(T, ['/Volumes/BnBfaSSD/USER/Valk/Genetics/Sleep/PLS/', filename]);
     end
    
end
  
% ROIs output for brainmap

rois = [88 192 91] % HCP genetic rois

for i = 1:length(rois)
    schaefer200 = load_nii('/Volumes/BnB_TEMP/Sofie/Masks/Schaefer2018_200Parcels_17Networks_order_FSLMNI152_2mm.nii')
    schaef200   = schaefer200.img(:);
    y = round(schaef200(schaef200>0),3);
    
    seed = (schaef200==rois(i));
    my_nii = schaefer200;
    my_nii.img = reshape(seed,size(my_nii.img));
    save_nii(my_nii,[RPATH 'Volume_', num2str(rois(i)) ,'.nii']);
end


%% supplement
% CTX - amount related
for k = 1
    [tst,idx] = unique(res.Family_ID);
    tmp3 = zeros(1206,1);
    tmp3(idx) = 1;
    variabl  = unres.PSQI_AmtSleep;
    z_variabl = zscore(variabl);
    keep    = mintersect(find(mean(HCP200_CT,2)>0), find(tmp3>=0))
    vark   = (variabl(keep));
    agek   = res.Age_in_Yrs(keep);
    sexk   = cellstr(unres.Gender(keep));
    vark(isinf(vark)) = 0;
    gb        = mean(global_thick_hcp(keep,:),2);
    sex_num   = grp2idx(unres.Gender(keep));
    ctrl      = (unres.PSQI_Score(keep));
    depre     = unres.PSQI_SleepMeds(keep);
    
    f = figure,
    hist(vark),
    exportfigbo(f,[RPATH 'F2.Amount_hist.HCP.all.png'],'png', 10)
    close(f)
    
    [r p] = partialcorri(HCP200_CT(keep,:),vark, [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman','rows','complete');
    min(p)
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05);
    max(h)
    
    sleep_amountRall = r;
    sleep_amountPall = p;
    sleep_amountFDRall= h;
    
    for plot_pvalue = 1
        p_valN = zeros(size(p));
        p_valP = p_valN;
        p_valN(p<(0.05/200)) = 1;
        p_val = h;
        
        ct_thick = zeros(1,20484);
        
        sleep = zeros(1,64984);
        for i = 1:200
            sleep(:,find(dlabel_HCP200_17.cdata==i)) = (h(i).*r(i));
        end
        
        ctx2 = cbrewer('div','Spectral',10);
        f = figure,
        BoSurfStatViewData(sleep,m.template, '')
        SurfStatColLim([-0.1 0.1])
        colormap(flipud(ctx2))
        exportfigbo(f,[RPATH 'F2.Amount_FDR.HCPall.png'],'png', 10)
        close(f)
    end
    
    for plot_pvalue_trends = 1
        p_valN = zeros(size(p));
        p_valP = p_valN;
        p_valN(p<=(0.01)) = 1;
        p_val = p_valN + p_valP;
        
        sleep = zeros(1,64984);
        for i = 1:200
            sleep(:,find(dlabel_HCP200_17.cdata==i)) = (p_val(i)).*(r(i));
        end
        
        ctx2 = cbrewer('div','Spectral',10);
        f = figure,
        BoSurfStatViewData(sleep,m.template, '')
        SurfStatColLim([-0.1 0.1])
        colormap(flipud(ctx2))
        exportfigbo(f,[RPATH 'F2.Amount_01.HCPall.png'],'png', 10)
        close(f)
    end
    
    
end
% CTX - score unrelated
for k = 1
    [tst,idx] = unique(res.Family_ID);
    tmp3 = zeros(1206,1);
    tmp3(idx) = 1;
    variabl  = unres.PSQI_Score;
    z_variabl = zscore(variabl);
    keep    = mintersect(find(mean(HCP200_CT,2)>0), find(tmp3>=0))
    vark   = (variabl(keep));
    agek   = res.Age_in_Yrs(keep);
    sexk   = cellstr(unres.Gender(keep));
    vark(isinf(vark)) = 0;
    gb        = mean(global_thick_hcp(keep,:),2);
    sex_num   = grp2idx(unres.Gender(keep));
    ctrl      = unres.PSQI_Score(keep);
    depre     = unres.PSQI_SleepMeds(keep);
    
    f = figure,
    hist(vark),
    exportfigbo(f,[RPATH 'F2.Score_HCPall_hist.png'],'png', 10)
    close(f)
    
    [r p] = partialcorri(HCP200_CT(keep,:),vark, [gb, sex_num, agek.*agek, sex_num.*agek],'type','spearman','rows','complete');
    min(p)
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05);
    max(h)
    
    sleep_scoreRall = r;
    sleep_scorePall = p;
    sleep_scoreFDRall = h;
    
    for plot_pvalue = 1
        p_valN = zeros(size(p));
        p_valP = p_valN;
        p_valN(p<(0.05/200)) = 1;
        p_val = h;
        
        ct_thick = zeros(1,20484);
        
        sleep = zeros(1,64984);
        for i = 1:200
            sleep(:,find(dlabel_HCP200_17.cdata==i)) = (h(i).*r(i));
        end
        
        ctx2 = cbrewer('div','Spectral',10);
        f = figure,
        BoSurfStatViewData(sleep,m.template, '')
        SurfStatColLim([-0.1 0.1])
        colormap(flipud(ctx2))
        exportfigbo(f,[RPATH 'F2.Score_FDR.HCPall.png'],'png', 10)
        close(f)
    end
    
    for plot_pvalue_trends = 1
        p_valN = zeros(size(p));
        p_valP = p_valN;
        p_valN(p<=(0.01)) = 1;
        p_val = p_valN + p_valP;
        
        sleep = zeros(1,64984);
        for i = 1:200
            sleep(:,find(dlabel_HCP200_17.cdata==i)) = (p_val(i)).*(r(i));
        end
        
        ctx2 = cbrewer('div','Spectral',10);
        f = figure,
        BoSurfStatViewData(sleep,m.template, '')
        SurfStatColLim([-0.1 0.1])
        colormap(flipud(ctx2))
        exportfigbo(f,[RPATH 'F2.Score_FDR.HCPp.png'],'png', 10)
        close(f)
    end
    [a, b, c,d,e] = Spearman(unres.PSQI_AmtSleep(keep),unres.PSQI_Score(keep),1,0.05)
    
end

% controlling for all in F4
for control_all_HCP = 1
   
    [tst,idx] = unique(res.Family_ID);
    tmp3 = zeros(1206,1);
    tmp3(idx) = 1;
    variabl   = unres.PSQI_AmtSleep;
    ctrl_var1 = res.DSM_Depr_Pct;
    ctrl_var2 = res.BMI;
    ctrl_var3 = unres.CogTotalComp_Unadj;
    keep        = mintersect(find(mean(HCP200_CT,2)>0), find(tmp3>0), find(~isnan(ctrl_var1)), find(~isnan(ctrl_var2)), find(~isnan(ctrl_var3)))
    vark        = (variabl(keep));
    agek        = res.Age_in_Yrs(keep);
    sexk        = cellstr(unres.Gender(keep));
    vark(isinf(vark)) = 0;
    gb        = mean(global_thick_hcp(keep,:),2);
    sex_num   = grp2idx(unres.Gender(keep));
    
    [r p] = partialcorri(HCP200_CT(keep,:),vark, [gb, ctrl_var1(keep), ctrl_var2(keep), ctrl_var3(keep), agek, sex_num, agek.*agek, sex_num.*agek],'type','spearman');
    min(p)
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.06);
    max(h)
    [a, b, c, d, e] =   Spearman(r, HCP200_grares,1,0.05)
     for plot_pvalue = 1
        p_valN = zeros(size(p));
        p_valP = p_valN;
        p_valN(p<(0.05/200)) = 1;
        p_val = h;
        
        ct_thick = zeros(1,20484);
        
        sleep = zeros(1,64984);
        for i = 1:200
            sleep(:,find(dlabel_HCP200_17.cdata==i)) = (h(i).*r(i));
        end
        
        ctx2 = cbrewer('div','Spectral',10);
        f = figure,
        BoSurfStatViewData(sleep,m.template, '')
        SurfStatColLim([-0.1 0.1])
        colormap(flipud(ctx2))
        exportfigbo(f,[RPATH 'Amount_FDR.png'],'png', 10)
        close(f)
    end
    
    for plot_pvalue_trends = 1
        p_valN = zeros(size(p));
        p_valP = p_valN;
        p_valN(p<=(0.01)) = 1;
        p_val = p_valN + p_valP;
        
        sleep = zeros(1,64984);
        for i = 1:200
            sleep(:,find(dlabel_HCP200_17.cdata==i)) = (p_val(i)).*(r(i));
        end
        
        ctx2 = cbrewer('div','Spectral',10);
        f = figure,
        BoSurfStatViewData(sleep,m.template, '')
        SurfStatColLim([-0.1 0.1])
        colormap(flipud(ctx2))
        exportfigbo(f,[RPATH 'Amount_01.png'],'png', 10)
        close(f)
    end
end

for contrl_all_eNKI = 1
    variabl   = eNKI_beh.sleep_amount;
    ctrl_var1 = eNKI_beh.depression;
    ctrl_var2 = eNKI_beh.bmi;
    ctrl_var3 = eNKI_beh.iq;
    
    keep     = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl)), find(~isnan(ctrl_var1)), find(~isnan(ctrl_var2)), find(~isnan(ctrl_var3)));%,find(eNKI_beh.age>20),find(eNKI_beh.age<40))
    vark   = (variabl(keep));
    agek   = eNKI_beh.age(keep);
    sexk   = eNKI_beh.gender(keep);
    gb     = mean(eNKI200_CT(:,keep))';
    sex_num = grp2idx(eNKI_beh.gender(keep));

    [r p]  = partialcorri(eNKI200_CT(:,keep)',vark,[gb, ctrl_var1(keep), ctrl_var2(keep), ctrl_var3(keep), sex_num, agek, agek.*agek, sex_num.*agek],'type','spearman');
    min(p)
    
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05);
    max(h)
    
    for plot_pvalue = 1
        p_valN = zeros(size(p));
        p_valP = p_valN;
        p_valN(p<(0.05/200)) = 1;
        p_val = h;
        
        heri_ct = zeros(1,20484);
        for i = 1:100
            heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
        end
        for i = 1:100
            heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
        end
        
        ctx2 = cbrewer('div','Spectral',10);
        f = figure,
        BoSurfStatViewData(heri_ct,SN,'')
        colormap(flipud(ctx2))
        SurfStatColLim([-0.1 0.1])
        exportfigbo(f,[RPATH, 'F2pheno', num2str(k), 'neo.FDR.eNKIs.png'],'png', 10)
        close(f)
    end
    
    for plot_pvalue = 1
        p_valN = zeros(size(p));
        p_valP = p_valN
        p_valN(p<(0.01)) = 1;
        p_val = p_valN + p_valP;
        
        heri_ct = zeros(1,20484);
        for i = 1:100
            heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
        end
        for i = 1:100
            heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
        end
        
        ctx2 = cbrewer('div','Spectral',10);
        f = figure,
        BoSurfStatViewData(heri_ct,SN,'')
        colormap(flipud(ctx2))
        SurfStatColLim([-0.1 0.1])
        exportfigbo(f,[RPATH, 'F2pheno', num2str(k), 'neo.p01.eNKIs.png'],'png', 10)
        close(f)
    end
end

% sleep duration post-hoc HCP

for sleepduration_long_short = 1
   %% R1: CTX - amount related short sleep
    for k = 1
         [tst,idx] = unique(res.Family_ID);
        tmp3 = zeros(1206,1);
        tmp3(idx) = 1;    
        variabl  = unres.PSQI_AmtSleep;
        z_variabl = zscore(variabl);
        keep    = mintersect(find(mean(HCP200_CT,2)>0), find(tmp3>0), find(variabl<9))
        vark   = (variabl(keep));
        agek   = res.Age_in_Yrs(keep);
        sexk   = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = (unres.PSQI_Score(keep));
        depre     = unres.PSQI_SleepMeds(keep);

        [r p] = partialcorri(HCP200_CT(keep,:),vark, [gb, sex_num, agek, agek.*agek, sex_num.*agek],'type','spearman','rows','complete');
        min(p)
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.03);
        max(h)
      
       
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN;
            p_valN(p<(0.05/200)) = 1;
            p_val = h;
            
            ct_thick = zeros(1,20484);
              
            sleep = zeros(1,64984);
            for i = 1:200
                sleep(:,find(dlabel_HCP200_17.cdata==i)) = (h(i).*r(i));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(sleep,m.template, '')
            SurfStatColLim([-0.1 0.1])
            colormap(flipud(ctx2))
            exportfigbo(f,[RPATH 'Amount_FDR_short.png'],'png', 10)
            close(f)
        end
        
        for plot_pvalue_trends = 1
            p_valN = zeros(size(p));
            p_valP = p_valN;
            p_valN(p<=(0.01)) = 1;
            p_val = p_valN + p_valP;
            
            sleep = zeros(1,64984);
            for i = 1:200
                sleep(:,find(dlabel_HCP200_17.cdata==i)) = (p_val(i)).*(r(i));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(sleep,m.template, '')
            SurfStatColLim([-0.1 0.1])
            colormap(flipud(ctx2))
            exportfigbo(f,[RPATH 'Amount_01_short.png'],'png', 10)
            close(f)
        end
        
     
    end  
     
    %% R1: CTX - amount related long sleep
    for k = 1
         [tst,idx] = unique(res.Family_ID);
        tmp3 = zeros(1206,1);
        tmp3(idx) = 1;    
        variabl  = unres.PSQI_AmtSleep;
        z_variabl = zscore(variabl);
        keep    = mintersect(find(mean(HCP200_CT,2)>0), find(tmp3>0), find(variabl>=7))
        vark   = (variabl(keep));
        agek   = res.Age_in_Yrs(keep);
        sexk   = cellstr(unres.Gender(keep));
        vark(isinf(vark)) = 0;
        gb        = mean(global_thick_hcp(keep,:),2);
        sex_num   = grp2idx(unres.Gender(keep));
        ctrl      = (unres.PSQI_Score(keep));
        depre     = unres.PSQI_SleepMeds(keep);

        [r p] = partialcorri(HCP200_CT(keep,:),vark, [gb, sex_num, agek, agek.*agek, sex_num.*agek],'type','spearman','rows','complete');
        min(p)
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.03);
        max(h)
      
       
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN;
            p_valN(p<(0.05/200)) = 1;
            p_val = h;
            
            ct_thick = zeros(1,20484);
              
            sleep = zeros(1,64984);
            for i = 1:200
                sleep(:,find(dlabel_HCP200_17.cdata==i)) = (h(i).*r(i));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(sleep,m.template, '')
            SurfStatColLim([-0.1 0.1])
            colormap(flipud(ctx2))
            exportfigbo(f,[RPATH 'Amount_FDR_long.png'],'png', 10)
            close(f)
        end
        
        for plot_pvalue_trends = 1
            p_valN = zeros(size(p));
            p_valP = p_valN;
            p_valN(p<=(0.01)) = 1;
            p_val = p_valN + p_valP;
            
            sleep = zeros(1,64984);
            for i = 1:200
                sleep(:,find(dlabel_HCP200_17.cdata==i)) = (p_val(i)).*(r(i));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(sleep,template, '')
            SurfStatColLim([-0.1 0.1])
            colormap(flipud(ctx2))
            exportfigbo(f,[RPATH 'Amount_01_long.png'],'png', 10)
            close(f)
        end
   
    end  
     
end
   
% sleep duration post-hoc eNKI

for sleepduration_long_short = 1
        %% short
   for k = 1
        variabl = eNKI_beh.sleep_amount;
        keep    = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl)),find((variabl<9)));%,find(eNKI_beh.age>20),find(eNKI_beh.age<40))
        vark    = (variabl(keep));
        agek    = eNKI_beh.age(keep);
        sexk    = eNKI_beh.gender(keep);
        gb      = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));

        [r p]   = partialcorri(eNKI200_CT(:,keep)',vark,[gb, sex_num, agek, agek.*agek, sex_num.*agek],'type','spearman', 'rows','complete');
        min(p)
 
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.02);
        max(h)
       
           
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN;
            p_valN(p<(0.05/200)) = 1;
            p_val = h;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'F2pheno', num2str(k), 'neo.FDR.eNKIs.short.png'],'png', 10)
            close(f)
        end
        
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN
            p_valN(p<(0.01)) = 1;
            p_val = p_valN + p_valP;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'F2pheno', num2str(k), 'neo.p01.eNKIs.short.png'],'png', 10)
            close(f)
        end
   end
  
        %% long
   for k = 1
        variabl = eNKI_beh.sleep_amount;
        keep    = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl)),find((variabl>=7)));%,find(eNKI_beh.age>20),find(eNKI_beh.age<40))
        vark    = (variabl(keep));
        agek    = eNKI_beh.age(keep);
        sexk    = eNKI_beh.gender(keep);
        gb      = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));
        [r p]   = partialcorri(eNKI200_CT(:,keep)',vark,[gb, sex_num, agek, agek.*agek, sex_num.*agek],'type','spearman', 'rows','complete');
        min(p)
 
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.02);
        max(h)
       
           
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN;
            p_valN(p<(0.05/200)) = 1;
            p_val = h;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'F2pheno', num2str(k), 'neo.FDR.eNKIs.long.png'],'png', 10)
            close(f)
        end
        
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN
            p_valN(p<(0.01)) = 1;
            p_val = p_valN + p_valP;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'F2pheno', num2str(k), 'neo.p01.eNKIs.long.png'],'png', 10)
            close(f)
        end
   end
  
end

% age x brain sleep effects eNKI

for sleepduration_age = 1
        %% eNKI sleep amount youth
    for k = 1
        variabl = eNKI_beh.sleep_amount;
        keep    = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl)),find(eNKI_beh.age<18)) %,find(eNKI_beh.age<40))
        vark    = (variabl(keep));
        agek    = eNKI_beh.age(keep);
        sexk    = eNKI_beh.gender(keep);
        gb      = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));
        [r p]   = partialcorri(eNKI200_CT(:,keep)',vark,[gb, sex_num, agek, agek.*agek, sex_num.*agek],'type','spearman', 'rows','complete');
        min(p)
        
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.02);
        max(h)
          
        eNKI_amount_phenoCT_FDRy(:,k) = h;
        eNKI_amount_phenoCT_Py(:,k) = p;
        eNKI_amount_phenoCT_Ry(:,k) = r;
        
        [r p] = corr(temp,vark,'type','spearman')
    
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN;
            p_valN(p<(0.05/200)) = 1;
            p_val = h;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'F2pheno', num2str(k), 'neo.FDR.eNKIs.png'],'png', 10)
            close(f)
        end
        
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN
            p_valN(p<(0.05)) = 1;
            p_val = p_valN + p_valP;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'F2pheno', num2str(k), 'neo.p01.eNKIs.png'],'png', 10)
            close(f)
        end
    end
  
    %% eNKI sleep amount middle age
    for k = 1
        variabl = eNKI_beh.sleep_amount;
        keep    = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl)),find(eNKI_beh.age>=18),find(eNKI_beh.age<50))
        vark    = (variabl(keep));
        agek    = eNKI_beh.age(keep);
        sexk    = eNKI_beh.gender(keep);
        gb      = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));
        [r p]   = partialcorri(eNKI200_CT(:,keep)',vark,[gb, sex_num, agek, agek.*agek, sex_num.*agek],'type','spearman', 'rows','complete');
        min(p)
        
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05);
        max(h)
        
        eNKI_amount_phenoCT_FDRm(:,k) = h;
        eNKI_amount_phenoCT_Pm(:,k) = p;
        eNKI_amount_phenoCT_Rm(:,k) = r;
   
        [r p] = corr(temp,vark,'type','spearman')
    
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN;
            p_valN(p<(0.05/200)) = 1;
            p_val = h;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'F2pheno', num2str(k), 'neo.FDR.eNKIs.png'],'png', 10)
            close(f)
        end
        
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN
            p_valN(p<=(0.01)) = 1;
            p_val = p_valN + p_valP;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'eNKI.middleage.duration.png'],'png', 10)
            close(f)
        end
    end
  
    %% eNKI sleep amount old age
    for k = 1
        variabl = eNKI_beh.sleep_amount;
        keep    = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl)),find(eNKI_beh.age>=50))
        vark    = (variabl(keep));
        agek    = eNKI_beh.age(keep);
        sexk    = eNKI_beh.gender(keep);
        gb      = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));
        [r p]   = partialcorri(eNKI200_CT(:,keep)',vark,[gb, sex_num, agek, agek.*agek, sex_num.*agek],'type','spearman', 'rows','complete');
        min(p)
        
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05);
        max(h)
 
        eNKI_amount_phenoCT_FDRo(:,k) = h;
        eNKI_amount_phenoCT_Po(:,k) = p;
        eNKI_amount_phenoCT_Ro(:,k) = r;
   
        [r p] = corr(temp,vark,'type','spearman')
    
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN;
            p_valN(p<(0.05/200)) = 1;
            p_val = h;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'F2pheno', num2str(k), 'neo.FDR.eNKIs.png'],'png', 10)
            close(f)
        end
        
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN
            p_valN(p<=(0.01)) = 1;
            p_val = p_valN + p_valP;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (r(i,:)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (r(i+100,:)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'eNKI.oldage.duration.png'],'png', 10)
            close(f)
        end
    end
  
    %% eNKI sleep amount youth vs adult
    for k = 1
        variabl = eNKI_beh.sleep_amount;
        keep    = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl)),find(eNKI_beh.age<50)) %,find(eNKI_beh.age<40))
        vark    = (variabl(keep));
        agek    = eNKI_beh.age(keep);
        sexk    = eNKI_beh.gender(keep);
        gb      = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));
        SL     = cell(size(keep));
        SL(find(agek<18)) = {'X'}; 
        SL(find(agek>=18)) = {'Y'};
        SLt = term(SL);
        Sx     = cell(size(keep));
        Sx(find(sex_num==1)) = {'X'}; 
        Sx(find(sex_num==2)) = {'Y'};
        Sxt = term(Sx);
      
        M       = 1 + term(gb) + term(agek) + term(Sx) + term(agek)*term(agek) + term(agek)*term(Sx) + term(vark) + term(SL) + term(vark)*term(SL);
        slm     = SurfStatLinMod(eNKI200_CT(:,keep)',M);
        slm     = SurfStatT(slm, (SLt.X.*vark)-(SLt.Y.*vark))
        p       = 1-tcdf(slm.t,slm.df)
     
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN
            p_valN(p<=(0.01)) = 1;
            p_val = p_valN + p_valP;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (slm.t(i)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (slm.t(i+100)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'eNKI.SD.y-a.png'],'png', 10)
            close(f)
        end
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN
            p_valN(p>=(0.99)) = 1;
            p_val = p_valN + p_valP;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (slm.t(i)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (slm.t(i+100)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'eNKI.SD.a-y.png'],'png', 10)
            close(f)
        end
    end
    
    %% eNKI sleep amount adult vs old
    for k = 1
        variabl = eNKI_beh.sleep_amount;
        keep    = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl)),find(eNKI_beh.age>=18)) %,find(eNKI_beh.age<40))
        vark    = (variabl(keep));
        agek    = eNKI_beh.age(keep);
        sexk    = eNKI_beh.gender(keep);
        gb      = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));
        SL     = cell(size(keep));
        SL(find(agek<50)) = {'X'}; 
        SL(find(agek>=50)) = {'Y'};
        SLt = term(SL);
        Sx     = cell(size(keep));
        Sx(find(sex_num==1)) = {'X'}; 
        Sx(find(sex_num==2)) = {'Y'};
        Sxt = term(Sx);
      
        M       = 1 + term(gb) + term(agek) + term(Sx) + term(agek)*term(agek) + term(agek)*term(Sx) + term(vark) + term(SL) + term(vark)*term(SL);
        slm     = SurfStatLinMod(eNKI200_CT(:,keep)',M);
        slm     = SurfStatT(slm, (SLt.X.*vark)-(SLt.Y.*vark))
        p       = 1-tcdf(slm.t,slm.df)
     
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN
            p_valN(p<=(0.01)) = 1;
            p_val = p_valN + p_valP;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (slm.t(i)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (slm.t(i+100)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'eNKI.SD.y-a.png'],'png', 10)
            close(f)
        end
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN
            p_valN(p>=(0.99)) = 1;
            p_val = p_valN + p_valP;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (slm.t(i)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (slm.t(i+100)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'eNKI.SD.e-a.png'],'png', 10)
            close(f)
        end
    end
    
    %% eNKI sleep amount young vs old
    for k = 1
        variabl = eNKI_beh.sleep_amount;
        keep    = mintersect(find(mean(eNKI200_CT)>0)', find(~isnan(variabl)),find((eNKI_beh.age<18)+(eNKI_beh.age>=50))) %,find(eNKI_beh.age<40))
        vark    = (variabl(keep));
        agek    = eNKI_beh.age(keep);
        sexk    = eNKI_beh.gender(keep);
        gb      = mean(eNKI200_CT(:,keep))';
        sex_num = grp2idx(eNKI_beh.gender(keep));
        SL      = cell(size(keep));
        SL(find(agek<50)) = {'X'}; 
        SL(find(agek>=50)) = {'Y'};
        SLt = term(SL);
        Sx     = cell(size(keep));
        Sx(find(sex_num==1)) = {'X'}; 
        Sx(find(sex_num==2)) = {'Y'};
        Sxt = term(Sx);
      
        M       = 1 + term(gb) + term(agek) + term(Sx) + term(agek)*term(agek) + term(agek)*term(Sx) + term(vark) + term(SL) + term(vark)*term(SL);
        slm     = SurfStatLinMod(eNKI200_CT(:,keep)',M);
        slm     = SurfStatT(slm, (SLt.X.*vark)-(SLt.Y.*vark))
        p       = 1-tcdf(slm.t,slm.df)
     
        
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN
            p_valN(p<=(0.01)) = 1;
            p_val = p_valN + p_valP;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (slm.t(i)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (slm.t(i+100)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'eNKI.SD.y-e.png'],'png', 10)
            close(f)
        end
        for plot_pvalue = 1
            p_valN = zeros(size(p));
            p_valP = p_valN
            p_valN(p>=(0.99)) = 1;
            p_val = p_valN + p_valP;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1)) = (slm.t(i)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200_17==i+1001)) = (slm.t(i+100)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(ctx2))
            SurfStatColLim([-0.1 0.1])
            exportfigbo(f,[RPATH, 'eNKI.SD.e-y.png'],'png', 10)
            close(f)
        end
    end
end
