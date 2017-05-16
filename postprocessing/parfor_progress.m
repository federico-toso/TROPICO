function completion_percentage = parfor_progress (total_evaluations)

if exist('temp_progress_file.mat','file')~=2
    completed_eval = 1;
    save('temp_progress_file.mat','completed_eval');
    completion_percentage=completed_eval/total_evaluations;
else
    load('temp_progress_file.mat','completed_eval');
    if completed_eval>=total_evaluations-1
        delete('temp_progress_file.mat')
        completion_percentage=1;
    else
        completed_eval=completed_eval+1;
        save('temp_progress_file.mat','completed_eval');
        completion_percentage=completed_eval/total_evaluations;
    end
end