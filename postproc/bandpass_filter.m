%%% band-pass filtering script
%%%
%%%

% GET THE NAME OF THE EXPERIMENT AND ITS DIRECTORY (function)
%first_dir = ['E:' sep 'Research_Leeds' sep 'Linux' sep 'code' sep 'C++' sep 'experiments' sep 'structuralGP_OLD' sep 'output' sep]
first_dir=['D:\relocated_users\relocated_umba\Documents\OneDrive - University of Leeds\Research\2020_Research_ML_AI_Vassili\output']
[experiment, directory_name] = get_experiment(first_dir);
%if experiment selection is cancelled, directory_name should be zero
%and nothing should happen
if (directory_name == 0)
  disp('Quitted normally.')
  return
end
[NVAR, NFITCASES, R, CONSTR0, CONSTR1]=read_INPUT_FILE(directory_name); 

% filter data with a band-pass filter
x_norm=R(:,1);
fs=3.334000000000000e-03;
bandpass(x_norm,[30 50],fs)