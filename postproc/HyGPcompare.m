function varargout = list_selection(varargin)
    % LIST_SELECTION M-file for list_selection.fig
    %      LIST_SELECTION, by itself, creates a new LIST_SELECTION or raises the existing
    %      singleton*.
    %
    %      H = LIST_SELECTION returns the handle to a new LIST_SELECTION or the handle to
    %      the existing singleton*.
    %
    %      LIST_SELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in LIST_SELECTION.M with the given input arguments.
    %
    %      LIST_SELECTION('Property','Value',...) creates a new LIST_SELECTION or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before list_selection_OpeningFunction gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to list_selection_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Copyright 2002-2003 The MathWorks, Inc.

    % Edit the above text to modify the response to help list_selection

    % Last Modified by GUIDE v2.5 24-Feb-2013 19:28:14

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @list_selection_OpeningFcn, ...
                   'gui_OutputFcn',  @list_selection_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT


% --- Executes just before list_selection is made visible.
function list_selection_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to list_selection (see VARARGIN)

    % initialize global variables
    handles.directory_list = '';    
    handles.experiment_list = '';
    handles.last_directory = '';
    % define output for list_selection
    handles.output = hObject;%handles.experiment_list; %'ciao';%
    
    handles.run_no = -1;
    
    % Update handles structure
    guidata(hObject, handles);

    
    % UIWAIT makes list_selection wait for user response (see UIRESUME)
    %uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = list_selection_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;
   
    % now delete the figure
    %delete(handles.figure1);
    
    
% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
    % hObject    handle to listbox1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from listbox1


% listbox1 CONSTRUCTOR --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to listbox1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: listbox controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    set(hObject, 'String','');


% Button REMOVE
function pushbutton2_Callback(hObject, eventdata, handles)
    %disp('REMOVE clicked')
    %get the values for the selected file names
    doomed_item = get(handles.listbox1,'Value')
    
    %is there is nothing to delete, nothing happens
    if (isempty(doomed_item) == 1 || doomed_item(1) == 0 || isempty(handles.directory_list))
        return
    end
 
    %erases the contents of highlighted item in data array
    handles.directory_list(doomed_item) = [];
    handles.experiment_list(doomed_item) = [];
 
    %updates the gui, erasing the selected item from the listbox
    set(handles.listbox1,'String',handles.experiment_list);
 
    %moves the highlighted item to an appropiate value or else will get error
    if doomed_item(end) > length(handles.experiment_list)
        set(handles.listbox1,'Value',length(handles.experiment_list));
    end
 
    % Update handles structure
    guidata(hObject, handles);
    
    

% Button ADD
function pushbutton3_Callback(hObject, eventdata, handles)
    % call the function to select experiment and directory
    [experiment_name, directory_name] = get_experiment(handles.last_directory);
    % experiment_name
    % if experiment selection is cancelled, directory_name should be zero
    % and nothing should happen
    if (directory_name == 0)
        %disp('Operation aborted.')
        return   
    end
    
    % add the most recent data file selected to the cell containing
    % all the data file names
    handles.experiment_list{end+1} = experiment_name;
    handles.directory_list{end+1} = directory_name;
   
    %updates the gui to display all filenames in the listbox
    set(handles.listbox1,'String', handles.experiment_list);  
    
    %make sure first file is always selected so it doesn't go out of range
    %the GUI will break if this value is out of range
    set(handles.listbox1,'Value',1);
    
    % to verify uncomment 
    %experiment_list
    %handles.experiment_list{:}
    %handles.directory_list{:}
    
    % update last visited directory
    handles.last_directory = directory_name;
    
    % Update handles structure
    guidata(hObject, handles);
    
% Button TRAIN data
function pushbutton1_Callback(hObject, eventdata, handles)
    % run GPBoxplot_fun function to visualize results
    %disp('TRAIN data pressed')
    sought_file = 'archives_best.txt';
    if check_file_existence(handles.experiment_list, handles.directory_list, sought_file)
        best_boxplot_fun(handles.experiment_list, handles.directory_list,'train');
    end
    % plot boxplots for the entire archives
    %archive_boxplot_fun(handles.experiment_list, handles.directory_list);


% Button TEST data
function pushbutton5_Callback(hObject, eventdata, handles)
    % hObject    handle to pushbutton5 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    sought_file = 'archives_best_TEST.txt';
    if check_file_existence(handles.experiment_list, handles.directory_list, sought_file);
        best_boxplot_fun(handles.experiment_list, handles.directory_list,'test');
    end

% Button SAVE in LATEX format
function pushbutton6_Callback(hObject, eventdata, handles)
    % hObject    handle to pushbutton6 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.directory_list
    print2file_best_latex(handles.experiment_list, handles.directory_list);


% Button EXIT
function pushbutton4_Callback(hObject, eventdata, handles)
    % hObject    handle to pushbutton4 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    close all;


% Button SIZE
function pushbutton7_Callback(hObject, eventdata, handles)
    % hObject    handle to pushbutton7 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    %disp('Size pressed')
    plot_size_dynamics(handles.experiment_list, handles.directory_list);
    


% Button DEPTH
function pushbutton8_Callback(hObject, eventdata, handles)
    % hObject    handle to pushbutton8 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    plot_depth_dynamics(handles.experiment_list, handles.directory_list);
    


%Button Plot runtime
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    plot_runtime(handles.experiment_list, handles.directory_list)


% Button Adaptive genetic operators: Rates and Performances
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    %get the value of the selected experiment
    selected_exp_no = get(handles.listbox1,'Value')

    if (handles.run_no == -1)
        plot_average_genetic_op_performances(handles.experiment_list, handles.directory_list, selected_exp_no)
    else
        plot_genetic_op_performances(handles.experiment_list, handles.directory_list, selected_exp_no, handles.run_no)
    end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
    % set run_no
    disp(get(hObject,'String'))
    handles.run_no = get(hObject,'String')
    guidata(hObject, handles);
    
% --- Executes when selected object is changed in uipanel7.
function uipanel7_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel7 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
    switch get(eventdata.NewValue,'Tag') % Get Tag of selected object
        case 'radiobutton1'
            %execute this code when radiobutton1 is selected
            %1) deactivate edit region
            % 2) set run_no to -1
            handles.run_no = -1
        case 'radiobutton2'
            %execute this code when radiobutton2 is selected
            % 1) activate edit region
            % 2) set run_no equal to edited value
            handles.run_no = get(handles.edit1,'String');
        otherwise
        % Code for when there is no match.
    end
    %updates the handles structure
    guidata(hObject, handles);


% average performances and rates
% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton1.
function radiobutton1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    %disp('average')

% single run performances
% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton2.
function radiobutton2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
     %disp('single run')



    

% button "average error"
% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    r = plot_average_training_error(handles.experiment_list, handles.directory_list)
