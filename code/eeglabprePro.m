% this script cannot be used immediately

% EEGLAB history file generated on the 23-Oct-2019
% ------------------------------------------------
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename','crq raw selected.set','filepath','D:\\Lab\\Hong_Lab\\processing\\new\\crq\\rest\\');
% at the fast-reference here
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve',2,'study',0); 
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
EEG = pop_eegfiltnew(EEG, 'locutoff',49,'hicutoff',51,'revfilt',1,'plotfreqz',1);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname','crq rest REST selected filtered','savenew','D:\\Lab\\Hong_Lab\\processing\\new\\crq\\rest\\crq rest REST selected filtered.set','comments',strvcat('Re-referencing to REST',' ','50Hz notched'),'overwrite','on','gui','off'); 
EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',40,'plotfreqz',1);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'savemode','resave');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
EEG = eeg_checkset( EEG );
EEG = pop_runica(EEG, 'icatype', 'binica', 'extended',1);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'savemode','resave');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw;
