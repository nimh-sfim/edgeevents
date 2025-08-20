%% PROJECT CONFIG

PROJ_DIR = pwd ; % put base path to project her
cd(PROJ_DIR)

%% add to the path

projPathDirs = {
    'src'
    'data'
    'bin'
} ;

for idx=1:length(projPathDirs)
    addpath(genpath(strcat(PROJ_DIR,'/',projPathDirs{idx})))
end

clear projPathDirs

%% setup global vars

OUTSTR = 'run1' ;

%% make output directory vars

DATADIR = strcat(PROJ_DIR , '/data/') ;
DD.INTERM = strcat(DATADIR, '/interim/' ) ;
DD.PROC = strcat(DATADIR, '/processed/' ) ;
DD.RAW = strcat(DATADIR, '/raw/' ) ;
