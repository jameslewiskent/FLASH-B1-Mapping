% read 

%fname = 'D:\Data2\4DFLow\adam\raw\meas_MID01428_FID51521_785A_4DFlow_ePAT3_retro_2_2mm_pat3.dat';
fname = '/media/data6t/disk1/aaron/flow/adam/2081208_raw/meas_MID01428_FID51521_785A_4DFlow_ePAT3_retro_2_2mm_pat3.dat';

    twix = mapVBVD(fname,'removeOS','ignoreSeg',true);  %,'averagereps'
    nset = length(twix);
    
    if(nset>1)
        reftwix = twix(1:nset-1);
        twix = twix{nset};
    end
    
    twix.image.flagIgnoreSeg = true;
    sz = twix.image.dataSize;
    kdata = twix.image('');
    
    % Readout partial fourier 
    kc = twix.image.centerCol(1);
    intended_ro = (sz(1) - kc/2)*2+1;
    ro_add = intended_ro - sz(1);
    kdata(intended_ro,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0;
    kdata = circshift(kdata,[ro_add 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
    
    %% PE partial fourier 
    kc = twix.image.centerLin(1);
    intended_pe = kc*2;
    pe_add = intended_pe - sz(3);
    kdata(:,:,intended_pe,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0;
    %kdata = circshift(kdata,[0 0 pe_add 0 0 0 0 0 0 0 0 0 0 0 0 0]);
    %%
    RO = 1;
    CH = 2;
    PE1 = 3;
    PE2 = 4;
    SLC=5;
    AVE=6;
    PHS=7;
    ECO = 8;
    REP = 9;
    SET = 10;
    OTHER5 = 11;
     % twix dims are (1)ro,ch,pe1,pe2,(5),(6),(7)phs,..(10)set
    % reorder to bart dims ro(1),pe1(2),pe2(3),ch(4),maps(5),(6),(7),rep(8),(9),(10),(11 cardiac),(12
    % resp)
    bart_order = zeros(1,11);
    bart_order(1) = RO;
    bart_order(2) = PE1;
    bart_order(3) = PE2;
    bart_order(4) = CH;
    bart_order(5) = OTHER5; % maps
    bart_order(6) = SLC;
    bart_order(7) = AVE;
    bart_order(8) = REP;
    bart_order(9) = SET;
    bart_order(10) = PHS;
    bart_order(11) = ECO;
    %RO PE1 PE2 CH OTHER5 SLC AVE
    kdata = permute(kdata,bart_order);
    sz = size(kdata);
    fname_cfl = strrep(fname,'.dat','_cfl');
    writecfl(fname_cfl,kdata);
    %kdata = readcfl(fname_cfl);
    
    fname1_cfl = strrep(fname,'.dat','_cfl1');
    writecfl(fname1_cfl,kdata(:,:,:,:,:,:,:,:,1,1:25));
    
    %% write ipat reference data to disk too
    fname_patref_cfl = '';
    if(isfield(twix,'refscan'))
        pat_ref= twix.refscan('');
        pat_ref= permute(pat_ref,bart_order);

        sz_ref = size(pat_ref);

        sz_diff = sz(1:4)-sz_ref(1:4);

        pat_ref(sz(1),sz(2),sz(3),sz(4)) = 0;
        pat_ref = circshift(pat_ref,round(sz_diff/2));

        fname_patref_cfl = strrep(fname,'.dat','patref_cfl');
        writecfl(fname_patref_cfl,pat_ref);
        
        
    end
    
    
    %% Generate smaps
    this_dir = pwd;
    cd '/home/ahess/bart';
    startup;
    fname_smaps_clf = strrep(fname,'.dat','smaps_cfl');
    bart(['ecalib -m1 -k5:5:3 ', fname_patref_cfl, ' ',fname_smaps_clf] );   
    
    %% CG SENSE recon
    my_cmd = sprintf('pics -d 3 -r 0.1 %s %s',fname1_cfl, fname_smaps_clf);
    im_cg=squeeze(bart(my_cmd));  % CG
    