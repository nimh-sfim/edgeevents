edgesort
==============================

sort the edges by different feats.

Project Organization
--------------------

    .
    ├── AUTHORS.md
    ├── LICENSE
    ├── README.md
    ├── config
    │   ├── config_hcp_sch200_1.m
    │   ├── config_hcp_sch200_2.m
    │   ├── config_hcp_sch400_1.m
    │   └── config_template.m
    ├── data
    │   ├── external
    │   │   ├── hcp.gradients.dscalar.nii
    │   │   ├── hcp_grad_sch400-yeo17.pscalar.nii
    │   │   ├── hpc_grad_sch200-yeo17.pscalar.nii
    │   │   └── nodes_2_canon.mat
    │   ├── interim
    │   │   ├── hcp352_nusregts_FIX2phys_schaefer200 -> /Users/faskowitzji/joshstuff/data/hcp352_nusregts_FIX2phys_schaefer200/
    │   │   ├── hcp352_nusregts_FIX2phys_schaefer400 -> /Users/faskowitzji/joshstuff/data/hcp352_nusregts_FIX2phys_schaefer400/
    │   │   ├── sub_list.txt
    │   │   ├── sub_list_subset1.txt
    │   │   ├── sub_list_subset2.txt
    │   │   ├── ts_autocorr.mat
    │   │   └── ts_variability.mat
    │   ├── processed
    │   │   ├── CARD
    │   │   ├── REST1_LR
    │   │   ├── REST1_RL
    │   │   ├── spike_hist_clust.mat
    │   │   ├── spk_conn_avg_sch200.mat
    │   │   ├── spk_corrWphys_sch200.mat
    │   │   ├── spk_corrwmot_sch200.mat
    │   │   ├── spk_exceedsurrthr_sch200.mat
    │   │   ├── spk_hist_sch200.mat
    │   │   ├── spk_rsspeaks_sch200.mat
    │   │   ├── spk_rsspeaksclust_sch200.mat
    │   │   ├── spk_rsspeakssim_sch200.mat
    │   │   ├── spk_slope_sch200.mat
    │   │   ├── spk_systrans_sch200.mat
    │   │   ├── spk_thr_dat_sch200.mat
    │   │   ├── spk_var_sch200.mat
    │   │   ├── surrogate2_sch200_2.25_spk.mat
    │   │   ├── surrogate2_sch200_2.25_spkcount.mat
    │   │   ├── surrogate2_sch200_2.25_spkrss.mat
    │   │   ├── surrogate3_sch200_2.25_spk.mat
    │   │   ├── surrogate3_sch200_2.25_spkcount.mat
    │   │   ├── surrogate3_sch200_2.25_spkrss.mat
    │   │   ├── surrogate_sch200_2.25_boothist.mat
    │   │   ├── surrogate_sch200_2.25_spk.mat
    │   │   ├── surrogate_sch200_2.25_spkcount.mat
    │   │   └── surrogate_sch200_2.25_spkrss.mat
    │   └── raw
    ├── docs
    ├── reports
    │   ├── etc
    │   │   ├── sub_list.txt
    │   │   ├── sub_list_subset1.txt
    │   │   └── sub_list_subset2.txt
    │   └── figures
    │       ├── figA
    │       ├── figB
    │       ├── figC
    │       ├── figD
    │       ├── figE
    │       ├── figInfo
    │       ├── hub_retention_pic.png
    │       └── supp
    └── src
        ├── ext
        │   ├── IPN_ccc.m
        │   ├── IPN_ssd.m
        │   ├── Modified_MannKendall_test.m
        │   ├── brainSync.m
        │   ├── centroid_extraction_sphere.m
        │   ├── consensus_und2.m
        │   ├── f_CCC.m
        │   ├── generate_phase_surrogates.m
        │   ├── ktaub.m
        │   └── simulate_BOLD_timecourse_func_v3.m
        ├── func
        │   ├── count_spks.m
        │   ├── ets_2_ts.m
        │   ├── get_acf_hwhm.m
        │   ├── get_blocky.m
        │   ├── get_blocky_ets.m
        │   ├── get_contact_times.m
        │   ├── get_entdrops.m
        │   ├── get_ets.m
        │   ├── get_ets_node_inds.m
        │   ├── get_ets_states.m
        │   ├── get_hub_score_wei_und.m
        │   ├── get_rss.m
        │   ├── lins_ccc.m
        │   ├── load_hcp_alldata.m
        │   ├── load_hcp_card.m
        │   ├── load_hcp_fd.m
        │   ├── load_hcp_regressors.m
        │   ├── make_tinda.m
        │   ├── make_trans_dens.m
        │   ├── make_trans_dens_multistate.m
        │   ├── make_trans_prob.m
        │   ├── myent_simp.m
        │   ├── node_versatility.m
        │   ├── norm_bin_model.m
        │   ├── parc_plot_wcolorbar.m
        │   ├── quick_pspec.m
        │   ├── rmssd.m
        │   ├── spk_lenmat.m
        │   └── spk_pks_sort.m
        ├── proc
        ├── scripts
        │   ├── anlyz_rel_var.m
        │   ├── clust_spike_hist.m
        │   ├── corr_rss_w_card.m
        │   ├── corr_rss_w_mot.asv
        │   ├── corr_rss_w_mot.m
        │   ├── corr_rss_w_phys.m
        │   ├── longest_edges_analy.m
        │   ├── make_autocorr_dat.m
        │   ├── make_sim_spike_lens.m
        │   ├── make_sim_spike_lens2.m
        │   ├── make_sim_spike_lens3.m
        │   ├── make_spike_conn.m
        │   ├── make_spike_hist.asv
        │   ├── make_spike_hubs.m
        │   ├── make_spike_length_data.m
        │   ├── make_sys_trans_data.m
        │   ├── make_variability_dat.m
        │   ├── pick_initial_spk_thr.m
        │   ├── simplestuff.m
        │   ├── spike_rss_analysis.m
        │   ├── spike_rss_analysis2.m
        │   ├── spike_rss_analysis3.m
        │   └── spike_slope_analy.m
        └── viz
            ├── basic_info_plots.m
            └── viz_conn_glassbrain.m

    29 directories, 107 files

