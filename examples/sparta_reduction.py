import sphere.SPARTA as SPARTA

####################################################@
# full reduction
#

#%% init reduction
reduction = SPARTA.Reduction('/Users/avigan/data/sphere-test-target/SPARTA/', log_level='info')

#%% configuration
reduction.config['misc_plot'] = True
reduction.show_config()

#%% reduction
reduction.full_reduction()

####################################################
# manual reduction
#

#%% init reduction
reduction = SPARTA.Reduction('/Users/avigan/data/sphere-test-target/SPARTA/', log_level='info')

#%% sorting
reduction.sort_files()

#%% science processing
reduction.sph_sparta_dtts(plot=False)
reduction.sph_sparta_wfs_parameters()
reduction.sph_sparta_atmospheric_parameters()
reduction.sph_query_databases(timeout=2)
reduction.sph_sparta_plot()

#%% cleaning
reduction.sph_sparta_clean(delete_raw=False, delete_products=False)
