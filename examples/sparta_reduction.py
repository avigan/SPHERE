import sphere.SPARTA as SPARTA

####################################################@
# full reduction
#

#%% init reduction
reduction = SPARTA.Reduction('/Users/avigan/data/sphere-test-target/SPARTA/', log_level='info')

###############################################################################
# It is possible to provide a default JSON configuration file to set some (or #
# all) of the reduction parameters to a default value different from the ones #
# hard-coded in the sphere package. This is done with the keyword:            #
#   user_config='... path to the file ...'                                    #
# The increasing priority for setting reduction parameters is the following:  #
#   0- default values hard-coded in the sphere package                        #
#   1- values contained in the file pointed by the user_config keyword, if a  #
#      file path is provided and exists                                       #
#   2- values contained in a reduction_config.json file left in the reduction #
#      directory by a previous reduction                                      #
#   3- values manually set by the user (see examples below)                   #
###############################################################################

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
